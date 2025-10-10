// negative_cycle_batch.cpp
#include "negative_cycle_batch.h"

static inline long long key64_pair(int a, int b) {
    if (a > b) std::swap(a,b);
    return (static_cast<long long>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint32_t>(b);
}

// Active stream for triangle-stage callbacks (thread-local, safe for parallel root/node separation)
static thread_local NegativeCycleBatch* g_tbb_active = nullptr;

// Strong definitions override weak ones in triangle_bucket_batch.cpp
extern "C" {

void TBB_on_emit(int edge_id, double used_density) {
    if (g_tbb_active) g_tbb_active->on_emit_(edge_id, used_density);
}

void TBB_on_accept(int edge_id, double density) {
    if (g_tbb_active) g_tbb_active->on_accept_(edge_id, density);
}

int TBB_budget_override(int base) {
    return g_tbb_active ? g_tbb_active->override_budget_(base) : base;
}

} // extern "C"

NegativeCycleBatch SignedGraphForMIP::open_negative_cycle_stream(bool cover, bool use_triangle_order) const {
    return NegativeCycleBatch(*this, cover, use_triangle_order);
}

NegativeCycleBatch::NegativeCycleBatch(const SignedGraphForMIP& G, bool cover, bool use_triangle_order)
    : G_(G),
      cover_(cover),
      use_tri_order_(use_triangle_order) {
    if (use_tri_order_) {
        try { neg_tri_vert_ = G_.negative_triangle_count_per_vertex(); }
        catch(...) { neg_tri_vert_.assign(igraph_vcount(&G_.g), 0); }
    } else {
        neg_tri_vert_.clear();
    }
    build_initial_state_();
}

int  NegativeCycleBatch::total_cycles_emitted() const { return static_cast<int>(total_found_); }
int  NegativeCycleBatch::batches_emitted()      const { return batches_emitted_; }

void NegativeCycleBatch::build_mask_for_batch_() {
    if (!saved_weights_pos_init_) return;
    for (long pe = 0; pe < (long)pos2full_eid_.size(); ++pe) {
        igraph_integer_t fe = pos2full_eid_[pe];
        if (use_lp_weights_ && fe < (igraph_integer_t)lp_score_full_.size()) {
            const double base = base_pos_[(size_t)fe];
            const double hint = lp_score_full_[(size_t)fe];
            VECTOR(saved_weights_pos_)[pe] = std::max(1e-12, lp_alpha_ * base / (1.0 + lp_beta_ * hint));
        } else {
            VECTOR(saved_weights_pos_)[pe] = base_pos_[(size_t)fe];
        }
    }
    if (saved_weights_init_) {
        for (long i = 0; i < (long)ecount_; ++i) VECTOR(saved_weights_)[i] = 0.0;
    }
    std::fill(reuse_accum_.begin(), reuse_accum_.end(), 0.0);
    alt_path_bump_ = med_base_pos_; // keep your prior behavior
}

void NegativeCycleBatch::set_lp_scores_full_edges(const std::vector<double>& s, double alpha, double beta) {
    lp_score_full_ = s; use_lp_weights_ = true; lp_alpha_ = alpha; lp_beta_ = beta;
    if (saved_weights_pos_init_) {
        for (long pe = 0; pe < (long)pos2full_eid_.size(); ++pe) {
            igraph_integer_t fe = pos2full_eid_[pe];
            const double base = base_pos_[(size_t)fe];
            const double hint = (fe < (igraph_integer_t)s.size() ? s[(size_t)fe] : 0.0);
            VECTOR(saved_weights_pos_)[pe] = std::max(1e-12, lp_alpha_ * base / (1.0 + lp_beta_ * hint));
        }
    }
}

void NegativeCycleBatch::build_initial_state_() {
    vcount_ = igraph_vcount(&G_.g);
    ecount_ = igraph_ecount(&G_.g);

    base_pos_.assign((size_t)ecount_, 1.0);
    reuse_accum_.assign((size_t)ecount_, 0.0);
    neg_edges_.clear(); neg_edges_.reserve((size_t)ecount_);
    neg_deg_.assign((size_t)vcount_, 0);
    pos_deg_.assign((size_t)vcount_, 0);

    std::vector<double> pos_bases; pos_bases.reserve((size_t)ecount_);
    for (igraph_integer_t eid = 0; eid < ecount_; ++eid) {
        igraph_integer_t u, v; igraph_edge(&G_.g, eid, &u, &v);
        const double w = G_.switched_weights[eid];
        if (w > 0.0) {
            base_pos_[(size_t)eid] = std::max(1e-12, std::fabs(w));
            pos_bases.push_back(base_pos_[(size_t)eid]);
            ++pos_deg_[(size_t)u]; ++pos_deg_[(size_t)v];
        } else if (w < 0.0) {
            neg_edges_.emplace_back((int)u, (int)v);
            ++neg_deg_[(size_t)u]; ++neg_deg_[(size_t)v];
        } else {
            base_pos_[(size_t)eid] = 1e-12;
            pos_bases.push_back(1e-12);
            ++pos_deg_[(size_t)u]; ++pos_deg_[(size_t)v];
        }
    }

    tri_cap_per_v_.assign((size_t)vcount_, tri_cap_per_vertex_);
    for (int v = 0; v < (int)vcount_; ++v) {
        int bump = (pos_deg_[(size_t)v] + 63) / 64;
        tri_cap_per_v_[(size_t)v] = std::min(64, tri_cap_per_vertex_ + bump);
    }

    if (!pos_bases.empty()) {
        const size_t mid = pos_bases.size()/2;
        std::nth_element(pos_bases.begin(), pos_bases.begin()+mid, pos_bases.end());
        med_base_pos_ = std::max(1e-12, pos_bases[mid]);
    } else {
        med_base_pos_ = 1.0;
    }

    // Build positive-only graph and mappings
    std::vector<igraph_integer_t> pos_edges; pos_edges.reserve((size_t)ecount_ * 2);
    full2pos_eid_.assign((size_t)ecount_, -1);
    pos2full_eid_.clear(); pos2full_eid_.reserve((size_t)ecount_);

    for (igraph_integer_t eid = 0; eid < ecount_; ++eid) {
        if (!edge_is_pos(eid)) continue;
        igraph_integer_t u, v; igraph_edge(&G_.g, eid, &u, &v);
        full2pos_eid_[(size_t)eid] = (igraph_integer_t)pos2full_eid_.size();
        pos2full_eid_.push_back(eid);
        pos_edges.push_back(u); pos_edges.push_back(v);
    }

    if (g_pos_built_) igraph_destroy(&g_pos_);
    igraph_vector_int_t edges_vec; igraph_vector_int_init(&edges_vec, (long)pos_edges.size());
    for (long i = 0; i < (long)pos_edges.size(); ++i) VECTOR(edges_vec)[i] = (igraph_integer_t)pos_edges[i];
    igraph_create(&g_pos_, &edges_vec, vcount_, IGRAPH_UNDIRECTED);
    igraph_vector_int_destroy(&edges_vec);
    g_pos_built_ = true;

    // Initialize working masks
    if (saved_weights_pos_init_) igraph_vector_destroy(&saved_weights_pos_);
    igraph_vector_init(&saved_weights_pos_, (long)pos2full_eid_.size());
    for (long pe = 0; pe < (long)pos2full_eid_.size(); ++pe) {
        igraph_integer_t fe = pos2full_eid_[pe];
        VECTOR(saved_weights_pos_)[pe] = base_pos_[(size_t)fe];
    }
    saved_weights_pos_init_ = true;

    if (!saved_weights_init_) {
        igraph_vector_init(&saved_weights_, (long)ecount_);
        for (long i = 0; i < (long)ecount_; ++i) VECTOR(saved_weights_)[i] = 0.0;
        saved_weights_init_ = true;
    }

    tri_used_per_vertex_.assign((size_t)vcount_, 0);
    neg_edge_covered_.assign((size_t)ecount_, 0);

    // Light stats
    int negE = (int)neg_edges_.size();
    double avgPosDeg = 0.0; int posV = 0; int hist[8] = {0};
    for (int v = 0; v < (int)vcount_; ++v) {
        const int d = pos_deg_[(size_t)v];
        if (d > 0) { avgPosDeg += d; ++posV; }
        int b = (d<=2?0:d<=4?1:d<=8?2:d<=16?3:d<=32?4:d<=64?5:d<=128?6:7);
        hist[b]++;
    }
    if (posV>0) avgPosDeg /= (double)posV;
    K_tri_per_neg_ = std::max(3, std::min(8, (int)std::sqrt(std::max(1, (int)avgPosDeg))));
    std::cout << "[TRI-DEG] negE=" << negE
              << ", V=" << vcount_ << ", E=" << ecount_
              << ", E[minDeg+]=" << std::fixed << std::setprecision(2) << avgPosDeg
              << ", hist=[0-2]=" << hist[0]
              << " [3-4]=" << hist[1]
              << " [5-8]=" << hist[2]
              << " [9-16]=" << hist[3]
              << " [17-32]=" << hist[4]
              << " [33-64]=" << hist[5]
              << " [65-128]=" << hist[6]
              << " [129+]=" << hist[7]
              << std::endl;
}

void NegativeCycleBatch::build_pos_adj_and_index_(
    TriangleBucketBatch::PosAdj& pos_adj,
    TriangleBucketBatch::EdgeIndex& edge_index) const
{
    pos_adj.assign((size_t)vcount_, {});
    // Map all edges to full-eid (so neg anchors resolve)
    for (igraph_integer_t eid = 0; eid < ecount_; ++eid) {
        igraph_integer_t u, v; igraph_edge(&G_.g, eid, &u, &v);
        edge_index.emplace(key64_pair((int)u,(int)v), (int)eid);
        // only positive edges go to adjacency
        if (edge_is_pos(eid)) {
            pos_adj[(size_t)u].push_back((int)v);
            pos_adj[(size_t)v].push_back((int)u);
        }
    }
}

int NegativeCycleBatch::override_budget_(int base) const {
    // Hook for annealing B_tri; keep passthrough for now.
    return base;
}

void NegativeCycleBatch::on_emit_(int full_eid, double used_density) {
    // WHY: steer Dijkstra away from overused positive edges within this batch
    if (full_eid < 0 || full_eid >= (int)full2pos_eid_.size()) return;
    igraph_integer_t pe = full2pos_eid_[(size_t)full_eid];
    if (pe < 0) return;
    // Tiny bump proportional to median and used_density
    const double bump = 0.05 * med_base_pos_ * std::max(0.0, used_density);
    VECTOR(saved_weights_pos_)[pe] = std::max(1e-12, VECTOR(saved_weights_pos_)[pe] + bump);
}

void NegativeCycleBatch::on_accept_(int full_eid, double /*density*/) {
    // WHY: persistent cross-batch drift on edges that were used in accepted triangles
    bump_cross_batch_(full_eid, /*|C|=*/3);
}

bool NegativeCycleBatch::run_triangle_first_batch_(
    std::vector<NegativeCycle>& out,
    std::vector<int>& covered_neg_eids)
{
    covered_neg_eids.clear();

    // Build TriangleBucketBatch inputs
    TriangleBucketBatch::PosAdj pos_adj;
    TriangleBucketBatch::EdgeIndex edge_index;
    build_pos_adj_and_index_(pos_adj, edge_index);

    // Neg edges under current switching
    std::vector<std::pair<int,int>> neg_edges_pairs; neg_edges_pairs.reserve(neg_edges_.size());
    for (const auto& e : neg_edges_) neg_edges_pairs.emplace_back(e.first, e.second);

    // Params
    TriangleBucketBatch::Params P;
    P.K_tri_per_neg  = std::max(1, K_tri_per_neg_);
    P.cap_per_vertex = std::max(1, tri_cap_per_vertex_);
    // generous default; let TBB_budget_override clamp it
    P.B_tri          = std::max(1, (int)neg_edges_pairs.size() * P.K_tri_per_neg);

    TriangleBucketBatch tbb(neg_edges_pairs, pos_adj, edge_index, P);

    // Scorer: primary = sum(1/ω'(uw), 1/ω'(wv)) using current working weights; secondary free
    auto scorer = [&](TriangleBucketBatch::Candidate& c) {
        auto pe_uw = (c.pos_eid_uw >= 0 && c.pos_eid_uw < (int)full2pos_eid_.size()) ? full2pos_eid_[(size_t)c.pos_eid_uw] : -1;
        auto pe_wv = (c.pos_eid_wv >= 0 && c.pos_eid_wv < (int)full2pos_eid_.size()) ? full2pos_eid_[(size_t)c.pos_eid_wv] : -1;
        double wuw = (pe_uw >= 0 ? VECTOR(saved_weights_pos_)[pe_uw] : 1e12);
        double wwv = (pe_wv >= 0 ? VECTOR(saved_weights_pos_)[pe_wv] : 1e12);
        c.score_primary   = (1.0 / std::max(1e-12, wuw)) + (1.0 / std::max(1e-12, wwv));
        c.score_secondary = 0.0;
        c.viol = 0.0; c.phi = 0.0;
    };

    // Build buckets and select with active bridges
    tbb.build_buckets(scorer);
    NegativeCycleBatch* prev = g_tbb_active;
    g_tbb_active = this;
    const auto& selected = tbb.select(covered_neg_eids);
    g_tbb_active = prev;

    // Emit accepted triangles as NegativeCycle (two positive edges path)
    for (const auto& c : selected) {
        std::vector<Edge> cyc_path;
        cyc_path.emplace_back(c.u, c.w);
        cyc_path.emplace_back(c.w, c.v);
        out.emplace_back(Edge{c.u, c.v}, std::move(cyc_path));
        ++total_found_;

        // Update per-vertex caps & coverage
        ++tri_used_per_vertex_[c.u];
        ++tri_used_per_vertex_[c.v];
        ++tri_used_per_vertex_[c.w];
        igraph_integer_t neg_eid; igraph_get_eid(&G_.g, &neg_eid, c.u, c.v, /*directed=*/0, /*error=*/1);
        if (neg_eid >= 0 && neg_eid < ecount_) neg_edge_covered_[(size_t)neg_eid] = 1;
    }

    return !selected.empty();
}

bool NegativeCycleBatch::next(std::vector<NegativeCycle>& out) {
    out.clear();
    if (finished_) return false;

    using clock = std::chrono::steady_clock;
    long long ms_dijkstra = 0, ms_check = 0, ms_emit = 0, ms_misc = 0;
    size_t neg_edges_scanned = 0, cycles_emitted_now = 0, triangles_emitted_now = 0;
    long long path_nodes_scanned = 0, pos_edges_on_paths = 0;

    build_mask_for_batch_();

    const size_t before_total = total_found_;

    // === Triangle-first (bucketed) ===
    std::vector<NegativeCycle> triangles;
    std::vector<int> covered_neg_eids;   // full-eids of anchors for which a bucket is nonempty (per TBB spec)
    if (run_triangle_first_batch_(triangles, covered_neg_eids)) {
        triangles_emitted_now += triangles.size();
        for (auto& C : triangles) out.push_back(std::move(C));
    }

    // Build uncovered list for SP (remove all nonempty bucket anchors)
    std::unordered_set<int> covered_set(covered_neg_eids.begin(), covered_neg_eids.end());
    std::vector<Edge> neg_edges_uncov; neg_edges_uncov.reserve(neg_edges_.size());
    for (const auto& e_neg : neg_edges_) {
        igraph_integer_t eid; igraph_get_eid(&G_.g, &eid, e_neg.first, e_neg.second, 0, 1);
        if (covered_set.find((int)eid) == covered_set.end()) neg_edges_uncov.push_back(e_neg);
    }

    std::vector<Edge> new_disconnected; new_disconnected.reserve(neg_edges_uncov.size());

    // === SP-based generation on uncovered anchors ===
    for (const auto& e_neg : neg_edges_uncov) {
        const int uu = e_neg.first, vv = e_neg.second;

        int accepted_here = 0;
        std::unordered_set<uint64_t> seen_paths;
        for (int rep = 0; rep < K_sp_per_neg_; ++rep) {
            igraph_vector_int_t path; igraph_vector_int_init(&path, 0);
            const auto t_dij0 = clock::now();
            igraph_get_shortest_path_dijkstra(&g_pos_, &path, nullptr, uu, vv, &saved_weights_pos_, IGRAPH_ALL);
            if (igraph_vector_int_size(&path) < 2) {
                igraph_vector_int_destroy(&path);
                new_disconnected.push_back(e_neg);
                break;
            }
            ms_dijkstra += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t_dij0).count();
            ++neg_edges_scanned;

            std::vector<Edge>                cycle_path_tmp;
            std::vector<igraph_integer_t>    path_pos_eids;
            std::vector<igraph_integer_t>    path_full_eids;
            const int nodes_on_path = igraph_vector_int_size(&path);
            cycle_path_tmp.reserve(nodes_on_path);
            path_pos_eids.reserve(nodes_on_path);
            path_full_eids.reserve(nodes_on_path);
            for (int i = 1; i < nodes_on_path; ++i) {
                const int a = VECTOR(path)[i - 1], b = VECTOR(path)[i];
                cycle_path_tmp.emplace_back(a, b);
                igraph_integer_t peid;
                if (igraph_get_eid(&g_pos_, &peid, a, b, 0, 0) == IGRAPH_SUCCESS) {
                    path_pos_eids.push_back(peid);
                    igraph_integer_t feid = pos2full_eid_[peid];
                    path_full_eids.push_back(feid);
                }
            }

            const auto t_chk0 = clock::now();
            bool valid = true;
            uint64_t sig = 1469598103934665603ull; // FNV-1a
            for (auto peid : path_pos_eids) { sig ^= (uint64_t)peid; sig *= 1099511628211ull; }
            auto ins = seen_paths.insert(sig);
            if (!ins.second) valid = false;
            ms_check += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t_chk0).count();

            if (valid) {
                const auto t_emit0 = clock::now();

                path_nodes_scanned += nodes_on_path;
                for (size_t i = 0; i < path_pos_eids.size(); ++i) ++pos_edges_on_paths;

                // If triangle-length path, enforce triangle per-vertex caps (legacy)
                if (nodes_on_path == 3) {
                    igraph_integer_t pe1 = path_pos_eids[0], pe2 = path_pos_eids[1];
                    igraph_integer_t e1  = pos2full_eid_[pe1], e2 = pos2full_eid_[pe2];
                    igraph_integer_t u1, v1, u2, v2;
                    igraph_edge(&G_.g, e1, &u1, &v1);
                    igraph_edge(&G_.g, e2, &u2, &v2);
                    igraph_integer_t w = (u1 == u2 || u1 == v2) ? u1 : v1;
                    if (w == uu || w == vv) w = (u1 == uu || u1 == vv) ? v1 : u1;

                    if (tri_used_per_vertex_[(int)uu] >= tri_cap_per_v_[(int)uu] ||
                        tri_used_per_vertex_[(int)vv] >= tri_cap_per_v_[(int)vv] ||
                        tri_used_per_vertex_[(int)w ] >= tri_cap_per_v_[(int)w ]) {
                        // exceed cap -> treat as valid path but don’t increment tri caps
                    } else {
                        ++tri_used_per_vertex_[(int)uu];
                        ++tri_used_per_vertex_[(int)vv];
                        ++tri_used_per_vertex_[(int)w];
                        ++triangles_emitted_now;
                    }
                }

                // Cross-batch bump along path
                const int L = nodes_on_path;
                for (auto feid : path_full_eids) bump_cross_batch_(feid, L);

                out.emplace_back(e_neg, std::move(cycle_path_tmp));
                ++total_found_; ++cycles_emitted_now; ++accepted_here;

                igraph_integer_t neg_eid;
                igraph_get_eid(&G_.g, &neg_eid, uu, vv, 0, 1);
                if (neg_eid >= 0 && neg_eid < ecount_) neg_edge_covered_[(size_t)neg_eid] = 1;

                // Soft mask bumps to encourage alternative path next rep
                for (auto peid : path_pos_eids) VECTOR(saved_weights_pos_)[peid] += alt_path_bump_;
                ms_emit += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t_emit0).count();
            } else {
                // Duplicate path -> bigger bump to break ties
                for (auto peid : path_pos_eids) VECTOR(saved_weights_pos_)[peid] += 2.0 * alt_path_bump_;
            }

            igraph_vector_int_destroy(&path);
            if (accepted_here == 0) break;
        }
    }

    ++batches_emitted_;

    const bool made_progress = (total_found_ > before_total);
    if (!made_progress) {
        finished_ = true;
    } else if (cover_) {
        // Keep only negatives that failed to produce *any* path under SP (disconnected)
        neg_edges_.swap(new_disconnected);
        if (neg_edges_.empty()) finished_ = true;
    } else {
        finished_ = true;
    }

    std::cout << "[TRI-PROFILE] negE_scanned=" << neg_edges_scanned
              << ", cycles_out=" << cycles_emitted_now
              << ", tri_out=" << triangles_emitted_now
              << ", sum_path_nodes=" << path_nodes_scanned
              << ", pos_edges_on_paths=" << pos_edges_on_paths
              << ", ms_dijkstra=" << ms_dijkstra
              << ", ms_check=" << ms_check
              << ", ms_emit=" << ms_emit
              << ", ms_misc=" << ms_misc
              << ", batches=" << batches_emitted_
              << ", batch_size=" << out.size()
              << std::endl;

    return !out.empty();
}
