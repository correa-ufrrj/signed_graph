// negative_cycle_batch.cpp
#include "negative_cycle_batch.h"
#include <unordered_set>
#include <algorithm>
#include <chrono>

// Use centralized TBB bridge (strongly defined in separation_pipeline.cpp)
extern "C" {
void TBB_set_active(void* ctx,
                    void (*emit)(void*, int, double),
                    void (*accept)(void*, int, double),
                    int  (*budget)(void*, int));
void TBB_clear_active();
}

// ---------- C-callable wrappers that dispatch into this instance ----------
// NOTE: external linkage; the class declares them as friends.
void ncb_emit(void* ctx, int eid, double used_density) {
    static_cast<NegativeCycleBatch*>(ctx)->on_emit_(eid, used_density);
}
void ncb_accept(void* ctx, int eid, double density) {
    static_cast<NegativeCycleBatch*>(ctx)->on_accept_(eid, density);
}
int ncb_budget(void* ctx, int base) {
    return static_cast<NegativeCycleBatch*>(ctx)->override_budget_(base);
}

// Simple scope guard for TBB hookup
struct TBBHookScope {
    TBBHookScope(void* ctx,
                 void (*emit)(void*, int, double),
                 void (*accept)(void*, int, double),
                 int  (*budget)(void*, int)) {
        TBB_set_active(ctx, emit, accept, budget);
    }
    ~TBBHookScope() {
        TBB_clear_active();
    }
};

// helper
static inline long long key64_pair(int a, int b) {
    if (a > b) std::swap(a,b);
    return (static_cast<long long>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint32_t>(b);
}

// -------------------------------------------------------------------------

NegativeCycleBatch::NegativeCycleBatch(const SignedGraphForMIP& G,
                                       bool cover,
                                       bool use_triangle_order)
    : NegativeCycleBatch(G, cover, use_triangle_order, Params{}) {}

NegativeCycleBatch::NegativeCycleBatch(const SignedGraphForMIP& G,
                                       bool cover,
                                       bool /*use_triangle_order*/,
                                       Params p)
    : P_(p), G_(G), cover_(cover)
{
    tri_cap_per_vertex_ = P_.tri_cap_per_vertex;
    build_initial_state_();
}

int  NegativeCycleBatch::total_cycles_emitted() const { return static_cast<int>(total_found_); }
int  NegativeCycleBatch::batches_emitted()      const { return batches_emitted_; }

void NegativeCycleBatch::build_mask_for_batch_() {
    if (!saved_weights_pos_init_) return;
    for (long pe = 0; pe < (long)pos2full_eid_.size(); ++pe) {
        igraph_integer_t fe = pos2full_eid_[pe];
        VECTOR(saved_weights_pos_)[pe] = base_pos_[(size_t)fe];
    }
    if (saved_weights_init_) {
        for (long i = 0; i < (long)ecount_; ++i) VECTOR(saved_weights_)[i] = 0.0;
    }
    std::fill(reuse_accum_.begin(), reuse_accum_.end(), 0.0);
    used_in_batch_pos_.assign(pos2full_eid_.size(), 0.0);
    alt_path_bump_ = std::max(1e-12, P_.alt_path_bump_scale * med_base_pos_);
}

void NegativeCycleBatch::build_initial_state_() {
    vcount_ = G_.vertex_count();
    ecount_ = G_.edge_count();

    base_pos_.assign((size_t)ecount_, 1.0);
    reuse_accum_.assign((size_t)ecount_, 0.0);
    neg_edges_.clear(); neg_edges_.reserve((size_t)ecount_);
    neg_deg_.assign((size_t)vcount_, 0);
    pos_deg_.assign((size_t)vcount_, 0);

    std::vector<double> pos_bases; pos_bases.reserve((size_t)ecount_);
    for (igraph_integer_t eid = 0; eid < ecount_; ++eid) {
        const auto se = G_.signs_view()[eid];
        igraph_integer_t u = se.points.first, v = se.points.second;
        const double w = G_.get_switched_weight()[eid];
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
        const auto se = G_.signs_view()[eid];
        igraph_integer_t u = se.points.first, v = se.points.second;
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
}

void NegativeCycleBatch::build_pos_adj_and_index_(
    TriangleBucketBatch::PosAdj& pos_adj,
    TriangleBucketBatch::EdgeIndex& edge_index) const
{
    pos_adj.assign((size_t)vcount_, {});
    // Map all edges to full-eid; only positive edges go to adjacency
    for (igraph_integer_t eid = 0; eid < ecount_; ++eid) {
        const auto se = G_.signs_view()[eid];
        igraph_integer_t u = se.points.first, v = se.points.second;
        edge_index.emplace(key64_pair((int)u,(int)v), (int)eid);
        if (edge_is_pos(eid)) {
            pos_adj[(size_t)u].push_back((int)v);
            pos_adj[(size_t)v].push_back((int)u);
        }
    }
}

int NegativeCycleBatch::override_budget_(int base) const {
    return base; // passthrough (anneal later if needed)
}

void NegativeCycleBatch::on_emit_(int full_eid, double used_density) {
    // Steer Dijkstra away from overused positive edges within this batch
    if (full_eid < 0 || full_eid >= (int)full2pos_eid_.size()) return;
    igraph_integer_t pe = full2pos_eid_[(size_t)full_eid];
    if (pe < 0) return;
    // Tiny bump proportional to median and used_density
    const double bump = 0.05 * med_base_pos_ * std::max(0.0, used_density);
    VECTOR(saved_weights_pos_)[pe] = std::max(1e-12, VECTOR(saved_weights_pos_)[pe] + bump);
    // Track usage density for persistent (ω/H) updates at commit time
    used_in_batch_pos_[(size_t)pe] += std::max(0.0, used_density);
}

void NegativeCycleBatch::on_accept_(int full_eid, double /*density*/) {
    // Persistent cross-batch drift on edges that were used in accepted triangles
    bump_cross_batch_(full_eid, /*|C|=*/3);
}

bool NegativeCycleBatch::run_triangle_first_batch_(
    std::vector<NegativeCycle>& out,
    std::vector<int>& covered_neg_eids)
{
    const auto edge_idx = G_.edge_index();
    covered_neg_eids.clear();

    // Build TriangleBucketBatch inputs
    TriangleBucketBatch::PosAdj pos_adj;
    TriangleBucketBatch::EdgeIndex edge_index;
    build_pos_adj_and_index_(pos_adj, edge_index);

    // Neg edges under current switching
    std::vector<std::pair<int,int>> neg_edges_pairs; neg_edges_pairs.reserve(neg_edges_.size());
    for (const auto& e : neg_edges_) neg_edges_pairs.emplace_back(e.first, e.second);

    // Params (use P_.tbb; ensure sane defaults)
    TriangleBucketBatch::Params P = P_.tbb;
    if (P.K_tri_per_neg <= 0) P.K_tri_per_neg = 8;
    if (P.cap_per_vertex <= 0) P.cap_per_vertex = tri_cap_per_vertex_;
    if (P.B_tri <= 0) P.B_tri = std::max(1, (int)neg_edges_pairs.size() * P.K_tri_per_neg);

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
    {
        TBBHookScope scope(static_cast<void*>(this), &ncb_emit, &ncb_accept, &ncb_budget);
        const auto& selected = tbb.select(covered_neg_eids);

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
            igraph_integer_t neg_eid = (igraph_integer_t) edge_idx[Edge{c.u, c.v}];
            if (neg_eid >= 0 && neg_eid < ecount_) neg_edge_covered_[(size_t)neg_eid] = 1;
        }
    }

    return !out.empty();
}

bool NegativeCycleBatch::next(std::vector<NegativeCycle>& out) {
    const auto edge_idx = G_.edge_index();
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
    std::vector<int> covered_neg_eids;   // full-eids of anchors with a nonempty bucket
    if (run_triangle_first_batch_(triangles, covered_neg_eids)) {
        triangles_emitted_now += triangles.size();
        for (auto& C : triangles) out.push_back(std::move(C));
    }

    // Build uncovered list for SP (remove all nonempty bucket anchors)
    std::unordered_set<int> covered_set(covered_neg_eids.begin(), covered_neg_eids.end());
    std::vector<Edge> neg_edges_uncov; neg_edges_uncov.reserve(neg_edges_.size());
    for (const auto& e_neg : neg_edges_) {
        igraph_integer_t eid = (igraph_integer_t) edge_idx[Edge{e_neg.first, e_neg.second}];
        if (covered_set.find((int)eid) == covered_set.end()) neg_edges_uncov.push_back(e_neg);
    }

    std::vector<Edge> new_disconnected; new_disconnected.reserve(neg_edges_uncov.size());

    // === SP-based generation on uncovered anchors ===
    for (const auto& e_neg : neg_edges_uncov) {
        const int uu = e_neg.first, vv = e_neg.second;

        int accepted_here = 0;
        std::unordered_set<uint64_t> seen_paths;
        for (int rep = 0; rep < std::max(1, P_.K_sp_per_neg); ++rep) {
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
            std::vector<int> path_nodes; path_nodes.reserve(nodes_on_path);
            cycle_path_tmp.reserve(nodes_on_path);
            path_pos_eids.reserve(nodes_on_path);
            path_full_eids.reserve(nodes_on_path);
            for (int i = 0; i < nodes_on_path; ++i) path_nodes.push_back(VECTOR(path)[i]);
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
            // Count usage on this shortest path for persistent repulsion
            for (auto peid : path_pos_eids) used_in_batch_pos_[(size_t)peid] += 1.0;

            const auto t_chk0 = clock::now();
            bool valid = true;
            uint64_t sig = 1469598103934665603ull; // FNV-1a signature
            for (auto peid : path_pos_eids) { sig ^= (uint64_t)peid; sig *= 1099511628211ull; }
            if (!seen_paths.insert(sig).second) valid = false;
            ms_check += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t_chk0).count();

            if (valid) {
                const auto t_emit0 = clock::now();

                path_nodes_scanned += nodes_on_path;
                for (size_t i = 0; i < path_pos_eids.size(); ++i) ++pos_edges_on_paths;

                // If triangle-length path, enforce triangle per-vertex caps
                if (nodes_on_path == 3) {
                    igraph_integer_t pe1 = path_pos_eids[0], pe2 = path_pos_eids[1];
                    igraph_integer_t e1  = pos2full_eid_[pe1], e2 = pos2full_eid_[pe2];
                    const auto se1 = G_.signs_view()[e1];
                    const auto se2 = G_.signs_view()[e2];
                    igraph_integer_t u1 = se1.points.first, v1 = se1.points.second;
                    igraph_integer_t u2 = se2.points.first, v2 = se2.points.second;
                    igraph_integer_t w = (u1 == u2 || u1 == v2) ? u1 : v1;
                    if (w == uu || w == vv) w = (u1 == uu || u1 == vv) ? v1 : u1;

                    if (tri_used_per_vertex_[(int)uu] < tri_cap_per_vertex_ &&
                        tri_used_per_vertex_[(int)vv] < tri_cap_per_vertex_ &&
                        tri_used_per_vertex_[(int)w ] < tri_cap_per_vertex_) {
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

                igraph_integer_t neg_eid = (igraph_integer_t) edge_idx[Edge{uu, vv}];
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
        // Keep only negatives that failed to produce any path under SP (disconnected)
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

void NegativeCycleBatch::accumulate_pos_usage_to_full(std::vector<double>& dst, double scale) const {
    if ((size_t)dst.size() < (size_t)ecount_) dst.resize((size_t)ecount_, 0.0);
    for (size_t pe = 0; pe < pos2full_eid_.size(); ++pe) {
        igraph_integer_t fe = pos2full_eid_[pe];
        if (fe >= 0) dst[(size_t)fe] += scale * used_in_batch_pos_[pe];
    }
}
