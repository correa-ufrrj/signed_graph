// negative_cycle_batch.cpp
#include "negative_cycle_batch.h"
#include "separation_pipeline_tls.h"
#include <unordered_set>
#include <algorithm>
#include <chrono>
#include <limits>
#include <cstdint>

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

// NEW utility used by the compat constructors
std::vector<Edge>
NegativeCycleBatch::collect_all_negatives_(const SignedGraphForMIP& G) {
    std::vector<Edge> out;
    out.reserve((size_t)G.edge_count()/3);
    const auto& sw = G.get_switched_weight();
    const auto& sv = G.signs_view();
    const int m = G.edge_count();
    for (int eid = 0; eid < m; ++eid) {
        if (sw[eid] < 0.0) {
            const int u = (int)sv[eid].points.first;
            const int v = (int)sv[eid].points.second;
            out.emplace_back(u, v);
        }
    }
    return out;
}

// REFACTORED main constructor
NegativeCycleBatch::NegativeCycleBatch(const SignedGraphForMIP& G,
                                       const std::vector<Edge>& neg_edges_uncov,
                                       bool cover,
                                       Params p)
    : P_(p), G_(G), cover_(cover)
{
    tri_cap_per_vertex_ = P_.tri_cap_per_vertex;
    // seed the anchors directly from caller
    neg_edges_.assign(neg_edges_uncov.begin(), neg_edges_uncov.end());
    build_initial_state_();   // now only builds POS graph & weights; does NOT enumerate triangles
}

// REFACTORED: this no longer fills neg_edges_ from the graph.
// It builds positive graph, base weights, and then computes neg degrees from the provided anchors.
void NegativeCycleBatch::build_initial_state_() {
    vcount_ = G_.vertex_count();
    ecount_ = G_.edge_count();

    base_pos_.assign((size_t)ecount_, 1.0);
    reuse_accum_.assign((size_t)ecount_, 0.0);
    neg_deg_.assign((size_t)vcount_, 0);
    pos_deg_.assign((size_t)vcount_, 0);

    const auto& signs_view = G_.signs_view();
    const auto& sw         = G_.get_switched_weight();

    std::vector<double> pos_bases; pos_bases.reserve((size_t)ecount_);
    for (igraph_integer_t eid = 0; eid < ecount_; ++eid) {
        const auto se = signs_view[eid];
        const double w = sw[eid];
        if (w > 0.0) {
            base_pos_[(size_t)eid] = std::max(1e-12, std::fabs(w));
            pos_bases.push_back(base_pos_[(size_t)eid]);
            ++pos_deg_[(size_t)se.points.first];
            ++pos_deg_[(size_t)se.points.second];
        } else if (w == 0.0) {
            base_pos_[(size_t)eid] = 1e-12;
            pos_bases.push_back(1e-12);
            ++pos_deg_[(size_t)se.points.first];
            ++pos_deg_[(size_t)se.points.second];
        }
        // NOTE: w < 0.0 handled by caller via neg_edges_ (we do NOT push here).
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
        const auto se = signs_view[eid];
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

    // Now compute negative degrees from the provided anchors
    for (const auto& e : neg_edges_) {
        ++neg_deg_[(size_t)e.first];
        ++neg_deg_[(size_t)e.second];
    }
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

void NegativeCycleBatch::build_pos_adj_and_index_(
    TriangleBucketBatch::PosAdj& pos_adj,
    TriangleBucketBatch::EdgeIndex& edge_index) const
{
    pos_adj.assign((size_t)vcount_, {});
    const auto& signs_view = G_.signs_view();
    const auto& sw = G_.get_switched_weight();

    // Map all edges to full-eid; only positive edges go to adjacency
    for (igraph_integer_t eid = 0; eid < ecount_; ++eid) {
        const auto se = signs_view[eid];
        igraph_integer_t u = se.points.first, v = se.points.second;
        edge_index.emplace(key64_pair((int)u,(int)v), (int)eid);
        if (sw[eid] > 0.0) {
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

bool NegativeCycleBatch::next(std::vector<NegativeCycle>& out) {
    const auto edge_idx = G_.edge_index();
    out.clear();
    if (finished_) return false;

    using clock = std::chrono::steady_clock;
    long long ms_dijkstra = 0, ms_check = 0, ms_emit = 0, ms_misc = 0;
    size_t neg_edges_scanned = 0, cycles_emitted_now = 0;
    size_t triangles_emitted_now = 0;  // kept for logging; will stay 0
    long long path_nodes_scanned = 0, pos_edges_on_paths = 0;

    build_mask_for_batch_();

    const size_t before_total = total_found_;

    // === NO triangle-first step here anymore ===
    // The uncovered set is exactly what we were given at construction.
    std::vector<Edge> neg_edges_uncov = neg_edges_;
    std::vector<Edge> new_disconnected; new_disconnected.reserve(neg_edges_uncov.size());

    // === SP-based generation on uncovered anchors (bucketed two-pass) ===
	{
	    using clock = std::chrono::steady_clock;
	
	    const int sp_budget = (P_.B_sp > 0 ? P_.B_sp : std::numeric_limits<int>::max());
	    const int sp_cap    = std::max(1, P_.sp_cap_per_vertex);
	
	    // Candidate per uncovered negative edge (bucket)
	    struct SPCand {
	        std::vector<Edge>           pos_path;      // positive edges on the path
	        std::vector<int>            nodes;         // nodes on the path (to check caps)
	        std::vector<igraph_integer_t> pos_peids;   // pos-eids along the path
	        std::vector<igraph_integer_t> full_eids;   // full-eids along the path
	        double                      score1 = 0.0;  // primary (alpha * sum 1/omega')
	        double                      score2 = 0.0;  // secondary (viol/|C|) – unavailable here
	        double                      cost   = 0.0;  // current path length (for tiebreaks)
	        int                         L      = 0;    // number of nodes on path
	    };
	
	    // Map buckets by FULL eid of the negative anchor
	    std::unordered_map<int, std::vector<SPCand>> buckets;
	    buckets.reserve(neg_edges_uncov.size());

		// Helper: score_primary = inverse harmonic mean of ω′ along the path (α not needed for ordering)
		auto score_primary = [&](const std::vector<igraph_integer_t>& pos_peids) -> double {
		    if (pos_peids.empty()) return 0.0;
		    const double eps = 1e-12; // or thread through config if available
		    double sum_inv = 0.0;
		    int    k = 0;
		    for (auto peid : pos_peids) {
		        double w = std::max(eps, VECTOR(saved_weights_pos_)[peid]);
		        sum_inv += 1.0 / w;
		        ++k;
		    }
		    return (k > 0) ? (sum_inv / (double)k) : 0.0; // = 1 / Hmean(ω′)
		};
	
	    // Helper: path cost with current working weights
	    auto path_cost = [&](const std::vector<igraph_integer_t>& pos_peids) -> double {
	        double c = 0.0;
	        for (auto peid : pos_peids) c += std::max(1e-12, VECTOR(saved_weights_pos_)[peid]);
	        return c;
	    };
	
	    // --- Build buckets -------------------------------------------------
	    for (const auto& e_neg : neg_edges_uncov) {
	        const int uu = e_neg.first, vv = e_neg.second;
	        const int neg_full_eid = (int)edge_idx[Edge{uu, vv}];
	        double best_len = std::numeric_limits<double>::infinity();
	
	        // deterministic repeat up to K_sp_per_neg
	        for (int rep = 0; rep < std::max(1, P_.K_sp_per_neg); ++rep) {
	            igraph_vector_int_t path; igraph_vector_int_init(&path, 0);
	            const auto t_dij0 = clock::now();
	            igraph_get_shortest_path_dijkstra(&g_pos_, &path, nullptr, uu, vv, &saved_weights_pos_, IGRAPH_ALL);
	            ms_dijkstra += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t_dij0).count();
	            ++neg_edges_scanned;
	
	            const int nodes_on_path = igraph_vector_int_size(&path);
	            if (nodes_on_path < 2) {
	                igraph_vector_int_destroy(&path);
	                new_disconnected.push_back(e_neg);
	                break;
	            }
	
	            // materialize path edges & ids
	            std::vector<Edge> pos_path; pos_path.reserve(nodes_on_path);
	            std::vector<int>  nodes;    nodes.reserve(nodes_on_path);
	            std::vector<igraph_integer_t> pos_peids; pos_peids.reserve(nodes_on_path);
	            std::vector<igraph_integer_t> full_eids; full_eids.reserve(nodes_on_path);
	
	            for (int i = 0; i < nodes_on_path; ++i) nodes.push_back(VECTOR(path)[i]);
	            for (int i = 1; i < nodes_on_path; ++i) {
	                const int a = VECTOR(path)[i - 1], b = VECTOR(path)[i];
	                pos_path.emplace_back(a, b);
	                igraph_integer_t peid;
	                if (igraph_get_eid(&g_pos_, &peid, a, b, 0, 0) == IGRAPH_SUCCESS) {
	                    pos_peids.push_back(peid);
	                    full_eids.push_back(pos2full_eid_[peid]);
	                }
	            }
	
	            const double len_now = path_cost(pos_peids);
	
	            // only keep strictly improving bumped length for this anchor
	            if (len_now + 1e-12 < best_len) {
	                best_len = len_now;
	
	                SPCand c;
	                c.pos_path  = std::move(pos_path);
	                c.nodes     = std::move(nodes);
	                c.pos_peids = pos_peids;
	                c.full_eids = full_eids;
	                c.L         = nodes_on_path;
	                c.cost      = len_now;
	                c.score1    = score_primary(c.pos_peids); // φ, viol unavailable here
	                c.score2    = 0.0;
	
	                buckets[neg_full_eid].push_back(std::move(c));
	
	                // within-batch density + emit bump ∝ 1/|C|
	                const double dens = 1.0 / std::max(1, nodes_on_path);
	                for (auto peid : pos_peids) {
	                    used_in_batch_pos_[(size_t)peid] += dens;
	                    // approximate beta_emit: soft bump
	                    VECTOR(saved_weights_pos_)[peid] = std::max(1e-12, VECTOR(saved_weights_pos_)[peid] + alt_path_bump_ * dens);
	                }
	            } else {
	                // discourage same path; larger bump
	                for (auto peid : pos_peids) {
	                    VECTOR(saved_weights_pos_)[peid] = std::max(1e-12, VECTOR(saved_weights_pos_)[peid] + 2.0 * alt_path_bump_);
	                }
	            }
	            igraph_vector_int_destroy(&path);
	        }
	
	        // sort & truncate bucket
	        auto& buck = buckets[neg_full_eid];
	        std::sort(buck.begin(), buck.end(), [](const SPCand& A, const SPCand& B){
	            if (A.score1 != B.score1) return A.score1 > B.score1; // higher better
	            if (A.score2 != B.score2) return A.score2 > B.score2; // higher better
	            return A.cost < B.cost; // shorter tie-break
	        });
	        if ((int)buck.size() > std::max(1, P_.K_sp_per_neg)) buck.resize(std::max(1, P_.K_sp_per_neg));
	    }
	
	    // --- Two-pass selection -------------------------------------------
	    std::vector<int> sp_used_per_vertex((size_t)vcount_, 0);
	    auto try_accept = [&](int neg_full_eid, const SPCand& c) -> bool {
	        // Vertex-cap test
	        for (int p : c.nodes) if (sp_used_per_vertex[(size_t)p] >= sp_cap) return false;
	
	        // Update SP per-vertex caps
	        for (int p : c.nodes) ++sp_used_per_vertex[(size_t)p];
	
	        // Triangle cap accounting for 2-pos-edge cycles
	        if (c.L == 3) {
	            // infer w and count using tri caps
	            const int uu = G_.signs_view()[neg_full_eid].points.first;
	            const int vv = G_.signs_view()[neg_full_eid].points.second;
	            int w = -1;
	            if (!c.nodes.empty()) {
	                // path nodes are [u, w, v]
	                if ((int)c.nodes.size() == 3) w = c.nodes[1];
	            }
	            if (w >= 0) {
	                if (tri_used_per_vertex_[(int)uu] < tri_cap_per_vertex_) ++tri_used_per_vertex_[(int)uu];
	                if (tri_used_per_vertex_[(int)vv] < tri_cap_per_vertex_) ++tri_used_per_vertex_[(int)vv];
	                if (tri_used_per_vertex_[(int)w ] < tri_cap_per_vertex_) ++tri_used_per_vertex_[(int)w];
	                ++triangles_emitted_now;
	            }
	        }
	
	        // Cross-batch drift along full edges (|C| = L)
	        for (auto feid : c.full_eids) bump_cross_batch_(feid, c.L);
	
	        // Within-batch density mask (selected)
	        const double dens = 1.0 / std::max(1, c.L);
	        for (auto peid : c.pos_peids) {
	            used_in_batch_pos_[(size_t)peid] += dens;
	        }
	
	        // Materialize the cycle
	        const auto& se = G_.signs_view()[neg_full_eid];
			std::vector<Edge> path = c.pos_path; // copy then move (candidate is const)
			out.emplace_back(Edge{(int)se.points.first,(int)se.points.second}, std::move(path));
	        ++total_found_; ++cycles_emitted_now;
	
	        // Mark coverage for this negative edge
	        neg_edge_covered_[(size_t)neg_full_eid] = 1;
	        return true;
	    };
	
	    int accepted_sp = 0;
	
	    // Pass 1: take at most one per bucket (preserve coverage), respect budget
	    for (auto& kv : buckets) {
	        if (accepted_sp >= sp_budget) break;
	        const int neg_full_eid = kv.first;
	        auto& buck = kv.second;
	        for (const auto& c : buck) {
	            if (accepted_sp >= sp_budget) break;
	            if (try_accept(neg_full_eid, c)) { ++accepted_sp; break; }
	        }
	    }
	
	    // Pass 2: fill under budget across remaining candidates
	    if (accepted_sp < sp_budget) {
	        // Build a flat list of remaining candidates, ordered by score
	        std::vector<std::pair<int, const SPCand*>> rem; rem.reserve(buckets.size()*2);
	        for (auto& kv : buckets) {
	            const int neg_full_eid = kv.first;
	            auto& buck = kv.second;
	            // skip the first candidate (potentially used in pass 1); add the rest
	            for (size_t i = 0; i < buck.size(); ++i) rem.emplace_back(neg_full_eid, &buck[i]);
	        }
	        std::sort(rem.begin(), rem.end(), [](auto& A, auto& B){
	            if (A.second->score1 != B.second->score1) return A.second->score1 > B.second->score1;
	            if (A.second->score2 != B.second->score2) return A.second->score2 > B.second->score2;
	            return A.second->cost < B.second->cost;
	        });
	
	        for (auto& pr : rem) {
	            if (accepted_sp >= sp_budget) break;
	            if (try_accept(pr.first, *pr.second)) ++accepted_sp;
	        }
	    }
	
	    g_sp_cycles_accepted = accepted_sp;
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
