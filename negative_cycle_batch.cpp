// negative_cycle_batch.cpp
#include "negative_cycle_batch.h"
#include "separation_pipeline_tls.h"
#include <unordered_set>
#include <algorithm>
#include <numeric>
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

// helper
static inline long long key64_pair(int a, int b) {
    if (a > b) std::swap(a,b);
    return (static_cast<long long>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint32_t>(b);
}

// -------------------------------------------------------------------------

NegativeCycleBatch::NegativeCycleBatch(const SignedGraphForMIP& G,
                                       const std::vector<Edge>& neg_edges_uncov,
                                       Params p)
    : P_(p), G_(G)
{
    tri_cap_per_vertex_ = P_.tri_cap_per_vertex;
    // seed the anchors directly from caller
    neg_edges_.assign(neg_edges_uncov.begin(), neg_edges_uncov.end());
    build_initial_state_();   // now only builds POS graph & weights; does NOT enumerate triangles
}

// REFACTORED: this no longer fills neg_edges_ from the graph.
// It builds positive graph, base weights, and then computes neg degrees from the provided anchors.
// Delegating convenience constructors (compat overloads)
NegativeCycleBatch::NegativeCycleBatch(const SignedGraphForMIP& G)
    : NegativeCycleBatch(G, G.get_negative_edges(), Params{}) {}

NegativeCycleBatch::NegativeCycleBatch(const SignedGraphForMIP& G,
                                       Params p)
    : P_(p), G_(G) {
    tri_cap_per_vertex_ = P_.tri_cap_per_vertex;
    neg_edges_ = G.get_negative_edges();
    build_initial_state_();
}

void NegativeCycleBatch::build_initial_state_() {
    vcount_ = G_.vertex_count();
    ecount_ = G_.edge_count();

    base_pos_.assign((size_t)ecount_, 1.0);
    reuse_accum_.assign((size_t)ecount_, 0.0);
    neg_deg_.assign((size_t)vcount_, 0);
    pos_deg_.assign((size_t)vcount_, 0);

    const auto& signs_view = G_.signs_view();
    const auto  pos_eids = G_.get_positive_eids();           // vector<int>

    // Build positive-only graph and mappings
    std::vector<igraph_integer_t> pos_edges; pos_edges.reserve(pos_eids.size() * 2);
    full2pos_eid_.assign((size_t)ecount_, -1);
    pos2full_eid_.clear(); pos2full_eid_.reserve(pos_eids.size());

    for (auto eid : pos_eids) {
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
        const igraph_integer_t fe = pos2full_eid_[pe];
        double val = base_pos_[(size_t)fe];
        if (P_.guide_len_full) {
            const double g = (*(P_.guide_len_full))[(size_t)fe];
            val = std::max(P_.guide_len_eps, g);
        }
        VECTOR(saved_weights_pos_)[pe] = val;
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
}

int NegativeCycleBatch::override_budget_(int base) const {
    return base; // passthrough (anneal later if needed)
}

void NegativeCycleBatch::on_emit_(int full_eid, double used_density) {
    // Steer Dijkstra away from overused positive edges within this batch
    if (full_eid < 0 || full_eid >= (int)full2pos_eid_.size()) return;
    igraph_integer_t pe = full2pos_eid_[(size_t)full_eid];
    if (pe < 0) return;

    used_density = std::max(0.0, used_density);
    // accumulate usage for commit-stage mirrors
    used_in_batch_pos_[(size_t)pe] += used_density;

    // ω′_local ← ω′_local + β_emit · used_density  (pipeline contract)
    double w = VECTOR(saved_weights_pos_)[pe] + P_.beta_emit * used_density;
    if (w < P_.omega_eps) w = P_.omega_eps;
    if (P_.omega_max > 0.0 && w > P_.omega_max) w = P_.omega_max;
    VECTOR(saved_weights_pos_)[pe] = w;
}

void NegativeCycleBatch::on_accept_(int full_eid, double /*density*/) {
    // Persistent cross-batch drift on edges that were used in accepted triangles
    bump_cross_batch_(full_eid, /*|C|=*/3);
}

// Inside NegativeCycleBatch (private section)
double NegativeCycleBatch::best_twohop_upper_bound_(int u, int v) const {
    double ub = std::numeric_limits<double>::infinity();

    // 1-hop (direct)
    igraph_integer_t e_uv;
    if (igraph_get_eid(&g_pos_, &e_uv, u, v, 0, 0) == IGRAPH_SUCCESS) {
        ub = std::min(ub, std::max(1e-12, (double) VECTOR(saved_weights_pos_)[e_uv]));
    }

    // 2-hop via triangles u-w-v
    igraph_vector_int_t Nu; igraph_vector_int_init(&Nu, 0);
    igraph_neighbors(&g_pos_, &Nu, u, IGRAPH_ALL);
    const igraph_integer_t nu = igraph_vector_int_size(&Nu);
    for (igraph_integer_t i = 0; i < nu; ++i) {
        const int w = VECTOR(Nu)[i];
        igraph_integer_t e_uw, e_wv;
        if (igraph_get_eid(&g_pos_, &e_uw, u, w, 0, 0) != IGRAPH_SUCCESS) continue;
        if (igraph_get_eid(&g_pos_, &e_wv, w, v, 0, 0) != IGRAPH_SUCCESS) continue;
        const double c =
            std::max(1e-12, (double) VECTOR(saved_weights_pos_)[e_uw]) +
            std::max(1e-12, (double) VECTOR(saved_weights_pos_)[e_wv]);
        if (c < ub) ub = c;
    }
    igraph_vector_int_destroy(&Nu);
    return ub;
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
        // Silence igraph "Couldn't reach some vertices" warnings just for this batch.
        struct IgraphWarningSilencer {
            igraph_warning_handler_t* prev = nullptr;
            static void noop(const char*, const char*, int) {}
            IgraphWarningSilencer() : prev(igraph_set_warning_handler(noop)) {}
            ~IgraphWarningSilencer() { igraph_set_warning_handler(prev); }
        } _silence_igraph_warnings_guard;
	    for (const auto& e_neg : neg_edges_uncov) {
	        const int uu = e_neg.first, vv = e_neg.second;
	        const int neg_full_eid = (int)edge_idx[Edge{uu, vv}];
	        double best_len = std::numeric_limits<double>::infinity();

            // Reserve bucket capacity up-front for this negative edge
            auto& buck_res = buckets[neg_full_eid];
            buck_res.reserve(std::max(1, P_.K_sp_per_neg));

			// Reuse a single igraph_vector_int_t across K repeats to avoid init/destroy churn
            igraph_vector_int_t path; igraph_vector_int_init(&path, 0);
	        // deterministic repeat up to K_sp_per_neg
	        for (int rep = 0; rep < std::max(1, P_.K_sp_per_neg); ++rep) {
				// Pre-Dijkstra prune for repeats using 1–2 hop upper bound on current ω′
				if (rep > 0 && std::isfinite(best_len)) {
				    const double ub2 = best_twohop_upper_bound_(uu, vv);
				    // If even the fastest 1–2 hop route cannot beat incumbent best_len, skip this repeat
				    if (!(ub2 + 1e-12 < best_len)) {
				        continue;
				    }
				}
				
				
				// 1) Collect peids for all triangles u-w-v (positive graph)
//				std::vector<igraph_integer_t> tri_peids;
//				tri_peids.reserve(64);
//				
//				const auto W = G_.common_pos_neighbors((size_t)uu, (size_t)vv); // Bitmap
//				for (size_t w_sz : W) {
//				    const int w = (int)w_sz;
//				    igraph_integer_t p1, p2;
//				    if (igraph_get_eid(&g_pos_, &p1, uu, w, /*directed*/0, /*error*/0) == IGRAPH_SUCCESS)
//				        tri_peids.push_back(p1);
//				    if (igraph_get_eid(&g_pos_, &p2, w, vv, /*directed*/0, /*error*/0) == IGRAPH_SUCCESS)
//				        tri_peids.push_back(p2);
//				}
				
				// 2) Snapshot old weights and bump to a big value
//				struct Saved { igraph_integer_t peid; double w; };
//				std::vector<Saved> saved;
//				saved.reserve(tri_peids.size());
//				constexpr double BIG = 1e9;  // large finite weight
//				
//				for (auto peid : tri_peids) {
//				    double &w = VECTOR(saved_weights_pos_)[peid];
//				    saved.push_back({peid, w});
//				    w *= 1.2;
//				}
				
				// 3) Run Dijkstra with triangle edges penalized
				igraph_vector_int_clear(&path);
				const auto t0 = clock::now();
				igraph_get_shortest_path_dijkstra(&g_pos_, &path, nullptr, uu, vv, &saved_weights_pos_, IGRAPH_ALL);
				ms_dijkstra += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t0).count();
	            ++neg_edges_scanned;
				
				// 4) Restore weights (keep your usual per-path bumps afterwards as you already do)
//				for (const auto& s : saved) VECTOR(saved_weights_pos_)[s.peid] = s.w;
	
	            const int nodes_on_path = igraph_vector_int_size(&path);
	            if (nodes_on_path < 2) {
					std::cout << "\n-\n\n";
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
	
	            path_nodes_scanned += nodes_on_path;
				pos_edges_on_paths += (int)pos_peids.size();
	            const double len_now = path_cost(pos_peids);
	
	            // only keep strictly improving bumped length for this anchor
	            if (len_now + 1e-12 < 1.0 * best_len) {
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
	
	                buck_res.push_back(std::move(c));
	
	                // within-batch density + emit bump ∝ 1/|C|
					const double dens = 1.0 / std::max(1, nodes_on_path);
					for (size_t j = 0; j < pos_peids.size(); ++j) {
					    const int full_eid = (int)full_eids[j];
					    on_emit_(full_eid, dens); // updates used_in_batch_pos_ and bumps ω′ with β_emit and clamps
					}
	            } else {
	                // discourage same path; larger bump
					for (auto peid : pos_peids) {
					    double w = VECTOR(saved_weights_pos_)[peid] + 2.0 * P_.beta_emit;
					    if (w < P_.omega_eps) w = P_.omega_eps;
					    if (P_.omega_max > 0.0 && w > P_.omega_max) w = P_.omega_max;
					    VECTOR(saved_weights_pos_)[peid] = w;
					}
	            }
	        }
            igraph_vector_int_destroy(&path);
	
	        // sort & truncate bucket
	        std::sort(buck_res.begin(), buck_res.end(), [](const SPCand& A, const SPCand& B){
	            if (A.score1 != B.score1) return A.score1 > B.score1; // higher better
	            if (A.score2 != B.score2) return A.score2 > B.score2; // higher better
	            return A.cost < B.cost; // shorter tie-break
	        });
	        if ((int)buck_res.size() > std::max(1, P_.K_sp_per_neg)) buck_res.resize(std::max(1, P_.K_sp_per_neg));
	    }
	
	    // --- Two-pass selection -------------------------------------------
	    const auto& sv = G_.signs_view();
	    std::vector<int> sp_used_per_vertex((size_t)vcount_, 0);
	    auto try_accept = [&](int neg_full_eid, const SPCand& c) -> bool {
	        // Vertex-cap test
	        for (int p : c.nodes) if (sp_used_per_vertex[(size_t)p] >= sp_cap) return false;
	
	        // Update SP per-vertex caps
	        for (int p : c.nodes) ++sp_used_per_vertex[(size_t)p];
	
	        // Triangle cap accounting for 2-pos-edge cycles
	        if (c.L == 3) {
	            // infer w and count using tri caps
	            const int uu = sv[neg_full_eid].points.first;
	            const int vv = sv[neg_full_eid].points.second;
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
	        const auto& se = sv[neg_full_eid];
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
	        struct Item { int neg_full_eid; const SPCand* c; };
	        auto cmp = [](const Item& A, const Item& B){
	            if (A.c->score1 != B.c->score1) return A.c->score1 < B.c->score1; // max-heap by score1
	            if (A.c->score2 != B.c->score2) return A.c->score2 < B.c->score2; // then by score2
	            if (A.c->cost   != B.c->cost  ) return A.c->cost   > B.c->cost;   // shorter first
	            return A.c < B.c; // deterministic tie-break
	        };
	
	        std::vector<Item> heap; heap.reserve(buckets.size()*2);
	        for (auto& kv : buckets) {
	            const int neg_full_eid = kv.first;
	            auto& buck = kv.second;
	            // Skip the first candidate (potentially used in Pass 1); add the rest
	            for (size_t i = 1; i < buck.size(); ++i) {
	                heap.push_back(Item{neg_full_eid, &buck[i]});
	            }
	        }
	
	        std::make_heap(heap.begin(), heap.end(), cmp);
	        while (accepted_sp < sp_budget && !heap.empty()) {
	            std::pop_heap(heap.begin(), heap.end(), cmp);
	            Item it = heap.back(); heap.pop_back();
	            if (try_accept(it.neg_full_eid, *it.c)) ++accepted_sp;
	        }
	    }
	
	    g_sp_cycles_accepted = accepted_sp;
	}
	
	++batches_emitted_;
    finished_ = true;

    std::cout << "[TRI_CYC-PROFILE] negE_scanned=" << neg_edges_scanned
    		  << ", K_sp_per_neg=" << P_.K_sp_per_neg
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
