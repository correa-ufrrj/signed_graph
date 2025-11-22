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
    : P_(p), G_(G), SP(G_.shortest_path_graph())
{
    tri_cap_per_vertex_ = P_.tri_cap_per_vertex;
    // seed the anchors directly from caller
    neg_edges_.assign(neg_edges_uncov.begin(), neg_edges_uncov.end());
    build_initial_state_();   // now only builds POS graph & weights; does NOT enumerate triangles
}

// It builds positive graph, base weights, and then computes neg degrees from the provided anchors.
// Delegating convenience constructors (compat overloads)
NegativeCycleBatch::NegativeCycleBatch(const SignedGraphForMIP& G)
    : NegativeCycleBatch(G, G.get_edges(-1), Params{}) {}

NegativeCycleBatch::NegativeCycleBatch(const SignedGraphForMIP& G,
                                       Params p)
    : P_(p), G_(G), SP(G_.shortest_path_graph()) {
    tri_cap_per_vertex_ = P_.tri_cap_per_vertex;
    neg_edges_ = G.get_edges(-1);
    build_initial_state_();
}

void NegativeCycleBatch::build_initial_state_() {
    vcount_ = G_.vertex_count();
    ecount_ = G_.edge_count();

    base_pos_.assign((size_t)ecount_, 1.0);
//	SP.notify_signs_changed();

    tri_used_per_vertex_.assign((size_t)vcount_, 0);
    neg_edge_covered_.assign((size_t)ecount_, 0);
}

int  NegativeCycleBatch::total_cycles_emitted() const { return static_cast<int>(total_found_); }
int  NegativeCycleBatch::batches_emitted()      const { return batches_emitted_; }

void NegativeCycleBatch::build_mask_for_batch_() {
    // Re-seed SP working weights from persistent base (positive edges only)
    SP.reseed_weights(base_pos_);
    // Reset within-batch density mask (full-edge indexing)
    used_in_batch_pos_.assign((size_t)G_.edge_count(), 0.0);
}

int NegativeCycleBatch::override_budget_(int base) const {
    return base; // passthrough (anneal later if needed)
}

// TBB bridge: within-batch density by FULL eid (kept minimal)
void NegativeCycleBatch::on_emit_(int full_eid, double used_density) {
    if ((size_t)used_in_batch_pos_.size() < (size_t)ecount_) used_in_batch_pos_.assign((size_t)ecount_, 0.0);
    if (full_eid >= 0 && full_eid < (int)used_in_batch_pos_.size()) {
        used_in_batch_pos_[(size_t)full_eid] += std::max(0.0, used_density);
	    SP.bump_weight(full_eid, P_.beta_emit * used_density);
    }
}

void NegativeCycleBatch::on_emit_(const ShortestPathGraph::Path& p) {
	// within-batch density + emit bump ∝ 1/|C|
    const double used_density = std::max(0.0, 1.0 / std::max(1, p.length()));
    // accumulate usage (full-edge indexing)
	SP.for_each_eid(p, [&](int fe){ used_in_batch_pos_[(size_t)fe] += used_density; });
    // ω′_local bump on SP by endpoints
    SP.bump_weights(p, P_.beta_emit * used_density);
}

void NegativeCycleBatch::on_accept_(int full_eid, double /*density*/) {
    // Persistent cross-batch drift on edges that were used in accepted triangles
    bump_cross_batch_(full_eid, /*|C|=*/3);
}

// negative_cycle_batch.cpp  (inside NegativeCycleBatch)
double NegativeCycleBatch::best_twohop_upper_bound_(int u, int v) const {
    double ub = std::numeric_limits<double>::infinity();
    // 1-hop (direct) if POS(u,v) exists; guard to avoid throwing on absent edge
    try {
        ub = std::min(ub, std::max(P_.guide_len_eps, SP.weight(u, v)));
    } catch (...) {
        /* no direct POS edge */
    }
    // 2-hop via triangles u-w-v over POSITIVE neighbors of u
    const auto Nu = G_.common_neighbors((size_t)u, (size_t)v, +1); // GraphCore::Bitmap from the owner
    for (size_t w_sz : Nu) {
        const int w = static_cast<int>(w_sz);
        const double c = std::max(1e-12, SP.weight(u, w)) + std::max(1e-12, SP.weight(w, v));
        if (c < ub) ub = c;
    }
    return ub;
}

bool NegativeCycleBatch::next(std::vector<NegativeCycle>& out) {
    const auto edge_idx = G_.edge_index();
    out.clear();
    if (finished_) return false;

    using clock = std::chrono::steady_clock;
    DBG_SP_DECL(
      long long ms_dijkstra = 0, ms_check = 0, ms_emit = 0, ms_misc = 0;
      size_t neg_edges_scanned = 0, cycles_emitted_now = 0;
      size_t triangles_emitted_now = 0;
      long long path_nodes_scanned = 0, pos_edges_on_paths = 0;
      const size_t before_total = total_found_;
    );

    build_mask_for_batch_();

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
	        ShortestPathGraph::Path		path;
	        double                      score1 = 0.0;  // primary (alpha * sum 1/omega')
	        double                      score2 = 0.0;  // secondary (viol/|C|) – unavailable here
	    };
	
	    // Map buckets by FULL eid of the negative anchor
	    std::unordered_map<int, std::vector<SPCand>> buckets;
	    buckets.reserve(neg_edges_uncov.size());

        // Helper: score_primary = inverse harmonic mean of ω′ along the path
        auto score_primary = [&](const ShortestPathGraph::Path& p) -> double {
            if (p.length() <= 1) return 0.0;
            double sum_inv = 0.0; int k = 0;
            SP.for_each_weight(p, [&](double w){
                sum_inv += 1.0 / std::max(P_.guide_len_eps, w);
                ++k;
            });
            return (k > 0) ? (sum_inv / (double)k) : 0.0; // = 1 / Hmean(ω′)
        };
        
	    // --- Build buckets in order of neg edges -------------------------------------------------
	    for (const auto& e_neg : neg_edges_uncov) {
	        const int uu = e_neg.first, vv = e_neg.second;
	        const int neg_full_eid = (int)edge_idx[Edge{uu, vv}];
	        double best_len = std::numeric_limits<double>::infinity();

            // Reserve bucket capacity up-front for this negative edge
            auto& buck_res = buckets[neg_full_eid];
            buck_res.reserve(std::max(1, P_.K_sp_per_neg));

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
				
                // 3) Run Dijkstra (through SP)
                DBG_SP_DECL(const auto t0 = clock::now(););
                auto p = SP.dijkstra(uu, vv);
                DBG_SP({
                  ms_dijkstra += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t0).count();
                  ++neg_edges_scanned;
                });
                if (!p.reachable() || p.length() < 2) {
					DBG_SP(std::cout << "-\n\n";);
	                break;
	            }
	            const double len_now = p.cost();
	
	            // only keep strictly improving bumped length for this anchor
	            if (len_now + 1e-12 < 1.0 * best_len) {
	                best_len = len_now;
	
	                SPCand c;
	                c.score1    = score_primary(p); // φ, viol unavailable here
	                c.score2    = 0.0;
	                c.path		= p;
	
	                buck_res.push_back(std::move(c));
	
	                // within-batch density + emit bump ∝ 1/|C|
					on_emit_(p);
	            } else {
	                // discourage same path; larger bump
	                SP.bump_weights(p, 2.0 * P_.beta_emit);
	            }
	        }
	
	        // sort & truncate bucket
	        std::sort(buck_res.begin(), buck_res.end(), [](const SPCand& A, const SPCand& B){
	            if (A.score1 != B.score1) return A.score1 > B.score1; // higher better
	            if (A.score2 != B.score2) return A.score2 > B.score2; // higher better
	            return A.path.cost() < B.path.cost(); // shorter tie-break
	        });
	        if ((int)buck_res.size() > std::max(1, P_.K_sp_per_neg)) buck_res.resize(std::max(1, P_.K_sp_per_neg));
	    }
	
	    // --- Two-pass selection -------------------------------------------
	    const auto& sv = G_.signs_view();
	    std::vector<int> sp_used_per_vertex((size_t)vcount_, 0);
	    auto try_accept = [&](int neg_full_eid, const SPCand& c) -> bool {
	        // Vertex-cap test
	        for (int p : c.path.nodes()) if (sp_used_per_vertex[(size_t)p] >= sp_cap) return false;
	
	        // Update SP per-vertex caps
	        for (int p : c.path.nodes()) ++sp_used_per_vertex[(size_t)p];
	
	        // Triangle cap accounting for 2-pos-edge cycles
	        if (c.path.length() == 3) {
	            // infer w and count using tri caps
	            const int uu = sv[neg_full_eid].points.first;
	            const int vv = sv[neg_full_eid].points.second;
	            int w = -1;
	            w = c.path.nodes()[1];
	            if (w >= 0) {
	                if (tri_used_per_vertex_[(int)uu] < tri_cap_per_vertex_) ++tri_used_per_vertex_[(int)uu];
	                if (tri_used_per_vertex_[(int)vv] < tri_cap_per_vertex_) ++tri_used_per_vertex_[(int)vv];
	                if (tri_used_per_vertex_[(int)w ] < tri_cap_per_vertex_) ++tri_used_per_vertex_[(int)w];
	                DBG_SP(++triangles_emitted_now;);
	            }
	        }

	        // Within-batch density mask (selected)
	        const double dens = 1.0 / std::max(1, c.path.length());
	        SP.for_each_eid(c.path, [&](int feid){
	            bump_cross_batch_(feid, c.path.length());
	            used_in_batch_pos_[(size_t)feid] += dens;
	        });
	
	        // Materialize the cycle
	        const auto& se = sv[neg_full_eid];
			out.emplace_back(Edge{(int)se.points.first,(int)se.points.second}, c.path.edges());
	        ++total_found_;
	        DBG_SP(++cycles_emitted_now;);
	
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
	            if (A.c->path.cost() != B.c->path.cost()) return A.c->path.cost() > B.c->path.cost(); // shorter first
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

    DBG_SP(std::cout << "[TRI_CYC-PROFILE] negE_scanned=" << neg_edges_scanned
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
              << std::endl;);

    return !out.empty();
}

void NegativeCycleBatch::accumulate_pos_usage_to_full(std::vector<double>& dst, double scale) const {
    if ((size_t)dst.size() < (size_t)ecount_) dst.resize((size_t)ecount_, 0.0);
    for (size_t fe = 0; fe < (size_t)ecount_; ++fe) {
        const double d = used_in_batch_pos_[fe];
        if (d > 0.0) dst[fe] += scale * d;
    }
}
