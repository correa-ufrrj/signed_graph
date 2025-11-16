// File: signed_graph_mip.cpp
#include "signed_graph_mip.h"
#include "negative_cycle_batch.h"

#include <algorithm>
#include <numeric>
#include <limits>
#include <cmath>

void ShortestPathGraph::build_graph_() {
	const int vcount = G_.vertex_count();
	const int ecount = G_.edge_count();

	// Count positive edges and materialize undirected edge list once.
	// We deliberately use the owner's fast bitmap neighbors to avoid duplicates.
	// We insert (u,v) only when v < u to create each undirected edge once.
	int pos_count = G_.edge_count(+1);
	igraph_vector_int_t edges_vec; igraph_vector_int_init(&edges_vec, pos_count << 1);

    // Build positive-only graph and mappings
    full2pos_eid_.assign((size_t)ecount, -1);
    pos2full_eid_.clear(); pos2full_eid_.reserve(pos_count);
    
	int i = 0;
    int pe = 0;
     for (int u = 0; u < vcount; ++u) {
         const auto& bm = G_.neighbors_bm(static_cast<size_t>(u), +1);
         // enumeration in increasing order of neighbors
         for (size_t vs : bm) {
             const int v = static_cast<int>(vs);
             if (v < u) {
				VECTOR(edges_vec)[i++] = u;
				VECTOR(edges_vec)[i++] = v;
				const int fe = G_.edge_index(u, v);
				pos2full_eid_.push_back(fe);
				full2pos_eid_[(size_t)fe] = pe++;
             }
             else break;
         }
     }

    // (Re)create the positive-only graph
    if (built_) igraph_destroy(&g_pos_);
    igraph_create(&g_pos_, &edges_vec, vcount, IGRAPH_UNDIRECTED);
    igraph_vector_int_destroy(&edges_vec);
    built_ = true;

    if (w_pos_init_) igraph_vector_destroy(&w_pos_);
    igraph_vector_init(&w_pos_, pos_count);
    for (long pe = 0; pe < pos_count; ++pe) {
        set_weight(pe, 1.0);
    }
    w_pos_init_ = true;
}

ShortestPathGraph::Path ShortestPathGraph::dijkstra(int s, int t) const {
    Path out;

    // ignore predecessor vector; we only need node sequence
    igraph_get_shortest_path_dijkstra(&g_pos_, &sp_nodes_, /*prev=*/nullptr,
                                      s, t, &w_pos_, IGRAPH_ALL);

    const int L = igraph_vector_int_size(&sp_nodes_);
    if (L >= 2) {
        out.reachable_ = true;
        out.nodes_.reserve(L);
        out.pos_peids_.reserve(L-1);
        for (int i = 0; i < L; ++i) out.nodes_.push_back(VECTOR(sp_nodes_)[i]);

        // derive pos_eids and total cost
        out.cost_ = 0.0;
        for (int i = 1; i < L; ++i) {
            int a = out.nodes_[i-1], b = out.nodes_[i];
            int eid = full2pos_eid_[(size_t)G_.edge_index(a, b)];
            out.pos_peids_.push_back(eid);
            out.cost_ += std::max(w_eps, weight(eid));
        }
    }
    return out;
}

// ─────────────────────────────────────────────────────────────
// Constructors / destructor (aligned with new header)
// ─────────────────────────────────────────────────────────────

SignedGraphForMIP::SignedGraphForMIP(const std::string& file_path)
    : SignedGraph(file_path)
{
    recompute_polarity_degrees_();
}

SignedGraphForMIP::SignedGraphForMIP(SignedGraph&& base)
    : SignedGraph(std::move(base))
{
    recompute_polarity_degrees_();
}

SignedGraphForMIP::SignedGraphForMIP(const SignedGraphForMIP* const other)
    : SignedGraph(static_cast<const SignedGraph* const>(other))
{
    if (other && (int)other->dtilde_pos_.size() == vertex_count()) {
        dtilde_pos_ = other->dtilde_pos_;
        dtilde_neg_ = other->dtilde_neg_;
    } else {
        recompute_polarity_degrees_();
    }
}

SignedGraphForMIP::SignedGraphForMIP(const SignedGraphForMIP* const other,
                                     std::vector<double> new_sigma)
    : SignedGraph(static_cast<const SignedGraph* const>(other), std::move(new_sigma))
{
    recompute_polarity_degrees_();
}

SignedGraphForMIP::SignedGraphForMIP(const SignedGraphForMIP* const other,
                                     std::vector<double> new_sigma,
                                     const std::vector<double>& xhat)
    : SignedGraph(static_cast<const SignedGraph* const>(other),
                  std::move(new_sigma),
                  [&]{
                      std::vector<double> s_new(other->vertex_count());
                      for (int u = 0; u < other->vertex_count(); ++u) {
                          double xu = std::clamp(xhat[u], 0.0, 1.0);
                          s_new[u] = 2.0 * xu - 1.0;
                      }
                      return s_new;
                  }())
{
    recompute_polarity_degrees_();
}

SignedGraphForMIP::~SignedGraphForMIP() = default;

size_t SignedGraphForMIP::count_crossing_neg_edges_(const GraphCore::Bitmap& anchor) const {
    size_t total = 0;
    // C̄ = V \ A
    GraphCore::Bitmap comp = anchor.complement();

    // Sum_v∈C̄ | N⁻(v) ∩ A |
    for (size_t vs : comp) {
        const int v = static_cast<int>(vs);
        GraphCore::Bitmap tmp = neighbors_bm((size_t)v, -1); // copy (we will mutate tmp)
        tmp.iand(anchor);
        total += tmp.count();
    }
    return total;
}


void SignedGraphForMIP::on_neg_flip_batch_(const GraphCore::Bitmap& anchor) {
	GraphCore::on_neg_flip_batch_(anchor);
	
    // C̄ = V \ A
    GraphCore::Bitmap comp = anchor.complement();

    // For each vertex v in the complement, convert all N⁻(v) ∩ A to positive
    for (size_t vs : comp) {
        const int v = static_cast<int>(vs);

        // av = set of anchor endpoints u with (u,v) currently negative
        GraphCore::Bitmap av = neighbors_bm((size_t)v, -1); // copy
        av.iand(anchor);

        const size_t k = av.count();
        if (k == 0) continue;

        // Accumulate polarity-degree deltas for v (and for each u in av)
        // Before flip: p_uv < 0; After flip: p'_uv = -p_uv > 0
        // Δ for v:   d̃⁻ -= Σ(-p_uv), d̃⁺ += Σ(-p_uv)
        // Δ for u∈av: same as above, and also move v between u's bitmaps + degrees.
        for (size_t us : av) {
            const int u = static_cast<int>(us);

            const double p = polarity_of(u, v); // negative by construction
            const double mag = -p;              // positive magnitude

            // Update d̃ for both endpoints
            dtilde_neg_[v] -= mag; dtilde_pos_[v] += mag;
            dtilde_neg_[u] -= mag; dtilde_pos_[u] += mag;
        }
    }
}

void SignedGraphForMIP::on_vertex_flip_(int u, double from, double to)
{
    // 1) Structural (flip-only) maintenance from base
    SignedGraph::on_vertex_flip_(u, from, to);

    // 2) Fractional polarity deltas (always)
    for (int v : neighbors_bm(u)) {
        const double p_old = polarity_if(u, v, from);
        const double p_new = polarity_if(u, v, to);

        const double dpos = std::max(0.0,  p_new) - std::max(0.0,  p_old);
        const double dneg = std::max(0.0, -p_new) - std::max(0.0, -p_old);

        // Update both endpoints’ split sums
        dtilde_pos_[u] += dpos; dtilde_neg_[u] += dneg;
        dtilde_pos_[v] += dpos; dtilde_neg_[v] += dneg;
    }
}

// ─────────────────────────────────────────────────────────────
// Integer projection / edge polarities / reseeding
// ─────────────────────────────────────────────────────────────

void SignedGraphForMIP::apply_integer_projection() {
    // Round s to {−1,+1} and rebuild bitmaps/degrees in a single O(m) pass.
    for (int u = 0; u < vertex_count(); ++u) {
        s_[(size_t)u] = round_pm1_(s_[(size_t)u]);
    }
}

SignedGraphForMIP SignedGraphForMIP::integer_projection() const {
    // Round s_ to ±1 and keep sigma_ intact; effective signs become integral
    SignedGraphForMIP out = *this;
    out.apply_integer_projection();
    return out;
}

const std::vector<double> SignedGraphForMIP::edge_polarities() const {
    const int m = edge_count();
    std::vector<double> vals((size_t)m, 0.0);
    for (int eid = 0; eid < m; ++eid) {
        igraph_integer_t u,v; igraph_edge(g_.get(), eid, &u, &v);
        vals[(size_t)eid] = polarity_of((int)u, (int)v,(int)eid);
    }
    return vals;
}

// ─────────────────────────────────────────────────────────────
// Local switching alignment (aligned with guarded/batched API)
// ─────────────────────────────────────────────────────────────

void SignedGraphForMIP::align_switching(int u, double xfix01) {
    const double xu = get_rounded_x(u);
    std::cout << "[ALIGN] x[v*]=" << xu << " target=" << xfix01 << "\n";
    if (std::fabs(xu - xfix01) > 1e-9) {
        s_.flip_all_no_bitmaps();   // global flip keeps edge signs unchanged
    }
    std::cout << "[ALIGNED] x[v*]=" << get_x(u) << " target=" << xfix01 << "\n";
}

// ─────────────────────────────────────────────────────────────
// Re-seed from x̂ (apply via per-vertex alignment to keep logic simple)
// ─────────────────────────────────────────────────────────────

void SignedGraphForMIP::reseed_switching(const std::vector<double>& xhat) {
    if ((int)xhat.size() != vertex_count())
        throw std::runtime_error("xhat size mismatch");

    for (int u = 0; u < vertex_count(); ++u) {
        s_[(size_t)u] = 2.0 * xhat[(size_t)u] - 1.0;
    }
}
// ─────────────────────────────────────────────────────────────
// Greedy switching wrappers (mutable flavor)
// ─────────────────────────────────────────────────────────────

void SignedGraphForMIP::apply_greedy_switching()
{
    const igraph_t* G = g_.get();
    const int n = vertex_count();
    const double EPS = 1e-12;

    // Incumbent bookkeeping
    int best_mminus = edge_count(-1);
    int cur_mminus  = best_mminus + 1;
    std::vector<double> best_s = s_.readonly();

    // --- heap key (Step A) ---
    struct Key { double d; double sal; int u; int ver; };
    auto cmp = [](const Key& a, const Key& b){
        if (a.d   != b.d  ) return a.d   > b.d;   // min-heap by d
        if (a.sal != b.sal) return a.sal < b.sal; // prefer higher salience
        return a.u > b.u;
    };

    bool advanced = true;
    int r = 0;

    while (advanced) {
        advanced = false;

        // ─── Step A: fractional greedy vertex flips ──────────────────────────
        {
            std::priority_queue<Key, std::vector<Key>, decltype(cmp)> pq(cmp);
            std::vector<int> ver(n, 0);

            for (int u = 0; u < n; ++u) {
                const double du = net_polarity(u);
                if (du < -EPS) pq.push({du, vertex_salience(u), u, ver[u]});
            }

            while (!pq.empty() && pq.top().d < -EPS) {
//				std::cout << "[GREEDY-B] top net polarity=" << pq.top().d
//				          << "\n";
                const int u    = pq.top().u;
                const int veru = pq.top().ver;
                pq.pop();
                if (ver[u] != veru) continue; // stale key

                // Flip u: this triggers on_switch_sign_changed_ and updates
                // bitmaps, degree buckets, and (dtilde_pos_/dtilde_neg_) in O(deg(u)).
                single_switching(u);
                advanced = true;

                // Refresh keys for u and all its neighbors (both signs)
                ++ver[u];
                for (size_t v : neighbors_bm((size_t)u)) {
                    ++ver[(int)v];
                    if (net_polarity((int)v) < -EPS) {
                        pq.push({net_polarity((int)v), vertex_salience((int)v), (int)v, ver[(int)v]});
                    }
                }
            }

            if (advanced) {
                cur_mminus = edge_count(-1);
                if (cur_mminus <= best_mminus) {
                    best_mminus = cur_mminus;
                    best_s = s_.readonly();
                }
                continue; // Re-enter Step A while it keeps improving
            }
        }

        // ─── Step B: boundary-descent clique over N⁺ (if Step A stalled) ─────
        ++r;
        // (B.1) Find best seed pair (u,w) with w∈N⁺(u) minimizing B({u,w})
        double best_pair_cost = -EPS; // std::numeric_limits<double>::infinity();
        int u0 = -1, w0 = -1;

        for (int u = 0; u < n; ++u) {
            int scanned = 0;
            for (size_t vs : neighbors_bm((size_t)u, +1)) { // only strictly positive neighbors
                int v = (int)vs;
                const double cost = net_polarity(u) + net_polarity(v) - 2.0 * polarity_of(u, v);
                if (cost < best_pair_cost) { best_pair_cost = cost; u0 = u; w0 = v;
                    std::cout << "[GREEDY-B] pair=(" << u0 << "," << w0 << ")"
                	  << " best_pair_cost=" << best_pair_cost
		              << " |N⁺(u)∩N⁺(v)|=" << common_neighbors((size_t)u0, (size_t)w0, +1).count()
		              << " |N⁻(u)∩N⁻(v)|=" << common_neighbors((size_t)u0, (size_t)w0, -1).count()
		              << " |N(u)∩N(v)|=" << common_neighbors((size_t)u0, (size_t)w0).count()
		              << "\n";
		        }
            }
        }

		if (u0 >= 0) {
		    const double du  = net_polarity(u0);
		    const double dv  = net_polarity(w0);
		    const double puv = polarity_of(u0, w0);                 // s[u0]*σ(u0,w0)*s[w0] ∈ [-1,1]
		    const int degp_u = (int)neighbors_bm((size_t)u0, +1).count();
		    const int degp_v = (int)neighbors_bm((size_t)w0, +1).count();
		    auto W0 = common_neighbors((size_t)u0, (size_t)w0, +1);   // N⁺(u0) ∩ N⁺(w0)
		    const size_t frontier0 = W0.count();
		
		    // B({u0,w0}) should equal du + dv − 2·puv (for sanity)
		    std::cout << "[GREEDY-B] r=" << r
		              << " pair=(" << u0 << "," << w0 << ")"
		              << " B=" << (du + dv - 2.0 * puv)
		              << " (du=" << du << ", dv=" << dv << ", 2*p_uv=" << (2.0 * puv) << ")"
		              << " |N⁺(u)|=" << degp_u
		              << " |N⁺(v)|=" << degp_v
		              << " |N⁺(u)∩N⁺(v)|=" << frontier0
		              << " cur_minus=" << cur_mminus
		              << " best_minus=" << best_mminus
		              << "\n";
		} else {
		    std::cout << "[GREEDY-B] r=" << r
		              << " pair=<none>"
		              << " cur_minus=" << cur_mminus
		              << " best_minus=" << best_mminus
		              << "\n";
		}

        if (u0 >= 0 && best_pair_cost < -EPS) {
            std::vector<int> Q; Q.reserve(32); Q.push_back(u0); Q.push_back(w0);
            double B = best_pair_cost; // boundary cost

            // frontier W = N⁺(u0) ∩ N⁺(w0)
            auto W = common_neighbors((size_t)u0, (size_t)w0, +1);

            // acc[z] = Σ_{q∈Q} p_{zq} for current Q
            std::vector<double> acc(n, 0.0);
            for (size_t x : W) {
                int ix = (int)x;
                acc[ix] = polarity_of(u0, ix) + polarity_of(w0, ix);
            }

            while (!W.empty()) {
                // choose z* minimizing ΔB(z|Q) = d̃(z) − 2·acc[z]
                double best_delta = std::numeric_limits<double>::infinity();
                int zstar = -1;
                for (size_t x : W) {
                    int ix = (int)x;
                    double delta = net_polarity(ix) - 2.0 * acc[ix];
                    if (delta < best_delta) { best_delta = delta; zstar = ix; }
                }
                // Stop if we cannot strictly decrease boundary
                if (zstar < 0 || !(best_delta < -EPS) || !(B + best_delta < B - EPS)) break;

                Q.push_back(zstar);
                B += best_delta;

                // W ← W ∩ N⁺(z*) and update acc on survivors
                W.iand(neighbors_bm((size_t)zstar, +1));
                for (size_t x : W) acc[(int)x] += polarity_of(zstar, (int)x);
            }

            std::cout << "[GREEDY-B] |Q|=" << Q.size() << "\n";

            if ((int)Q.size() >= 2) {
                for (int v : Q) {
                    // Flip v: keeps bitmaps and (dtilde_pos_/neg_) in sync
                    single_switching(v);
                }
                advanced = true;

                cur_mminus = edge_count(-1);
                if (cur_mminus <= best_mminus) {
                    best_mminus = cur_mminus;
                    best_s = s_.readonly();
                }
            }
        }
        
		// ─── Step C: reconnect POS subgraph if disconnected (no per-vertex switching) ───
		// Preconditions: Step A stalled AND Step B did not advance this iteration.
		{
		    const auto comps = connected_components(+1);
		    if (comps.size() > 1) {
		        // Pick anchor = argmax over components of (# negative edges crossing it)
		        size_t best_i = 0;
		        size_t best_cross = 0;
		        for (size_t i = 0; i < comps.size(); ++i) {
		            const size_t cross = count_crossing_neg_edges_(comps[i]); // C.2 helper
		            if (cross > best_cross) { best_cross = cross; best_i = i; }
		        }
		
		        if (best_cross > 0) {
		            neg_flip_batch_(comps[best_i]); // C.3 helper
		            advanced = true;
	                cur_mminus = edge_count(-1);
                    best_mminus = cur_mminus;
                    best_s = s_.readonly();
		            continue; // Re-enter the while(advanced) loop
		        }
		    }
		} // ─── End of Step C ───
    } // while (advanced)

    // Commit the best incumbent switching back into s_
    // (assign_ is guarded and will flip bitmaps and update caches only where signs change)
    for (int u = 0; u < vertex_count(); ++u) {
        s_[u] = best_s[u];
    }
}

SignedGraphForMIP SignedGraphForMIP::greedy_switching() const {
    SignedGraphForMIP sg(this);      // copy
    sg.apply_greedy_switching(); // exposed via using-declaration
    return sg;
}

    
ShortestPathGraph SignedGraphForMIP::shortest_path_graph() const { return ShortestPathGraph(*this); }

// ─────────────────────────────────────────────────────────────
// Negative-cycle finder stream (unchanged API)
// ─────────────────────────────────────────────────────────────

NegativeCycleBatch SignedGraphForMIP::open_negative_cycle_stream() const {
    return NegativeCycleBatch(*this);
}

//std::vector<NegativeCycle>
//SignedGraphForMIP::find_switched_lower_bound(bool cover) {
//    std::vector<NegativeCycle> out;
//    auto stream = open_negative_cycle_stream(cover);
//    std::vector<NegativeCycle> batch;
//    while (stream.next(batch)) {
//        out.insert(out.end(),
//                   std::make_move_iterator(batch.begin()),
//                   std::make_move_iterator(batch.end()));
//    }
//    return out;
//}
//
//std::vector<std::vector<NegativeCycle>>
//SignedGraphForMIP::find_switched_lower_bound_grouped(bool cover) const {
//    std::vector<std::vector<NegativeCycle>> groups;
//    auto stream = open_negative_cycle_stream(cover);
//    std::vector<NegativeCycle> batch;
//    while (stream.next(batch)) {
//        groups.emplace_back();
//        groups.back().swap(batch);
//    }
//    return groups;
//}
