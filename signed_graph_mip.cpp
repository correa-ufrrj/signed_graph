// File: signed_graph_mip.cpp
#include "signed_graph_mip.h"
#include "negative_cycle_batch.h"

#include <algorithm>
#include <numeric>
#include <limits>
#include <cmath>

// ─────────────────────────────────────────────────────────────
// Constructors / destructor (aligned with new header)
// ─────────────────────────────────────────────────────────────

SignedGraphForMIP::SignedGraphForMIP(const std::string& file_path)
    : SignedGraph(file_path) {}

SignedGraphForMIP::SignedGraphForMIP(SignedGraph&& base)
    : SignedGraph(std::move(base)) {}

SignedGraphForMIP::SignedGraphForMIP(const SignedGraphForMIP* const other)
    : SignedGraph(static_cast<const SignedGraph* const>(other)) {}

SignedGraphForMIP::SignedGraphForMIP(const SignedGraphForMIP* const other,
                                     std::vector<double> new_sigma)
    : SignedGraph(static_cast<const SignedGraph* const>(other), std::move(new_sigma)) {}

SignedGraphForMIP::SignedGraphForMIP(const SignedGraphForMIP* const other,
                                     std::vector<double> new_sigma,
                                     const std::vector<double>& xhat)
    : SignedGraph(static_cast<const SignedGraph* const>(other),
                  std::move(new_sigma),
                  [&](){
                      std::vector<double> s_new; s_new.resize(other->vertex_count());
                      for (int u = 0; u < other->vertex_count(); ++u) {
                          double xu = xhat[u]; if (xu < 0.0) xu = 0.0; else if (xu > 1.0) xu = 1.0;
                          s_new[u] = 2.0 * xu - 1.0;
                      }
                      return s_new;
                  }()) {}

SignedGraphForMIP::~SignedGraphForMIP() = default;

// Called BEFORE s[u] is assigned 'to'. Base updates (+/−) rows and integer counts,
// then we adjust polarity-degree splits d̃⁺/d̃⁻ by local deltas (O(deg(u))).
void SignedGraphForMIP::on_switch_sign_changed_(int u, double from, double to,
                                                igraph_vector_int_t* incident)
{
    // First, perform the structural flip (degrees + bitmaps)
    SignedGraph::on_switch_sign_changed_(u, from, to, incident);

    const igraph_t* G = g_.get();

    igraph_vector_int_t local;
    bool owned = false;
    if (!incident) { igraph_vector_int_init(&local, 0); incident = &local; owned = true; }
    if (igraph_vector_int_size(incident) == 0) {
        igraph_incident(G, incident, u, IGRAPH_ALL);
    }

    const int deg = (int)igraph_vector_int_size(incident);

    // Only s[u] changed: p'_{uv} − p_{uv} = (to − from) * σ_uv * s[v]
    for (int k = 0; k < deg; ++k) {
        const int eid = VECTOR(*incident)[k];

        igraph_integer_t a, b; igraph_edge(G, eid, &a, &b);
        const int v = ((int)a == u) ? (int)b : (int)a;

		const double p_old = polarity_of(u, v, eid); // BEFORE change
		const double p_new = polarity_if(u, v, eid, to); // AFTER change

		// Vertex u: you can either
		//   (i) swap dtilde_pos[u] / dtilde_neg[u] up-front like you did for counts,
		//       then add the above deltas symmetrically for u, or
		//   (ii) just apply the same delta logic to u as well.
		
		const double pos_delta = std::max(0.0, p_new) - std::max(0.0, p_old);
		const double neg_delta = std::max(0.0,-p_new) - std::max(0.0,-p_old);

        // Update u’s split sums
        dtilde_pos_[u] += pos_delta;
        dtilde_neg_[u] += neg_delta;

        // And v’s split sums (each edge contributes to both endpoints)
		dtilde_pos_[v] += pos_delta;
		dtilde_neg_[v] += neg_delta;
    }

    if (owned) igraph_vector_int_destroy(&local);
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

void SignedGraphForMIP::apply_greedy_switching(const GreedyKickOptions& opts)
{
    const igraph_t* G = g_.get();
    const int n = vertex_count();
    const double EPS = 1e-12;

    // Live accessor for fractional net-degree d̃(u) = sum_v p_uv
    auto dtilde_of = [&](int u) -> double {
        return dtilde_pos_[u] - dtilde_neg_[u];
    };

    // Incumbent bookkeeping
    int best_mminus = negative_edge_count();
    int cur_mminus  = best_mminus + 1;
    std::vector<double> best_s = s_.readonly();

    // --- heap key (Step A) ---
    struct Key { double d; double sal; int u; int ver; };
    auto cmp = [](const Key& a, const Key& b){
        if (a.d   != b.d  ) return a.d   > b.d;   // min-heap by d
        if (a.sal != b.sal) return a.sal < b.sal; // prefer higher salience
        return a.u > b.u;
    };

    const int R = std::max(1, opts.R_max);
    bool advanced = true;
    int r = 0;

    while (advanced) {
        advanced = false;

        // ─── Step A: fractional greedy vertex flips ──────────────────────────
        {
            std::priority_queue<Key, std::vector<Key>, decltype(cmp)> pq(cmp);
            std::vector<int> ver(n, 0);

            for (int u = 0; u < n; ++u) {
                const double du = dtilde_of(u);
                if (du < -EPS) pq.push({du, vertex_salience(u), u, ver[u]});
            }

            while (!pq.empty() && pq.top().d < -EPS) {
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
                for (size_t v : neighbors_bm((size_t)u, +1)) {
                    ++ver[(int)v];
                    if (dtilde_of((int)v) < -EPS) {
                        pq.push({dtilde_of((int)v), vertex_salience((int)v), (int)v, ver[(int)v]});
                    }
                }
                for (size_t v : neighbors_bm((size_t)u, -1)) {
                    ++ver[(int)v];
                    if (dtilde_of((int)v) < -EPS) {
                        pq.push({dtilde_of((int)v), vertex_salience((int)v), (int)v, ver[(int)v]});
                    }
                }
            }

            if (advanced) {
                cur_mminus = negative_edge_count();
                if (cur_mminus <= best_mminus) {
                    best_mminus = cur_mminus;
                    best_s = s_.readonly();
                }
                continue; // Re-enter Step A while it keeps improving
            }
        }

        // ─── Step B: boundary-descent clique over N⁺ (if Step A stalled) ─────
        ++r;
        if (r <= R) {
            const auto EIDX = edge_index(); // readonly eid map

            auto p_of = [&](int a, int b) -> double {
                const int eid = EIDX[Edge{a,b}];
                return polarity_of(a, b, eid); // s[a]*σ_ab*s[b] in [-1,1]
            };

            // (B.1) Find best seed pair (u,w) with w∈N⁺(u) minimizing B({u,w})
            double best_pair_cost = std::numeric_limits<double>::infinity();
            int u0 = -1, w0 = -1;

            for (int u = 0; u < n; ++u) {
                // p_max(u) = max_v p_uv over all neighbors (any sign)
                double pmax = -std::numeric_limits<double>::infinity();
                for (size_t vv : neighbors_bm((size_t)u)) {
                    int v = (int)vv;
                    double pv = p_of(u, v);
                    if (pv > pmax) pmax = pv;
                }
                // gate: d̃(u) < pmax ⇒ there exists a positive pair that can reduce boundary
                if (!(dtilde_of(u) < pmax - EPS)) continue;

                int scanned = 0;
                for (size_t vs : neighbors_bm((size_t)u, +1)) { // only strictly positive neighbors
                    if (opts.neighbor_cap > 0 && scanned++ >= opts.neighbor_cap) break;
                    int v = (int)vs;
                    const double cost = dtilde_of(u) + dtilde_of(v) - 2.0 * p_of(u, v);
                    if (cost < best_pair_cost) { best_pair_cost = cost; u0 = u; w0 = v; }
                }
            }

            std::cout << "[GREEDY-B] r=" << r
                      << " |Z0|=0"
                      << " cur_minus=" << cur_mminus
                      << " best_minus=" << best_mminus
                      << " neg edges=" << negative_edge_count()
                      << " project=" << integer_projection().negative_edge_count()
                      << "\n";

            if (u0 >= 0 && best_pair_cost < -EPS) {
                std::vector<int> Q; Q.reserve(32); Q.push_back(u0); Q.push_back(w0);
                double B = best_pair_cost; // boundary cost

                // frontier W = N⁺(u0) ∩ N⁺(w0)
                auto W = common_pos_neighbors((size_t)u0, (size_t)w0);

                // acc[z] = Σ_{q∈Q} p_{zq} for current Q
                std::vector<double> acc(n, 0.0);
                for (size_t x : W) {
                    int ix = (int)x;
                    acc[ix] = p_of(u0, ix) + p_of(w0, ix);
                }

                while (!W.empty()) {
                    // choose z* minimizing ΔB(z|Q) = d̃(z) − 2·acc[z]
                    double best_delta = std::numeric_limits<double>::infinity();
                    int zstar = -1;
                    for (size_t x : W) {
                        int ix = (int)x;
                        double delta = dtilde_of(ix) - 2.0 * acc[ix];
                        if (delta < best_delta) { best_delta = delta; zstar = ix; }
                    }
                    // Stop if we cannot strictly decrease boundary
                    if (zstar < 0 || !(best_delta < -EPS) || !(B + best_delta < B - EPS)) break;

                    Q.push_back(zstar);
                    B += best_delta;

                    // W ← W ∩ N⁺(z*) and update acc on survivors
                    W.iand(neighbors_bm((size_t)zstar, +1));
                    for (size_t x : W) acc[(int)x] += p_of(zstar, (int)x);
                }

                std::cout << "[GREEDY-B] |Q|=" << Q.size() << "\n";

                if ((int)Q.size() >= 2) {
                    for (int v : Q) {
                        // Flip v: keeps bitmaps and (dtilde_pos_/neg_) in sync
                        single_switching(v);
                    }
                    advanced = true;

                    cur_mminus = negative_edge_count();
                    if (cur_mminus <= best_mminus) {
                        best_mminus = cur_mminus;
                        best_s = s_.readonly();
                    }
                }
            }
        }

        // No Step C (removed)
    }

    // Commit the best incumbent switching back into s_
    // (assign_ is guarded and will flip bitmaps and update caches only where signs change)
    for (int u = 0; u < vertex_count(); ++u) {
        s_[u] = best_s[u];
    }
}

SignedGraphForMIP SignedGraphForMIP::greedy_switching(const GreedyKickOptions& opts) const {
    SignedGraphForMIP sg(this);      // copy
    sg.apply_greedy_switching(opts); // exposed via using-declaration
    return sg;
}

SignedGraphForMIP SignedGraphForMIP::greedy_switching() const {
    GreedyKickOptions opts; return greedy_switching(opts);
}

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
