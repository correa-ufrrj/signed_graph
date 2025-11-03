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

// ─────────────────────────────────────────────────────────────
// Integer projection (reuse base SignedGraph::integer_projection)
// ─────────────────────────────────────────────────────────────

void SignedGraphForMIP::apply_integer_projection() {
    round_pm1_inplace(s_);
    compute_degrees();
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
        vals[(size_t)eid] = s_[(int)u] * static_cast<double>(sigma_[(int)eid]) * s_[(int)v];
    }
    return vals;
}

// ─────────────────────────────────────────────────────────────
// Local switching alignment (keep existing logic)
// ─────────────────────────────────────────────────────────────

void SignedGraphForMIP::align_switching(int u, double xfix01) {
    const double xu = get_rounded_x(u);
    std::cout << "[ALIGN] x[v*]=" << xu << " target=" << xfix01 << "\n";
    if (std::fabs(xu - xfix01) > 1e-9) {
        for (double& sv : s_) sv = -sv;   // global flip keeps edge signs unchanged
    }
    std::cout << "[ALIGNED] x[v*]=" << get_x(u) << " target=" << xfix01 << "\n";
}

// ─────────────────────────────────────────────────────────────
// Greedy switching wrappers (mutable flavor)
// ─────────────────────────────────────────────────────────────

void SignedGraphForMIP::apply_greedy_switching(const GreedyKickOptions& opts)
{
    const igraph_t* G = g_.get();
    const int n = vertex_count();
    const int m = edge_count();
    const double EPS = 1e-12;

    // initialize fractional net degrees d~(u) = sum_{v∈N(u)} s[u]*sigma_uv*s[v]
    std::vector<double> dtilde(n, 0.0);
    for (int eid = 0; eid < m; ++eid) {
        igraph_integer_t u, v; igraph_edge(G, eid, &u, &v);
        const double p = s_[(int)u] * static_cast<double>(sigma_[eid]) * s_[(int)v];
        dtilde[(int)u] += p; dtilde[(int)v] += p;
    }

    // incumbent
    int best_mminus = negative_edge_count();
    int cur_mminus  = best_mminus + 1;
    std::vector<double> best_s = s_;

    // --- chooser for step B (weighted by salience). If frac_y==nullptr, fall back to vertex_salience ---
    auto pick_u_star_in_Q = [&](const std::vector<int>& Q) -> int {
        if (Q.empty()) return -1;
        igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
        double best_score = -std::numeric_limits<double>::infinity();
        int    best_u     = -1;
        for (int u : Q) {
            double score = 0.0;
            igraph_incident(G, &inc, u, IGRAPH_ALL);
            for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
                const int eid = VECTOR(inc)[ii];
                const int v   = IGRAPH_OTHER(G, eid, u);
                if (!is_pos_edge(u, v)) continue; // strictly positive only
                score += edge_salience(u, v);
            }
            igraph_vector_int_clear(&inc);
            if (score > best_score || (score == best_score && u < best_u)) {
                best_score = score; best_u = u;
            }
        }
        igraph_vector_int_destroy(&inc);
        return best_u;
    };

    struct Key { double d; double sal; int u; int ver; };
    auto cmp = [](const Key& a, const Key& b){
        if (a.d   != b.d  ) return a.d   > b.d;   // min-heap
        if (a.sal != b.sal) return a.sal < b.sal; // prefer higher salience
        return a.u > b.u;
    };

    const int R = std::max(1, opts.R_max);
    bool advanced = true;
    int r = 0;
    while (advanced) {
        advanced = false;

        // A) Fractional greedy flips (min-heap on d~, tie-break by salience)
        {
            std::priority_queue<Key, std::vector<Key>, decltype(cmp)> pq(cmp);
            std::vector<int> ver(n, 0);
            for (int u = 0; u < n; ++u) if (is_neg_(dtilde[u])) pq.push({dtilde[u], vertex_salience(u), u, ver[u]});

            igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
            while (!pq.empty() && pq.top().d < -EPS) {
                const int u = pq.top().u;
                const int veru = pq.top().ver;
                pq.pop();
                if (ver[u] != veru) continue; // stale key

                flip_and_refresh_net_degree(u, dtilde);
                advanced = true;

                // refresh keys for u and its neighbors
                igraph_incident(G, &inc, u, IGRAPH_ALL);
                ++ver[u];
                for (int k = 0; k < (int)igraph_vector_int_size(&inc); ++k) {
                    const int v = IGRAPH_OTHER(G, VECTOR(inc)[k], u);
                    ++ver[v];
                    if (dtilde[v] < -EPS) {
                        pq.push({dtilde[v], vertex_salience(v), v, ver[v]});
                    }
                }
            }
            igraph_vector_int_destroy(&inc);
            if (advanced || (!pq.empty() && is_neg_(pq.top().d))) {
                cur_mminus = negative_edge_count();
                if (cur_mminus <= best_mminus) { best_mminus = cur_mminus; best_s = s_; }
            }
        }

        // B) Zero-clique step (strictly positive edges; weighted by salience)
        advanced = false;
        ++r;
        if (r <= R) {
            std::vector<int> Z0; Z0.reserve(n);
            for (int u = 0; u < n; ++u) if (std::fabs(dtilde[u]) <= EPS) Z0.push_back(u);

            std::cout << "[GREEDY-B] r=" << r
            		  << " |Z0|=" << Z0.size()
                      << " cur_minus=" << cur_mminus
                      << " best_minus=" << best_mminus
                      << " neg edges=" << negative_edge_count()
                      << " project=" << integer_projection().negative_edge_count()
                      << "\n";

            if (!Z0.empty()) {
                auto Q = maximal_salience_clique_strict_pos(Z0);

                std::cout << "[GREEDY-B] |Q|=" << Q.size() << "\n";

                if ((int)Q.size() >= 2) {
                    for (int v: Q) {
                        flip_and_refresh_net_degree(v, dtilde);
                    }
                    advanced = true;
                }
            }
        }
        if (advanced) continue;

        // C) Integer detour (only if s is nonintegral; length ≤ L_max; gate d~(u*) ≤ Δ)
        if (r <= R && has_fractional_switching()) {
            // 1) integer projection and a short integer greedy pass there
            SignedGraphForMIP SGint = integer_projection();
            GreedyKickOptions int_opts = opts; int_opts.R_max = 0;
            SGint.apply_greedy_switching(int_opts); // integer s_, so this naturally skips its own Step C

            // 2) compute d_int on SGint from its degree buffers (binary net degree)
            std::vector<double> d_int = SGint.net_degrees();

            // 3) short replay sequence S from a weighted Z0-int clique (strictly positive edges in SGint)
            std::vector<int> S; S.reserve(R);
            for (int t = 0; t < R; ++t) {
                std::vector<int> Z0i; Z0i.reserve(n);
                for (int u = 0; u < n; ++u) if (std::fabs(d_int[u]) <= EPS) Z0i.push_back(u);
                if (Z0i.empty()) break;
                auto Qi = SGint.maximal_salience_clique_strict_pos(Z0i);
                if ((int)Qi.size() < 2) break;
                const int ustar = Qi.front();

                // gate on CURRENT fractional net degree
                if (dtilde[ustar] > (double)opts.Delta) break;

                // flip u* in SGint and refresh d_int locally
                igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
                S.push_back(ustar);
                SGint.single_switching(ustar, &inc);
                d_int[ustar] = SGint.d_pos[ustar] - SGint.d_neg[ustar];
                for (int k = 0; k < (int)igraph_vector_int_size(&inc); ++k) {
                    const int v  = IGRAPH_OTHER(G, VECTOR(inc)[k], ustar);
                    d_int[v] = SGint.d_pos[v] - SGint.d_neg[v];
                }
                igraph_vector_int_destroy(&inc);
            }

            // 4) Replay S on the working state; accept iff any d~ becomes negative
            bool created_neg = false; std::vector<int> replayed; replayed.reserve(S.size());
            for (int u : S) {
                flip_and_refresh_net_degree(u, dtilde);
                replayed.push_back(u);
                if (dtilde[u] < -EPS) { created_neg = true; break; }
                igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
                igraph_incident(G, &inc, u, IGRAPH_ALL);
                for (int k = 0; k < (int)igraph_vector_int_size(&inc); ++k) {
                    const int v = IGRAPH_OTHER(G, VECTOR(inc)[k], u);
                    if (dtilde[v] < -EPS) { created_neg = true; break; }
                }
                igraph_vector_int_destroy(&inc);
                if (created_neg) break;
            }

            if (created_neg) {
                advanced = true;
            } else {
                for (int i = (int)replayed.size() - 1; i >= 0; --i) flip_and_refresh_net_degree(replayed[i], dtilde);
            }
        }
    }

    // commit incumbent
    s_ = std::move(best_s);
    compute_degrees();
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
// Re-seed from x̂ (apply via per-vertex alignment to keep logic simple)
// ─────────────────────────────────────────────────────────────

void SignedGraphForMIP::reseed_switching(const std::vector<double>& xhat) {
    if ((int)xhat.size() != vertex_count())
        throw std::runtime_error("xhat size mismatch");

    // Rebuild s_ directly from x̂ (no per-vertex alignment; preserve original logic)
    for (int u = 0; u < vertex_count(); ++u) {
        double xu = xhat[u];
        if (xu < 0.0) xu = 0.0;
        else if (xu > 1.0) xu = 1.0;
        s_[u] = 2.0 * xu - 1.0;
    }
    // Refresh aggregates once; no salience arrays exist anymore
    compute_degrees();
}

// ─────────────────────────────────────────────────────────────
// Negative-cycle finder stream (unchanged API)
// ─────────────────────────────────────────────────────────────

NegativeCycleBatch SignedGraphForMIP::open_negative_cycle_stream(bool cover,
                                                                 bool use_triangle_order) const {
    return NegativeCycleBatch(*this, cover, use_triangle_order);
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
