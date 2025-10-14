// File: signed_graph_mip.cpp
#include "signed_graph_mip.h"
#include "negative_cycle_batch.h"

#include <algorithm>
#include <numeric>
#include <limits>
#include <cmath>

// ─────────────────────────────────────────────────────────────
// Constructors / destructor
// ─────────────────────────────────────────────────────────────

SignedGraphForMIP::SignedGraphForMIP(const SignedGraph* const other)
    : SignedGraph(other),
      frac_weights(other->edge_count(), 0.0),
      mask_weights(other->edge_count(), 0.0) {
    // Build local edge -> eid map (canonicalized Edge keys)
    const auto& edge_map = this->edge_index();
    for (const auto& kv : edge_map) {
        const Edge& e = kv.first;
        igraph_integer_t eid = kv.second;
        edge_to_eid[e] = eid;
    }
}

SignedGraphForMIP::SignedGraphForMIP(const SignedGraph* const other, std::vector<double> new_weights)
    : SignedGraph(other, std::move(new_weights)),
      frac_weights(other->edge_count(), 0.0),
      mask_weights(other->edge_count(), 0.0) {
    const auto& edge_map = this->edge_index();
    for (const auto& kv : edge_map) {
        const Edge& e = kv.first;
        igraph_integer_t eid = kv.second;
        edge_to_eid[e] = eid;
    }
}

SignedGraphForMIP::SignedGraphForMIP(const std::string& file_path)
    : SignedGraph(file_path) {
    const igraph_integer_t m = edge_count();
    // Build local edge -> eid map (canonicalized Edge keys)
    for (igraph_integer_t eid = 0; eid < m; ++eid) {
        igraph_integer_t u, v;
        igraph_edge(&g, eid, &u, &v);
        edge_to_eid[Edge{static_cast<int>(u), static_cast<int>(v)}] = eid;
    }
    // Initialize working arrays
    frac_weights.assign((size_t)m, 0.0);
    mask_weights.assign((size_t)m, 0.0);
    salience_full_.clear();
}

SignedGraphForMIP::~SignedGraphForMIP() = default;

// ─────────────────────────────────────────────────────────────
// Fractional guidance → salience (edge-aligned)
// ─────────────────────────────────────────────────────────────

bool SignedGraphForMIP::weighting_from_fractional(const std::vector<double>& x,
                                                  const std::vector<double>& y) {
    const igraph_integer_t m = edge_count();
    if ((int)frac_weights.size() != m)   frac_weights.assign((size_t)m, 0.0);
    if ((int)mask_weights.size() != m)   mask_weights.assign((size_t)m, 0.0);
    // keep salience_full_ empty; edge_salience_view() will return mask_weights if empty

    const auto signs = signs_view();

    bool changed = false;
    for (igraph_integer_t eid = 0; eid < m; ++eid) {
        igraph_integer_t u, v; igraph_edge(&g, eid, &u, &v);
        const int s = signs[eid].sign; // ±1

        // τ(x,y) = 4y_uv − 2x_u − 2x_v + 1 ∈ [-1,1]
        const double tau = 4.0 * y[(size_t)eid]
                         - 2.0 * ((size_t)u < x.size() ? x[(size_t)u] : 0.0)
                         - 2.0 * ((size_t)v < x.size() ? x[(size_t)v] : 0.0)
                         + 1.0;

        // fractional "weight" in [0,2], then map to [0,1]
        const double frac = 1.0 - tau * static_cast<double>(s); // ∈ [0,2]
        frac_weights[(size_t)eid] = frac;

        const double t01 = 0.5 * frac; // ∈ [0,1]
        // salience: 1 at 0.5, fades to 0 at 0 or 1
        double sal = 1.0 - std::min(1.0, 2.0 * std::fabs(t01 - 0.5));
        if (sal < 0.0) sal = 0.0;
        if (sal > 1.0) sal = 1.0;
        mask_weights[(size_t)eid] = sal;

        // crude change detector (optional, used by callers heuristically)
        const double y_hat = ( (size_t)u < x.size() && (size_t)v < x.size() )
                           ? (x[(size_t)u] * x[(size_t)v])
                           : 0.0;
        if (std::fabs(y[(size_t)eid] - y_hat) > 1e-5) changed = true;
    }
    return changed;
}

// ─────────────────────────────────────────────────────────────
// Switching helpers
// ─────────────────────────────────────────────────────────────

const std::vector<int> SignedGraphForMIP::greedy_switching() {
    SignedGraph::GreedyKickOptions opts; // pure-integer pass
    return this->greedy_switching_base(/*cmp_fn=*/nullptr, opts);
}

std::optional<std::shared_ptr<const std::vector<int>>>
SignedGraphForMIP::fractional_greedy_switching(const SignedGraph::GreedyKickOptions& user_opts) {
    auto opts = user_opts; // copy (caller may have already set frac_x/frac_y)
    auto sp = std::make_shared<const std::vector<int>>(
        this->greedy_switching_base(/*cmp_fn=*/nullptr, opts));
    return sp;
}

std::optional<std::shared_ptr<const std::vector<int>>>
SignedGraphForMIP::fractional_greedy_switching() {
    SignedGraph::GreedyKickOptions opts;
    auto sp = std::make_shared<const std::vector<int>>(
        this->greedy_switching_base(/*cmp_fn=*/nullptr, opts));
    return sp;
}

std::optional<std::shared_ptr<const std::vector<int>>>
SignedGraphForMIP::fractional_greedy_switching(const std::vector<double>& x,
                                               const std::vector<double>& y) {
    SignedGraph::GreedyKickOptions opts;
    opts.frac_x = &x;
    opts.frac_y = &y;
    auto sp = std::make_shared<const std::vector<int>>(
        this->greedy_switching_base(/*cmp_fn=*/nullptr, opts));
    return sp;
}

// ─────────────────────────────────────────────────────────────
// Negative-cycle finder wrappers
// ─────────────────────────────────────────────────────────────

NegativeCycleBatch SignedGraphForMIP::open_negative_cycle_stream(bool cover,
                                                                 bool use_triangle_order) const {
    return NegativeCycleBatch(*this, cover, use_triangle_order);
}

std::vector<NegativeCycle>
SignedGraphForMIP::find_switched_lower_bound(bool cover) {
    std::vector<NegativeCycle> out;
    auto stream = open_negative_cycle_stream(cover);
    std::vector<NegativeCycle> batch;
    while (stream.next(batch)) {
        out.insert(out.end(),
                   std::make_move_iterator(batch.begin()),
                   std::make_move_iterator(batch.end()));
    }
    return out;
}

std::vector<std::vector<NegativeCycle>>
SignedGraphForMIP::find_switched_lower_bound_grouped(bool cover) const {
    std::vector<std::vector<NegativeCycle>> groups;
    auto stream = open_negative_cycle_stream(cover);
    std::vector<NegativeCycle> batch;
    while (stream.next(batch)) {
        groups.emplace_back();
        groups.back().swap(batch);
    }
    return groups;
}
