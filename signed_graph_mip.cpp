// File: signed_graph_mip.cpp
#include "signed_graph_mip.h"
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <set>
#include <boost/heap/pairing_heap.hpp>
#include <numeric>

// Correct signature for igraph ≥ 0.10.x
void silent_warning_handler(const char* reason, const char* file, int line) {
    // Suppress all warnings
}

// SignedGraphForMIP default clone constructor
SignedGraphForMIP::SignedGraphForMIP(const SignedGraph* const other)
  : SignedGraph(other), frac_weights(other->edge_count(), 0.0),
    mask_weights(other->edge_count(), 0.0) {
    const auto& edge_map = this->edge_index();
    for (const auto& [edge, eid] : edge_map) {
        edge_to_eid[edge] = eid;
    }
}

// SignedGraphForMIP weight-based clone constructor
SignedGraphForMIP::SignedGraphForMIP(const SignedGraph* const other, std::vector<double> new_weights)
  : SignedGraph(other, std::move(new_weights)), frac_weights(other->edge_count(), 0.0),
    mask_weights(other->edge_count(), 0.0) {
    const auto& edge_map = this->edge_index();
    for (const auto& [edge, eid] : edge_map) {
        edge_to_eid[edge] = eid;
    }
}

SignedGraphForMIP::SignedGraphForMIP(const std::string& file_path)
    : SignedGraph(file_path)
{
    for (int eid = 0; eid < edge_count(); ++eid) {
        igraph_integer_t u, v;
        igraph_edge(&g, eid, &u, &v);
        edge_to_eid[{static_cast<int>(u), static_cast<int>(v)}] = eid;
    }

    frac_weights.resize(weights.size());
	mask_weights.resize(weights.size());
	std::transform(weights.begin(), weights.end(), frac_weights.begin(),
	               [](double v){ return static_cast<double>(v); });
	std::fill(mask_weights.begin(), mask_weights.end(), 0.0);
}

SignedGraphForMIP::~SignedGraphForMIP() {
}

bool SignedGraphForMIP::weighting_from_fractional(const std::vector<double>& x,
                                                  const std::vector<double>& y) {
    igraph_integer_t ecount = edge_count();
    bool has_changed = false;
    auto signs = signs_view();

    if (mask_weights.size() != static_cast<size_t>(ecount))
        mask_weights.assign(ecount, 0.0);

    for (igraph_integer_t eid = 0; eid < ecount; ++eid) {
        igraph_integer_t from, to; igraph_edge(&g, eid, &from, &to);
        const int old_w = signs[eid].sign; // ±1 (or ±0.0 treated consistently)
        // frac_weights in [0,2] since tau∈[-1,1]; map to [0,1] via 0.5×
        const double tau = 4.0 * y[eid] - 2.0 * x[from] - 2.0 * x[to] + 1.0; // ∈[-1,1]
        frac_weights[eid] = 1.0 - tau * old_w;            // ∈[0,2]
        const double t01 = 0.5 * frac_weights[eid];       // ∈[0,1]
        // salience: 1 at 0.5, fades to 0 at 0 or 1
        mask_weights[eid] = 1.0 - std::min(1.0, 2.0 * std::abs(t01 - 0.5));

        has_changed |= std::abs(y[eid] - x[from] * x[to]) > 1e-5;
    }
    return has_changed;
}

const std::vector<int> SignedGraphForMIP::greedy_switching() {
    SignedGraph::GreedyKickOptions opts;  // pure integer greedy: no fractional guidance
    return this->greedy_switching_base(/*cmp_fn=*/nullptr, opts);
}

// ---- fractional_greedy_switching overloads (SignedGraphForMIP) ----
std::optional<std::shared_ptr<const std::vector<int>>>
SignedGraphForMIP::fractional_greedy_switching(const SignedGraph::GreedyKickOptions& user_opts) {
    SignedGraph::GreedyKickOptions opts = user_opts; // local copy
    // If you have LP vectors here, set:
    //   opts.frac_x = &lp_x;
    //   opts.frac_y = &lp_y;  // optional
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

std::optional<std::shared_ptr<const std::vector<int>>>
SignedGraphForMIP::fractional_greedy_switching() {
    SignedGraph::GreedyKickOptions opts;
    auto sp = std::make_shared<const std::vector<int>>(
        this->greedy_switching_base(/*cmp_fn=*/nullptr, opts));
    return sp;
}

// --- SignedGraphForMIP wrappers over the stream ---
std::vector<NegativeCycle>
SignedGraphForMIP::find_switched_lower_bound(bool cover) {
    std::vector<NegativeCycle> flat;
    auto stream = open_negative_cycle_stream(cover);
    std::vector<NegativeCycle> batch;
    while (stream.next(batch)) {
        flat.insert(flat.end(),
                    std::make_move_iterator(batch.begin()),
                    std::make_move_iterator(batch.end()));
    }
    return flat;
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
