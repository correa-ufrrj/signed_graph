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

SignedGraphForMIP::SignedGraphForMIP(const SignedGraphForMIP* const other)
    : SignedGraph(other),
	   reference_s_(other->reference_s_),
	   reference_weights_(other->reference_weights_) {
    is_switching_ = true;
}

SignedGraphForMIP::SignedGraphForMIP(const SignedGraphForMIP* const other, std::vector<double> new_weights)
    : SignedGraph(other, new_weights),
	   reference_s_(other->reference_s_),
	   reference_weights_(other->reference_weights_) {
    is_switching_ = true;
}

SignedGraphForMIP::SignedGraphForMIP(const SignedGraphForMIP* const other,
                                     std::vector<double> new_weights,
                                     std::vector<double> new_s)
: SignedGraph(other, std::move(new_weights), std::move(new_s)),
  reference_s_(other->reference_s_),
  reference_weights_(other->reference_weights_)
{
    is_switching_ = true;
}

SignedGraphForMIP::SignedGraphForMIP(const std::string& file_path)
    : SignedGraph(file_path),
      reference_s_(s_), reference_weights_(weights_) {
    is_switching_ = true;
}

SignedGraphForMIP::~SignedGraphForMIP() = default;

// ─────────────────────────────────────────────────────────────
// Fractional guidance → salience (edge-aligned)
// ─────────────────────────────────────────────────────────────

SignedGraphForMIP SignedGraphForMIP::integer_projection() const {
    const int n = vertex_count();
    const int m = edge_count();

    std::vector<double> s_int(n);
    for (int u = 0; u < n; ++u) {
        s_int[u] = (s_[u] >= 0.0) ? +1.0 : -1.0; // round to {±1}
    }

//	const auto s = signs_view();
    std::vector<double> sigma_int(m);
    for (int eid = 0; eid < m; ++eid) {
		igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
        const double w0 = reference_weights_[eid];         // ±1 from the base graph
        sigma_int[eid] = s_int[(int)u] * w0 * s_int[(int)v];
//		sigma_int[eid] = s[eid].sign; // σ_{s_int} = s_int ⊙ σ_s ⊙ s_int
    }

    SignedGraphForMIP sg(this, sigma_int, s_int);
    sg.is_switching_ = true;
    return sg;
}

void SignedGraphForMIP::align_switching(const int u, const double xfix01) {
    const double xu = get_x(u);           // in {0,1}
    std::cout << "[ALIGN] x[v*]=" << xu << " target=" << xfix01 << "\n";
    if (std::fabs(xu - xfix01) > 1e-9) {
        for (double& sv : s_) sv = -sv;   // global flip keeps edge signs unchanged
    }
    std::cout << "[ALIGNED] x[v*]=" << get_x(u) << " target=" << xfix01 << "\n";
}

// ─────────────────────────────────────────────────────────────
// Switching helpers
// ─────────────────────────────────────────────────────────────

SignedGraphForMIP SignedGraphForMIP::greedy_switching(const GreedyKickOptions& opts) {
    SignedGraphForMIP sg(this, weights_, s_);
    sg.is_switching_ = true;
    sg.apply_greedy_switching(opts);
    return sg;
}

SignedGraphForMIP SignedGraphForMIP::greedy_switching() {
    GreedyKickOptions opts; return greedy_switching(opts);
}

// Rebuild current state from a new x̂ (overwrites current s_, weights)
void SignedGraphForMIP::reseed_switching(const std::vector<double>& xhat) {
    if ((int)xhat.size() != vertex_count())
        throw std::runtime_error("xhat size mismatch");

	weights_ = reference_weights_;
    s_ = reference_s_;

    // s = 2x̂ - 1
    std::vector<double> s(vertex_count());
    for (int u = 0; u < vertex_count(); ++u) {
		s[u] = 2.0 * xhat[u] - 1.0;
		if (xhat[u] >= 0.499 && xhat[u] <= 0.501) std::cout << "."; else std::cout << u;
    }
    std::cout  << "\n";
    
    std::cout << "[SWITCH-RESEED] before sum s_=" << std::accumulate(s_.begin(), s_.end(), 0.0)
    		  << " before sum s=" << std::accumulate(s.begin(), s.end(), 0.0)
    		  << " before sum xhat=" << std::accumulate(xhat.begin(), xhat.end(), 0.0)
    		  << "\n";
    apply_compose_switching(s);
    
    std::cout << "[SWITCH-RESEED] after sum s_=" << std::accumulate(s_.begin(), s_.end(), 0.0)
    		  << "\n";
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
