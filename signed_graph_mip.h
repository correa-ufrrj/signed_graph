// signed_graph_mip.h
#pragma once

#include "signed_graph.h"
#include <optional>
#include <vector>
#include <chrono>

// Forward-declare the stream class
class NegativeCycleBatch;

class SignedGraphForMIP : public SignedGraph {
private:
	const std::vector<double> reference_s_;        // size |V|, entries ±1
	const std::vector<double> reference_weights_; // this instance’s edge weights (±1 here)

public:
    friend class ::NegativeCycleBatch; // stream needs internals
	using SignedGraph::SignedGraph;
  
    SignedGraphForMIP(const std::string& file_path);
    SignedGraphForMIP(const SignedGraphForMIP* const other);
    SignedGraphForMIP(const SignedGraphForMIP* const other, std::vector<double> new_weights);
    SignedGraphForMIP(const SignedGraphForMIP* const other, std::vector<double> new_weights,
    					const std::vector<double> new_s);
    ~SignedGraphForMIP();

	SignedGraphForMIP integer_projection() const;
	void align_switching(const int u, const	double xu);
    SignedGraphForMIP greedy_switching(const GreedyKickOptions& opts);
    SignedGraphForMIP greedy_switching();
	// Rebuild current state from a new x̂ (overwrites current s_, weights)
	void reseed_switching(const std::vector<double>& xhat);
	inline double get_x(int u) { return (s_[u]+1.0)/2.0; }
	inline double get_y(int u, int v) { return get_x(u)*get_x(v); }
	inline double get_y(Edge uv) { return get_y(uv.first, uv.second); }
	
    // Public “finder” APIs (now wrappers over the stream)
    std::vector<NegativeCycle> find_switched_lower_bound(bool cover = false);
    std::vector<std::vector<NegativeCycle>> find_switched_lower_bound_grouped(bool cover) const;

    // Stream factory
    NegativeCycleBatch open_negative_cycle_stream(bool cover, bool use_triangle_order = false) const;
};
