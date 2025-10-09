// signed_graph_mip.h
#pragma once

#include "signed_graph.h"
#include <optional>
#include <vector>
#include <chrono>

// Forward-declare the stream class
class NegativeCycleBatchStream;

class SignedGraphForMIP : public SignedGraph {
private:
    mutable std::unordered_map<Edge, igraph_integer_t, EdgeHash> edge_to_eid;
    std::vector<double> frac_weights;
    std::vector<double> mask_weights;     // |frac - 0.5| (edge-aligned)
    std::vector<double> salience_full_;   // 1 - min(1, 2*|frac-0.5|) in [0,1]

public:
    friend class ::NegativeCycleBatchStream; // stream needs internals

	const std::vector<double>& get_mask_weights() const { return mask_weights; }

    SignedGraphForMIP(const std::string& file_path);
    SignedGraphForMIP(const SignedGraph* const other);
    SignedGraphForMIP(const SignedGraph* const other, std::vector<double> new_weights);
    ~SignedGraphForMIP();
    
    // Prefer normalized [0,1] salience; fallback to mask_weights if empty
    const std::vector<double>& edge_salience_view() const {
        return salience_full_.empty() ? mask_weights : salience_full_;
    };

    bool weighting_from_fractional(const std::vector<double>& x, const std::vector<double>& y);
    const std::vector<int> greedy_switching();
    std::optional<std::shared_ptr<const std::vector<int>>> fractional_greedy_switching(const SignedGraph::GreedyKickOptions& opts);
    std::optional<std::shared_ptr<const std::vector<int>>> fractional_greedy_switching();

    // Public “finder” APIs (now wrappers over the stream)
    std::vector<NegativeCycle> find_switched_lower_bound(bool cover = false);
    std::vector<std::vector<NegativeCycle>> find_switched_lower_bound_grouped(bool cover) const;

    // Stream factory
    NegativeCycleBatchStream open_negative_cycle_stream(bool cover, bool use_triangle_order = false) const;
};
