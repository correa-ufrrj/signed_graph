#pragma once

#include <signed_graph_mip.h>

// --- On-demand negative-cycle stream (declaration only) ---------------------
class NegativeCycleBatchStream {
public:
    explicit NegativeCycleBatchStream(const SignedGraphForMIP& G, bool cover, bool use_triangle_order = false);
    bool next(std::vector<NegativeCycle>& out);

    int  total_cycles_emitted() const;
    int  batches_emitted()      const;
	void set_lp_scores_full_edges(const std::vector<double>& s, double alpha=1.0, double beta=0.3);
	

private:
    const SignedGraphForMIP& G_;
    igraph_integer_t vcount_;
    igraph_integer_t ecount_;
    std::vector<double> base_pos_;
    std::vector<double> neg_hard_;
    std::vector<double> reuse_accum_;
    double              med_base_pos_;
    std::vector<Edge>   neg_edges_;
    std::vector<Edge>   disconnected_;
    igraph_vector_t     saved_weights_;
    bool                saved_weights_init_;
    bool                cover_;
    bool                finished_;
    size_t              total_found_;
    int                 batches_emitted_;
    // triangle-aware ordering (enabled only when requested)
    bool               use_tri_order_ = false;
    std::vector<int>   neg_tri_vert_;  // per-vertex negative triangle counts
    inline int edge_tri_score_(igraph_integer_t eid) const;
    // Cross-batch edge penalty = scale + 1/|C|.
    // Default scale = 1.0 for initial experiments.
    double cross_batch_penalty_scale_ = 3.0;
    // Positive-only graph + mappings
    igraph_t g_pos_{};
    bool g_pos_built_ = false;
    std::vector<igraph_integer_t> full2pos_eid_; // size = ecount_, -1 for negative edges
    std::vector<igraph_integer_t> pos2full_eid_; // size = ecount_pos
    // Triangle-per-vertex soft cap (in-batch); set small for effectiveness
    int tri_cap_per_vertex_ = 6;
    std::vector<int> tri_used_per_vertex_; // sized to vcount_ in build_initial_state_
    std::vector<int> neg_deg_; // per-vertex degree over negative edges
    std::vector<char> neg_edge_covered_; // 0/1 per negative edge id
    std::vector<int> pos_deg_; // per-vertex positive degree (stored for bucketing)

    // Weights for positive-only graph (parallel to pos edges)
    igraph_vector_t saved_weights_pos_;
    bool saved_weights_pos_init_ = false;

    // New knobs for selection strategy
    int    K_tri_per_neg_ = 3;       // K=1 => strict disjointness within batch
    double overlap_penalty_gamma_ = 0.25; // penalize reused positive edges proportionally
	std::vector<double> lp_score_full_;
	bool use_lp_weights_ = false;
    double lp_alpha_ = 1.0, lp_beta_ = 0.3;       // persist LP blending across batches

    std::vector<int> tri_cap_per_v_;              // adaptive per-vertex triangle caps (member)

    inline void bump_cross_batch_(int eid, int cycle_len) {
        if (cycle_len <= 0) return;
        if (eid >= 0 && eid < static_cast<int>(base_pos_.size())) {
            base_pos_[eid] += cross_batch_penalty_scale_ + 1.0 / static_cast<double>(cycle_len);
        }
    }

    inline bool edge_is_pos(igraph_integer_t eid) const;
    void build_initial_state_();
    void build_mask_for_batch_();
};
