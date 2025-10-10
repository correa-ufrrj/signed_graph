// negative_cycle_batch.h
#pragma once

#include <signed_graph_mip.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "triangle_bucket_batch.h"

class NegativeCycleBatch {
public:
    explicit NegativeCycleBatch(const SignedGraphForMIP& G, bool cover, bool use_triangle_order = false);
    bool next(std::vector<NegativeCycle>& out);

    int  total_cycles_emitted() const;
    int  batches_emitted()      const;
    void set_lp_scores_full_edges(const std::vector<double>& s, double alpha=1.0, double beta=0.3);

    // Grant the global C hooks friend access to call private internals
    friend void TBB_on_emit(int, double);
    friend void TBB_on_accept(int, double);
    friend int  TBB_budget_override(int);

private:
    const SignedGraphForMIP& G_;
    igraph_integer_t vcount_{0};
    igraph_integer_t ecount_{0};

    // Persistent base weights on full graph edges (positive edges only used)
    std::vector<double> base_pos_;
    std::vector<double> neg_hard_;
    std::vector<double> reuse_accum_;
    double              med_base_pos_{1.0};

    // Negative edges (anchors) under current switching; may shrink when cover_==true
    std::vector<Edge>   neg_edges_;
    std::vector<Edge>   disconnected_;

    // Saved working masks (full graph and pos-only)
    igraph_vector_t     saved_weights_{};
    bool                saved_weights_init_{false};

    bool                cover_{false};
    bool                finished_{false};
    size_t              total_found_{0};
    int                 batches_emitted_{0};

    // triangle-aware ordering (kept for compatibility; triangles now run via TBB)
    bool               use_tri_order_ = false;
    std::vector<int>   neg_tri_vert_;  // per-vertex negative triangle counts
    inline int edge_tri_score_(igraph_integer_t /*eid*/) const { return 0; } // legacy no-op

    // Cross-batch: bump base_pos_ when a cycle/triangle is accepted
    double cross_batch_penalty_scale_ = 3.0;

    // Positive-only graph + mappings
    igraph_t g_pos_{};
    bool g_pos_built_ = false;
    std::vector<igraph_integer_t> full2pos_eid_; // size = ecount_, -1 for negative edges
    std::vector<igraph_integer_t> pos2full_eid_; // size = ecount_pos

    // Triangle-per-vertex cap (uniform for TBB)
    int tri_cap_per_vertex_ = 6;
    std::vector<int> tri_used_per_vertex_; // sized to vcount_ in build_initial_state_
    std::vector<int> neg_deg_;             // per-vertex degree over negative edges
    std::vector<char> neg_edge_covered_;   // 0/1 per negative full-eid
    std::vector<int> pos_deg_;             // per-vertex positive degree

    // Working weights (pos-only) for Dijkstra
    igraph_vector_t saved_weights_pos_{};
    bool saved_weights_pos_init_ = false;

    // Scoring/LP blending (carried across batches)
    int    K_tri_per_neg_ = 3;
    double overlap_penalty_gamma_ = 0.25;
    std::vector<double> lp_score_full_;
    bool   use_lp_weights_ = false;
    double lp_alpha_ = 1.0, lp_beta_ = 0.3;

    std::vector<int> tri_cap_per_v_; // adaptive (legacy), still used for SP acceptance

    inline void bump_cross_batch_(int eid, int cycle_len) {
        if (cycle_len <= 0) return;
        if (eid >= 0 && eid < static_cast<int>(base_pos_.size())) {
            base_pos_[(size_t)eid] += cross_batch_penalty_scale_ + 1.0 / static_cast<double>(cycle_len);
        }
    }

    inline bool edge_is_pos(igraph_integer_t eid) const {
        return G_.switched_weights[eid] > 0.0;
    }

    void build_initial_state_();
    void build_mask_for_batch_();

    // ---------- Triangle-first batch (via TriangleBucketBatch) ----------
    bool run_triangle_first_batch_(std::vector<NegativeCycle>& out, std::vector<int>& covered_neg_eids);

    void build_pos_adj_and_index_(
        TriangleBucketBatch::PosAdj& pos_adj,
        TriangleBucketBatch::EdgeIndex& edge_index) const;

    // ---------- Bridges for within-batch/cross-batch updates ------------
    void on_emit_(int full_eid, double used_density);   // within-batch ω'_t bump
    void on_accept_(int full_eid, double density);      // cross-batch ω drift
    int  override_budget_(int base) const;              // annealing hook (kept simple)

    // allow k alternate shortest paths per neg edge (soft mask)
    int    K_sp_per_neg_ = 3;
    double alt_path_bump_ = 1.0; // scaled from med_base_pos_ at build
};

// --- global bridges implemented in .cpp ---
// (TriangleBucketBatch links against these)
// triangle_bucket_batch.h
extern "C" {
void TBB_on_emit(int /*edge_id*/, double /*used_density*/);
void TBB_on_accept(int /*edge_id*/, double /*density*/);
int  TBB_budget_override(int base);
}
