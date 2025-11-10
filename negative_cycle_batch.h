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

// Forward declarations (C linkage) for TBB hooks used by TriangleBucketBatch.
// These are defined with strong linkage in negative_cycle_batch.cpp.
extern "C" {
void TBB_on_emit(int, double);
void TBB_on_accept(int, double);
int  TBB_budget_override(int);
}

class NegativeCycleBatch {
public:
    // Minimal, decoupled parameter pack for SP stage (mirrors TBB pattern)
    struct Params {
        int    B_sp               = 0;   // per-batch SP budget
        int    K_sp_per_neg       = 3;    // alt paths per negative anchor
        int    sp_cap_per_vertex  = 8;    // per-vertex cap for accepted SP cycles
        int    tri_cap_per_vertex = 6;    // cap reused when NCB enforces per-vertex limits in triangle-first accept
        double cross_batch_penalty_scale = 3.0; // drift base_pos_ when cycles are accepted
        TriangleBucketBatch::Params tbb{}; // pass-through for the internal triangle stage

        // NEW: optional per-full-edge positive lengths (e.g., ω′ or a function of salience).
        // Size must be m = |E|; only entries for positive edges are used.
        // If provided, these values are copied into saved_weights_pos_ and
        // used directly by Dijkstra for path finding.
        const std::vector<double>* guide_len_full = nullptr;

	    // new: use same knobs as pipeline for emit bumps and clamps
	    double beta_emit  = 0.5;   // ← SeparationConfig.weights.beta_emit
	    double omega_eps  = 1e-12; // ← SeparationConfig.weights.omega_eps (min)
	    double omega_max  = 0.0;   // 0 = no cap; else clamp to this

        // Min length for numerical stability
        double guide_len_eps = 1e-12;
    };

// NegativeCycleBatch.h
    template <class Sink>
    void flush_emitted_to(Sink&& sink) {
        // used_in_batch_pos_ is indexed by pos-edge id (peid) in g_pos_
        const int Epos = (int)used_in_batch_pos_.size();
        for (int peid = 0; peid < Epos; ++peid) {
            double d = used_in_batch_pos_[peid];
            if (d <= 0.0) continue;
            const int full_eid = pos2full_eid_[peid];
            sink(full_eid, d);
            used_in_batch_pos_[peid] = 0.0; // clear for next batch
        }
    }

	NegativeCycleBatch(const SignedGraphForMIP& G,
	                   const std::vector<Edge>& neg_edges_uncov,
	                   Params p);
	// negative_cycle_batch.h  (only the constructor declarations shown here)
	explicit NegativeCycleBatch(const SignedGraphForMIP& G);
	
	explicit NegativeCycleBatch(const SignedGraphForMIP& G,
	                            Params p);
    bool next(std::vector<NegativeCycle>& out);

    void set_params(const Params& p) { P_ = p; }
    const Params& params() const     { return P_; }

    int  total_cycles_emitted() const;
    int  batches_emitted()      const;

    // Stateless de-dup & reheat moved to the pipeline driver; NCB no longer owns them.
    // Read-only view of per-positive-edge usage density accumulated in this batch.
    const std::vector<double>& pos_usage_density() const { return used_in_batch_pos_; }
    // Accumulate the pos-graph usage density into a full-edge vector (dst is resized if needed).
    void accumulate_pos_usage_to_full(std::vector<double>& dst, double scale = 1.0) const;

    // Grant the global C hooks friend access to call private internals
    friend void ::TBB_on_emit(int, double);
    friend void ::TBB_on_accept(int, double);
    friend int  ::TBB_budget_override(int);
	friend void ncb_emit(void*, int, double);
	friend void ncb_accept(void*, int, double);
	friend int  ncb_budget(void*, int);

private:
    Params P_{};
    const SignedGraphForMIP& G_;
    igraph_integer_t vcount_{0};
    igraph_integer_t ecount_{0};

    // Persistent base weights on full graph edges (positive edges only used)
    std::vector<double> base_pos_;
    std::vector<double> neg_hard_;
    std::vector<double> reuse_accum_;

    // Negative edges (anchors) under current switching; may shrink when cover_==true
    std::vector<Edge>   neg_edges_;

    // Saved working masks (full graph and pos-only)
    igraph_vector_t     saved_weights_{};
    bool                saved_weights_init_{false};

    bool                finished_{false};
    size_t              total_found_{0};
    int                 batches_emitted_{0};

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

    // Within-batch usage density on positive edges (pos-graph indexing)
    std::vector<double> used_in_batch_pos_;

    // Scoring/LP blending (carried across batches)
    
    inline void bump_cross_batch_(int eid, int cycle_len) {
        if (cycle_len <= 0) return;
        if (eid >= 0 && eid < static_cast<int>(base_pos_.size())) {
            base_pos_[(size_t)eid] += P_.cross_batch_penalty_scale + 1.0 / static_cast<double>(cycle_len);
        }
    }
	double best_twohop_upper_bound_(int u, int v) const;

    void build_initial_state_();
    void build_mask_for_batch_();

    // ---------- Bridges for within-batch/cross-batch updates ------------
    void on_emit_(int full_eid, double used_density);   // within-batch ω'_t bump
    void on_accept_(int full_eid, double density);      // cross-batch ω drift
    int  override_budget_(int base) const;              // annealing hook (kept simple)
};