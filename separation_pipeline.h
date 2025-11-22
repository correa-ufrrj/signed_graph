// separation_pipeline.h
#pragma once

#include "fm_debug.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include <cstdint>
#include <string>
#include <utility>
#include <cmath>

#include "signed_graph_mip.h"      // SignedGraphForMIP
#include "cycle_key.h"             // fmkey::CycleKey
#include "separation_config.h"     // SeparationConfig + adapters
#include "triangle_bucket_batch.h"
#include "negative_cycle_batch.h"   // for SP stage
#include "separation_pipeline_tls.h"

// Extern hooks used by TriangleBucketBatch. We provide strong definitions
// in separation_pipeline.cpp and declare them as friends inside the class so
// they can call the pipeline's private hook methods.
extern "C" {
void TBB_on_emit(int edge_id, double used_density);
void TBB_on_accept(int edge_id, double density);
int  TBB_budget_override(int base);
}

// This header carves out the persistent/per-round state and declares the round
// driver shell. There is NO behavior change yet. Methods are stubs and will be
// filled in during steps 2–8 of the incremental plan.


struct SeparationPersistent {
    // Edge-aligned arrays over the FULL graph (size = |E|)
    std::vector<double> omega;         // persistent base weights ω ≥ eps (init 1.0)
    std::vector<double> H;             // repulsion EMA ledger

    // Stateless de-dup (switching-agnostic keys)
    std::unordered_set<fmkey::CycleKey, fmkey::CycleKeyHash, fmkey::CycleKeyEq> in_model_keys;  // already added
    std::unordered_set<fmkey::CycleKey, fmkey::CycleKeyHash, fmkey::CycleKeyEq> recent_keys;    // seen recently

    // Annealing state (triangle budget across rounds)
    int    B_tri_cur        = 0;   // effective budget used by the next triangle stage
    double bar_v_triangle   = 0.0; // running mean/EMA of last-round triangle violation

    // Convenience: initialize sizes once when graph is known
    void init_sizes_if_needed(const SignedGraphForMIP& G) {
        const int m = G.edge_count();
        if ((int)omega.size()      != m) omega.assign(m, 1.0);
        if ((int)H.size()          != m) H.assign(m, 0.0);
    }
};

struct SeparationRound {
    // Edge-aligned arrays over FULL graph (size = |E|)
    std::vector<double> selected_count_full;  // per-round density increments (∑ 1/|C| for edges in accepted cuts)
    std::vector<double> sal_full;             // salience over FULL edges (normalized [0,1])

    // Edge-aligned arrays over the FULL graph (size = |E|)
    std::vector<double> used_in_batch_full;    // within-batch “used” counters on FULL eids

    // Bookkeeping (diagnostics only in step 1)
    int emitted_triangles = 0;
    int emitted_cycles    = 0;

    void clear() {
        emitted_triangles = emitted_cycles = 0;
        selected_count_full.clear();
        used_in_batch_full.clear();
    }
};

// store the switching used for this round
struct RoundSwitching {
    std::shared_ptr<const std::vector<int>> s;  // partition {+1,-1}
    bool from_fractional = false;               // true if came from fractional_greedy_switching
};

class TriangleCyclePipeline {
public:

    TriangleCyclePipeline(SeparationPersistent& persistent,
                          SeparationConfig cfg = {});

    // STEP 1: shell only; returns an empty set. Later steps will fill it.
    struct Result {
        std::vector<fmkey::CycleKey> accepted_keys; // placeholder; material rows will be built by the caller
        int triangles_accepted = 0;
        int cycles_accepted    = 0;
    };

    Result run_round(const SignedGraphForMIP& G_);

    const SeparationRound& round() const { return round_; }
    SeparationRound&       round()       { return round_; }

private:
    // Hooks called by TBB free functions (declared friend at end of class)
    void on_emit_(int full_eid, double used_density);
    void on_accept_(int full_eid, double density);
    int  override_budget_(int base) const;

    // Storage of triangle stage outputs for this round
    std::vector<TriangleBucketBatch::Candidate> tri_selected_;
    GraphCore::Bitmap uncovered_neg_eids_;  // bits set = NEG edges still uncovered

    // Per-round exported keys (triangles now; SP in later step). Cleared at run start.
    std::vector<fmkey::CycleKey> accepted_keys_;

    // Allow TBB C hooks to call our private methods
    friend void ::TBB_on_emit(int, double);
    friend void ::TBB_on_accept(int, double);
    friend int  ::TBB_budget_override(int);
    // Core references
    SeparationPersistent&   S_;
    SeparationConfig        C_;

    // Per-round scratch
    SeparationRound          round_;
    
    // debug: keep a few FULL eids that were actually used/accepted this round
	std::vector<int> dbg_sample_feids_;

    // === Step-2+ stubs (no-ops in step 1) ===
    void seed_sp_weights_(const SignedGraphForMIP& G);
    void triangle_stage_(const SignedGraphForMIP& G_);
    void sp_stage_(const SignedGraphForMIP& G_);
    void commit_stage_(const SignedGraphForMIP& G_);
    void resize_round_full_arrays_(int m_full);    
};