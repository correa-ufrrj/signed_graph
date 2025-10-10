// separation_pipeline.h
#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include <cstdint>
#include <string>
#include <utility>
#include <cmath>

#include "signed_graph_mip.h"      // SignedGraphForMIP
#include "reheat_pool.h"           // ReheatPool
#include "cycle_key.h"             // fmkey::CycleKey

// Extern hooks used by TriangleBucketBatch. We provide strong definitions
// in separation_pipeline.cpp and declare them as friends inside the class so
// they can call the pipeline's private hook methods.
extern "C" {
void TBB_on_emit(int edge_id, double used_density);
void TBB_on_accept(int edge_id, double density);
int  TBB_budget_override(int base);
}


// STEP 1 NOTES
// ------------
// This header carves out the persistent/per-round state and declares the round
// driver shell. There is NO behavior change yet. Methods are stubs and will be
// filled in during steps 2–8 of the incremental plan.

struct SeparationConfig {
    // Ranking knobs
    double alpha        = 0.3;   // weight inversion multiplier for 1/ω′ term
    double theta        = 0.5;   // salience bonus multiplier

    // Working & persistent weight knobs
    double lambda_hist  = 0.20;  // history → log(1+H) multiplier
    double lambda_LP    = 0.0;   // LP guidance (kept 0.0 for step 2)
    double beta_emit    = 0.02;  // within-batch bump magnitude
    double beta_sel     = 0.05;  // cross-batch drift magnitude
    double omega_eps    = 1e-8;  // lower clamp for ω and ω′
    double omega_max    = 64.0;  // upper clamp for ω

    // Budgets / caps (defaults; may be specialized per phase)
    int    B_tri        = 0;     // if 0: auto from #candidates (triangle stage)
    int    K_tri_per_neg= 8;     // per-bucket prefilter cap (triangles)
    int    tri_cap_per_vertex = 6;

    int    B_sp         = 0;     // if 0: auto (SP stage)
    int    K_sp_per_neg = 3;     // per-bucket cap (SP)
    int    sp_cap_per_vertex  = 8;

    // Annealing schedule for B_tri (fractional phase)
    int    B_tri_min    = 1;
    double gamma_min    = 0.70;
    double gamma_max    = 0.92;
    double v0           = 0.10;  // target violation level
    double tau          = 0.10;  // softness

    // Recent-key TTL (rounds)
    int    recent_ttl   = 3;
};

struct SeparationPersistent {
    // Edge-aligned arrays over the FULL graph (size = |E|)
    std::vector<double> omega;         // persistent base weights ω ≥ eps (init 1.0)
    std::vector<double> pool_count;    // density counts across ACCEPTED cuts
    std::vector<double> H;             // repulsion EMA ledger

    // Stateless de-dup (switching-agnostic keys)
    std::unordered_set<fmkey::CycleKey, fmkey::CycleKeyHash, fmkey::CycleKeyEq> in_model_keys;  // already added
    std::unordered_set<fmkey::CycleKey, fmkey::CycleKeyHash, fmkey::CycleKeyEq> recent_keys;    // seen recently

    // Reheat storage (switching-agnostic)
    ReheatPool P_reheat;

    // Annealing state (triangle budget across rounds)
    int    B_tri_cur        = 0;   // effective budget used by the next triangle stage
    double bar_v_triangle   = 0.0; // running mean/EMA of last-round triangle violation

    // EMA parameters for H
    double ema_delta  = 0.85;   // decay (0<delta<1)
    double ema_rho    = 0.02;   // reheated-accepted weight
    double ema_kappa  = 0.005;  // emitted-but-rejected weight

    // Convenience: initialize sizes once when graph is known
    void init_sizes_if_needed(const SignedGraphForMIP& G) {
        const int m = G.edge_count();
        if ((int)omega.size()      != m) omega.assign(m, 1.0);
        if ((int)pool_count.size() != m) pool_count.assign(m, 0.0);
        if ((int)H.size()          != m) H.assign(m, 0.0);
    }
};

struct SeparationRound {
    // Edge-aligned arrays over FULL graph (size = |E|)
    std::vector<double> selected_count_full;  // per-round density increments (∑ 1/|C| for edges in accepted cuts)
    std::vector<double> sal_full;             // salience over FULL edges (normalized [0,1])

    // Edge-aligned arrays over POSITIVE subgraph only (size = |E⁺| for current switching)
    std::vector<double> omega_prime_pos;      // working weights ω′ on E⁺
    std::vector<double> used_in_batch_pos;    // within-batch density counters on E⁺

    // Bookkeeping (diagnostics only in step 1)
    int emitted_triangles = 0;
    int emitted_cycles    = 0;

    void clear() {
        emitted_triangles = emitted_cycles = 0;
        selected_count_full.clear();
        sal_full.clear();
        omega_prime_pos.clear();
        used_in_batch_pos.clear();
    }
};

struct RoundGraphView {
    // maps
    std::vector<int> full2pos;    // |E| → {0..m_pos-1} or -1 if not in E⁺
    std::vector<int> pos2full;    // |E⁺| → full eid

    // masks on FULL edge ids under current switching
    std::vector<char> pos_mask;   // |E|, 1 if edge ∈ E⁺_{σ_s}, else 0
    std::vector<char> neg_mask;   // |E|, 1 if edge ∈ E⁻_{σ_s}, else 0

    int m_pos = 0;                // |E⁺| for fast access
    void clear() {
        full2pos.clear(); pos2full.clear();
        pos_mask.clear(); neg_mask.clear();
        m_pos = 0;
    }
};

// Optional: store the switching used for this round (filled in step 2)
struct RoundSwitching {
    std::shared_ptr<const std::vector<int>> s;  // partition {+1,-1}
    bool from_fractional = false;               // true if came from fractional_greedy_switching
};

class TriangleCyclePipeline {
public:
    enum class Phase { Build, Fractional };

    TriangleCyclePipeline(SignedGraphForMIP& G,
                          SeparationPersistent& persistent,
                          SeparationConfig cfg = {});

    // STEP 1: shell only; returns an empty set. Later steps will fill it.
    struct Result {
        std::vector<fmkey::CycleKey> accepted_keys; // placeholder; material rows will be built by the caller
        int triangles_accepted = 0;
        int cycles_accepted    = 0;
    };

    Result run_round(const std::vector<double>* x_hat,
                     const std::vector<double>* y_hat,
                     Phase phase);

    const SeparationRound& round() const { return round_; }
    SeparationRound&       round()       { return round_; }

    // Precompute full↔pos id maps for current switching; called at round start
    void rebuild_pos_maps();

private:

    // Lightweight record of an accepted triangle this round (decoupled from TBB types)
    struct TriAcc {
        int u = -1, v = -1, w = -1;      // vertices (u,w,v) with anchor (u,v)
        int neg_eid = -1;                // full edge id of (u,v)
        int pos_eid_uw = -1;             // full edge id of (u,w)
        int pos_eid_wv = -1;             // full edge id of (w,v)
        double score_primary = 0.0;
        double score_secondary = 0.0;
        double viol = 0.0;               // not used in step 3 (kept for step 5+)
        double phi  = 0.0;               // avg salience over triangle edges
    };

    // Hooks called by TBB free functions (declared friend at end of class)
    void on_emit_(int full_eid, double used_density);
    void on_accept_(int full_eid, double density);
    int  override_budget_(int base) const;

    // Storage of triangle stage outputs for this round
    std::vector<TriAcc> tri_selected_;
    std::vector<int>    covered_neg_eids_;   // anchors that received at least one candidate

    // Allow TBB C hooks to call our private methods
    friend void ::TBB_on_emit(int, double);
    friend void ::TBB_on_accept(int, double);
    friend int  ::TBB_budget_override(int);
    // Core references
    SignedGraphForMIP& G_;
    SeparationPersistent&    S_;
    SeparationConfig         C_;

    // Per-round scratch
    SeparationRound          round_;

    // Full↔pos mapping under *current* switching (rebuilt each round)
    RoundGraphView Gt_;

    // === Step-2+ stubs (no-ops in step 1) ===
    void build_salience_(const std::vector<double>* x_hat,
                         const std::vector<double>* y_hat);
    void build_omega_prime_pos_();
    void pre_enumeration_reheat_();
    void triangle_stage_();
    void sp_stage_();
    void commit_stage_();
    void resize_round_pos_arrays_(int m_pos);
    void resize_round_full_arrays_(int m_full);    
};