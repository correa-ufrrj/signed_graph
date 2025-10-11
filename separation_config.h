// separation_config.h
#pragma once

#include "triangle_bucket_batch.h" // for TriangleBucketBatch::Params adapter
#include "negative_cycle_batch.h"  // for NegativeCycleBatch::Params adapter

// Unified separation configuration with grouped sub-structs.
// Keeps a single source of truth while preserving readability.
struct SeparationConfig {
    // ---------------------- Ranking & LP guidance ----------------------
    struct Ranking {
        double alpha       = 0.3;   // triangle primary score weight (1/ω′)
        double theta       = 0.5;   // triangle salience blend
        double lambda_hist = 0.20;  // historical repulsion blend (H)
        double lambda_LP   = 0.25;  // LP salience blend (step 7+)
    } ranking;

    // ------------------------- Weight dynamics -------------------------
    struct Weights {
        double beta_emit = 0.02;    // within-batch ω′ bump per emit (used_density-scaled)
        double beta_sel  = 0.05;    // cross-batch ω  drift per accepted edge
        double omega_eps = 1e-8;    // lower clamp for ω, ω′
        double omega_max = 64.0;    // upper clamp for ω
    } weights;

    // -------------------------- Budgets / caps -------------------------
    struct Budgets {
        // Triangles (TBB)
        int B_tri            = 64;  // global per-batch budget
        int K_tri_per_neg    = 8;   // per-bucket prefilter
        int tri_cap_per_vertex = 6; // per-vertex cap during selection
        // Shortest-path cycles (NCB)
        int B_sp             = 64;  // cycles budget (if used by driver)
        int K_sp_per_neg     = 3;   // alt paths per negative anchor
        int sp_cap_per_vertex= 8;   // per-vertex cap for SP acceptance
    } budget;

    // ----------------------- Annealing for B_tri -----------------------
    struct AnnealTri {
        int    B_min   = 1;
        double gamma_min = 0.70;
        double gamma_max = 0.92;
        double v0        = 0.10;  // target violation
        double tau       = 0.10;  // smoothing
    } anneal_tri;

    // ---------------------- EMA (H) & recent TTL -----------------------
    struct EMA {
        double delta = 0.85;   // decay (0<delta<1)
        double rho   = 0.02;   // accepted-weight
        double kappa = 0.005;  // emitted-but-rejected weight
        int    recent_ttl = 3; // stateless dedup TTL rounds
    } ema;

    // ===================== Thin adapters for engines ====================
    // TriangleBucketBatch adapter
    TriangleBucketBatch::Params to_tbb_params() const {
        TriangleBucketBatch::Params p;
        p.K_tri_per_neg  = budget.K_tri_per_neg;
        p.B_tri          = budget.B_tri;
        p.cap_per_vertex = budget.tri_cap_per_vertex;
        return p;
    }

    // Adapter for NegativeCycleBatch parameter pack (decoupled engine)
NegativeCycleBatch::Params to_ncb_params() const {
    NegativeCycleBatch::Params p;
    p.B_sp               = budget.B_sp;
    p.K_sp_per_neg       = budget.K_sp_per_neg;
    p.sp_cap_per_vertex  = budget.sp_cap_per_vertex;
    p.tri_cap_per_vertex = budget.tri_cap_per_vertex;
    p.alt_path_bump_scale = 1.0;              // default: keep current behavior
    p.cross_batch_penalty_scale = 3.0;        // default drift
    p.tbb = to_tbb_params();                  // pass through triangle knobs for NCB's internal TBB
    return p;
}
};

// Optional: presets
inline SeparationConfig build_profile() {
    SeparationConfig c; 
    c.ranking.lambda_LP = 0.0; 
    c.budget.B_tri = 64; c.budget.B_sp = 64; 
    return c;
}
inline SeparationConfig fractional_profile() {
    SeparationConfig c; 
    c.ranking.lambda_LP = 0.25; 
    return c;
}
