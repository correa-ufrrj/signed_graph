// separation_config.h
#pragma once

#include "triangle_bucket_batch.h" // for TriangleBucketBatch::Params adapter
#include "negative_cycle_batch.h"  // for NegativeCycleBatch::Params adapter

// Unified separation configuration with grouped sub-structs.
// Keeps a single source of truth while preserving readability.
struct SeparationConfig {
    // ---------------------- Ranking & LP guidance ----------------------
    struct Ranking {
        double alpha       = 1.2;   // triangle primary score weight (1/ω′)
        double theta       = 0.5;  // triangle salience blend
        double lambda_hist = 0.4;  // historical repulsion blend (H)
        double lambda_LP   = 0.45;  // LP salience blend (step 7+)
        double viol_tol    = 1e-6;
        double inv_power   = 1.0;
    } ranking;

    // ------------------------- Weight dynamics -------------------------
    struct Weights {
        double beta_emit = 10;    // within-batch ω′ bump per emit (used_density-scaled)
        double beta_sel  = 0.3;    // cross-batch ω  drift per accepted edge
        double omega_eps = 1e-6;    // lower clamp for ω, ω′
        double omega_max = 32.0;    // upper clamp for ω
    } weights;

    // -------------------------- Budgets / caps -------------------------
    struct Budgets {
        // Triangles (TBB)
        int B_tri            = 1024;  // global per-batch budget
        int K_tri_per_neg    = 1;   // per-bucket prefilter
        int tri_cap_per_vertex = 1; // per-vertex cap during selection
        // Shortest-path cycles (NCB)
        int B_sp             = 1024;  // cycles budget (if used by driver)
        int K_sp_per_neg     = 1;   // alt paths per negative anchor
        int sp_cap_per_vertex= 6;   // per-vertex cap for SP acceptance
    } budget;

    // ----------------------- Annealing for B_tri -----------------------
    struct AnnealTri {
	    int    B_min     = 1;     // never drop below 1 via the multiplicative rule
	    double gamma_min = 0.90;  // faster decay when v is high (but still gentle)
	    double gamma_max = 0.985; // very gentle decay when v is low
	    double v0        = 0.25;  // target proxy-violation
	    double tau       = 0.15;  // slower EMA
    } anneal_tri;

    // ---------------------- EMA (H) & recent TTL -----------------------
    struct EMA {
        double delta = 0.992;   // decay (0<delta<1)
        double rho   = 0.3;     // accepted-weight
        double kappa = 0.15;    // emitted-but-rejected weight
        int    recent_ttl = 3;  // stateless dedup TTL rounds
    } ema;

    // -------------------------- Modes / gating -------------------------
    struct Modes {
        bool use_triangles  = true;  // driver should set from FrustrationModel::TRIANGLE_CUTS
        bool use_negcycles  = true;  // driver should set from FrustrationModel::NEGATIVE_CYCLE_CUTS
    } modes;

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
	    p.cross_batch_penalty_scale = 3.0;        // default drift
	    p.tbb = to_tbb_params();                  // pass through triangle knobs for NCB's internal TBB

	    // unify with triangle contract
	    p.beta_emit  = weights.beta_emit;
	    p.omega_eps  = weights.omega_eps;
	    p.omega_max  = weights.omega_max;
	    return p;
	}
};

// Optional: presets
inline SeparationConfig build_profile() {
    SeparationConfig c; 
    c.ranking.lambda_LP = 0.0; 
    c.budget.B_tri = 64;
    c.budget.B_sp = 64; 
    return c;
}
inline SeparationConfig fractional_profile() {
    SeparationConfig c; 
    c.ranking.lambda_LP = 0.25; 
    return c;
}
