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
        double theta       = 0.05;   // triangle salience blend
        double lambda_hist = 0.10;  // historical repulsion blend (H)
        double lambda_LP   = 1.0;  // LP salience blend (step 7+)
    } ranking;

    // ------------------------- Weight dynamics -------------------------
    struct Weights {
        double beta_emit = 2.8;    // within-batch ω′ bump per emit (used_density-scaled)
        double beta_sel  = 0.05;    // cross-batch ω  drift per accepted edge
        double omega_eps = 1e-6;    // lower clamp for ω, ω′
        double omega_max = 32.0;    // upper clamp for ω
    } weights;

    // -------------------------- Budgets / caps -------------------------
    struct Budgets {
        // Triangles (TBB)
        int B_tri            = 512;  // global per-batch budget
        int K_tri_per_neg    = 4;   // per-bucket prefilter
        int tri_cap_per_vertex = 8; // per-vertex cap during selection
        // Shortest-path cycles (NCB)
        int B_sp             = 512;  // cycles budget (if used by driver)
        int K_sp_per_neg     = 3;   // alt paths per negative anchor
        int sp_cap_per_vertex= 8;   // per-vertex cap for SP acceptance
    } budget;

    // ----------------------- Annealing for B_tri -----------------------
    struct AnnealTri {
	    int    B_min     = 1;     // never drop below 1 via the multiplicative rule
	    double gamma_min = 0.94;  // faster decay when v is high (but still gentle)
	    double gamma_max = 0.995; // very gentle decay when v is low
	    double v0        = 0.30;  // target proxy-violation
	    double tau       = 0.10;  // slower EMA
    } anneal_tri;

    // ---------------------- EMA (H) & recent TTL -----------------------
    struct EMA {
        double delta = 0.85;   // decay (0<delta<1)
        double rho   = 0.02;   // accepted-weight
        double kappa = 0.005;  // emitted-but-rejected weight
        int    recent_ttl = 3; // stateless dedup TTL rounds
    } ema;
    
    // ---------------------- Reheat pool -----------------------
	struct ReheatPool {
	  int   reheat_stage_cap = 512;      // how many non-violated cycles we stage per round
	  size_t reheat_pool_cap = 10000;    // global cap to keep memory bounded
	  int   reheat_default_ttl = 3;      // TTL given to newly staged items
	  double reheat_ema_alpha = 0.5;     // EMA smoothing for violation
	  double reheat_ema_min_keep = 1e-4; // purge extremely cold items when trimming
	} reheat;


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
    c.budget.B_tri = 64;
    c.budget.B_sp = 64; 
    return c;
}
inline SeparationConfig fractional_profile() {
    SeparationConfig c; 
    c.ranking.lambda_LP = 0.25; 
    return c;
}
