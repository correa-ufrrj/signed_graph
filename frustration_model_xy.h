// File: frustration_model_xy.h
#pragma once

#include "frustration_model.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <memory>
#include <utility>
#include "cycle_key.h"
#include "reheat_pool.h"
#include "separation_config.h"
#include "separation_pipeline.h"
// Explicit registration mode for separation callbacks
enum class SeparationMode {
    TrianglesOnly,  // TriangleBucketBatch only
    CyclesOnly,     // NegativeCycleBatch only
    Pipeline        // Triangle → SP round driver
};

// ----------------------------------------------------------------------------
// FrustrationModelXY — API ONLY (implementations deferred)
// ----------------------------------------------------------------------------
// Responsibilities relevant to separation steps 5–7:
//  - Owns the separation persistent state (ω, H, pool_count, de-dup sets, reheat pool)
//  - Hosts a TriangleCyclePipeline driver (round orchestrator)
//  - Exposes callback classes that delegate to the driver (triangle-only / SP-only)
//  - Provides utilities to translate accepted CycleKeys → CPLEX cuts (materialization)
//  - Leaves all method bodies to the .cpp (not included here to save canvas space)
// ----------------------------------------------------------------------------
class FrustrationModelXY : public FrustrationModel {
public:
    // Construct with the working graph and cut flags.
    explicit FrustrationModelXY(SignedGraphForMIP& g, int cut_flags = 0);

    // Standard model lifecycle
    void build() override;
    void solve() override;
    void export_solution(const std::string& file_prefix, bool with_svg) const;

    // ----------------------------- Separation API ----------------------------
    // Configure the separation knobs (copied into the driver at creation time).
    void configure_separation(const SeparationConfig& cfg);

    // Access the current configuration/persistent state used by the driver.
    const SeparationConfig&      separation_config()   const { return sep_cfg_; }
    SeparationConfig&            separation_config()         { return sep_cfg_; }
    const SeparationPersistent&  separation_state()    const { return sep_state_; }
    SeparationPersistent&        separation_state()          { return sep_state_; }

    // Driver access (created lazily on first use).
    TriangleCyclePipeline&       driver();
    const TriangleCyclePipeline& driver() const;

    // Register exactly one separation callback, per the selected mode.
//  - TrianglesOnly : NegativeTriangleCutGenerator (TBB-only)
//  - CyclesOnly    : NegativeCycleCutGenerator (SP-only)
//  - Pipeline      : TriangleCycleCutGenerator (triangle→SP→commit)
void attach_separation_callbacks(IloCplex& cplex,
                                 SeparationMode mode = SeparationMode::Pipeline);

    // Phase policy for callbacks (default: Fractional phase enabled)
    void set_fractional_phase_enabled(bool on) { fractional_phase_enabled_ = on; }
    bool fractional_phase_enabled() const { return fractional_phase_enabled_; }

    // Convenience: run driver using current LP solution (x̂,ŷ) from CPLEX
    TriangleCyclePipeline::Result run_triangle_only_from_lp();
    TriangleCyclePipeline::Result run_cycle_only_from_lp();
    TriangleCyclePipeline::Result run_pipeline_round_from_lp();

    // Helper to extract current LP solution to plain vectors
    void snapshot_lp_solution(std::vector<double>& xhat, std::vector<double>& yhat) const;

    // Clear & reinit the persistent separation state (ω←1, H←0, pool_count←0, de-dup cleared).
    void reset_separation_state();

    // Materialize a set of cycle keys into CPLEX cuts (standard inequality family).
    // The returned vector carries (IloRange, debug_label) pairs.
    std::vector<std::pair<IloRange, std::string>>
    build_cycle_cuts_from_keys(IloEnv& env,
                               const std::vector<fmkey::CycleKey>& keys) const;

    // Convenience: build triangle cuts that were buffered during the previous round
    // (used by the conditional cycle generator in some strategies).
    std::vector<std::pair<IloRange, std::string>>
    generate_pending_triangle_cuts(IloEnv& env, const TriangleInequalities& t) override;

    // Lightweight stats/telemetry for debugging the separation driver.
    void print_separation_stats(std::ostream& os) const;

private:
    // ========================== Variables & bookkeeping ======================
    std::vector<IloBoolVar> x;   // vertex binaries
    std::vector<IloNumVar>  y;   // edge binaries (undirected)

    // Whether the last round added triangle cuts (used by ConditionalCycleCutGenerator)
    bool 
    // Global policy for callbacks: Fractional vs Build phase
    bool fractional_phase_enabled_ = true;

    // (Moved) Reheat + de-dup state now lives in SeparationPersistent (sep_state_):
//   - sep_state_.P_reheat
//   - sep_state_.in_model_keys
//   - sep_state_.recent_keys
// Configuration for reheat/dedup can be hosted in SeparationConfig if needed.

    // ===== Step 6 persistent updates (ω, H, pool_count) live here =====
    SeparationConfig   sep_cfg_{};   // knobs (LP guidance/annealing also live here)
    SeparationPersistent sep_state_{}; // ω, H, pool_count, de-dup sets and reheat pool (mirrored)

    // Triangle/Cycle pipeline driver (created on demand)
    std::unique_ptr<TriangleCyclePipeline> driver_;

    // Ensure the driver exists and is synchronized with current config/state
    void ensure_driver_();

    // Seed buckets from reheated items that are currently violated.
    // (Implementation in .cpp — uses current LP solution x̂/ŷ.)
    void reheat_seed_buckets_(
        /*in/out*/ std::unordered_map<int, std::vector<int>>& by_neg_eid,
        const std::vector<char>& neg_edge_mask,
        const IloNumArray& xhat,
        const IloNumArray& yhat);

    // After selection, decay TTLs and prune the pool; update in_model/recent sets as needed.
    void reheat_after_selection_(
        const std::vector<int>& chosen_cids,
        const std::vector<std::vector<int>>& cand_pos_edges,
        const IloNumArray& xhat,
        const IloNumArray& yhat);

    // ============================ Cut builders ===============================
	IloRange generate_cycle_cut_standard(IloEnv& env, const std::vector<Edge>& all_edges);
	std::vector<std::pair<IloRange, std::string>> generate_cycle_cuts(IloEnv& env, const std::vector<Edge>& all_edges) override;
	std::vector<std::pair<IloRange, std::string>> generate_positive_triangle_cuts(IloEnv& env, const std::vector<Edge>& all_edges);

public:
    // ============================== Callbacks ================================
    class NegativeTriangleCutGenerator : public IloCplex::UserCutCallbackI {
    public:
        FrustrationModelXY& owner;

        NegativeTriangleCutGenerator(IloEnv env,
                                     FrustrationModelXY& owner);
        ~NegativeTriangleCutGenerator();

        // Delegates to the pipeline (triangle stage only);
        // collects accepted keys and adds materialized cuts to the model.
        void main() override;

        IloCplex::CallbackI* duplicateCallback() const override {
            return new (getEnv()) NegativeTriangleCutGenerator(getEnv(), owner);
        }
    };

    class NegativeCycleCutGenerator : public IloCplex::UserCutCallbackI {
    public:
        FrustrationModelXY& owner;

        NegativeCycleCutGenerator(IloEnv env,
                                  FrustrationModelXY& owner);
		~NegativeCycleCutGenerator();

        // Delegates to the pipeline (SP stage on uncovered anchors);
        // collects accepted keys and adds materialized cuts to the model.
        void main() override;

        IloCplex::CallbackI* duplicateCallback() const override {
            return (new (getEnv()) NegativeCycleCutGenerator(getEnv(), owner));
        }
    };

    // New: full pipeline driver callback (triangle→SP→commit)
    class TriangleCycleCutGenerator : public IloCplex::UserCutCallbackI {
    public:
        FrustrationModelXY& owner;

        TriangleCycleCutGenerator(IloEnv env, FrustrationModelXY& owner)
            : IloCplex::UserCutCallbackI(env), owner(owner) {}

        void main() override;

        IloCplex::CallbackI* duplicateCallback() const override {
            return new (getEnv()) TriangleCycleCutGenerator(getEnv(), owner);
        }
    };

    class SwitchingHeuristicCallback : public IloCplex::HeuristicCallbackI {
    private:
        FrustrationModelXY& owner;

    public:
        SwitchingHeuristicCallback(IloEnv env,
                                   FrustrationModelXY& owner)
            : IloCplex::HeuristicCallbackI(env), owner(owner) {}

        void main() override;

        IloCplex::CallbackI* duplicateCallback() const override {
            return (new (getEnv()) SwitchingHeuristicCallback(getEnv(), owner));
        }
    };
};