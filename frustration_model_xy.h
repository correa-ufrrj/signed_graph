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
#include "separation_config.h"
#include "separation_pipeline.h"
#include <chrono>
#include <limits>
#include <cmath>

// ----------------------------------------------------------------------------
// FrustrationModelXY — API
// ----------------------------------------------------------------------------
class FrustrationModelXY : public FrustrationModel {
public:
    explicit FrustrationModelXY(SignedGraphForMIP& g, int cut_flags = 0);

    // Standard model lifecycle
    void build() override;
    void solve() override;
    void export_solution(const std::string& file_prefix, bool with_svg) const;

    // ----------------------------- Separation API ----------------------------
    void configure_separation(const SeparationConfig& cfg);

    const SeparationConfig&      separation_config()   const { return sep_cfg_; }
    SeparationConfig&            separation_config()         { return sep_cfg_; }
    const SeparationPersistent&  separation_state()    const { return sep_state_; }
    SeparationPersistent&        separation_state()          { return sep_state_; }

    // Driver access (created lazily on first use).
    TriangleCyclePipeline&       driver();
    const TriangleCyclePipeline& driver() const;

    // Phase policy for callback (default: Fractional phase enabled)
    void set_fractional_phase_enabled(bool on) { fractional_phase_enabled_ = on; }
    bool fractional_phase_enabled() const { return fractional_phase_enabled_; }

    // Helper to extract current LP solution to plain vectors
    void snapshot_lp_solution(const LpAccessor& acc) const;

    // Clear & reinit persistent separation state
    void reset_separation_state();

    // Materialize CycleKeys → CPLEX cuts (standard inequality family).
    std::vector<std::pair<IloRange, std::string>>
    build_cycle_cuts_from_keys(IloEnv& env,
                               const std::vector<fmkey::CycleKey>& keys) const;

    void print_separation_stats(std::ostream& os) const;

private:
    // ========================== Variables & bookkeeping ======================
    IloNumVarArray x;   // vertex binaries
    IloNumVarArray y;   // edge binaries (undirected)

    mutable IloNumArray        xval;
    mutable std::vector<double> xhat;

    int    sym_vstar = -1;
    double sym_xfix  = 0.0;

    // Global policy for callback: Fractional vs Build phase
    bool fractional_phase_enabled_ = true;

    // ===== persistent separation state (ω, H, pool_count) =====
    SeparationConfig      sep_cfg_{};
    SeparationPersistent  sep_state_{};

    // Triangle/Cycle pipeline driver (created on demand)
    std::unique_ptr<TriangleCyclePipeline> driver_;

    // Ensure the driver exists and is synchronized with current config/state
    void ensure_driver_();

    // ============================ Cut builders ===============================
    std::vector<std::pair<IloRange, std::string>>
    generate_cycle_cuts(IloEnv& env, const std::vector<Edge>& all_edges) const override;

public:
    // ============================== Callbacks ================================
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
