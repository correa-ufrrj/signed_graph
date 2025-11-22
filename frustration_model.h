// File: frustration_model.h
#pragma once

#include "signed_graph_mip.h"
#include <ilcplex/ilocplex.h>
#include <string>
#include <iostream>
#include <chrono>
#include <numeric>
#include <vector>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <array>
#include <limits>
#include <cmath>
#include <iomanip>

struct TriangleInequalities {
    std::vector<Edge> triangle;
    std::vector<int> pending_cut_ids; // -1 for full triangle cut, vertex id for vertex-based cuts
};

// 1) Accessor interface
struct LpAccessor {
    virtual ~LpAccessor() = default;
    virtual void getValues(IloNumArray& out, const IloNumVarArray& vars) const = 0;
};

// 2) Model-side accessor (safe outside callbacks)
struct ModelAccessor final : LpAccessor {
    explicit ModelAccessor(IloCplex& cpx);
    void getValues(IloNumArray& out, const IloNumVarArray& vars) const override;
private:
    IloCplex& cpx_;
};

// 3) Callback-side accessor (safe inside legacy callbacks)
struct UserCutCallbackAccessor final : LpAccessor {
    explicit UserCutCallbackAccessor(IloCplex::UserCutCallbackI& cb);
    void getValues(IloNumArray& out, const IloNumVarArray& vars) const override;
private:
    IloCplex::UserCutCallbackI& cb_;
};

struct HeuristicCallbackAccessor final : LpAccessor {
    explicit HeuristicCallbackAccessor(IloCplex::HeuristicCallbackI& cb);
    void getValues(IloNumArray& out, const IloNumVarArray& vars) const override;
private:
    IloCplex::HeuristicCallbackI& cb_;
};

class FrustrationModel {
public:
    static constexpr int NO_CUTS = 0;
    static constexpr int NEGATIVE_CYCLE_CUTS = 1;
    static constexpr int NET_DEGREE_CUTS = 2;
    static constexpr int TRIANGLE_CUTS = 4;

    // ───────────────────────── Round Control Policy (base) ─────────────────
    struct RoundStats {
        int    rows_added      = 0;
        int    cand_seen       = 0;
        double time_s          = 0.0;
        double best_integer    = std::numeric_limits<double>::infinity();
        double best_bound      = -std::numeric_limits<double>::infinity();
        double prev_best_bound = -std::numeric_limits<double>::infinity();
        int tri_rows = 0;
        int cyc_rows = 0;
        int std_rows = 0;
        int rev_rows = 0;
        double accept_rate_hint = 0.0;
    };

    struct RoundPolicyConfig {
        int    max_rounds_at_root     = 200;
        int    min_rows_per_round     = 32;
        double min_accept_rate        = 0.03;
        double min_gap_gain_rel       = 1e-4;
        double min_gap_gain_speed     = 5e-5; // rel-gap per second
        int    weak_patience          = 5;
        double min_tri_cycle_ratio    = 0.25;
        double min_std_rev_ratio      = 0.30;
        double gap_close_rel          = 3e-3;

        // Where to apply the separation pipeline:
        //  - RootOnly : only at the root node of the B&B tree
        //  - AllNodes : at every node
        enum class Scope { RootOnly, AllNodes };
        Scope pipeline_scope = Scope::AllNodes;
	};

    enum class RoundDecision { Continue, StopCutting };

    class RoundPolicy {
    public:
        RoundPolicy() : cfg_() {}
        explicit RoundPolicy(const RoundPolicyConfig& cfg) : cfg_(cfg) {}
        void reset_node() {
            rounds_ = 0; weak_streak_ = 0;
            last_gap_rel_ = std::numeric_limits<double>::infinity();
        }
        RoundDecision decide(const double gap_rel, const RoundStats& rs) {
            ++rounds_;
            const double gap_impr = std::max(0.0, last_gap_rel_ - gap_rel);
            const double speed    = (rs.time_s > 1e-12) ? (gap_impr / rs.time_s) : 0.0;
            const double denom    = std::max(1, rs.tri_rows + rs.cyc_rows);
            const double tri_ratio= static_cast<double>(rs.tri_rows) / denom;
            const double std_ratio= static_cast<double>(rs.std_rows) / std::max(1, rs.std_rows + rs.rev_rows);
            const double acc_rate = (rs.accept_rate_hint > 0.0)
                ? rs.accept_rate_hint
                : (rs.cand_seen > 0 ? (static_cast<double>(rs.rows_added)/rs.cand_seen) : 0.0);

            const bool good_gap    = (gap_impr >= cfg_.min_gap_gain_rel);
            const bool good_speed  = (speed   >= cfg_.min_gap_gain_speed);
            const bool enough_rows = (rs.rows_added >= cfg_.min_rows_per_round);
            const bool decent_acc  = (acc_rate >= cfg_.min_accept_rate);
            const bool comp_ok     = (tri_ratio >= cfg_.min_tri_cycle_ratio) && (std_ratio >= cfg_.min_std_rev_ratio);
            const bool nearly_done = (gap_rel   <= cfg_.gap_close_rel);

            const bool weak = !(good_gap && good_speed && enough_rows && decent_acc);
            weak_streak_ = weak ? (weak_streak_ + 1) : 0;

            std::cout.setf(std::ios::fixed); std::cout << std::setprecision(6);
            std::cout << "[ROUND-POLICY] k=" << rounds_
                      << " rows=" << rs.rows_added
                      << " acc="  << acc_rate
                      << " gap_rel=" << gap_rel
                      << " dgap=" << gap_impr
                      << " speed=" << speed
                      << " tri%=" << tri_ratio
                      << " std%=" << std_ratio
                      << " weak=" << weak
                      << " weak_streak=" << weak_streak_
                      << (nearly_done ? " near_close=1" : "")
                      << (comp_ok ? " comp_ok=1" : " comp_ok=0")
                      << "\n";

            if (rounds_ >= cfg_.max_rounds_at_root) return RoundDecision::StopCutting;
            if (nearly_done)                         return RoundDecision::StopCutting;
            const int patience = comp_ok ? cfg_.weak_patience : std::max(1, cfg_.weak_patience - 2);
            if (weak_streak_ >= patience)           return RoundDecision::StopCutting;

            last_gap_rel_ = gap_rel;
            return RoundDecision::Continue;
        }
        const RoundPolicyConfig& config() const { return cfg_; }
        RoundPolicyConfig&       config()       { return cfg_; }
        int rounds() const { return rounds_; }
        int weak_streak() const { return weak_streak_; }
    private:
        RoundPolicyConfig cfg_{};
        int rounds_ = 0;
        int weak_streak_ = 0;
        double last_gap_rel_ = std::numeric_limits<double>::infinity();
    };

    // Base accessors used by all derived models/callbacks
    RoundPolicy&       round_policy()       { return policy_; }
    const RoundPolicy& round_policy() const { return policy_; }
    bool cuts_enabled() const { return !stop_cutting_root_; }
    void disable_cuts() { stop_cutting_root_ = true; }
    void enable_cuts()  { stop_cutting_root_ = false; }

    FrustrationModel(SignedGraphForMIP& g, int cut_flags = 0);
    virtual ~FrustrationModel() {
        env.end();
    };

    virtual void build() = 0;
    virtual void solve();
    virtual double get_frustration_index() const;
    void print_solution() const;
    std::string active_cut_names() const;
    virtual void export_solution(const std::string& file_prefix, bool with_svg) const = 0;

protected:
    // Round-policy state (shared by all derived models)
    RoundPolicy policy_{};
    bool stop_cutting_root_ = false;
    // When StopCutting happens in AllNodes scope, remember the node id so we
    // can keep the rest of THIS node silent and only recover on the next one.
    long long last_disabled_node_id_ = std::numeric_limits<long long>::min();

    SignedGraphForMIP& graph;
    int m_minus;
    std::vector<TriangleInequalities> uncut_triangles; // pending triangle inequalities
    IloEnv env;
    IloModel model;
    IloCplex cplex;
    IloObjective objective;
    int use_cut_generator;
    long int lower_bound;
    double f_index;
    int net_degree_cut_count = 0;
    int standard_cycle_cuts_build = 0;
    int reversed_cycle_cuts_build = 0;
    int standard_cycle_cuts_cutgen = 0;
    int reversed_cycle_cuts_cutgen = 0;
    int injected_heuristic_solutions = 0;
    std::chrono::steady_clock::time_point start_time;
    std::chrono::steady_clock::time_point end_time;
    std::shared_ptr<const std::vector<int>> switching_solution = nullptr;
    SignedGraph::SignedEdgesView signs;
    std::vector<int> signs0;   // size = m, entries in {+1,-1}; frozen at construction
    EdgeIndexesView edge_index;
    std::unordered_map<int, Edge> edge_reverse;

    double frustration_index(double obj_val) const;

    void initialize_uncut_triangles();
    virtual std::vector<std::pair<IloRange, std::string>> generate_cycle_cuts(IloEnv& env, const std::vector<Edge>& all_edges) const = 0;
    std::vector<std::pair<IloRange, std::string>> generate_triangle_cuts(IloEnv& env);
    virtual void export_solution(const std::string& file_prefix, bool with_svg, std::vector<int> partition) const;

    // Helpers to interpret flags consistently across models/callbacks
    inline bool pipeline_enabled() const {
        return ((use_cut_generator & TRIANGLE_CUTS) != 0) &&
               ((use_cut_generator & NEGATIVE_CYCLE_CUTS) != 0);
    }
    inline bool legacy_negcycles_enabled() const {
        return (use_cut_generator & NEGATIVE_CYCLE_CUTS) != 0;
    }
    inline bool legacy_triangles_enabled() const {
        return (use_cut_generator & TRIANGLE_CUTS) != 0;
    }

    // Helpers for scope-aware gating
    inline void mark_disabled_at_node(long long node_id) {
        stop_cutting_root_   = true;
        last_disabled_node_id_ = node_id;
    }
    inline long long last_disabled_node_id() const noexcept {
        return last_disabled_node_id_;
    }
    
    // Assemble RoundStats from a CPLEX callback, evaluate the policy, and
	// optionally disable further cuts at the root if the policy says to stop.
	// - cb: any IloCplex::CallbackI (UserCut, Lazy, Heuristic, etc.)
	// - rows_added: total rows added in this round
	// - cand_seen: candidates inspected (pass rows_added if unknown)
	// - time_s: wall time of the round (seconds)
	// - tri_rows, cyc_rows, std_rows, rev_rows: composition tallies
	// - accept_rate_hint: optional precomputed acceptance rate; pass 0.0 if unknown
	// - auto_disable: if true, will call disable_cuts() when StopCutting
	// Pass best_integer/best_bound/rel_gap explicitly (taken from the concrete callback).
	RoundDecision apply_round_policy(int rows_added,
	                                int cand_seen,
	                                double time_s,
	                                int tri_rows,
	                                int cyc_rows,
	                                int std_rows,
	                                int rev_rows,
	                                double best_integer,
	                                double best_bound,
	                                double rel_gap,
	                                double accept_rate_hint = 0.0,
	                                bool auto_disable = true);
};
