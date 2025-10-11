// File: frustration_model_xy.cpp

#include "frustration_model_xy.h"
#include "frustration_model.h"
#include "cycle_key.h"
#include "reheat_pool.h"
#include "separation_pipeline.h"
#include "separation_config.h"

#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <cmath>
#include <optional>
#include <iomanip>
#include <stdexcept>

// ────────────────────────────────────────────────────────────────────────────
// Small utilities
// ────────────────────────────────────────────────────────────────────────────

static inline int sgn_from_sign_entry(const SignedEdge& se) {
    // SignedEdge::sign is expected to be +1 for positive, -1 for negative
    return (se.sign >= 0 ? +1 : -1);
}

static inline IloRange make_cycle_cut_from_meta(
    IloEnv& env,
    const std::vector<int>& y_idx,
    const std::vector<Edge>& edges,
    const std::vector<int>& eff_sign,
    int rhs,
    const std::vector<IloNumVar>& y,
    const std::vector<IloBoolVar>& x)
{
    IloExpr expr(env);
    for (size_t j = 0; j < edges.size(); ++j) {
        const Edge& e = edges[j];
        const int s   = eff_sign[j];
        expr += (y[y_idx[j]] - 0.5 * x[e.first] - 0.5 * x[e.second]) * s;
    }
    IloRange cut = (expr <= rhs);
    expr.end();
    return cut;
}

// Added: helper used by generate_cycle_cut_standard
static IloRange make_cycle_cut_from_prelookups(
    IloEnv& env,
    const std::vector<Edge>& all_edges,
    const std::vector<int>& idx,
    const std::vector<int>& sgn,
    const std::vector<int>* flip /*nullable*/,
    const std::vector<IloNumVar>& y,
    const std::vector<IloBoolVar>& x)
{
    IloExpr expr(env);
    int rhs = 0;
    for (size_t j = 0; j < all_edges.size(); ++j) {
        const Edge& e = all_edges[j];
        int eff = sgn[j] * (flip ? (*flip)[j] : 1);
        expr += (y[idx[j]] - 0.5 * x[e.first] - 0.5 * x[e.second]) * eff;
        if (eff < 0) ++rhs;
    }
    IloRange cut = (expr <= std::floor(rhs / 2.0));
    expr.end();
    return cut;
}

// ────────────────────────────────────────────────────────────────────────────
// FrustrationModelXY basics
// ────────────────────────────────────────────────────────────────────────────

FrustrationModelXY::FrustrationModelXY(SignedGraphForMIP& g, int cut_flags)
    : FrustrationModel(g, cut_flags) {}

// The original build() from your project remains as-is (unchanged here).
// If you have a custom build() in a different TU already, keep that one.
// Below is the same build() version you shared earlier (trimmed only for logs).
void FrustrationModelXY::build() {
    using clock = std::chrono::steady_clock;
    const auto TB_ALL0 = clock::now();
    auto TB_mark = TB_ALL0;

    int n = graph.vertex_count();
    int m = graph.edge_count();

    x.resize(n);
    y.resize(m);

    const auto& d_plus  = graph.plus_degrees_view();
    const auto& d_minus = graph.minus_degrees_view();

    for (int i = 0; i < n; ++i) x[i] = IloBoolVar(env, ("x_" + std::to_string(i)).c_str());
    for (int i = 0; i < m; ++i) y[i] = IloNumVar(env, 0.0, 1.0, ILOFLOAT, ("y_" + std::to_string(i)).c_str());

    // Objective
    IloExpr obj_expr(env);
    for (int i = 0; i < n; ++i) {
        int net_deg = d_plus[i] - d_minus[i];
        obj_expr += net_deg * x[i];
    }
    for (const auto& [edge, weight] : weights)
        obj_expr += -2 * weight * y[edge_index[edge]];
    objective = IloMinimize(env, obj_expr);
    model.add(objective);
    obj_expr.end();

    // Base inequalities
    IloRangeArray base(env);
    for (const auto& [edge, sign] : signs) {
        int u = edge.first;
        int v = edge.second;
        auto idxe = edge_index[edge];
        if (sign > 0) {
            base.add(y[idxe] - x[u] <= 0.0);
            base.add(y[idxe] - x[v] <= 0.0);
        } else {
            base.add(y[idxe] - x[u] - x[v] >= -1.0);
        }
    }
    if (base.getSize() > 0) model.add(base);

    // Greedy switching (sets current switched signature)
    graph.restore_switched_sign();
    std::vector<int> s = graph.greedy_switching();
    graph.switching_from_partition(s);

    // Optionally seed some triangle/cycle inequalities in build() stage if you want
    // (keeping your original behavior is fine). Nothing required for Round 5.

    // Heuristic feasible start from switching s
    IloNumArray start_vals(env);
    IloNumVarArray xy_vars(env);
    auto s01 = [&](int si) -> IloNum { return (si == -1) ? 1.0 : 0.0; };

    for (int i = 0; i < n; ++i) {
        model.add(x[i]);
        xy_vars.add(x[i]);
        start_vals.add(s01(s[i]));
    }
    for (const auto& [edge, idxe] : edge_index) {
        int u = edge.first, v = edge.second;
        model.add(y[idxe]); xy_vars.add(y[idxe]);

        int xu = static_cast<int>(s01(s[u]));
        int xv = static_cast<int>(s01(s[v]));
        double w = weights[edge];
        double y0 = (w > 0.0) ? static_cast<double>(std::min(xu, xv))
                              : static_cast<double>(std::max(0, xu + xv - 1));
        start_vals.add(y0);
    }

    // Some priorities on x to help the solver
    for (int i = 0; i < n; ++i) {
        int priority = d_plus[i] + d_minus[i];
        cplex.setPriority(x[i], priority);
    }

    cplex.addMIPStart(xy_vars, start_vals);
    injected_heuristic_solutions++;

    // Symmetry break
    igraph_integer_t max_vertex = graph.max_degree_vertex();
    model.add(x[max_vertex] == (s[max_vertex] == -1 ? 1 : 0));

    // Finalize
    std::cout << "[BUILD] n=" << n << " m=" << m << " done in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_ALL0).count()
              << " ms\n";
}

void FrustrationModelXY::configure_separation(const SeparationConfig& cfg) {
    sep_cfg_ = cfg;
    // Recreate driver with new config on next access
    driver_.reset();
}

TriangleCyclePipeline& FrustrationModelXY::driver() {
    ensure_driver_();
    return *driver_;
}

const TriangleCyclePipeline& FrustrationModelXY::driver() const {
    const_cast<FrustrationModelXY*>(this)->ensure_driver_();
    return *driver_;
}

void FrustrationModelXY::ensure_driver_() {
    if (!driver_) {
        // Keep persistent state sized to current graph
        sep_state_.init_sizes_if_needed(graph);
        driver_ = std::make_unique<TriangleCyclePipeline>(graph, sep_state_, sep_cfg_);
    }
}

void FrustrationModelXY::solve() {
    // Attach exactly ONE set of separation callbacks via the new delegators
    attach_separation_callbacks(cplex, SeparationMode::Pipeline);
    // Keep your heuristic switching callback as before
    cplex.use(new (env) SwitchingHeuristicCallback(env, *this));
    // Solve using base implementation
    FrustrationModel::solve();
}

void FrustrationModelXY::export_solution(const std::string& file_prefix, bool with_svg) const {
    std::ofstream xfile(file_prefix + "_x.csv");
    std::ofstream yfile(file_prefix + "_y.csv");

    std::vector<int> partition;
    for (std::size_t i = 0; i < x.size(); ++i) {
        double val = cplex.getValue(x[i]);
        xfile << val << "\n";
        partition.push_back(static_cast<int>(std::round(val)));
    }
    for (std::size_t i = 0; i < y.size(); ++i)
        yfile << cplex.getValue(y[i]) << "\n";

    xfile.close();
    yfile.close();

    FrustrationModel::export_solution(file_prefix, with_svg, partition);
}

// ────────────────────────────────────────────────────────────────────────────
// Separation utilities
// ────────────────────────────────────────────────────────────────────────────

void FrustrationModelXY::attach_separation_callbacks(IloCplex& cpx, SeparationMode mode) {
    switch (mode) {
        case SeparationMode::TrianglesOnly:
            cpx.use(new (env) NegativeTriangleCutGenerator(env, *this));
            break;
        case SeparationMode::CyclesOnly:
            cpx.use(new (env) NegativeCycleCutGenerator(env, *this));
            break;
        case SeparationMode::Pipeline:
        default:
            cpx.use(new (env) TriangleCycleCutGenerator(env, *this)); // full round
            break;
    }
}

void FrustrationModelXY::snapshot_lp_solution(std::vector<double>& xhat, std::vector<double>& yhat) const {
    xhat.assign(x.size(), 0.0);
    yhat.assign(y.size(), 0.0);

    // Try bulk pull via arrays (faster in callbacks and solves)
    try {
        IloNumVarArray X(env), Y(env);
        IloNumArray    Xv(env), Yv(env);
        for (const auto& xv : x) X.add(xv);
        for (const auto& yv : y) Y.add(yv);
        cplex.getValues(Xv, X);
        cplex.getValues(Yv, Y);
        for (IloInt i = 0; i < Xv.getSize() && i < (IloInt)xhat.size(); ++i) xhat[(size_t)i] = Xv[i];
        for (IloInt i = 0; i < Yv.getSize() && i < (IloInt)yhat.size(); ++i) yhat[(size_t)i] = Yv[i];
    } catch (...) {
        // Fallback
        for (size_t i = 0; i < x.size(); ++i) xhat[i] = cplex.getValue(x[i]);
        for (size_t i = 0; i < y.size(); ++i) yhat[i] = cplex.getValue(y[i]);
    }
}

TriangleCyclePipeline::Result
FrustrationModelXY::run_triangle_only_from_lp() {
    std::vector<double> xhat, yhat; snapshot_lp_solution(xhat, yhat);
    // Triangle-only mode: we call the driver in Build phase (no LP guidance required here).
    // Later stages (SP) in the driver are harmless if no-op; keys will reflect triangle stage.
    return driver().run_round(&xhat, &yhat, TriangleCyclePipeline::Phase::Build);
}

TriangleCyclePipeline::Result
FrustrationModelXY::run_cycle_only_from_lp() {
    std::vector<double> xhat, yhat; snapshot_lp_solution(xhat, yhat);
    // SP-only mode: for now we still call run_round; in later steps the driver
    // will respect mode/budgets. Keys returned will include SP cycles once implemented.
    return driver().run_round(&xhat, &yhat, TriangleCyclePipeline::Phase::Build);
}

TriangleCyclePipeline::Result
FrustrationModelXY::run_pipeline_round_from_lp() {
    std::vector<double> xhat, yhat; snapshot_lp_solution(xhat, yhat);
    const bool use_frac = fractional_phase_enabled_;
    return driver().run_round(&xhat, &yhat,
                              use_frac ? TriangleCyclePipeline::Phase::Fractional
                                       : TriangleCyclePipeline::Phase::Build);
}

void FrustrationModelXY::reset_separation_state() {
    sep_state_.in_model_keys.clear();
    sep_state_.recent_keys.clear();
    sep_state_.P_reheat.clear();

    sep_state_.init_sizes_if_needed(graph);
    std::fill(sep_state_.omega.begin(), sep_state_.omega.end(), 1.0);
    std::fill(sep_state_.pool_count.begin(), sep_state_.pool_count.end(), 0.0);
    std::fill(sep_state_.H.begin(), sep_state_.H.end(), 0.0);
}

// Materialize standard cycle cuts from switching-invariant fmkey::CycleKey
std::vector<std::pair<IloRange, std::string>>
FrustrationModelXY::build_cycle_cuts_from_keys(IloEnv& env,
                                               const std::vector<fmkey::CycleKey>& keys) const
{
    std::vector<std::pair<IloRange, std::string>> out;
    out.reserve(keys.size());

    for (const auto& k : keys) {
        // Neg edge (a,b) + all positive edges from key
        std::vector<Edge> all_edges;
        all_edges.reserve(1 + k.pos.size());

        // neg (a,b)
        all_edges.emplace_back(k.a, k.b);

        // positives, already canonicalized as (min,max)
        for (const auto& p : k.pos) {
            all_edges.emplace_back(p.first, p.second);
        }

        // Build the *standard* inequality (no vertex-reversed variants here)
        IloRange rng = generate_cycle_cut_standard(env, all_edges);
        out.emplace_back(std::move(rng),
                         (k.pos.size() == 2 ? "triangle" : "cycle"));
    }
    return out;
}

IloRange FrustrationModelXY::generate_cycle_cut_standard(IloEnv& env, const std::vector<Edge>& all_edges) const {
    // Pre-lookup indices & signs once
    std::vector<int> idx; idx.reserve(all_edges.size());
    std::vector<int> sgn; sgn.reserve(all_edges.size());
    for (const auto& e : all_edges) {
        idx.push_back(edge_index[e]);
        sgn.push_back(signs[e].sign);
    }
    // Build via helper (no flips)
    return make_cycle_cut_from_prelookups(env, all_edges, idx, sgn, nullptr, y, x);
}

// Keep your original triangle helper (unchanged). Required by base class.
std::vector<std::pair<IloRange, std::string>>
FrustrationModelXY::generate_pending_triangle_cuts(IloEnv& env, const TriangleInequalities& t) {
    std::vector<std::pair<IloRange, std::string>> cuts;
    std::unordered_set<int> seen;
    for (const auto& e : t.triangle) {
        seen.insert(e.first);
        seen.insert(e.second);
    }

    auto add_term = [&](const Edge& e, IloExpr& expr, int& rhs, int sign_flip = 1) {
        const SignedEdge& edge_info = signs[e];
        int index = edge_index[e];
        expr += ((y[index] - 0.5 * x[e.first] - 0.5 * x[e.second]) * edge_info.sign * sign_flip);
        if (edge_info.sign * sign_flip < 0) rhs++;
    };

    for (int id : t.pending_cut_ids) {
        if (id == -1) {
            IloExpr base_expr(env);
            int base_rhs = 0;
            for (const auto& e : t.triangle)
                add_term(e, base_expr, base_rhs);
            cuts.emplace_back(base_expr <= std::floor(base_rhs / 2.0), "standard");
            base_expr.end();
        } else {
            IloExpr flipped_expr(env);
            int flipped_rhs = 0;
            for (const auto& e : t.triangle) {
                int flip = (e.first == id || e.second == id) ? -1 : 1;
                add_term(e, flipped_expr, flipped_rhs, flip);
            }
            cuts.emplace_back(flipped_expr <= std::floor(flipped_rhs / 2.0), "reversed");
            flipped_expr.end();
        }
    }

    return cuts;
}

void FrustrationModelXY::print_separation_stats(std::ostream& os) const {
    os << "[SEP] omega.size=" << sep_state_.omega.size()
       << " H.size=" << sep_state_.H.size()
       << " pool_count.size=" << sep_state_.pool_count.size()
       << " recent_keys=" << sep_state_.recent_keys.size()
       << " in_model_keys=" << sep_state_.in_model_keys.size()
       << std::endl;
}

// ────────────────────────────────────────────────────────────────────────────
// Separation callback delegators (Round 5)
// ────────────────────────────────────────────────────────────────────────────

FrustrationModelXY::NegativeTriangleCutGenerator::NegativeTriangleCutGenerator(IloEnv env, FrustrationModelXY& owner)
    : IloCplex::UserCutCallbackI(env), owner(owner) {}

FrustrationModelXY::NegativeTriangleCutGenerator::~NegativeTriangleCutGenerator() {}

// Full round (triangle → SP → commit)
void FrustrationModelXY::TriangleCycleCutGenerator::main() {
    IloEnv ENV = getEnv();
    auto R = owner.run_pipeline_round_from_lp();
    if (!R.accepted_keys.empty()) {
        auto cuts = owner.build_cycle_cuts_from_keys(ENV, R.accepted_keys);
        for (auto& [rng, label] : cuts) add(rng).setName(label.c_str());
    }
}

FrustrationModelXY::NegativeCycleCutGenerator::NegativeCycleCutGenerator(IloEnv env, FrustrationModelXY& owner)
    : IloCplex::UserCutCallbackI(env), owner(owner) {}

FrustrationModelXY::NegativeCycleCutGenerator::~NegativeCycleCutGenerator() {}

// Triangle-only driver
void FrustrationModelXY::NegativeTriangleCutGenerator::main() {
    IloEnv ENV = getEnv();  // make an lvalue to bind to IloEnv&
    auto R = owner.run_triangle_only_from_lp();
    if (!R.accepted_keys.empty()) {
        auto cuts = owner.build_cycle_cuts_from_keys(ENV, R.accepted_keys);
        for (auto& [rng, label] : cuts) add(rng).setName(label.c_str());
    }
}

// SP-only driver
void FrustrationModelXY::NegativeCycleCutGenerator::main() {
    IloEnv ENV = getEnv();
    auto R = owner.run_cycle_only_from_lp();
    if (!R.accepted_keys.empty()) {
        auto cuts = owner.build_cycle_cuts_from_keys(ENV, R.accepted_keys);
        for (auto& [rng, label] : cuts) add(rng).setName(label.c_str());
    }
}

// ────────────────────────────────────────────────────────────────────────────
// Heuristic callback (unchanged from your existing implementation)
// ────────────────────────────────────────────────────────────────────────────

void FrustrationModelXY::SwitchingHeuristicCallback::main() {
    // Pull LP solution
    std::vector<double> x_vals(owner.x.size()), y_vals(owner.y.size());
    for (int i = 0; i < (int)owner.x.size(); ++i) x_vals[i] = getValue(owner.x[i]);
    for (int i = 0; i < (int)owner.y.size(); ++i) y_vals[i] = getValue(owner.y[i]);

    // Reweight for the UB pass (local only!)
    owner.graph.weighting_from_fractional(x_vals, y_vals);

    // Root-only polish knobs measured on the reweighted view
    int zero_nd = 0;
    {
        const auto& d_plus_pol  = owner.graph.plus_degrees_view();
        const auto& d_minus_pol = owner.graph.minus_degrees_view();
        for (int u = 0; u < (int)owner.x.size(); ++u)
            if ((d_plus_pol[u] - d_minus_pol[u]) == 0) ++zero_nd;
    }
    const bool   is_root    = (getNnodes64() == 0);
    const double zero_ratio = owner.x.empty() ? 0.0 : double(zero_nd) / double(owner.x.size());
    const int    max_kicks  = (is_root && zero_ratio >= 0.30) ? 3 : 1;

    const SignedGraph::GreedyKickOptions UB_OPTS{
        /*neg_edge_threshold_abs*/  -1,
        /*neg_edge_threshold_frac*/ 0.03,
        /*max_kicks*/               max_kicks,
        /*use_weighted_degree*/     true,
        /*use_triangle_tiebreak*/   true,
        /*triangle_beta*/           (is_root ? 0.08 : 0.05),
        /*neighbor_cap*/            1024,
        /*triangle_cap_per_u*/      1024
    };

    // Try fractional heuristic
    auto frac_result = owner.graph.fractional_greedy_switching(UB_OPTS);

    // Always restore global state before evaluating/printing/applying
    owner.graph.restore_switched_sign();

    // Normalize to shared_ptr<const vector<int>>
    std::shared_ptr<const std::vector<int>> s_ptr;
    if (frac_result.has_value()) {
        s_ptr = frac_result.value();  // already shared_ptr
        std::cout << "[Heuristic] switching solution (fractional) size=" << s_ptr->size() << "\n";
    } else {
        auto s_fb = owner.graph.greedy_switching();  // plain vector<int>
        if (s_fb.empty()) {
            std::cout << "[Heuristic] No switching solution found.\n";
            return;
        }
        s_ptr = std::make_shared<const std::vector<int>>(std::move(s_fb));
        std::cout << "[Heuristic] switching solution (fallback) size=" << s_ptr->size() << "\n";
    }

    const std::vector<int>& s = *s_ptr;

    // Evaluate in the 0/1 model space (views AFTER restore)
    const auto& d_plus_eval  = owner.graph.plus_degrees_view();
    const auto& d_minus_eval = owner.graph.minus_degrees_view();
    auto s01 = [](int si) -> double { return (si == -1) ? 1.0 : 0.0; };

    double obj_val = 0.0;
    IloNumVarArray xy_vars(getEnv());
    IloNumArray    xy_vals(getEnv());

    // x part
    for (int i = 0; i < (int)owner.x.size(); ++i) {
        double xi = s01(s[i]);
        xy_vars.add(owner.x[i]);
        xy_vals.add(xi);
        obj_val += (d_plus_eval[i] - d_minus_eval[i]) * xi;
    }

    // y part (tight, constraint-consistent initialization)
    for (const auto& [edge, idxe] : owner.edge_index) {
        int u = edge.first, v = edge.second;
        int xu = (int)s01(s[u]);
        int xv = (int)s01(s[v]);

        double w = owner.weights[edge];

        double y0 = (w > 0.0)
                  ? (double)std::min(xu, xv)              // y ≤ x_u, y ≤ x_v
                  : (double)std::max(0, xu + xv - 1);     // y ≥ x_u + x_v − 1

        xy_vars.add(owner.y[idxe]);
        xy_vals.add(y0);

        obj_val += -2.0 * w * y0;
    }

    double incumbent = hasIncumbent() ? getIncumbentObjValue() : IloInfinity;
    if (obj_val < incumbent) {
        setSolution(xy_vars, xy_vals);
        std::cout << "[Heuristic] Proposed solution with objective: " << obj_val << "\n";
        owner.f_index = obj_val / owner.graph.edge_count();
        owner.injected_heuristic_solutions++;
    }
}
