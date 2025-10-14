// File: frustration_model_xy.cpp

#include "frustration_model_xy.h"
#include "frustration_model.h"
#include "cycle_key.h"
#include "separation_pipeline.h"
#include "separation_pipeline_tls.h"  // g_reheat_pool
#include "separation_config.h"

#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <cmath>
#include <optional>
#include <iomanip>
#include <stdexcept>
#include <fstream>

// ────────────────────────────────────────────────────────────────────────────
// Small utilities
// ────────────────────────────────────────────────────────────────────────────

// Helper used by generate_cycle_cut_standard
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

void FrustrationModelXY::build() {
    using clock = std::chrono::steady_clock;
    const auto T0 = clock::now();

    const int n = graph.vertex_count();
    const int m = graph.edge_count();

    x.resize(n);
    y.resize(m);

    const auto& d_plus  = graph.plus_degrees_view();
    const auto& d_minus = graph.minus_degrees_view();

    // ================== Vars ==================
    for (int i = 0; i < n; ++i)
        x[i] = IloBoolVar(env, ("x_" + std::to_string(i)).c_str());
    for (int i = 0; i < m; ++i)
        y[i] = IloNumVar(env, 0.0, 1.0, ILOFLOAT, ("y_" + std::to_string(i)).c_str());

    for (int i = 0; i < n; ++i) model.add(x[i]);
    for (int i = 0; i < m; ++i) model.add(y[i]);

    // ================== Objective ==================
    {
        IloExpr obj(env);
        for (int i = 0; i < n; ++i) obj += (d_plus[i] - d_minus[i]) * x[i];
        for (const auto& [e, w] : weights) obj += -2.0 * w * y[edge_index[e]];
        objective = IloMinimize(env, obj);
        model.add(objective);
        obj.end();
    }

    // ============ Base formulation (edge constraints) ============
    {
        IloRangeArray base(env);
        for (const auto& [e, sgn] : signs) {
            const int u = e.first, v = e.second;
            const int idx = edge_index[e];
            if (sgn > 0) {
                base.add(y[idx] - x[u] <= 0.0);
                base.add(y[idx] - x[v] <= 0.0);
            } else {
                base.add(y[idx] - x[u] - x[v] >= -1.0);
            }
        }
        if (base.getSize() > 0) model.add(base);
    }

    // ============ Greedy switching ============
    graph.restore_switched_sign();
    std::vector<int> s = graph.greedy_switching();      // {-1,+1}
    graph.switching_from_partition(s);                  // apply σ_s

    // =============== MIP start ===============
    {
        IloNumVarArray vars(env);
        IloNumArray vals(env);
        auto s01 = [&](int si){ return (si == -1) ? 1.0 : 0.0; };

        for (int i = 0; i < n; ++i) { vars.add(x[i]); vals.add(s01(s[i])); }
        for (const auto& [e, idx] : edge_index) {
            const double xu = s01(s[e.first]), xv = s01(s[e.second]);
            vars.add(y[idx]); vals.add(xu * xv);
        }

        for (int i = 0; i < n; ++i) cplex.setPriority(x[i], d_plus[i] + d_minus[i]);
        cplex.addMIPStart(vars, vals);
        injected_heuristic_solutions++;
    }

    // ========= Symmetry break (no leaks) ========
    {
        const int vstar = (int)graph.max_degree_vertex();
        model.add(x[vstar] == (s[vstar] == -1 ? 1 : 0));
    }
    
	// ===== Net-degree cuts (iff --netdeg) — exact spec =====
	if (use_cut_generator & NET_DEGREE_CUTS) {
	    int added = 0, eligible = 0;
	
	    for (int u = 0; u < n; ++u) {
	        const int dpos = d_plus[u];
	        const int dneg = d_minus[u];
	        const int dsig = dpos - dneg;                      // d_u(σ)
	        const int rhs  = (int)std::floor(0.5 * (double)dsig); // ⌊d_u(σ)/2⌋
	
	        // Only add when not dominated: ⌊d_u(σ)/2⌋ < d_u^-(σ)
	        if (rhs < dneg) {
	            ++eligible;
	
	            IloExpr E(env);
	            E += (double)dsig * x[u];
	
	            // Sum over edges incident to u, using σ(u,v) ∈ {+1,-1}
	            for (const auto& [e, sgn] : signs) {
	                if (e.first == u || e.second == u) {
	                    const int v   = (e.first == u) ? e.second : e.first;
	                    const int idx = edge_index[e];
	                    // -(2 y_uv - x_v) * σ(u,v)
	                    E -= (2.0 * y[idx] - x[v]) * (double)sgn;
	                }
	            }
	
	            model.add(E <= (double)rhs);
	            E.end();
	            ++added;
	        }
	    }
	
	    net_degree_cut_count += added;
	    std::cout << "[INIT-NETDEG] eligible=" << eligible
	              << " added=" << added << "\n";
	}

	// ===== Step 3–4: Weighted–SP seed of hard edge-disjoint cycles =====
	// Temporary working costs ω^seed on E^+_{σ_s}: start at 1; inflate by |V| on every accepted cycle.
	if (use_cut_generator & NEGATIVE_CYCLE_CUTS) {
	    const int N = n; // |V|
	    const int M = m; // |E|
	
	    // Switched weights decide sign under current σ_s
	    const auto& sw = graph.get_switched_weight();
	
	    // Build adjacency on the positive subgraph under σ_s with (neighbor, full_eid)
	    std::vector<std::vector<std::pair<int,int>>> adj((size_t)N);
	    for (int eid = 0; eid < M; ++eid) {
	        if (sw[eid] > 0.0) {
	            const Edge& e = edge_reverse.at(eid);
	            adj[(size_t)e.first ].emplace_back(e.second, eid);
	            adj[(size_t)e.second].emplace_back(e.first , eid);
	        }
	    }
	
	    // Working costs ω^seed : only meaningful on positive edges; keep an array size M for convenience.
	    std::vector<double> cost((size_t)M, 1e12);
	    for (int eid = 0; eid < M; ++eid) if (sw[eid] > 0.0) cost[(size_t)eid] = 1.0;
	
	    auto dijkstra = [&](int src, int dst,
	                        std::vector<int>& parent_v,
	                        std::vector<int>& parent_e) -> bool
	    {
	        const double INF = 1e100;
	        parent_v.assign(N, -1);
	        parent_e.assign(N, -1);
	        std::vector<double> dist((size_t)N, INF);
	
	        using QN = std::pair<double,int>;
	        std::priority_queue<QN, std::vector<QN>, std::greater<QN>> pq;
	
	        dist[(size_t)src] = 0.0;
	        pq.emplace(0.0, src);
	
	        while (!pq.empty()) {
	            auto [du, u] = pq.top(); pq.pop();
	            if (du != dist[(size_t)u]) continue;
	            if (u == dst) break;
	            for (const auto& [v, feid] : adj[(size_t)u]) {
	                const double w = cost[(size_t)feid]; // ω^seed(feid)
	                const double nd = du + w;
	                if (nd + 1e-15 < dist[(size_t)v]) {
	                    dist[(size_t)v] = nd;
	                    parent_v[(size_t)v] = u;
	                    parent_e[(size_t)v] = feid; // full eid of positive edge (u,v)
	                    pq.emplace(nd, v);
	                }
	            }
	        }
	        return dist[(size_t)dst] < INF/2;
	    };
	
	    // Collect all negative edges in current switching
	    std::vector<Edge> neg_edges; neg_edges.reserve(M);
	    for (int eid = 0; eid < M; ++eid) if (sw[eid] < 0.0) neg_edges.push_back(edge_reverse.at(eid));
	
	    int tri_seed_added = 0, cyc_seed_added = 0, seen = 0;
	
	    std::vector<int> par_v, par_e, order_nodes;
	    order_nodes.reserve((size_t)N);
	
	    for (const Edge& neg : neg_edges) {
	        const int u = neg.first, v = neg.second;
	
	        if (!dijkstra(u, v, par_v, par_e)) continue; // disconnected: no cycle from this anchor
	
	        // Reconstruct positive path nodes u -> ... -> v
	        order_nodes.clear();
	        for (int cur = v; cur != -1; cur = par_v[(size_t)cur]) order_nodes.push_back(cur);
	        std::reverse(order_nodes.begin(), order_nodes.end());
	        if ((int)order_nodes.size() < 2) continue;
	
	        // Inflate ω^seed by |V| along the positive edges of the path to induce edge-disjointness
	        std::vector<int> pos_full_eids; pos_full_eids.reserve(order_nodes.size());
	        for (int i = 1; i < (int)order_nodes.size(); ++i) {
	            const int a = order_nodes[i-1], b = order_nodes[i];
	            const int fe = par_e[(size_t)b]; // full eid of (a,b) on the shortest path
	            if (fe >= 0) {
	                pos_full_eids.push_back(fe);
	                cost[(size_t)fe] = std::max(1e-12, cost[(size_t)fe] + (double)N); // ω^seed(e) += |V|
	            }
	        }
	
	        // Decide whether to add this cycle according to input flags
	        const bool is_triangle = ((int)pos_full_eids.size() == 2);
	        const bool want_tri    = (use_cut_generator & TRIANGLE_CUTS) != 0;
	        const bool want_cyc    = (use_cut_generator & NEGATIVE_CYCLE_CUTS) != 0;
	
	        if ( (is_triangle && !want_tri) || (!is_triangle && !want_cyc) ) {
	            ++seen; // counted for seed cap even if we didn't add a cut
	            continue;
	        }
	
	        // Build edge sequence [neg_edge, pos_path...] for the generic cycle cut generator
	        std::vector<Edge> cycle_edges;
	        cycle_edges.reserve(1 + pos_full_eids.size());
	        cycle_edges.push_back(neg);
	        for (int fe : pos_full_eids) cycle_edges.push_back(edge_reverse.at(fe));
	
	        auto cuts = generate_cycle_cuts(env, cycle_edges); // generic (works for triangles and longer cycles)
	        for (auto& kv : cuts) model.add(kv.first);
	
	        if (is_triangle) tri_seed_added += (int)cuts.size();
	        else             cyc_seed_added += (int)cuts.size();
	
	        ++seen;
	    }
	
	    std::cout << "[INIT-SEED] triangles_added=" << tri_seed_added
	              << " negcycles_added=" << cyc_seed_added
	              << " scanned=" << seen << "/" << neg_edges.size() << "\n";
	}

    // ===== One build-guided pipeline round (no LP hints) just to advance state =====
    // Reuse the switching we just computed; x̂ = s01(s), ŷ = x̂_u x̂_v.
//    try {
//        ensure_driver_();
//        driver().reuse_current_switching_once();
//
//        std::vector<double> xhat(n, 0.0), yhat(m, 0.0);
//        auto s01v = [&](int si){ return (si == -1) ? 1.0 : 0.0; };
//        for (int i = 0; i < n; ++i) xhat[i] = s01v(s[i]);
//        for (const auto& [e, idx] : edge_index)
//            yhat[idx] = xhat[e.first] * xhat[e.second];
//
//        auto R = driver().run_round(&xhat, &yhat, TriangleCyclePipeline::Phase::Build);
//	
//	    // Filter keys by input flags: triangles vs. longer cycles
//	    const bool want_tri = (use_cut_generator & TRIANGLE_CUTS) != 0;
//	    const bool want_cyc = (use_cut_generator & NEGATIVE_CYCLE_CUTS) != 0;
//	
//	    std::size_t tri_keys = 0, cyc_keys = 0;
//	    int cuts_added = 0;
//	
//	    for (const auto& k : R.accepted_keys) {
//	        const bool is_triangle = (k.pos.size() == 2);
//	        if ((is_triangle && !want_tri) || (!is_triangle && !want_cyc)) {
//	            continue; // respect flags
//	        }
//	
//	        // Materialize the cycle as a list of canonical edges [neg, pos...]
//	        std::vector<Edge> all_edges;
//	        all_edges.reserve(1 + k.pos.size());
//	        all_edges.push_back(k.neg);
//	        all_edges.insert(all_edges.end(), k.pos.begin(), k.pos.end());
//	
//	        // Build the corresponding linear cuts and add to the model
//	        auto cuts = generate_cycle_cuts(env, all_edges);
//	        for (auto& kv : cuts) model.add(kv.first);
//	        cuts_added += static_cast<int>(cuts.size());
//	
//	        if (is_triangle) ++tri_keys; else ++cyc_keys;
//	    }
//	
//	    std::cout << "[INIT-INJECT] tri_sel=" << R.triangles_accepted
//	              << " sp_sel=" << R.cycles_accepted
//	              << " keys_total=" << R.accepted_keys.size()
//	              << " tri_keys=" << tri_keys
//	              << " cyc_keys=" << cyc_keys
//	              << " cuts_added=" << cuts_added << "\n";
//    } catch (const std::exception& e) {
//        std::cout << "[INIT-INJECT] error: " << e.what() << "\n";
//    } catch (...) {
//        std::cout << "[INIT-INJECT] error: unknown exception\n";
//    }

    // ============== Done ==============
    std::cout << "[BUILD] n=" << n << " m=" << m << " done in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - T0).count()
              << " ms\n";
}

void FrustrationModelXY::configure_separation(const SeparationConfig& cfg) {
    sep_cfg_ = cfg;
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
        sep_state_.init_sizes_if_needed(graph);
        // Sync EMA knobs into persistent (mirrors SeparationConfig)
        sep_state_.ema_delta = sep_cfg_.ema.delta;
        sep_state_.ema_rho   = sep_cfg_.ema.rho;
        sep_state_.ema_kappa = sep_cfg_.ema.kappa;
        driver_ = std::make_unique<TriangleCyclePipeline>(graph, sep_state_, sep_cfg_);
    }
}

void FrustrationModelXY::solve() {
    // One-shot probe of the configured separation knobs (before callbacks)
    const auto& CFG = sep_cfg_;
    std::cout << "[SEP-CONFIG] B_tri=" << CFG.budget.B_tri
              << ", alpha=" << CFG.ranking.alpha
              << ", theta=" << CFG.ranking.theta
              << ", lambda_hist=" << CFG.ranking.lambda_hist
              << ", lambda_LP=" << CFG.ranking.lambda_LP
              << ", omega_eps=" << CFG.weights.omega_eps
              << ", omega_max=" << CFG.weights.omega_max
              << "\n";
    std::cout << "[SEP-CONFIG] anneal{tau,v0,gamma_min,gamma_max,B_min}="
              << "{" << CFG.anneal_tri.tau
              << "," << CFG.anneal_tri.v0
              << "," << CFG.anneal_tri.gamma_min
              << "," << CFG.anneal_tri.gamma_max
              << "," << CFG.anneal_tri.B_min << "}" << "\n";

    // Register full pipeline separation callback
    cplex.use(new (env) TriangleCycleCutGenerator(env, *this));
    // Keep heuristic
    cplex.use(new (env) SwitchingHeuristicCallback(env, *this));

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

    FrustrationModel::export_solution(file_prefix, with_svg, partition);
}

// ────────────────────────────────────────────────────────────────────────────
// Separation utilities
// ────────────────────────────────────────────────────────────────────────────

void FrustrationModelXY::snapshot_lp_solution(const LpAccessor& acc, std::vector<double>& xhat, std::vector<double>& yhat) const {
    xhat.assign(x.size(), 0.0);
    yhat.assign(y.size(), 0.0);

    try {
        IloNumVarArray X(env), Y(env);
        IloNumArray    Xv(env), Yv(env);
        for (const auto& xv : x) X.add(xv);
        for (const auto& yv : y) Y.add(yv);
        acc.getValues(Xv, X);
        acc.getValues(Yv, Y);
        for (IloInt i = 0; i < Xv.getSize() && i < (IloInt)xhat.size(); ++i) xhat[(size_t)i] = Xv[i];
        for (IloInt i = 0; i < Yv.getSize() && i < (IloInt)yhat.size(); ++i) yhat[(size_t)i] = Yv[i];
    } catch (...) {
        for (size_t i = 0; i < x.size(); ++i) xhat[i] = cplex.getValue(x[i]);
        for (size_t i = 0; i < y.size(); ++i) yhat[i] = cplex.getValue(y[i]);
    }
}


TriangleCyclePipeline::Result
FrustrationModelXY::run_pipeline_round_from_lp(const LpAccessor& acc) {
    std::vector<double> xhat, yhat; snapshot_lp_solution(acc, xhat, yhat);
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
        all_edges.emplace_back(k.neg.first, k.neg.second);

        // positives, already canonicalized as (min,max)
        for (const auto& p : k.pos) all_edges.emplace_back(p.first, p.second);

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
        idx.push_back(edge_index[e]);     // <- use operator[] on view
        sgn.push_back(signs[e].sign);     // <- use operator[] on view
    }
    // Build via helper (no flips)
    return make_cycle_cut_from_prelookups(env, all_edges, idx, sgn, nullptr, y, x);
}

// Base virtual in FrustrationModel is pure, so provide a minimal override
std::vector<std::pair<IloRange, std::string>>
FrustrationModelXY::generate_cycle_cuts(IloEnv& env, const std::vector<Edge>& all_edges) {
    std::vector<std::pair<IloRange, std::string>> cuts;
    cuts.emplace_back(generate_cycle_cut_standard(env, all_edges), "standard");
    return cuts;
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
        const SignedEdge& edge_info = signs[e];     // <- use operator[] on view
        int index = edge_index[e];                  // <- use operator[] on view
        expr += ((y[index] - 0.5 * x[e.first] - 0.5 * x[e.second]) * edge_info.sign * sign_flip);
        if (edge_info.sign * sign_flip < 0) rhs++;
    };

    for (int id : t.pending_cut_ids) {
        if (id == -1) {
            IloExpr base_expr(env);
            int base_rhs = 0;
            for (const auto& e : t.triangle) add_term(e, base_expr, base_rhs);
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
// Separation callback delegator (pipeline)
// ────────────────────────────────────────────────────────────────────────────
void FrustrationModelXY::TriangleCycleCutGenerator::main() {
    IloEnv ENV = getEnv();
    try {
        // Callback-safe accessor (must not use model-side API in legacy callbacks)
        UserCutCallbackAccessor acc(*this);

        // One driver round at the current LP solution
        auto R = owner.run_pipeline_round_from_lp(acc);

        // Tuning/knobs
        const double VIOL_TOL = 1e-6;
        const auto&  rcfg     = owner.sep_cfg_.reheat;   // SeparationConfig::ReheatPool
        const double ema_alpha = rcfg.reheat_ema_alpha;

        int cuts_added = 0;
        int reheat_added_now = 0;
        int reheat_requeued  = 0;

        // Stage up to reheat_stage_cap non-violated keys (we commit after the loop)
        std::vector<std::pair<fmkey::CycleKey, ReheatItem>> staged;
        staged.reserve(static_cast<size_t>(rcfg.reheat_stage_cap));

        // Process each accepted key once
        for (const auto& key : R.accepted_keys) {
            // Build all ranges this cycle induces (could be >1 per cycle)
            std::vector<fmkey::CycleKey> singleton{key};
            auto cuts = owner.build_cycle_cuts_from_keys(ENV, singleton);

            bool any_violated = false;
            double max_viol = 0.0;

            // Evaluate each inequality at the current LP point
            for (auto& kv : cuts) {
                IloRange  rng   = kv.first;
                const char* lbl = kv.second.c_str();

                const double lb  = rng.getLB();
                const double ub  = rng.getUB();
                const double act = getValue(rng.getExpr());

                // Generic violation amount (handles <=, >=, ==)
                double viol_amt = 0.0;
                if (lb != -IloInfinity) viol_amt = std::max(viol_amt, lb - act);
                if (ub !=  IloInfinity) viol_amt = std::max(viol_amt, act - ub);

                const bool violated = (viol_amt > VIOL_TOL);
                max_viol = std::max(max_viol, violated ? viol_amt : 0.0);

                if (violated) {
                    add(rng).setName(lbl);
                    ++cuts_added;
                    any_violated = true;
                }

                // Concert copies the range into the pool on add(); safe to end handle
                rng.end();
            }

            // If nothing was violated, stage this key for reheat (respect per-round cap)
            if (!any_violated && static_cast<int>(staged.size()) < rcfg.reheat_stage_cap) {
                ReheatItem item;
                item.ttl       = rcfg.reheat_default_ttl;
                item.ema_viol  = max_viol;   // 0 if no viol; we still seed EMA
                item.last_viol = max_viol;
                // item.cyc_vertices left empty (optional); we can reconstruct later from key
                staged.emplace_back(key, std::move(item));
            } else if (any_violated) {
                // If this came from the pool, count as "added from reheat"
                if (g_reheat_pool.contains(key)) ++reheat_added_now;
                // Ensure it’s not lingering in the pool anymore
                if (g_reheat_pool.contains(key)) {
                    // erase by key
                    if (auto it = g_reheat_pool.find_mutable(key)) {
                        // iterator erase API available on ReheatPool
                        for (auto it2 = g_reheat_pool.begin(); it2 != g_reheat_pool.end(); ) {
                            if (fmkey::CycleKeyEq{}(it2->first, key)) { it2 = g_reheat_pool.erase(it2); }
                            else { ++it2; }
                        }
                    }
                }
            } else {
                // was non-violated but staging cap already reached → requeue signal (best-effort)
                if (g_reheat_pool.contains(key)) ++reheat_requeued;
            }
        }

        // === Commit staged (non-violated) items to the pool, then enforce global cap ===
        for (auto& kv : staged) {
            const auto& key  = kv.first;
            auto&       item = kv.second;

            if (g_reheat_pool.contains(key)) {
                // Refresh TTL and EMA for already-present key
                g_reheat_pool.update_ema(key, item.last_viol, ema_alpha);
                if (auto* p = g_reheat_pool.find_mutable(key)) p->ttl = rcfg.reheat_default_ttl;
                ++reheat_requeued;
            } else {
                g_reheat_pool.upsert(key, item);
                ++reheat_added_now;
            }
        }

        // Hard clip to global cap (no TTL decrement here—keep decay in your scheduled stage)
        if (g_reheat_pool.size() > rcfg.reheat_pool_cap) {
            std::size_t overshoot = g_reheat_pool.size() - rcfg.reheat_pool_cap;
            for (auto it = g_reheat_pool.begin(); it != g_reheat_pool.end() && overshoot > 0; ) {
                it = g_reheat_pool.erase(it);
                --overshoot;
            }
        }

        // Telemetry
        std::cout << "[CALLBACK] tri_selected=" << R.triangles_accepted
                  << " sp_cycles="  << R.cycles_accepted
                  << " cuts_added="  << cuts_added
                  << " reheat_staged=" << staged.size()
                  << " reheat_added="  << reheat_added_now
                  << " reheat_requeued=" << reheat_requeued
                  << " pool_size_now=" << g_reheat_pool.size()
                  << "\n";

        // Budget summary
        const int tri_sel = R.triangles_accepted;
        const int sp_cyc  = R.cycles_accepted;
        const int base_B  = owner.separation_config().budget.B_tri;
        const int next_B  = owner.separation_state().B_tri_cur;
        const int used_B  = (next_B > 0 ? std::min(base_B, next_B) : base_B);

        std::cout << "[CALLBACK] tri_selected=" << tri_sel
                  << " sp_cycles=" << sp_cyc
                  << " cuts_added=" << cuts_added << "\n";
        std::cout << "[CALLBACK] budget_used(B_tri)=" << used_B
                  << "  anneal_next=" << next_B << "\n";
    } catch (const IloCplex::Exception& e) {
        std::cout << "[CALLBACK] exception: " << e << "\n";
    }
}

// ────────────────────────────────────────────────────────────────────────────
 // Heuristic callback (unchanged)
// ────────────────────────────────────────────────────────────────────────────

void FrustrationModelXY::SwitchingHeuristicCallback::main() {
    std::vector<double> x_vals(owner.x.size()), y_vals(owner.y.size());
    for (int i = 0; i < (int)owner.x.size(); ++i) x_vals[i] = getValue(owner.x[i]);
    for (int i = 0; i < (int)owner.y.size(); ++i) y_vals[i] = getValue(owner.y[i]);

    owner.graph.weighting_from_fractional(x_vals, y_vals);

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

    auto frac_result = owner.graph.fractional_greedy_switching(UB_OPTS);
    owner.graph.restore_switched_sign();

    std::shared_ptr<const std::vector<int>> s_ptr;
    if (frac_result.has_value()) {
        s_ptr = frac_result.value();
        std::cout << "[Heuristic] switching solution (fractional) size=" << s_ptr->size() << "\n";
    } else {
        auto s_fb = owner.graph.greedy_switching();
        if (s_fb.empty()) {
            std::cout << "[Heuristic] No switching solution found.\n";
            return;
        }
        s_ptr = std::make_shared<const std::vector<int>>(std::move(s_fb));
        std::cout << "[Heuristic] switching solution (fallback) size=" << s_ptr->size() << "\n";
    }

    const std::vector<int>& s = *s_ptr;

    const auto& d_plus_eval  = owner.graph.plus_degrees_view();
    const auto& d_minus_eval = owner.graph.minus_degrees_view();
    auto s01 = [](int si) -> double { return (si == -1) ? 1.0 : 0.0; };

    double obj_val = 0.0;
    IloNumVarArray xy_vars(getEnv());
    IloNumArray    xy_vals(getEnv());

    for (int i = 0; i < (int)owner.x.size(); ++i) {
        double xi = s01(s[i]);
        xy_vars.add(owner.x[i]);
        xy_vals.add(xi);
        obj_val += (d_plus_eval[i] - d_minus_eval[i]) * xi;
    }

    for (const auto& [edge, idxe] : owner.edge_index) {
        int u = edge.first, v = edge.second;
        int xu = (int)s01(s[u]);
        int xv = (int)s01(s[v]);

        double w = owner.weights[edge];  // <- use operator[] on view

        double y0 = (w > 0.0)
                  ? (double)std::min(xu, xv)
                  : (double)std::max(0, xu + xv - 1);

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
