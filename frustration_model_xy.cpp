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
    const IloNumVarArray& y,
    const IloNumVarArray& x)
{
    IloExpr expr(env);
    int rhs = 0;
    for (size_t j = 0; j < all_edges.size(); ++j) {
        const Edge& e = all_edges[j];
        int eff = sgn[j];
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
    : FrustrationModel(g, cut_flags), x(env), y(env), xval(env) {}

void FrustrationModelXY::build() {
    using clock = std::chrono::steady_clock;
    const auto T0 = clock::now();

    const int n = graph.vertex_count();
    const int m = graph.edge_count();

    const auto& d_plus  = graph.get_pos_degrees();
    const auto& d_minus = graph.get_neg_degrees();

    // ================== Vars ==================
    for (int i = 0; i < n; ++i)
        x.add(IloBoolVar(env, ("x_" + std::to_string(i)).c_str()));
    for (int i = 0; i < m; ++i)
        y.add(IloNumVar(env, 0.0, 1.0, ILOFLOAT, ("y_" + std::to_string(i)).c_str()));

	model.add(x);
	model.add(y);
    xhat.assign(x.getSize(), 0.0);
	
    // ================== Objective ==================
    {
        IloExpr obj(env);
        for (int i = 0; i < n; ++i) obj += (d_plus[i] - d_minus[i]) * x[i];
        for (const auto& [e, sgn, _] : signs) obj += -2.0 * sgn * y[edge_index[e]];
        objective = IloMinimize(env, obj);
        model.add(objective);
        obj.end();
    }

    // ============ Base formulation (edge constraints) ============
    {
        IloRangeArray base(env);
        for (const auto& [e, sgn, _] : signs) {
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
	            for (const auto& [e, sgn, _] : signs) {
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

    // ============ Greedy switching ============
    SignedGraphForMIP sg_sw = graph.greedy_switching();       // switched graph (by value)

	// ===== Step 3–4: Weighted–SP seed of hard edge-disjoint cycles =====
	// Temporary working costs ω^seed on E^+_{σ_s}: start at 1; inflate by |V| on every accepted cycle.
	if (use_cut_generator & NEGATIVE_CYCLE_CUTS) {
	    const int N = n; // |V|
	    const int M = m; // |E|
	
	    // Build adjacency on the positive subgraph under σ_s with (neighbor, full_eid)
	    std::vector<std::vector<std::pair<int,int>>> adj((size_t)N);
	    for (int eid = 0; eid < M; ++eid) {
	        if (sg_sw.is_pos_edge(eid)) {
	            const Edge& e = edge_reverse.at(eid); // same topology ⇒ same eid↔edge mapping
	            adj[(size_t)e.first ].emplace_back(e.second, eid);
	            adj[(size_t)e.second].emplace_back(e.first , eid);
	        }
	    }
	
	    // Working costs ω^seed : only meaningful on positive edges; keep an array size M for convenience.
	    std::vector<double> cost((size_t)M, 1e12);
	    for (int eid = 0; eid < M; ++eid) if (sg_sw.is_pos_edge(eid)) cost[(size_t)eid] = 1.0;
	
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
	
	    // Enumerate negative edges from the switched graph (no full scan needed)
	    int tri_seed_added = 0, cyc_seed_added = 0, seen = 0;
	    std::vector<int> par_v, par_e, order_nodes;
	    order_nodes.reserve((size_t)N);
	
	    for (auto [neg, neg_eid] : sg_sw.negative_edges_view()) {
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
	
	        // Build edge sequence [neg_edge, pos_path...] for the generic cycle cut generator
	        std::vector<Edge> cycle_edges;
	        cycle_edges.reserve(1 + pos_full_eids.size());
	        cycle_edges.push_back(neg);
	        for (int fe : pos_full_eids) cycle_edges.push_back(edge_reverse.at(fe));
	
	        auto cuts = generate_cycle_cuts(env, cycle_edges); // generic (works for triangles and longer cycles)
	        for (auto& kv : cuts) model.add(kv.first);
	
	        const bool is_triangle = ((int)cycle_edges.size() == 3);
	        if (is_triangle) tri_seed_added += (int)cuts.size();
	        else             cyc_seed_added += (int)cuts.size();
	
	        ++seen;
	    }
	
	    std::cout << "[INIT-SEED] triangles_added=" << tri_seed_added
	              << " negcycles_added=" << cyc_seed_added
	              << " scanned=" << seen << "\n";
	}

    // =============== MIP start ===============
    {
        IloNumVarArray vars(env);
        IloNumArray vals(env);

		// Build candidate start from sg_sw
		for (int u = 0; u < n; ++u) {
		    // s[u] is ±1 here; turn into {0,1}
		    const double xu = sg_sw.get_x(u);
		    vars.add(x[u]); vals.add(xu);
		}
		for (const auto& [e, idx] : edge_index) {
		    const double yprod = sg_sw.get_y(e);     // product semantics
            vars.add(y[idx]); vals.add(yprod);
		}

        for (int i = 0; i < n; ++i) cplex.setPriority(x[i], d_plus[i]); // + d_minus[i]);
        cplex.addMIPStart(vars, vals);
        injected_heuristic_solutions++;
    }

    // ========= Symmetry break (no leaks) ========
    {
		const int vstar = (int)graph.max_pos_degree_vertex();
		const double xfix = (sg_sw.get_x(vstar) >= 0.5 ? 1.0 : 0.0);
		model.add(x[vstar] == xfix);

		sym_vstar = vstar;
		sym_xfix  = xfix;
    }

    // ============== Done ==============
    std::cout << "[BUILD] n=" << n << " m=" << m << " done in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - T0).count()
              << " ms\n";
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
    if (use_cut_generator & (NEGATIVE_CYCLE_CUTS | TRIANGLE_CUTS))
    	cplex.use(new (env) TriangleCycleCutGenerator(env, *this));
    // Keep heuristic
    cplex.use(new (env) SwitchingHeuristicCallback(env, *this));

    FrustrationModel::solve();
}

// ────────────────────────────────────────────────────────────────────────────
// Separation utilities
// ────────────────────────────────────────────────────────────────────────────

void FrustrationModelXY::snapshot_lp_solution(const LpAccessor& acc) const {
    acc.getValues(xval, x);
    for (IloInt i = 0; i < (IloInt)x.getSize(); ++i) xhat[(size_t)i] = xval[i];
}

// Materialize standard cycle cuts from switching-invariant fmkey::CycleKey
std::vector<std::pair<IloRange, std::string>>
FrustrationModelXY::build_cycle_cuts_from_keys(IloEnv& env,
                                               const std::vector<fmkey::CycleKey>& keys) const
{
    std::vector<std::pair<IloRange, std::string>> out;
    out.reserve(keys.size()); // lower bound; actual count may be larger

    for (const auto& k : keys) {
        // Neg edge (a,b) + all positive edges from key (ordered edges form a cycle)
        std::vector<Edge> all_edges;
        all_edges.reserve(1 + k.pos.size());

        // neg (a,b)
        all_edges.emplace_back(k.neg.first, k.neg.second);

        // positives (expected to be in the cycle’s order)
        for (const auto& p : k.pos) all_edges.emplace_back(p.first, p.second);

        // Generate ALL inequalities associated to this cycle
        auto cuts = generate_cycle_cuts(env, all_edges);

        // Append (move) into output
        out.insert(out.end(),
                   std::make_move_iterator(cuts.begin()),
                   std::make_move_iterator(cuts.end()));
    }

    return out;
}

// Base virtual in FrustrationModel is pure, so provide a minimal override
std::vector<std::pair<IloRange, std::string>>
FrustrationModelXY::generate_cycle_cuts(IloEnv& env, const std::vector<Edge>& all_edges) const
{
    // --- Precompute edge ids and base (LIVE) signs in the order of all_edges
    std::vector<int> idx;    idx.reserve(all_edges.size());
    std::vector<int> sgn0;   sgn0.reserve(all_edges.size());
    for (const auto& e : all_edges) {
        const int eid = edge_index[e];      // full-edge id for e
        idx.push_back(eid);
        sgn0.push_back(signs0[eid]);      // LIVE sign (post-switching)
    }

    // --- Collect unique vertices in deterministic order and build incidence map
    std::vector<int> verts; verts.reserve(all_edges.size());
    std::unordered_map<int, std::vector<int>> incident; // v -> edge indices in all_edges
    incident.reserve(all_edges.size());

    auto touch_vertex = [&](int v){
        if (incident.find(v) == incident.end()) {
            incident.emplace(v, std::vector<int>{});
            verts.push_back(v);
        }
    };
    for (int i = 0; i < (int)all_edges.size(); ++i) {
        const auto& e = all_edges[(size_t)i];
        touch_vertex(e.first);
        incident[e.first].reserve(2); incident[e.first].push_back(i);
    }
    for (int i = 0; i < (int)all_edges.size(); ++i) {
        const auto& e = all_edges[(size_t)i];
        incident[e.second].push_back(i);
    }

    const int n = (int)verts.size();
    const int kmax = n / 2;

    // Enumerate all subsets S with |S| ≤ floor(n/2) via DFS (sequential flips).
    std::vector<std::pair<IloRange, std::string>> cuts;
    cuts.reserve(1u << std::min(kmax, 12));  // heuristic pre-reserve

    std::vector<int> sgn = sgn0; // mutable working signs
    std::string cycle_type = (all_edges.size() == 3) ? "triangle" : "cycle";

    std::function<void(int,int)> dfs = [&](int start, int chosen){
        // Emit inequality for the current subset S (of size = chosen)
        if (chosen <= kmax) {
            cuts.emplace_back(make_cycle_cut_from_prelookups(env, all_edges, idx, sgn, y, x),
            (chosen == 0) ? std::string("std-")+cycle_type : std::string("rev-")+cycle_type);
        } else {
            return; // guard (shouldn’t happen due to the for-loop pruning below)
        }
        // Try adding more vertices to S, respecting |S| ≤ kmax
        for (int i = start; i < n; ++i) {
            if (chosen + 1 > kmax) break;
            const int v = verts[i];

            // --- Flip signs for all edges incident to v (sequentially)
            for (int ei : incident[v]) sgn[(size_t)ei] *= -1;

            dfs(i + 1, chosen + 1);

            // --- Undo flips (backtrack)
            for (int ei : incident[v]) sgn[(size_t)ei] *= -1;
        }
    };

    dfs(0, 0);

    return cuts;
}

// ────────────────────────────────────────────────────────────────────────────
// Separation callback delegator (pipeline)
// ────────────────────────────────────────────────────────────────────────────
void FrustrationModelXY::TriangleCycleCutGenerator::main() {
    IloEnv ENV = getEnv();
    try {
        // Callback-safe accessor (must not use model-side API in legacy callbacks)
        UserCutCallbackAccessor acc(*this);

        // Snapshot LP and reseed switching from x̂ for coherence with heuristic callback
        owner.snapshot_lp_solution(acc);
        owner.graph.reseed_switching(owner.xhat);

		const SignedGraphForMIP::GreedyKickOptions UB_OPTS = [&]{
		    SignedGraphForMIP::GreedyKickOptions o;
		    o.neg_edge_threshold_abs        = -1;
		    o.neg_edge_threshold_frac       = 0.03;
		    o.max_kicks                     = 3;
		    o.use_weighted_degree           = true;
		    o.use_triangle_tiebreak         = true;
		    o.triangle_beta                 = 0.08;
		    o.neighbor_cap                  = 1024;
		    o.triangle_cap_per_u            = 1024;
		    o.relax_to_all_pos_if_Z0_empty  = true;
		    o.delta_m_minus_cap             = 512;
		    o.delta_m_minus_penalty         = 0.0;
		    o.R_max                         = 10;
		    o.Delta                         = 1;
		    return o;
		}();

		owner.graph.apply_greedy_switching(UB_OPTS);

        // One driver round at the current LP solution
        auto R = owner.driver().run_round(owner.graph);

        // Knobs from SeparationConfig
        const double VIOL_TOL = owner.separation_config().ranking.viol_tol;

        int cuts_added = 0;
        struct CutTypeTally {
		  int std_triangle = 0;
		  int std_cycle    = 0;
		  int rev_triangle = 0;
		  int rev_cycle    = 0;
		} g_cut_tally;

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
                if (violated) {
                    add(rng).setName(lbl);
                    ++cuts_added;
					if      (kv.second == "std-triangle") ++g_cut_tally.std_triangle;
					else if (kv.second == "std-cycle")    ++g_cut_tally.std_cycle;
					else if (kv.second == "rev-triangle") ++g_cut_tally.rev_triangle;
					else if (kv.second == "rev-cycle")    ++g_cut_tally.rev_cycle;
                    any_violated = true;
                    if (viol_amt > max_viol) max_viol = viol_amt;
                }

                // Concert copies the range into the pool on add(); safe to end handle
                rng.end();
            }
        }

        // Telemetry
        std::cout << "[CALLBACK] tri_selected=" << R.triangles_accepted
                  << " sp_cycles="  << R.cycles_accepted
                  << " cuts_added="  << cuts_added
                  << "\n";
		std::cout << "[CALLBACK] cuts_by_type"
		          << " std-triangle=" << g_cut_tally.std_triangle
		          << " std-cycle="    << g_cut_tally.std_cycle
		          << " rev-triangle=" << g_cut_tally.rev_triangle
		          << " rev-cycle="    << g_cut_tally.rev_cycle
		          << "\n";

        // Budget summary
        const int base_B  = owner.separation_config().budget.B_tri;
        const int next_B  = owner.separation_state().B_tri_cur;
        const int used_B  = (next_B > 0 ? std::min(base_B, next_B) : base_B);

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
    // 1) read LP (for reseed)
	HeuristicCallbackAccessor acc(*this);
	owner.snapshot_lp_solution(acc);

    // 2) reseed σ̃ from x̂
	owner.graph.reseed_switching(owner.xhat);

    // 3) crude “root node” policy and zero net-degree ratio
    const auto d_pol = owner.graph.net_degrees();
    const int n = (int)d_pol.size();
    int zero_nd = 0;
    for (int u = 0; u < n; ++u)
        if (std::fabs(d_pol[u]) <= 1e-12) ++zero_nd;

    const bool is_root = (getNnodes64() == 0);
    const double zero_ratio = double(zero_nd) / double(owner.x.getSize());
    const int max_kicks = (is_root && zero_ratio >= 0.30) ? 3 : 1;

    // 4) UB options
    SignedGraphForMIP::GreedyKickOptions UB_OPTS;
    UB_OPTS.neg_edge_threshold_abs  = -1;
    UB_OPTS.neg_edge_threshold_frac = 0.03;
    UB_OPTS.max_kicks               = max_kicks;
    UB_OPTS.use_weighted_degree     = true;
    UB_OPTS.use_triangle_tiebreak   = true;
    UB_OPTS.triangle_beta           = (is_root ? 0.08 : 0.05);
    UB_OPTS.neighbor_cap            = 1024;
    UB_OPTS.triangle_cap_per_u      = 1024;
    UB_OPTS.relax_to_all_pos_if_Z0_empty = true;
    UB_OPTS.delta_m_minus_cap       = 512;
    UB_OPTS.delta_m_minus_penalty   = 0.0;
    UB_OPTS.R_max = 10; UB_OPTS.Delta = 0;

    // 5) run greedy switching
    SignedGraphForMIP sg = owner.graph.greedy_switching(UB_OPTS);
    sg.apply_integer_projection();
    UB_OPTS.R_max = 0;
    sg.apply_greedy_switching(UB_OPTS);
   	sg.align_switching(owner.sym_vstar, owner.sym_xfix);

    IloNumVarArray xy_vars(getEnv());
    IloNumArray    xy_vals(getEnv());

	for (int u = 0; u < n; ++u) {
	    // s[u] is ±1 here; turn into {0,1}
	    const double xu = sg.get_rounded_x(u);
	    xy_vars.add(owner.x[u]); xy_vals.add(xu);
	}
	for (const auto& [e, idx] : owner.edge_index) {
	    const double yprod = sg.get_rounded_y(e);
        xy_vars.add(owner.y[idx]); xy_vals.add(yprod);
	}

	const auto neg = sg.negative_edge_count();
	const auto projNeg = sg.integer_projection().negative_edge_count();
	std::cout << "[HEURISTIC-SOL] num neg edges=" << neg
	          << " projNeg=" << projNeg
	          << " (should correlate with frustration)" << std::endl;
	
	auto eval_objective_for_assignment = [&](
	    const IloNumVarArray& vars,
	    const IloNumArray& vals)
	{
	    // Build fast lookup: var id -> value
	    std::unordered_map<IloInt, IloNum> val_by_id;
	    val_by_id.reserve(static_cast<std::size_t>(vars.getSize()));
	    for (IloInt i = 0; i < vars.getSize(); ++i)
	        val_by_id.emplace(vars[i].getId(), vals[i]);
	
	    const IloEnv env = getEnv();
	    IloExpr e(env);
	    e += owner.objective.getExpr();  // take a copy to iterate
	
	    IloNum z = e.getConstant();
	
	    // Linear terms
	    for (IloExpr::LinearIterator it = e.getLinearIterator(); it.ok(); ++it) {
	        const auto p = val_by_id.find(it.getVar().getId());
	        if (p != val_by_id.end()) z += it.getCoef() * p->second;
	    }
	
	    e.end();
	    return z;  // same numeric value as CPLEX reports (no sense flip needed)
	};
	
	// Build once (reuse across checks if you like)
	auto build_val_map = [&](const IloNumVarArray& vars, const IloNumArray& vals){
	    std::unordered_map<IloInt, double> val_by_id;
	    val_by_id.reserve(static_cast<size_t>(vars.getSize()));
	    for (IloInt i = 0; i < vars.getSize(); ++i)
	        val_by_id.emplace(vars[i].getId(), static_cast<double>(vals[i]));
	    return val_by_id;
	};
	
	auto max_violations = [&](const IloNumVarArray& vars, const IloNumArray& vals){
	    const auto val_by_id = build_val_map(vars, vals);
	    auto val = [&](IloNumVar v) -> double {
	        auto it = val_by_id.find(v.getId());
	        return (it == val_by_id.end()) ? 0.0 : it->second;
	    };
	    double base_max = 0.0, netdeg_max = 0.0;
	
	    // Base constraints
	    for (const auto& [e, sgn, _] : sg.signs_view()) {
	        int u=e.first, v=e.second, idx=owner.edge_index[e];
	        double xu = val(owner.x[u]), xv = val(owner.x[v]), yuv = val(owner.y[idx]);
	
	        if (sgn>0){
	            base_max = std::max(base_max, yuv - xu);
	            base_max = std::max(base_max, yuv - xv);
	        }else{
	            base_max = std::max(base_max, (xu + xv - 1.0) - yuv); // y >= xu+xv-1
	        }
	    }
	
	    // Net-degree cuts (same formula you used in build)
	    for (int u=0; u<sg.vertex_count(); ++u) {
	        const int dsig = sg.net_degree(u);
	        if (dsig < 0) netdeg_max++;
	    }
	
	    std::cout << "[CHECK] base_max_violation=" << base_max
	              << " netdeg_max_violation=" << netdeg_max << "\n";
	};

    // 7) post solution if improved
    const auto obj_val = eval_objective_for_assignment(xy_vars, xy_vals);
    const double incumbent = hasIncumbent() ? getIncumbentObjValue() : IloInfinity;
	
	std::cout << "[HEURISTIC-SOL] num neg edges=" << sg.negative_edge_count()
			  << " obj_val=" << obj_val
			  << " incumbent=" << incumbent
			  << "\n";
	max_violations(xy_vars, xy_vals);
    if (obj_val < incumbent) {
        setSolution(xy_vars, xy_vals, obj_val);
        owner.f_index = obj_val / owner.graph.edge_count();
        owner.injected_heuristic_solutions++;
        std::cout << "[Heuristic] injected switching-based solution; obj=" << obj_val << "\n";
    }
}

void FrustrationModelXY::configure_separation(const SeparationConfig& cfg) {
	std::cout << "[CONFIGURE SEPARATION] use_triangles=" << cfg.modes.use_triangles << "\n";
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
        sep_state_.ema_kappa = sep_cfg_.ema.kappa;
        driver_ = std::make_unique<TriangleCyclePipeline>(sep_state_, sep_cfg_);
    }
}

void FrustrationModelXY::reset_separation_state() {
    sep_state_.in_model_keys.clear();
    sep_state_.recent_keys.clear();

    sep_state_.init_sizes_if_needed(graph);
    std::fill(sep_state_.omega.begin(), sep_state_.omega.end(), 1.0);
    std::fill(sep_state_.H.begin(), sep_state_.H.end(), 0.0);
}

void FrustrationModelXY::print_separation_stats(std::ostream& os) const {
    os << "[SEP] omega.size=" << sep_state_.omega.size()
       << " H.size=" << sep_state_.H.size()
       << " recent_keys=" << sep_state_.recent_keys.size()
       << " in_model_keys=" << sep_state_.in_model_keys.size()
       << std::endl;
}

void FrustrationModelXY::export_solution(const std::string& file_prefix, bool with_svg) const {
    std::ofstream xfile(file_prefix + "_x.csv");
    std::ofstream yfile(file_prefix + "_y.csv");

    std::vector<int> partition;
    for (std::size_t i = 0; i < x.getSize(); ++i) {
        double val = cplex.getValue(x[i]);
        xfile << val << "\n";
        partition.push_back(static_cast<int>(std::round(val)));
    }
    for (std::size_t i = 0; i < y.getSize(); ++i)
        yfile << cplex.getValue(y[i]) << "\n";

    FrustrationModel::export_solution(file_prefix, with_svg, partition);
}
