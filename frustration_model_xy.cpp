// File: frustration_model_xy.cpp

#include "frustration_model_xy.h"
#include "frustration_model.h"
#include "cycle_key.h"
#include "reheat_pool.h"
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <cmath>

// ============================================================
// Helpers: key construction, reheat seeding, EMA weights
// (placed early to avoid needing to touch the end of file)
// ============================================================

// Canonical CycleKey from an ordered vertex cycle.
// Strategy:
//  - map each consecutive pair (vi, v_{i+1}) including closure to y-index
//  - choose rotation + direction (forward/reverse) with lexicographically smallest y-index sequence
//  - set rhs=0 for 1-neg cycles (our pipeline) and record 'reversed' chosen
fmkey::CycleKey
FrustrationModelXY::NegativeCycleCutGenerator::make_key_from_vertices_(
    const std::vector<int>& cyc_vertices) const
{
    fmkey::CycleKey key;
    if (cyc_vertices.size() < 3) return key;

    // Build edge index ring (undirected)
    std::vector<int> ring;
    ring.reserve(cyc_vertices.size());
    const int L = (int)cyc_vertices.size();
    for (int i = 0; i < L; ++i) {
        int u = cyc_vertices[i];
        int v = cyc_vertices[(i + 1) % L];
        auto uv = std::minmax(u, v);
        auto it = owner.edge_index.find({uv.first, uv.second});
        if (it == owner.edge_index.end()) {
            // missing edge; leave key empty (will be ignored)
            key.y_idx.clear();
            return key;
        }
        ring.push_back(it->second);
    }

    // All rotations of forward and reverse; pick lexicographically smallest
    auto best = ring;
    bool best_rev = false;
    // forward rotations
    for (int s = 1; s < L; ++s) {
        std::vector<int> cand(L);
        for (int i = 0; i < L; ++i) cand[i] = ring[(i + s) % L];
        if (cand < best) best = std::move(cand);
    }
    // reverse base
    std::vector<int> ring_rev(ring.rbegin(), ring.rend());
    auto best_rev_seq = ring_rev;
    for (int s = 1; s < L; ++s) {
        std::vector<int> cand(L);
        for (int i = 0; i < L; ++i) cand[i] = ring_rev[(i + s) % L];
        if (cand < best_rev_seq) best_rev_seq = std::move(cand);
    }
    if (best_rev_seq < best) {
        best = std::move(best_rev_seq);
        best_rev = true;
    }

    key.y_idx = std::move(best);
    key.rhs = 0;           // 1-neg cycles → floor(1/2)=0
    key.reversed = best_rev;
    return key;
}

// Build working weights on E⁺ using EMA repulsion (hist_edge_ema_)
// and LP salience (ŷ), then clamp to [ε, w_cap_].
void
FrustrationModelXY::NegativeCycleCutGenerator::build_work_weights_(
    const IloNumArray& /*x_hat*/, const IloNumArray& y_hat)
{
    const int m = (int)owner.edge_index.size();
    if ((int)hist_edge_ema_.size() != m) hist_edge_ema_.assign(m, 0.0);
    if ((int)work_w_pos_.size() != m)    work_w_pos_.assign(m, 1.0);

    auto sal = [&](double yh)->double {
        double d = std::abs(yh - 0.5);
        double s = 1.0 - std::min(1.0, 2.0 * d);
        return std::max(0.0, std::min(1.0, s));
    };

    for (const auto& kv : owner.edge_index) {
        int eid = kv.second;
        double base = 1.0; // persistent ω; currently flat (kept minimal here)
        double rep  = lambda_hist_ * std::log1p(std::max(0.0, hist_edge_ema_[eid]));
        double lp   = (eid < y_hat.getSize() ? lambda_lp_ * sal(y_hat[eid]) : 0.0);
        double w    = base + rep - lp;
        if (w < w_eps_) w = w_eps_;
        if (w_cap_ > 0.0) w = std::min(w, w_cap_);
        work_w_pos_[eid] = w;
    }
}

// Reheat: if an item now has positive slack, anchor it to its unique
// currently negative edge (if unique) and mark that anchor as "covered".
// We DO NOT emit rows here; we only seed anchors/buckets.
void
FrustrationModelXY::NegativeCycleCutGenerator::reheat_seed_buckets_(
    std::unordered_set<int>& reheated_anchors,
    const std::vector<char>& neg_edge_mask,
    const IloNumArray& /*x_hat*/, const IloNumArray& /*y_hat*/)
{
    if (reheat_pool_.empty()) return;

    for (auto it = reheat_pool_.begin(); it != reheat_pool_.end(); ) {
        const ReheatItem& item = it->second;
        const auto& cyc = item.cyc_vertices;
        if (cyc.size() < 3) { it = reheat_pool_.erase(it); continue; }

        std::vector<int> neg_eids;
        const int L = (int)cyc.size();
        for (int i = 0; i < L; ++i) {
            int u = cyc[i], v = cyc[(i+1)%L];
            auto uv = std::minmax(u, v);
            auto eit = owner.edge_index.find({uv.first, uv.second});
            if (eit == owner.edge_index.end()) continue;
            int eid = eit->second;
            if (eid >= 0 && eid < (int)neg_edge_mask.size() && neg_edge_mask[eid]) {
                neg_eids.push_back(eid);
            }
        }

        if (neg_eids.size() == 1) {
            reheated_anchors.insert(neg_eids[0]);
            ++it;
        } else {
            ++it; // keep; might become 1-neg later
        }
    }
}

// After selection: update EMA from chosen cycles (edge densities) and
// decay TTL/prune reheat entries that remain non-violated.
void
FrustrationModelXY::NegativeCycleCutGenerator::reheat_after_selection_(
    const std::vector<int>& /*chosen_cids*/,
    const std::vector<std::vector<int>>& cand_pos_edges,
    const IloNumArray& /*x_hat*/, const IloNumArray& /*y_hat*/)
{
    for (const auto& pos_edges : cand_pos_edges) {
        const double inc = (pos_edges.empty() ? 0.0 : 1.0 / (double)pos_edges.size());
        for (int eid : pos_edges) {
            if (eid >= 0 && eid < (int)hist_edge_ema_.size()) {
                hist_edge_ema_[eid] = ema_decay_ * hist_edge_ema_[eid] + inc;
            }
        }
    }
    for (auto it = reheat_pool_.begin(); it != reheat_pool_.end(); ) {
        if (it->second.ttl > 0) {
            --(it->second.ttl);
            ++it;
        } else {
            it = reheat_pool_.erase(it);
        }
    }
}
#include <algorithm>
#include <stdexcept>

FrustrationModelXY::FrustrationModelXY(SignedGraphForMIP& g, int cut_flags)
    : FrustrationModel(g, cut_flags) {}

void FrustrationModelXY::build() {
    using clock = std::chrono::steady_clock;
    const auto TB_ALL0 = clock::now();
    auto TB_mark = TB_ALL0;
    int n = graph.vertex_count();
    std::cout << "[BUILD-PHASE] start model vars" << std::endl;
    int m = graph.edge_count();
    x.resize(n);
    y.resize(m);
    std::cout << "[BUILD-PHASE] vars done in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_mark).count() << "ms" << std::endl;
    TB_mark = clock::now();

    const auto& d_plus = graph.plus_degrees_view();
    const auto& d_minus = graph.minus_degrees_view();

    for (int i = 0; i < n; ++i)
        x[i] = IloBoolVar(env, ("x_" + std::to_string(i)).c_str());
    for (int i = 0; i < m; ++i)
        y[i] = IloNumVar(env, 0.0, 1.0, ILOFLOAT, ("y_" + std::to_string(i)).c_str());

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
    std::cout << "[BUILD-PHASE] objective built in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_mark).count() << "ms" << std::endl;
    TB_mark = clock::now();

    // Basic inequalities (bulk add)
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
    std::cout << "[BUILD-PHASE] base ineq in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_mark).count() << "ms" << std::endl;
    TB_mark = clock::now();

    std::cout << "[FrustrationModel] Greedy swtching is being called." << std::endl;

	// Upper bound
	graph.restore_switched_sign();
	std::vector<int> s = graph.greedy_switching();
	
	// NEW: apply the switching so the graph's positive/negative edges match s
	graph.switching_from_partition(s);

    std::cout << "[BUILD-PHASE] greedy switching in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_mark).count() << "ms" << std::endl;
    TB_mark = clock::now();
    
    std::cout << "[FrustrationModel] Greedy swtching found." << std::endl;

    // Improving inequalities, if required
    if (use_cut_generator != NO_CUTS) {
        const auto TB_CUTS0 = TB_mark;
        long int nc = 0;
		initialize_uncut_triangles();
        // Net degree inequalities
        if (use_cut_generator & NET_DEGREE_CUTS) {
            const auto TB_ND0 = clock::now();
            IloRangeArray nd(env);
            for (int u = 0; u < graph.vertex_count(); ++u) {
                int d_u = d_plus[u] - d_minus[u];
                if ((d_u & 1) == 0) continue; // only odd
                double rhs = std::floor(d_u / 2.0);
                if (rhs >= d_minus[u]) continue; // not tightening
                IloExpr constraint(env);
                constraint += d_u * x[u];
                for (int v : graph.neighbors(u)) {
                    Edge uv = {v, u};
                    auto w = weights[uv];
                    int index = edge_index[uv];
                    constraint += w * x[v];
                    constraint -= 2 * w * y[index];
                }
                nd.add(constraint <= rhs); constraint.end();
                ++net_degree_cut_count;
            }
            if (nd.getSize() > 0) model.add(nd);
            std::cout << "[BUILD-PHASE] net-degree in "
                      << std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_ND0).count() << "ms" << std::endl;
            TB_mark = clock::now();
        }    if (use_cut_generator & NEGATIVE_CYCLE_CUTS) {
            const auto TB_NC0 = clock::now();
            // Seed the initial model with triangles only (bulk add, first triangle batch only)
            size_t std_added = 0, rev_added = 0, tri_seeded = 0;
            auto stream = graph.open_negative_cycle_stream(/*cover=*/n < 20000, /*use_triangle_order=*/true);
            std::vector<NegativeCycle> batch;
            const auto TB_NC_stream0 = clock::now();
            bool got = stream.next(batch); // fetch only the first (triangle-ordered) batch
            auto TB_NC_stream_ms = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_NC_stream0).count();
            long long TB_NC_gen_ms = 0;
            long long TB_NC_add_ms = 0;
            if (got) {
                IloRangeArray nc_arr(env);
                for (const auto& cyc : batch) {
                    const auto& pos = cyc.pos_edges();
                    // Keep only triangles: exactly 1 negative + 2 positive edges
                    if (pos.size() != 2) continue;
                    std::vector<Edge> all_edges; all_edges.reserve(3);
                    all_edges.push_back(cyc.neg_edge());
                    all_edges.insert(all_edges.end(), pos.begin(), pos.end());
                    const auto TB_NC_gen0 = clock::now();
                    IloRange rng = generate_cycle_cut_standard(env, all_edges); // standard only at build
                    TB_NC_gen_ms += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_NC_gen0).count();
                    ++tri_seeded;
                    nc_arr.add(rng); ++std_added;
                }
                const auto TB_NC_add0 = clock::now();
                if (nc_arr.getSize() > 0) model.add(nc_arr);
                TB_NC_add_ms += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_NC_add0).count();
            }
            std::cout << "[INIT-NC] batches=1, triangles=" << tri_seeded
                      << ", std_cuts=" << std_added
                      << ", rev_cuts=" << rev_added
                      << std::endl;
            std::cout << "[INIT-NC-TIME] stream=" << TB_NC_stream_ms << "ms gen=" << TB_NC_gen_ms << "ms add=" << TB_NC_add_ms << "ms" << std::endl;

            // Heuristic: if triangle seeding is far below the (original) negative edge count, fetch 1 more general batch
			int negE0 = 0; for (const auto& [e,sig] : signs) if (sig < 0) ++negE0;
            const double TRI_FRACTION_TAU = 0.20; // trigger if triangles < 20% of neg edges
            if (tri_seeded > 0 && negE0 > 0 && tri_seeded < (int)(TRI_FRACTION_TAU * (double)negE0)) {
                std::vector<NegativeCycle> extra_batch;
                auto stream2 = graph.open_negative_cycle_stream(/*cover=*/false, /*use_triangle_order=*/false);
                const auto TB_NC2_stream0 = clock::now();
                if (stream2.next(extra_batch) && !extra_batch.empty()) {
                    long long TB_NC2_gen_ms = 0, TB_NC2_add_ms = 0; size_t extra_added = 0;
                    IloRangeArray nc2(env);
                    const int MAX_EXTRA = 5000; // safety cap
                    int used = 0;
                    for (const auto& cyc : extra_batch) {
                        std::vector<Edge> all_edges; all_edges.reserve(1 + cyc.pos_edges().size());
                        all_edges.push_back(cyc.neg_edge());
                        const auto& pos = cyc.pos_edges();
                        all_edges.insert(all_edges.end(), pos.begin(), pos.end());
                        const auto TB_NC2_gen0 = clock::now();
                        IloRange rng = generate_cycle_cut_standard(env, all_edges);
                        TB_NC2_gen_ms += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_NC2_gen0).count();
                        nc2.add(rng); ++extra_added; if (++used >= MAX_EXTRA) break;
                    }
                    const auto TB_NC2_add0 = clock::now();
                    if (nc2.getSize() > 0) model.add(nc2);
                    TB_NC2_add_ms += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_NC2_add0).count();
                    auto TB_NC2_stream_ms = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_NC2_stream0).count();
                    std::cout << "[INIT-NC-EXTRA] added=" << extra_added << ", stream=" << TB_NC2_stream_ms
                              << "ms gen=" << TB_NC2_gen_ms << "ms add=" << TB_NC2_add_ms << "ms" << std::endl;
                    std_added += extra_added;
                }
            }
            std::cout << "[BUILD-PHASE] init-neg-cycles in "
                      << std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_NC0).count() << "ms" << std::endl;
            TB_mark = clock::now();
        }

        if (use_cut_generator & TRIANGLE_CUTS) {
            const auto TB_TR0 = clock::now();
            auto cuts = generate_triangle_cuts(env);
            for (auto& [cut_expr, type] : cuts) {
                model.add(cut_expr);
                if (type == "standard") {
                    ++standard_cycle_cuts_build;
                }
                else if (type == "reversed") ++reversed_cycle_cuts_build;
            }
            auto it = uncut_triangles.begin();
            while (it != uncut_triangles.end()) {
                auto& triangle = *it;
                auto id_it = triangle.pending_cut_ids.begin();
                while (id_it != triangle.pending_cut_ids.end()) {
                    id_it = triangle.pending_cut_ids.erase(id_it);
                }
                it = uncut_triangles.erase(it);
            }
            std::cout << "[BUILD-PHASE] triangles in "
                      << std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_TR0).count() << "ms" << std::endl;
            TB_mark = clock::now();
        }

        if (use_cut_generator & NEGATIVE_CYCLE_COVER_CUTS) {
            const auto TB_COV0 = clock::now();
	        auto cuts = generate_negative_cycle_cover_cuts(env);
    	    for (auto& [cut_expr, type] : cuts) {
	            model.add(cut_expr);
    	        if (type == "standard") {
					++standard_cycle_cuts_build;
	    	    }
    	        else if (type == "reversed") ++reversed_cycle_cuts_build;
            }
            std::cout << "[BUILD-PHASE] covers in "
                      << std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_COV0).count() << "ms" << std::endl;
        }
        TB_mark = clock::now();
        lower_bound = std::max(lower_bound, nc);
	}
    
	// Initial solution (feasible)
	IloNumArray start_vals(env);
	IloNumVarArray xy_vars(env);
	
	// map s ∈ {+1,-1} to x̂ ∈ {0,1}; choose 1 for s=-1, 0 for s=+1
	auto s01 = [&](int si) -> IloNum { return (si == -1) ? 1.0 : 0.0; };
	
	// x-start
	for (int i = 0; i < n; ++i) {
	    model.add(x[i]);                // ensure extracted
	    xy_vars.add(x[i]);
	    start_vals.add(s01(s[i]));      // in {0,1}
	}
	
	// y-start consistent with base constraints (tight but feasible)
	for (const auto& [edge, idxe] : edge_index) {
	    int u = edge.first, v = edge.second;
	    model.add(y[idxe]);             // ensure extracted
	    xy_vars.add(y[idxe]);
	
	    int xu = static_cast<int>(s01(s[u]));
	    int xv = static_cast<int>(s01(s[v]));
	    double w = weights[edge];       // sign of the edge
	
	    // For pos edges (w>0): y ≤ x_u and y ≤ x_v -> set y = min(x_u, x_v)
	    // For neg edges (w<0): y ≥ x_u + x_v - 1       -> set y = max(0, x_u + x_v - 1)
	    double y0 = (w > 0.0) ? static_cast<double>(std::min(xu, xv))
	                          : static_cast<double>(std::max(0, xu + xv - 1));
	    start_vals.add(y0);             // in [0,1] and satisfies constraints
	}
    for (int i = 0; i < n; ++i) {
        int priority = d_plus[i] + d_minus[i];
        cplex.setPriority(x[i], priority);
    }
    cplex.addMIPStart(xy_vars, start_vals);
    injected_heuristic_solutions++;
    
    // Fix one variable to reduce symmetries
	igraph_integer_t max_vertex = graph.max_degree_vertex();
	model.add(x[max_vertex] == (s[max_vertex] == -1 ? 1 : 0));   // 0/1
    std::cout << "[BUILD-PHASE] symmetry fix in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_mark).count() << "ms" << std::endl;

    std::cout << "[BUILD-PHASE] TOTAL build in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TB_ALL0).count() << "ms" << std::endl;
}

std::vector<std::pair<IloRange, std::string>> FrustrationModelXY::generate_pending_triangle_cuts(IloEnv& env, const TriangleInequalities& t) {
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

// Internal helper: build a (standard or flipped) cycle cut using pre-looked indices & signs.
// Place next to the existing helper (same name, different parameter types)
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
        const int eff_sign = sgn[j] * (flip ? (*flip)[j] : 1);
        expr += (y[idx[j]] - 0.5 * x[e.first] - 0.5 * x[e.second]) * eff_sign;
        if (eff_sign < 0) ++rhs;
    }
    IloRange cut = (expr <= std::floor(rhs / 2.0));
    expr.end();
    return cut;
}

// ---- Lightweight metadata to evaluate slacks without materializing IloExpr ----
struct CutMeta {
    std::vector<int> y_idx;          // indices in y
    std::vector<Edge> edges;         // endpoints per y edge (u,v)
    std::vector<int> eff_sign;       // +1 / -1 per term
    int rhs = 0;                     // floor(#neg / 2)
    bool reversed = false;           // reversed-vertex variant?
};

// Build standard (or vertex-reversed) metadata from prelookups
static CutMeta make_cutmeta_from_prelookups(
    const std::vector<Edge>& all_edges,
    const std::vector<int>& idx,
    const std::vector<int>& sgn,
    const int* flipped_vertex /*nullable*/)
{
    CutMeta cm;
    const int n = static_cast<int>(all_edges.size());
    cm.y_idx = idx;
    cm.edges = all_edges;
    cm.eff_sign.resize(n);
    int negs = 0;
    for (int j = 0; j < n; ++j) {
        const Edge& e = all_edges[j];
        int eff = sgn[j];
        if (flipped_vertex && (e.first == *flipped_vertex || e.second == *flipped_vertex)) eff = -eff;
        cm.eff_sign[j] = eff;
        if (eff < 0) ++negs;
    }
    cm.rhs = static_cast<int>(std::floor(negs / 2.0));
    cm.reversed = (flipped_vertex != nullptr);
    return cm;
}

// Materialize IloRange from metadata (used only for selected survivors)
static IloRange materialize_from_meta(IloEnv& env,
                                      const CutMeta& cm,
                                      const std::vector<IloNumVar>& y,
                                      const std::vector<IloBoolVar>& x)
{
    IloExpr expr(env);
    for (size_t j = 0; j < cm.edges.size(); ++j) {
        const Edge& e = cm.edges[j];
        const int s = cm.eff_sign[j];
        expr += (y[cm.y_idx[j]] - 0.5 * x[e.first] - 0.5 * x[e.second]) * s;
    }
    IloRange cut = (expr <= cm.rhs);
    expr.end();
    return cut;
}

IloRange FrustrationModelXY::generate_cycle_cut_standard(IloEnv& env, const std::vector<Edge>& all_edges) {
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

std::vector<std::pair<IloRange, std::string>> FrustrationModelXY::generate_cycle_cuts(IloEnv& env, const std::vector<Edge>& all_edges) {
    std::vector<std::pair<IloRange, std::string>> cuts;

    // Pre-lookup indices & signs once; collect unique vertices
    std::unordered_set<int> seen;
    std::vector<int> idx; idx.reserve(all_edges.size());
    std::vector<int> sgn; sgn.reserve(all_edges.size());
    for (const auto& e : all_edges) {
        seen.insert(e.first);
        seen.insert(e.second);
        idx.push_back(edge_index[e]);
        sgn.push_back(signs[e].sign);
    }

    // Base inequality (standard)
    {
        IloRange rng = make_cycle_cut_from_prelookups(env, all_edges, idx, sgn, nullptr, y, x);
        cuts.emplace_back(std::move(rng), "standard");
    }

    // Reversed vertex versions (flip incident edges' sign)
    for (int v : seen) {
        std::vector<int> flip(all_edges.size(), 1);
        for (size_t j = 0; j < all_edges.size(); ++j) {
            const Edge& e = all_edges[j];
            if (e.first == v || e.second == v) flip[j] = -1;
        }
        IloRange rng = make_cycle_cut_from_prelookups(env, all_edges, idx, sgn, &flip, y, x);
        cuts.emplace_back(std::move(rng), "reversed");
    }

    return cuts;
}

void FrustrationModelXY::solve() {
    if (use_cut_generator & (NEGATIVE_CYCLE_CUTS | NEGATIVE_CYCLE_COVER_CUTS)) {
        std::cout << "[Solver] Cut generator enabled." << std::endl;
        cplex.use(new (env) NegativeTriangleCutGenerator(env, *this));
        cplex.use(new (env) ConditionalCycleCutGenerator(env, *this));
    } else {
        std::cout << "[Solver] Cut generator disabled." << std::endl;
    }
    cplex.use(new (env) SwitchingHeuristicCallback(env, *this));

    FrustrationModel::solve();
}

FrustrationModelXY::NegativeTriangleCutGenerator::NegativeTriangleCutGenerator(IloEnv env, FrustrationModelXY& owner)
    : IloCplex::UserCutCallbackI(env), owner(owner)
{
}

FrustrationModelXY::NegativeTriangleCutGenerator::~NegativeTriangleCutGenerator()
{
}

void FrustrationModelXY::NegativeTriangleCutGenerator::main() {
	if (getNnodes() > 0) {
        owner.uncut_triangles.clear();
        return;
    } // Apply only at the root node

    IloEnv env = getEnv();

    static int invocation_count = 0;
    static double total_slack = 0.0;
    static int total_cuts = 0;
    static double total_unadded_slack = 0.0;
    static int total_unadded_cuts = 0;
    static bool printed_header = false;
    ++invocation_count;

    long int nc = 0;
    double slack_sum = 0.0;
    int unadded_cuts = 0;
    double unadded_slack_sum = 0.0;

    auto& triangles = owner.uncut_triangles;
    auto it = triangles.begin();
    while (it != triangles.end()) {
        bool added_any = false;
        auto& triangle_data = *it;

        auto cuts = owner.generate_pending_triangle_cuts(env, triangle_data);

        auto cuts_it = cuts.begin();
        auto id_it = triangle_data.pending_cut_ids.begin();
        while (cuts_it != cuts.end() && id_it != triangle_data.pending_cut_ids.end()) {
            auto& [cut_expr, name] = *cuts_it;
            IloNum slack = getValue(cut_expr.getExpr()) - cut_expr.getUB();
            if (slack > 2.5e-2) {
                add(cut_expr).setName(name.c_str());
                slack_sum += slack;
                ++nc;
                added_any = true;
                cuts_it = cuts.erase(cuts_it);
                id_it = triangle_data.pending_cut_ids.erase(id_it);
                owner.triangle_cut_added_last_round = true;
            } else {
                unadded_slack_sum += slack;
                ++unadded_cuts;
                cut_expr.end();
                ++cuts_it;
                ++id_it;
            }
        }

        if (triangle_data.pending_cut_ids.empty()) {
            it = triangles.erase(it);
        } else {
            ++it;
        }
    }

    total_slack += slack_sum;
    total_cuts += nc;
    total_unadded_slack += unadded_slack_sum;
    total_unadded_cuts += unadded_cuts;

    if (invocation_count == 5) {
        double mean_slack = (total_cuts > 0) ? (total_slack / total_cuts) : 0.0;
        double mean_unadded_slack = (total_unadded_cuts > 0) ? (total_unadded_slack / total_unadded_cuts) : 0.0;
        if (!printed_header) {
            std::cout << std::setw(15) << "Invocation"
                      << std::setw(15) << "Total Cuts"
                      << std::setw(20) << "Mean Slack"
                      << std::setw(20) << "Unadded Cuts"
                      << std::setw(30) << "Mean Slack (Unadded)"
                      << std::setw(25) << "Lower Bound"
                      << std::endl;
            printed_header = true;
        }
        std::cout << std::setw(15) << invocation_count
                  << std::setw(15) << total_cuts
                  << std::setw(20) << mean_slack
                  << std::setw(20) << total_unadded_cuts
                  << std::setw(30) << mean_unadded_slack
                  << std::setw(25) << getObjValue()
                  << std::endl;

        invocation_count = 0;
        total_slack = 0.0;
        total_cuts = 0;
        total_unadded_slack = 0.0;
        total_unadded_cuts = 0;
    }
}

FrustrationModelXY::NegativeCycleCutGenerator::NegativeCycleCutGenerator(
	IloEnv env,
    FrustrationModelXY& owner)
          : IloCplex::UserCutCallbackI(env), owner(owner)
{
}

FrustrationModelXY::NegativeCycleCutGenerator::~NegativeCycleCutGenerator() {
}

// Free helper to generate cut metas (no header changes)
static std::vector<CutMeta> generate_cycle_cut_metas_free(
        FrustrationModelXY& owner,
        const SignedGraphForMIP::EdgeIndexesView& edge_index,
        const SignedGraph::SignedEdgesView& signs,
        const std::vector<Edge>& all_edges,
        const std::vector<double>* x_vals_opt,
        const std::vector<int>* tri_counts_opt) {
    std::vector<CutMeta> metas;
    std::unordered_set<int> seen;
    std::vector<int> idx; idx.reserve(all_edges.size());
    std::vector<int> sgn; sgn.reserve(all_edges.size());
    for (const auto& e : all_edges) {
        seen.insert(e.first);
        seen.insert(e.second);
        idx.push_back(edge_index[e]);
        sgn.push_back(signs[e].sign);
    }
    metas.emplace_back(make_cutmeta_from_prelookups(all_edges, idx, sgn, nullptr));
    if (x_vals_opt) {
        const auto& xv = *x_vals_opt;
        const double TAU_X = 0.15;
        int top_v = -1; int top_tri = -1;
        if (tri_counts_opt && !tri_counts_opt->empty()) {
            for (int v : seen) { int t = (*tri_counts_opt)[v]; if (t > top_tri) { top_tri = t; top_v = v; } }
        }
        for (int v : seen) {
            const bool near_frac = (std::fabs(xv[v] - 0.5) <= TAU_X);
            const bool is_top = (v == top_v);
            if (!(near_frac || is_top)) continue;
            metas.emplace_back(make_cutmeta_from_prelookups(all_edges, idx, sgn, &v));
        }
    }
    return metas;
}

void FrustrationModelXY::NegativeCycleCutGenerator::main(){ // compact signature to force replacement
	IloEnv env = getEnv();
	
	using clock = std::chrono::steady_clock;
	const auto TALL0 = clock::now();
	long long T_harvest_ms = 0, T_gen_ms = 0, T_eval_ms = 0, T_select_ms = 0, T_add_ms = 0;
	
	static int invocation_count = 0;
	static double total_slack = 0.0;
	static int total_cuts = 0;
	static double total_unadded_slack = 0.0;
	static int total_unadded_cuts = 0;
	static bool printed_header = false;
	++invocation_count;
	
	#ifdef ILOCPLEX
	    long long nodes_proc = getNnodes64();
	#else
	    long long nodes_proc = 0;
	#endif
	const bool IS_ROOT = (nodes_proc == 0);
	
	IloInt rows_before = -1, rows_after = -1;
	IloInt user_before = -1, user_after = -1;   // user cut counters
	
	try { rows_before = owner.cplex.getNrows(); } catch(...) { rows_before = -1; }
	try { user_before = owner.cplex.getNcuts(IloCplex::CutUser); } catch(...) { user_before = -1; }
	// Note: CutLazyConstraint only increases in LazyConstraintCallback; keep for completeness.
	
	static const int   ROOT_K = 2;
	static const double ROOT_STALL_EPS = 2e-3;
	static const double ROOT_SAT_TAU  = 0.05;
	static double root_deltas[ROOT_K] = {0.0};
	static double root_ratios[ROOT_K] = {0.0};
	static int    root_pos = 0;
	static int    root_count = 0;
	static double root_lb_prev = 0.0; static bool root_lb_prev_set = false;
	
	if (IS_ROOT && root_count >= ROOT_K) {
	    double stall_avg = 0.0, sat_avg = 0.0;
	    for (int i=0;i<ROOT_K;++i) { stall_avg += std::fabs(root_deltas[i]); sat_avg += root_ratios[i]; }
	    stall_avg /= ROOT_K; sat_avg /= ROOT_K;
	    if (stall_avg < ROOT_STALL_EPS && sat_avg < ROOT_SAT_TAU) {
	        std::cout << "[ROOT-STOP] stall_avg=" << stall_avg << ", sat_avg=" << sat_avg << " -> branching" << std::endl;
	        return;
	    }
	}
	
	const double BASE_MIN_SLACK = 1.0e-4;
	
	// (1) Bulk-fetch values to avoid per-var callback overhead
	std::vector<double> x_vals(owner.x.size());
	std::vector<double> y_vals(owner.y.size());
	{
	    static bool vars_built = false;
	    static IloNumVarArray Xvars(getEnv());
	    static IloNumVarArray Yvars(getEnv());
	    static IloNumArray    Xvals(getEnv());
	    static IloNumArray    Yvals(getEnv());
	    if (!vars_built) {
	        Xvars.clear(); Yvars.clear();
	        for (int i = 0; i < (int)owner.x.size(); ++i) Xvars.add(owner.x[i]);
	        for (int i = 0; i < (int)owner.y.size(); ++i) Yvars.add(owner.y[i]);
	        vars_built = true;
	    }
	    Xvals.clear(); Yvals.clear();
	    try {
	        getValues(Xvals, Xvars);
	        getValues(Yvals, Yvars);
	        // Copy IloNumArray to std::vector (IloNumArray has no STL iterators)
	        int xN = Xvals.getSize();
	        if ((int)x_vals.size() != xN) x_vals.resize(xN);
	        for (int i = 0; i < xN; ++i) x_vals[i] = Xvals[i];
	        int yN = Yvals.getSize();
	        if ((int)y_vals.size() != yN) y_vals.resize(yN);
	        for (int i = 0; i < yN; ++i) y_vals[i] = Yvals[i];
	    } catch(...) {
	        // Fallback (should rarely execute)
	        for (int i = 0; i < (int)owner.x.size(); ++i) x_vals[i] = getValue(owner.x[i]);
	        for (int i = 0; i < (int)owner.y.size(); ++i) y_vals[i] = getValue(owner.y[i]);
	    }
	}
	
	// (3) Cache triangle counts and recompute lazily
	static std::vector<int> tri_counts_cache; static bool tri_valid = false; 
	const int TRI_PERIOD = 3; // recompute every 3 root invocations (or first time)
	if (!tri_valid || (IS_ROOT && (invocation_count % TRI_PERIOD == 1))) {
	    try { tri_counts_cache = owner.graph.negative_triangle_count_per_vertex(); tri_valid = true; }
	    catch(...) { tri_counts_cache.assign(owner.x.size(), 0); tri_valid = true; }
	}
	const std::vector<int>& tri_counts = tri_counts_cache;
	  
	// (2) Throttle fractional_greedy_switching: root or every 3rd invocation; lighter caps away from root
	static int fs_inv = 0; ++fs_inv;
	bool do_fs = IS_ROOT || (fs_inv % 3 == 0);
	bool has_changed = do_fs; 
	if (do_fs) {
		has_changed = owner.graph.weighting_from_fractional(x_vals, y_vals);
		if (has_changed) {
			SignedGraph::GreedyKickOptions CUTS_OPTS{
			    /*neg_edge_threshold_abs*/  -1,
			    /*neg_edge_threshold_frac*/ 0.05,
			    /*max_kicks*/               1,
			    /*use_weighted_degree*/     true,
			    /*use_triangle_tiebreak*/   IS_ROOT,
			    /*triangle_beta*/           0.05,
		        /*neighbor_cap*/            IS_ROOT ? 512 : 256,   // cheaper scans off-root
		        /*triangle_cap_per_u*/      IS_ROOT ? 768 : 384
			};
			CUTS_OPTS.kick_salience_bias = 0.5;                   // moderate nudge
			CUTS_OPTS.relax_to_all_pos_if_Z0_empty = true;        // allow KICK even if Z0=∅
			CUTS_OPTS.delta_m_minus_cap = std::max(256, (int)(0.008 * owner.graph.edge_count()));
			CUTS_OPTS.delta_m_minus_penalty = 0.0;                // keep scoring tight
			
			owner.graph.restore_switched_sign();
			auto s_opt = owner.graph.fractional_greedy_switching(CUTS_OPTS);
		    has_changed = s_opt.has_value();
    		if (has_changed) owner.switching_solution = s_opt.value();
	    }
	}
    if (!has_changed) owner.switching_solution = nullptr;

    auto stream = owner.graph.open_negative_cycle_stream(IS_ROOT /*cover*/);

    // ===== BEGIN: reheat + working weights wiring (separation only) =====
    // Rebuild salience-aware, nonnegative working weights on E^+ using (x_vals,y_vals).
    IloNumArray xhat2(env, (IloInt)owner.x.size());
    IloNumArray yhat2(env, (IloInt)owner.y.size());
    for (IloInt i = 0; i < xhat2.getSize(); ++i) xhat2[i] = (i < (IloInt)x_vals.size() ? x_vals[i] : 0.0);
    for (IloInt i = 0; i < yhat2.getSize(); ++i) yhat2[i] = (i < (IloInt)y_vals.size() ? y_vals[i] : 0.0);
    build_work_weights_(xhat2, yhat2);

    // Build current negative-edge mask under the (switched) signature.
    std::vector<char> neg_edge_mask(owner.edge_index.size(), 0);
    for (const auto& kv : owner.edge_index) {
        const Edge& e = kv.first;
        auto sit = owner.signs.find(e);
        if (sit != owner.signs.end() && sit->second.sign < 0) neg_edge_mask[kv.second] = 1;
    }

    // Seed buckets from the reheat pool: anchors to treat as already covered in this round.
    std::unordered_set<int> reheated_anchors;
    reheat_seed_buckets_(reheated_anchors, neg_edge_mask, xhat2, yhat2);
    // Optionally guard your per-negative-edge enumerators with:
    // if (reheated_anchors.count(neg_eid)) continue;
    // ===== END: reheat + working weights wiring =====

    auto mix64 = [](uint64_t x){
        x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27; x *= 0x94d049bb133111ebULL;
        x ^= x >> 31; return x;
    };
    auto sig_from_meta = [&](const CutMeta& cm)->uint64_t{
        uint64_t h = 0xcbf29ce484222325ULL;
        for (int id : cm.y_idx) h ^= mix64((uint64_t)id + 0x9e3779b97f4a7c15ULL);
        h ^= mix64((uint64_t)cm.rhs + 0x517cc1b727220a95ULL);
        h ^= mix64((uint64_t)cm.reversed);
        return h;
    };

    std::unordered_set<uint64_t> seen; seen.reserve(65536);

    std::vector<CutMeta> metas; metas.reserve(32768);
    std::vector<int> meta_batch; meta_batch.reserve(32768);
    std::vector<double> slacks; slacks.reserve(32768);
    std::vector<double> densities; densities.reserve(32768);  // (B) per-cut slack density = viol/|C|
    std::vector<int> candidates; candidates.reserve(16384);

    const double QUAL_EPS = 0.02;  // mean slack threshold
    const double QUAL_TAU = 0.05;  // candidate ratio threshold
    int bad_streak = 0;

    std::vector<NegativeCycle> batch_nc;
    size_t used_groups = 0;

    while (true) {
        const auto Th0 = clock::now();
        bool got = stream.next(batch_nc);
        T_harvest_ms += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - Th0).count();
        if (!got) break;
        const size_t batch_idx = used_groups++;

        const auto Tg0 = clock::now();
        double sum_cycle_len = 0.0; size_t metas_before = metas.size();
        for (const auto& cycle : batch_nc) {
            std::vector<Edge> all_edges = { cycle.neg_edge() };
            const auto& pos = cycle.pos_edges();
            all_edges.insert(all_edges.end(), pos.begin(), pos.end());
            sum_cycle_len += 1.0 + (double)pos.size();
            auto batch_metas = generate_cycle_cut_metas_free(owner, owner.edge_index, owner.signs, all_edges, &x_vals, &tri_counts);
            for (auto& cm : batch_metas) {
                uint64_t sig = sig_from_meta(cm);
                if (seen.insert(sig).second) { metas.emplace_back(std::move(cm)); meta_batch.push_back((int)batch_idx); }
            }
        }
        T_gen_ms += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - Tg0).count();

        const auto Te0 = clock::now();
        double slack_sum = 0.0; double density_sum = 0.0; int batch_added = 0; int batch_cand = 0;
        for (size_t i = metas_before; i < metas.size(); ++i) {
            const CutMeta& cm = metas[i];
            double lhs = 0.0;
            for (size_t j = 0; j < cm.y_idx.size(); ++j) {
                const Edge& e = cm.edges[j];
                const int s = cm.eff_sign[j];
                lhs += s * ( y_vals[ cm.y_idx[j] ] - 0.5 * x_vals[e.first] - 0.5 * x_vals[e.second] );
            }
            double viol = lhs - (double)cm.rhs;
            slacks.push_back(viol);
            // (B) density: normalize by cycle length to de-bias long cycles
            int clen = (int)cm.y_idx.size();
            double dens = viol / std::max(2, clen);
            densities.push_back(dens);
            density_sum += dens;
            ++batch_added; if (viol > 1e-4) { candidates.push_back((int)i); ++batch_cand; }
            slack_sum += viol;
        }
        T_eval_ms += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - Te0).count();

        double avg_len = batch_nc.empty() ? 0.0 : (sum_cycle_len / (double)batch_nc.size());
        std::cout << "[BATCH] idx=" << batch_idx
                  << ", cycles=" << batch_nc.size()
                  << ", avgLen=" << std::fixed << std::setprecision(2) << avg_len
                  << ", cuts_accum=" << metas.size()
                  << std::endl;

        double meanS = (batch_added ? (slack_sum / (double)batch_added) : 0.0);
        double meanD = (batch_added ? (density_sum / (double)batch_added) : 0.0);
        double r = (batch_added ? ((double)batch_cand / (double)batch_added) : 0.0);
        std::cout << "[BATCH-QUAL] idx=" << batch_idx
                  << ", cuts=" << batch_added
                  << ", meanSlack=" << std::fixed << std::setprecision(4) << meanS
                  << ", meanDens=" << std::fixed << std::setprecision(4) << meanD
                  << ", candRatio=" << std::fixed << std::setprecision(3) << r
                  << std::endl;
        if (meanS < QUAL_EPS && r < QUAL_TAU) ++bad_streak; else bad_streak = 0;
        if (bad_streak >= 2) break;
    }

    std::cout << "[GROUPS] used=" << used_groups << "/?" << ", raw_cuts=" << metas.size() << std::endl;

    // Branching hints
    std::vector<int> y_priority(owner.y.size(), 0);
    for (int idx : candidates) for (int e : metas[idx].y_idx) if (e >= 0 && e < (int)y_priority.size()) ++y_priority[e];

    int tri_max = 0; for (int t : tri_counts) tri_max = std::max(tri_max, t);
    std::vector<std::pair<double,int>> x_rank; x_rank.reserve(owner.x.size());
    for (int u = 0; u < (int)owner.x.size(); ++u) {
        double frac_sc = 1.0 - std::fabs(x_vals[u] - 0.5);
        double tri_sc  = (tri_max > 0 ? (double)tri_counts[u] / tri_max : 0.0);
        double sc = 0.6 * frac_sc + 0.4 * tri_sc; x_rank.emplace_back(sc, u);
    }
    std::sort(x_rank.begin(), x_rank.end(), [](const auto& A, const auto& B){ return A.first > B.first; });

    std::vector<std::pair<int,int>> y_rank; y_rank.reserve(owner.y.size());
    for (int e = 0; e < (int)owner.y.size(); ++e) y_rank.emplace_back(y_priority[e], e);
    std::sort(y_rank.begin(), y_rank.end(), [](const auto& A, const auto& B){ return A.first > B.first; });

    int topk = std::min(10, (int)owner.y.size());
    std::cout << "[BRANCH-HINT] top-x (frac⨁tri):";
    for (int i = 0; i < std::min(10, (int)x_rank.size()); ++i)
        std::cout << ' ' << x_rank[i].second << '(' << std::setprecision(3) << x_rank[i].first << ')';
    std::cout << "\n[BRANCH-HINT] top-y (coverage):";
    for (int i = 0; i < topk; ++i)
        std::cout << ' ' << y_rank[i].second << '(' << y_rank[i].first << ')';
    std::cout << std::endl;

    // Selection & addition phase
    // Compute global candidate density mean
    double cand_density_sum = 0.0;
    for (int id : candidates) cand_density_sum += (id < (int)densities.size() ? densities[id] : 0.0);
    double cand_density_mean = candidates.empty() ? 0.0 : cand_density_sum / (double)candidates.size();

    const double DENSITY_LOW_TAU = 0.15;   // (D) trigger for greedy-top-B
    const double RATIO_LOW_TAU   = 0.15;   // (D) trigger for greedy-top-B

    // Try to fetch iteration counter (ItCnt). Not all contexts support it; best-effort.
    long long iter_now = -1; static long long iter_prev = -1; long long d_iter = -1;
    try { iter_now = getNiterations64(); } catch(...) { iter_now = -1; }
    if (iter_prev >= 0 && iter_now >= 0) d_iter = std::max(0LL, iter_now - iter_prev);
    iter_prev = iter_now;

    // Adaptive budget B (smaller when density is low or iterations explode)
    int B = (int)candidates.size();
    bool use_greedy_topB = (cand_density_mean < DENSITY_LOW_TAU) || ((double)candidates.size() / std::max(1,(int)metas.size()) < RATIO_LOW_TAU);
    if (use_greedy_topB) {
        // Base budget: 12% of candidates, clamped to [400, 4500]
        int baseB = std::clamp((int)std::ceil(0.12 * (double)candidates.size()), owner.graph.vertex_count()/2, 4500);
        // If iteration count jumped a lot and densities are poor, shrink B further
        const long long ITER_EXPLODE = 150000; // heuristic
        if (d_iter >= 0 && d_iter > ITER_EXPLODE && cand_density_mean < 0.08) baseB = std::max(200, baseB / 2);
        B = std::min(baseB, (int)candidates.size());
    }

    // Ranking by benefit = density first, then raw slack (both descending), prefer shorter cycles implicitly via density
    std::vector<int> chosen;
    chosen.reserve(B);
    if (!use_greedy_topB) {
        // Add all candidates
        chosen = candidates;
    } else {
        struct Item { double dens; double viol; int id; };
        std::vector<Item> heap; heap.reserve(candidates.size());
        for (int id : candidates) {
            double dens = (id < (int)densities.size() ? densities[id] : 0.0);
            double vio  = (id < (int)slacks.size()    ? slacks[id]    : 0.0);
            heap.push_back({dens, vio, id});
        }
        std::partial_sort(heap.begin(), heap.begin() + B, heap.end(), [](const Item& a, const Item& b){
            if (a.dens != b.dens) return a.dens > b.dens; // higher density first
            return a.viol > b.viol;                        // tie-break by larger violation
        });
        for (int i = 0; i < B; ++i) chosen.push_back(heap[i].id);
        std::sort(chosen.begin(), chosen.end()); // ensure stable order for reproducibility
    }

    long int nc = 0; double slack_sum = 0.0; int unadded_cuts = 0; double unadded_slack_sum = 0.0;
    const auto Ta0 = clock::now();
    {
        // Use a set for fast membership during materialization
        std::vector<char> mark(metas.size(), 0);
        for (int id : chosen) if (id >= 0 && id < (int)mark.size()) mark[id] = 1;
        for (int i = 0; i < (int)metas.size(); ++i) {
            const double s = (i < (int)slacks.size() ? slacks[i] : 0.0);
            if (i < (int)mark.size() && mark[i]) {
                IloRange rng = materialize_from_meta(env, metas[i], owner.y, owner.x);
                add(rng);
                slack_sum += s; ++nc;
                if (metas[i].reversed) ++owner.reversed_cycle_cuts_cutgen; else ++owner.standard_cycle_cuts_cutgen;
            } else { unadded_slack_sum += s; ++unadded_cuts; }
        }
    }
    T_add_ms += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - Ta0).count();

    const long long T_total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - TALL0).count();
	try { rows_after = owner.cplex.getNrows(); } catch(...) { rows_after = -1; }
	try { user_after = owner.cplex.getNcuts(IloCplex::CutUser); } catch(...) { user_after = -1; }
	
	// Compute deltas safely
	auto safe_delta = [](IloInt a, IloInt b)->IloInt {
	    return (a >= 0 && b >= 0) ? (a - b) : -1;
	};
	const IloInt rows_delta = safe_delta(rows_after, rows_before);
	const IloInt user_delta = safe_delta(user_after, user_before);

    std::cout << "[SCAN] total=" << metas.size()
              << ", kept_candidates=" << candidates.size()
              << ", ratio=" << std::fixed << std::setprecision(3) << (metas.empty()?0.0:(double)candidates.size()/metas.size())
              << ", time_ms=" << T_eval_ms
              << (IS_ROOT ? ", root=1" : ", root=0")
              << ", meanDens*=" << std::setprecision(4) << cand_density_mean
              << ", ItΔ=" << d_iter
              << ", add_policy=" << (use_greedy_topB?"topB":"all")
              << ", B=" << B
              << std::endl;

    std::cout << "[PHASE] harvest=" << T_harvest_ms
              << "ms gen=" << T_gen_ms
              << "ms eval=" << T_eval_ms
              << "ms select=" << T_select_ms
              << "ms add=" << T_add_ms
              << "ms total=" << T_total_ms << "ms"
              << " rows_b=" << rows_before << " rows_a=" << rows_after
              << " rows_added=" << (rows_after>=0 && rows_before>=0 ? (rows_after-rows_before) : -1)
	          << " user_added=" << (user_delta >= 0 ? user_delta : (IloInt)nc)   // fallback to attempts
              << std::endl;

    std::cout << "[SELECT] cand=" << candidates.size()
              << ", kept_all_in_batches=" << chosen.size()
              << ", batches_used=" << used_groups
              << std::endl;

    total_slack += slack_sum; total_cuts += (int)nc;
    total_unadded_slack += unadded_slack_sum; total_unadded_cuts += unadded_cuts;

    if (IS_ROOT) {
        double curr_lb = getObjValue();
        double d = 0.0; if (root_lb_prev_set) d = curr_lb - root_lb_prev;
        root_lb_prev = curr_lb; root_lb_prev_set = true;
        root_deltas[root_pos] = d;
        root_ratios[root_pos] = (metas.empty()?0.0:(double)candidates.size()/metas.size());
        root_pos = (root_pos + 1) % ROOT_K;
        if (root_count < ROOT_K) ++root_count;
        std::cout << "[ROOT-TRACK] dLB=" << d << ", ratio=" << root_ratios[(root_pos+ROOT_K-1)%ROOT_K]
                  << ", count=" << root_count << "/" << ROOT_K << std::endl;
    }

    if (invocation_count == 5) {
        double mean_slack_added    = (total_cuts > 0) ? (total_slack / total_cuts) : 0.0;
        double mean_slack_unadded  = (total_unadded_cuts > 0) ? (total_unadded_slack / total_unadded_cuts) : 0.0;
        static bool   prev_lb_set = false; static double prev_lb = 0.0;
        const double curr_lb = getObjValue();
        double delta_lb = 0.0; if (prev_lb_set) delta_lb = curr_lb - prev_lb; prev_lb = curr_lb; prev_lb_set = true;

        if (!printed_header) {
            std::cout << std::setw(15) << "Invocation"
                      << std::setw(15) << "Total Cuts"
                      << std::setw(20) << "Mean Slack"
                      << std::setw(20) << "Unadded Cuts"
                      << std::setw(30) << "Mean Slack (Unadded)"
                      << std::setw(25) << "Lower Bound"
                      << std::setw(12) << "ΔLB"
                      << std::setw(12) << "Batches"
                      << std::setw(15) << "Budget"
                      << std::setw(10) << "%Keep"
                      << std::setw(12) << "Time(ms)"
                      << std::endl;
            printed_header = true;
        }
        int adaptive_budget = (int)chosen.size();
        std::cout << std::setw(15) << invocation_count
                  << std::setw(15) << total_cuts
                  << std::setw(20) << mean_slack_added
                  << std::setw(20) << total_unadded_cuts
                  << std::setw(30) << mean_slack_unadded
                  << std::setw(25) << curr_lb
                  << std::setw(12) << delta_lb
                  << std::setw(12) << used_groups
                  << std::setw(15) << adaptive_budget
                  << std::setw(10) << (metas.empty()?0:(100* (int)chosen.size() / (int)metas.size()))
                  << std::setw(12) << 0
                  << std::endl;

        invocation_count = 0;
        total_slack = 0.0; total_cuts = 0;
        total_unadded_slack = 0.0; total_unadded_cuts = 0;
    }
}

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

    // If fractional returned the trivial assignment (=> all x_i = 0), fall back immediately
    {
        const std::vector<int>& s_test = *s_ptr;
        const bool all_pos = std::all_of(s_test.begin(), s_test.end(),
                                         [](int si){ return si == +1; });
        if (all_pos) {
            auto s_fb = owner.graph.greedy_switching(); // base heuristic on original graph
            if (!s_fb.empty()) {
                s_ptr = std::make_shared<const std::vector<int>>(std::move(s_fb));
                std::cout << "[Heuristic] fractional was trivial; using base greedy instead.\n";
            }
        }
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

        // WeightsView supports operator[], not at()
        double w = owner.weights[edge];

        double y0 = (w > 0.0)
                  ? (double)std::min(xu, xv)              // y ≤ x_u, y ≤ x_v
                  : (double)std::max(0, xu + xv - 1);     // y ≥ x_u + x_v − 1

        xy_vars.add(owner.y[idxe]);
        xy_vals.add(y0);

        obj_val += -2.0 * w * y0;
    }

    double incumbent = hasIncumbent() ? getIncumbentObjValue() : IloInfinity;
    std::cout << "Incumbent: " << incumbent << "\n";
    std::cout << "Objective value: " << obj_val << "\n";

    if (obj_val < incumbent) {
        setSolution(xy_vars, xy_vals);
        std::cout << "[Heuristic] Proposed a solution with objective: " << obj_val << "\n";
        owner.f_index = obj_val / owner.graph.edge_count();
        owner.injected_heuristic_solutions++;
    }
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

