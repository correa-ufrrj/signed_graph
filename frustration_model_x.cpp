// File: frustration_model_x.cpp

#include "frustration_model_x.h"
#include <algorithm>
#include <stdexcept>

FrustrationModelX::FrustrationModelX(SignedGraphForMIP& g, int cut_flags)
    : FrustrationModel(g, cut_flags) {}

void FrustrationModelX::build() {
    int n = graph.vertex_count();
    x.resize(n);

    for (int i = 0; i < n; ++i)
        x[i] = IloBoolVar(env, ("x_" + std::to_string(i)).c_str());
    z = IloNumVar(env, "z");

    const auto& d_plus = graph.plus_degrees_view();
    const auto& d_minus = graph.minus_degrees_view();

    IloExpr obj_expr(env);
    obj_expr += z;

    objective = IloMinimize(env, obj_expr);
    model.add(objective);
    obj_expr.end();

    // Initial solution using greedy switching
    graph.restore_switched_sign();
    std::vector<int> s = graph.greedy_switching();

    // Improving inequalities, if required
    if (use_cut_generator != NO_CUTS) {
        long int nc = 0;
        // Net degree inequalities
//        if (use_cut_generator & NET_DEGREE_CUTS) {
//            for (int u = 0; u < graph.vertex_count(); ++u) {
//                int d_u = d_plus[u] - d_minus[u];
//                if(abs(d_u) % 2 == 0) continue;
//                IloExpr constraint(env);
//                constraint += d_u * x[u];
//
//                for (int v : graph.neighbors(u)) {
//                    Edge uv = {v, u};
//                    SignedEdge e = signs[uv];
//                    int sign = e.sign;
//                    constraint += sign * x[v];
//                    int index = edge_index[uv];
//                    constraint -= 2 * sign * y[index];
//                }
//               	model.add(constraint <= std::floor(d_u / 2.0));
//                constraint.end();
//                net_degree_cut_count++;
//            }
//        }
//
//	    if (use_cut_generator & NEGATIVE_CYCLE_CUTS) {
//	        auto cuts = generate_negative_cycle_cuts(env);
//    	    for (auto& [cut_expr, type] : cuts) {
//	            model.add(cut_expr);
//    	        if (type == "standard") {
//					++standard_cycle_cuts_build;
//	    	        nc++;
//	    	    }
//    	        else if (type == "reversed") ++reversed_cycle_cuts_build;
////	            cut_expr.end();
//	        }
//	    }

	    if (use_cut_generator & TRIANGLE_CUTS) {
	        auto cuts = generate_triangle_cuts(env);
    	    for (auto& [cut_expr, type] : cuts) {
    	    	std::cout << "Adding triangle cut: " << cut_expr << std::endl;
	            model.add(cut_expr);
    	        if (type == "standard") {
					++standard_cycle_cuts_build;
	    	    }
    	        else if (type == "reversed") ++reversed_cycle_cuts_build;
//	            cut_expr.end();
	        }
	    }

//	    if (use_cut_generator & NEGATIVE_CYCLE_COVER_CUTS) {
//	        auto cuts = generate_negative_cycle_cover_cuts(env);
//    	    for (auto& [cut_expr, type] : cuts) {
//	            model.add(cut_expr);
//    	        if (type == "standard") {
//					++standard_cycle_cuts_build;
//	    	    }
//    	        else if (type == "reversed") ++reversed_cycle_cuts_build;
////	            cut_expr.end();
//	        }
//	    }
        lower_bound = std::max(lower_bound, nc);
	}

    IloNumArray start_vals(env);
    IloNumVarArray vars(env);
    for (int i = 0; i < n; ++i) {
        model.add(x[i]);
        vars.add(x[i]);
        start_vals.add(static_cast<IloNum>(s[i]));
    }

    auto fe = graph.frustrated_edges(s);
    std::cout << "fe = " << fe << std::endl;
    model.add(z);
    vars.add(z);
    start_vals.add(static_cast<IloNum>(fe));

    auto [ordered, undecided] = generate_edge_partition_ineq(start_vals);
    auto lifted = lift_edge_partition_cuts(undecided, start_vals);

    ordered.insert(ordered.end(), lifted.begin(), lifted.end());
    IloRange range = build_edge_partition_cut(ordered);
	model.add(range);

    cplex.addMIPStart(vars, start_vals);
    injected_heuristic_solutions++;

    igraph_integer_t max_vertex = graph.max_degree_vertex();
    model.add(x[max_vertex] == s[max_vertex]);
}

std::vector<std::pair<IloRange, std::string>> FrustrationModelX::generate_pending_triangle_cuts(IloEnv& env, const TriangleInequalities& t) {
    std::vector<std::pair<IloRange, std::string>> result;
	return result;
}

std::vector<std::pair<IloRange, std::string>> FrustrationModelX::generate_cycle_cuts(IloEnv& env, const std::vector<Edge>& all_edges) {
//    std::set<int> vertices_set;
//    for (const auto& edge : all_edges) {
//        vertices_set.insert(edge.first);
//        vertices_set.insert(edge.second);
//    }
//
//    std::vector<int> vertices(vertices_set.begin(), vertices_set.end());
    std::vector<std::pair<IloRange, std::string>> result;

//    int total_combinations = 1 << vertices.size();
//
//    for (int mask = 0; mask < total_combinations; ++mask) {
//        IloNumArray x_vals(env, x.size());
//        for (std::size_t i = 0; i < vertices.size(); ++i) {
//            x_vals[vertices[i]] = (mask >> i) & 1;
//        }
//
//    	auto cuts = generate_cycle_ineq(all_edges, x_vals);
//        lift_edge_partition_cuts(cuts, x_vals);
//        for (auto& [range, _] : cuts) {
//            result.emplace_back(range, "standard");
//        }
//    }

    return result;
}

inline std::vector<std::pair<IloRange, std::vector<std::pair<int, int>>>> FrustrationModelX::generate_cycle_ineq(const std::vector<Edge>& all_edges, const IloNumArray& x_vals) {
    std::vector<std::pair<IloRange, std::vector<std::pair<int, int>>>> result;
    IloExpr sum_expr(env);
    std::vector<std::pair<int, int>> excluded_edges;

    double rhs = 1.0;
    for (const auto& edge : all_edges) {
        int sign = signs[edge].sign;
        int u = edge.first;
        int v = edge.second;
        double xu = x_vals[u];
        double xv = x_vals[v];

        if (sign > 0 && std::abs(xu - xv) > 1e-6) {
            sum_expr += (xu > xv) ? (x[u] - x[v]) : (x[v] - x[u]);
        } else if (sign < 0) {
        	rhs--;
			double val = xu + xv - 1.0;
			if (std::abs(val) > 1e-6) {
				sum_expr += (val > 0) ? (x[u] + x[v] - 1.0) : (1.0 - x[u] - x[v]);
			} else {
	            excluded_edges.emplace_back(u, v);
	        }
        } else {
            excluded_edges.emplace_back(u, v);
        }
    }

    result.emplace_back(sum_expr >= 1.0, excluded_edges);
    sum_expr.end();
    return result;
}

bool FrustrationModelX::vertex_partition_compare(int u, int v, const IloNumArray& x_vals) const {
    auto count_matching_neighbors = [&](int node) {
        int count = 0;
        for (int w : graph.neighbors(node)) {
            Edge edge = {node, w};
            int sign = signs[edge].sign;
            double x_node = x_vals[node];
            double x_w = x_vals[w];
            if ((sign > 0 && std::abs(x_node - x_w) < 1e-6) || (sign < 0 && std::abs(x_node + x_w - 1.0) < 1e-6)) {
                count++;
            } else {
            	count--;
            }
        }
        return count;
    };

    int du = count_matching_neighbors(u);
    int dv = count_matching_neighbors(v);
    return du > dv || (du == dv && u < v);
};

// FrustrationModelX::generate_edge_partition_ineq returns ordered and undecided pairs
inline std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>>
FrustrationModelX::generate_edge_partition_ineq(const IloNumArray& x_vals) {
    std::vector<std::pair<int, int>> ordered_pairs;
    std::vector<std::pair<int, int>> undecided_pairs;

    for (const auto& [edge, sign] : signs) {
        int u = edge.first;
        int v = edge.second;
        double xu = x_vals[u];
        double xv = x_vals[v];

        if (sign > 0 && std::abs(xu - xv) > 1e-6) {
            ordered_pairs.emplace_back((xu > xv) ? std::make_pair(u, v) : std::make_pair(v, u));
        } else if (sign < 0 && std::abs(xu + xv - 1.0) > 1e-6) {
            double val = xu + xv - 1.0;
            ordered_pairs.emplace_back((val > 0) ? std::make_pair(u, v) : std::make_pair(v, u));
        } else {
            undecided_pairs.emplace_back(u, v);
        }
    }

    return {ordered_pairs, undecided_pairs};
}

// lift_edge_partition_cuts uses vertex_partition_compare for undecided pairs
inline std::vector<std::pair<int, int>> FrustrationModelX::lift_edge_partition_cuts(
    const std::vector<std::pair<int, int>>& undecided,
    const IloNumArray& x_vals) {

    std::vector<std::pair<int, int>> lifted;
    for (const auto& [u, v] : undecided) {
        bool cmp = vertex_partition_compare(u, v, x_vals);
        lifted.emplace_back(cmp ? std::make_pair(u, v) : std::make_pair(v, u));
    }

    return lifted;
}

inline IloRange FrustrationModelX::build_edge_partition_cut(
    const std::vector<std::pair<int, int>>& ordered_pairs) {

    IloExpr expr(env);
    double lb = 0.0;

    for (const auto& [u, v] : ordered_pairs) {
        int sign = signs[{u, v}].sign;
        if (sign > 0) {
            expr += x[u] - x[v];
        } else {
            expr += x[u] + x[v];
            lb += 1.0;
        }
    }

    IloRange cut = (z - expr >= lb);
    expr.end();
    return cut;
}

inline void FrustrationModelX::EdgePartitionCutGenerator::main() {
    IloNumArray x_vals(getEnv());
    for (int i = 0; i < owner.x.size(); ++i)
        x_vals.add(static_cast<IloNum>(getValue(owner.x[i])));

    auto [ordered, undecided] = owner.generate_edge_partition_ineq(x_vals);
    auto lifted = owner.lift_edge_partition_cuts(undecided, x_vals);

    ordered.insert(ordered.end(), lifted.begin(), lifted.end());
    IloRange cut = owner.build_edge_partition_cut(ordered);

    IloNum slack = getValue(cut.getExpr()) - cut.getLB();
    if (slack < -5e-2) {
        add(cut);
    }
}

FrustrationModelX::EdgePartitionCutGenerator::EdgePartitionCutGenerator(IloEnv env, FrustrationModelX& owner)
	: IloCplex::UserCutCallbackI(env), owner(owner) {}

FrustrationModelX::EdgePartitionCutGenerator::~EdgePartitionCutGenerator() {}

FrustrationModelX::NegativeCycleCutGenerator::NegativeCycleCutGenerator(IloEnv env, FrustrationModelX& owner)
    : IloCplex::UserCutCallbackI(env), owner(owner) {}

FrustrationModelX::NegativeCycleCutGenerator::~NegativeCycleCutGenerator() {}

void FrustrationModelX::NegativeCycleCutGenerator::main() {
    // Placeholder: implement custom x-only user cut generation
}

FrustrationModelX::LazyEdgePartitionSeparator::LazyEdgePartitionSeparator(IloEnv env, FrustrationModelX& owner)
    : IloCplex::LazyConstraintCallbackI(env), owner(owner) {}

FrustrationModelX::LazyEdgePartitionSeparator::~LazyEdgePartitionSeparator() {}

void FrustrationModelX::LazyEdgePartitionSeparator::main() {
    IloNumArray x_vals(getEnv());
    for (int i = 0; i < owner.x.size(); ++i)
        x_vals.add(static_cast<IloNum>(getValue(owner.x[i])));

    auto [ordered, undecided] = owner.generate_edge_partition_ineq(x_vals);
    auto lifted = owner.lift_edge_partition_cuts(undecided, x_vals);

    ordered.insert(ordered.end(), lifted.begin(), lifted.end());
    IloRange cut = owner.build_edge_partition_cut(ordered);

    IloNum slack = getValue(cut.getExpr()) - cut.getLB();
    if (slack < -5e-2) {
        add(cut);
    }
}

void FrustrationModelX::SwitchingHeuristicCallback::main() {
    // Placeholder: implement heuristic switching solution proposal
}

inline void FrustrationModelX::solve() {
    // EdgePartitionCutGenerator is always added and has priority
    cplex.use(new (env) EdgePartitionCutGenerator(env, *this));

    // NegativeCycleCutGenerator is added only if requested, and has lower priority
    if (use_cut_generator & NEGATIVE_CYCLE_CUTS) {
        cplex.use(new (env) NegativeCycleCutGenerator(env, *this));
    }

    cplex.use(new (env) LazyEdgePartitionSeparator(env, *this));
    cplex.use(new (env) SwitchingHeuristicCallback(env, *this));

    FrustrationModel::solve();
}

void FrustrationModelX::export_solution(const std::string& file_prefix, bool with_svg) const {
    std::ofstream xfile(file_prefix + "_x.csv");
    std::vector<int> partition;
    for (std::size_t i = 0; i < x.size(); ++i) {
        double val = cplex.getValue(x[i]);
        xfile << val << "\n";
        partition.push_back(static_cast<int>(std::round(val)));
    }
    xfile.close();

    FrustrationModel::export_solution(file_prefix, with_svg, partition);
}
