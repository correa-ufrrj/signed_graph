// File: frustration_model.cpp
#include "frustration_model.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <unordered_set>

// Accessor definitions
ModelAccessor::ModelAccessor(IloCplex& cpx) : cpx_(cpx) {}
void ModelAccessor::getValues(IloNumArray& out, const IloNumVarArray& vars) const {
    cpx_.getValues(out, vars);
}

UserCutCallbackAccessor::UserCutCallbackAccessor(IloCplex::UserCutCallbackI& cb) : cb_(cb) {}
void UserCutCallbackAccessor::getValues(IloNumArray& out, const IloNumVarArray& vars) const {
    cb_.getValues(out, vars);
}

FrustrationModel::FrustrationModel(SignedGraphForMIP& g, int cut_flags)
    : graph(g)
    , model(env)
    , cplex(model)
    , use_cut_generator(cut_flags)
    , lower_bound(IloIntMin)
    , f_index(g.edge_count())
    , edge_index(g.edge_index())
    , signs(g.signs_view())
    , weights(g.weights_view()) {

    m_minus = 0;
    for (const auto& [_, sign] : signs) if (sign < 0) m_minus++;
    for (const auto& [e, idx] : edge_index) {
        edge_reverse[idx] = e;
    }
}

void FrustrationModel::initialize_uncut_triangles() {
    int n = graph.vertex_count();

    for (int u = 0; u < n; ++u) {
        const auto& neighbors_u = graph.neighbors(u);
        std::unordered_set<int> neighbor_set_u(neighbors_u.begin(), neighbors_u.end());

        for (int v : neighbors_u) {
            if (v <= u) continue;

            Edge uv{u, v};
            const auto& neighbors_v = graph.neighbors(v);
            for (int w : neighbors_v) {
                if (w <= v || w == u) continue;

                if (neighbor_set_u.count(w)) {
                    Edge vw{v, w}, uw{u, w};
                    int sign_product = signs[uv].sign * signs[vw].sign * signs[uw].sign;
                    if (sign_product < 0) {
                        TriangleInequalities triangle_entry;
                        triangle_entry.triangle = {uv, vw, uw};
                        triangle_entry.pending_cut_ids = {-1, uv.first, uv.second, vw.second};
                        uncut_triangles.push_back(std::move(triangle_entry));
                    }
                }
            }
        }
    }
}

std::vector<std::pair<IloRange, std::string>> FrustrationModel::generate_triangle_cuts(IloEnv& env) {
    std::vector<std::pair<IloRange, std::string>> triangle_cuts;
    auto it = uncut_triangles.begin();
    while (it != uncut_triangles.end()) {
        auto& triangle = *it;
        auto cuts = generate_pending_triangle_cuts(env, triangle);
        triangle_cuts.insert(triangle_cuts.end(), cuts.begin(), cuts.end());
        ++it;
    }
    return triangle_cuts;
}

std::vector<std::pair<IloRange, std::string>> FrustrationModel::generate_negative_cycle_cuts(IloEnv& env) {
    std::vector<std::pair<IloRange, std::string>> all_cuts;
    auto cycles = graph.find_switched_lower_bound();

    for (const auto& cycle : cycles) {
        // If triangle cuts are enabled, skip 2-pos cycles here to avoid duplicates.
        if (use_cut_generator & TRIANGLE_CUTS)
            if (cycle.pos_edges().size() == 2) continue;

        std::vector<Edge> all_edges = { cycle.neg_edge() };
        const auto& pos = cycle.pos_edges();
        all_edges.insert(all_edges.end(), pos.begin(), pos.end());

        auto cuts = generate_cycle_cuts(env, all_edges);
        all_cuts.insert(all_cuts.end(), cuts.begin(), cuts.end());
    }
    return all_cuts;
}

double FrustrationModel::frustration_index(double obj_val) const {
    return (obj_val + m_minus) / graph.edge_count();
}

double FrustrationModel::get_frustration_index() const {
    if (!cplex.isExtracted(model)) {
        throw std::runtime_error("Model not extracted; cannot retrieve objective value.");
    }
    IloNum obj_val = cplex.getObjValue();
    return frustration_index(obj_val);
}

void FrustrationModel::solve() {
    start_time = std::chrono::steady_clock::now();

    // Turn off most CPLEX cuts/heuristics to focus on our separation
    cplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);
    cplex.setParam(IloCplex::Param::MIP::Cuts::Cliques, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::ZeroHalfCut, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::LiftProj, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::MIRCut, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::GUBCovers, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, -1);

    cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);
    cplex.setParam(IloCplex::Param::MIP::Strategy::Probe, -1);
    cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur, -1);
    cplex.setParam(IloCplex::Param::MIP::Strategy::FPHeur, -1);

    cplex.setParam(IloCplex::Param::TimeLimit, 12600);

    if (!cplex.solve()) {
        std::cerr << "[Solver] Optimization failed." << std::endl;
    } else {
        f_index = get_frustration_index();
        std::cout << "[Solver] Optimization succeeded. Objective: " << f_index << std::endl;
    }

    end_time = std::chrono::steady_clock::now();
}

std::string FrustrationModel::active_cut_names() const {
    std::ostringstream oss;
    bool first = true;

    // When both flags are on, advertise the combined pipeline name
    if (pipeline_enabled()) {
        oss << "TriangleCyclePipeline";
        first = false;
    } else {
        if (use_cut_generator & NEGATIVE_CYCLE_CUTS) {
            oss << (first ? "" : ";") << "NegativeCycles";
            first = false;
        }
        if (use_cut_generator & TRIANGLE_CUTS) {
            oss << (first ? "" : ";") << "Triangles";
            first = false;
        }
    }

    if (use_cut_generator & NET_DEGREE_CUTS) {
        oss << (first ? "" : ";") << "NetDegree";
    }
    return oss.str();
}

void FrustrationModel::print_solution() const {
    std::cout << "Status: " << cplex.getStatus() << "\n";
    std::cout << "Lower bound (cycles): " << lower_bound << std::endl;
    std::cout << "Objective value: " << cplex.getObjValue() << std::endl;
    std::cout << "Frustration index: " << f_index << std::endl;
}

void FrustrationModel::export_solution(const std::string& file_prefix, bool with_svg, std::vector<int> partition) const {
    std::ofstream meta(file_prefix + "_summary.csv");

    int frustrated = graph.frustrated_edges(partition);
    auto neg_degrees = graph.minus_degrees_view();
    int num_neg_edges = std::accumulate(neg_degrees.begin(), neg_degrees.end(), 0) / 2;
    double runtime = std::chrono::duration<double>(end_time - start_time).count();

    meta << "objective,status,lower_bound,num_nodes,num_edges,num_neg_edges,use_cuts,runtime_sec,"
         << "frustrated_edges,net_degree_cuts,neg_cycle_build_std,neg_cycle_build_rev,"
         << "neg_cycle_cut_std,neg_cycle_cut_rev,nodes_explored,num_injected_heuristic\n";

    meta << std::fixed << std::setprecision(4);
    meta << f_index << "," << cplex.getStatus() << "," << lower_bound << "," << graph.vertex_count() << "," << graph.edge_count()
         << "," << num_neg_edges << "," << active_cut_names() << "," << runtime << "," << frustrated
         << "," << net_degree_cut_count
         << "," << standard_cycle_cuts_build
         << "," << reversed_cycle_cuts_build
         << "," << standard_cycle_cuts_cutgen
         << "," << reversed_cycle_cuts_cutgen
         << "," << cplex.getNnodes() << "," << injected_heuristic_solutions
         << "\n";

    meta.close();

    if (with_svg) {
        graph.save_partition_svg(partition, file_prefix + "_partition.svg", true);
    }
}
