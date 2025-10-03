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

// FrustrationModel.h (partial)
struct TriangleInequalities {
    std::vector<Edge> triangle;
    std::vector<int> pending_cut_ids; // -1 for full triangle cut, vertex id for vertex-based cuts
};

class FrustrationModel {
public:
    static constexpr int NO_CUTS = 0;
    static constexpr int NEGATIVE_CYCLE_CUTS = 1;
    static constexpr int NET_DEGREE_CUTS = 2;
    static constexpr int TRIANGLE_CUTS = 4;
    static constexpr int NEGATIVE_CYCLE_COVER_CUTS = 8;

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
    SignedGraphForMIP& graph;
    int m_minus;
    std::vector<TriangleInequalities> uncut_triangles; // Stores pending inequalities for 
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
    SignedGraph::WeightsView weights;
    SignedGraphForMIP::EdgeIndexesView edge_index;
	std::unordered_map<int, Edge> edge_reverse;

	double frustration_index(double obj_val) const;

    void initialize_uncut_triangles();
	virtual std::vector<std::pair<IloRange, std::string>> generate_cycle_cuts(IloEnv& env, const std::vector<Edge>& all_edges) = 0;
    virtual std::vector<std::pair<IloRange, std::string>> generate_pending_triangle_cuts(IloEnv& env, const TriangleInequalities& t) = 0;
    std::vector<std::pair<IloRange, std::string>> generate_negative_cycle_cuts(IloEnv& env);
    std::vector<std::pair<IloRange, std::string>> generate_negative_cycle_cover_cuts(IloEnv& env);
    std::vector<std::pair<IloRange, std::string>> generate_triangle_cuts(IloEnv& env);
    virtual void export_solution(const std::string& file_prefix, bool with_svg, std::vector<int> partition) const;
};
