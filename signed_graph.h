#pragma once

#define IGRAPH_ENABLE_LGPL
#include <igraph/igraph.h>
#include <igraph/igraph_vector.h>
#include <igraph/igraph_attributes.h>
#include <vector>
#include <string>
#include <utility>
#include <unordered_set>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <queue>
#include <set>
#include <unordered_map>
#include <functional>
#include <memory>
#include <cmath>

struct Edge {
    int first, second;
    Edge() : first(0), second(0) {}
    Edge(int a, int b) : first(std::min(a, b)), second(std::max(a, b)) {}

    bool operator==(const Edge& other) const {
        return first == other.first && second == other.second;
    }

    bool is_adjacent_to(const Edge& other) const {
        return first == other.first || first == other.second ||
               second == other.first || second == other.second;
    }
};

std::ostream& operator<<(std::ostream& os, const Edge& e);

struct SignedEdge {
    Edge points;
    int sign;

    bool operator==(const SignedEdge& other) const {
        return (points.first == other.points.first && points.second == other.points.second) ||
               (points.first == other.points.second && points.second == other.points.first);
    }
};

struct EdgeHash {
    size_t operator()(const Edge& edge) const {
        int a = std::min(edge.first, edge.second);
        int b = std::max(edge.first, edge.second);
        return std::hash<int>()(a) ^ std::hash<int>()(b << 1);
    }
};

std::ostream& operator<<(std::ostream& os, const SignedEdge& e);

class NegativeCycle {
public:
    NegativeCycle(const Edge& neg_edge, std::vector<Edge>&& pos_edges);

    const Edge& neg_edge() const;
    const std::vector<Edge>& pos_edges() const;

private:
    Edge neg_edge_;
    std::vector<Edge> pos_edges_;
};

std::ostream& operator<<(std::ostream& os, const NegativeCycle& nc);

class SignedGraph {
    friend class SignedGraphForMIP;
protected:
    igraph_t g;
    std::vector<double> weights;          // original ±1 signs as doubles (+1.0, -1.0)
    std::vector<double> switched_weights; // current switched signs (±1.0, possibly ±0.0 if ever used)

    SignedGraph(const SignedGraph* const other);
    SignedGraph(const SignedGraph* const other, std::vector<double> new_weights);

public:
    struct GreedyKickOptions {
        int    neg_edge_threshold_abs = -1;
        double neg_edge_threshold_frac = -1;
        int    max_kicks = 0;
        bool   use_weighted_degree = true;
        bool   use_triangle_tiebreak = false;
        double triangle_beta = 0.05;
        int    neighbor_cap = 1024;
        int    triangle_cap_per_u = 2048;
        const std::vector<double>* edge_salience = nullptr;
        double kick_salience_bias = 0.5;
        bool   relax_to_all_pos_if_Z0_empty = true;
        int    delta_m_minus_cap = 512;
        double delta_m_minus_penalty = 0.0;

        // New heuristic knobs
        int R_max = 2;     // round cap
        int K_max = 1;     // max consecutive zero-clique flips per round
        int L_max = 2;     // max integer zero-clique flips in the detour
        int Delta = 0;     // gate on current working signature (deg(u) <= Delta)

        // Fractional inputs (optional). If null, integer-only mode.
        const std::vector<double>* frac_x = nullptr;
        const std::vector<double>* frac_y = nullptr;
    };

protected:
    std::vector<int> greedy_switching_base(const std::function<int(int, int)>& cmp_fn,
                                           const GreedyKickOptions& opts);

private:
    using MapType = std::unordered_map<Edge, int, EdgeHash>;
    std::vector<double> d_plus;
    std::vector<double> d_minus;
    std::vector<double> switched_d_plus;
    std::vector<double> switched_d_minus;
    MapType _edge_index;

    void compute_degrees();
    static bool strictly_positive(double w);
    static double vertex_salience(double x);
    static double edge_salience_from_y(double y);
    static double tau_from_x(double xu, double xv);

    std::vector<int> maximal_salience_clique_strict_pos(
        const std::vector<int>& z0_vertices,
        const std::vector<char>& inZ0,
        const std::vector<double>& switched_frac_edge_sign,
        const std::vector<double>* frac_x) const;

    int pick_u_star_in_Q(
        const std::vector<int>& Q,
        const std::vector<double>& switched_frac_edge_sign,
        const std::vector<double>* frac_y) const;

    // NOTE: refactored — base_edge_sign is the current working signature (±1 per edge)
    void integer_projection_from_x(
        const std::vector<double>& x,
        const std::vector<double>& base_edge_sign,
        std::vector<int>& s_int,
        std::vector<double>& sigma_int_edge_sign) const;

    void integer_greedy_pass_minheap(
        std::vector<double>& sigma_int_edge_sign,
        std::vector<int>& s_int,
        std::vector<int>& d_int) const;

    void single_switching(int u, igraph_vector_int_t* incident);

public:
    std::vector<int> negative_triangle_count_per_vertex() const;
    int negative_triangle_count_of_vertex(int u) const;

    class SignedEdgesView {
    private:
        const igraph_t* g;
        const std::vector<double>* weights;

    public:
        explicit SignedEdgesView(const igraph_t* graph, const std::vector<double>* weights)
                : g(graph), weights(weights) {}

        class const_iterator {
        private:
            const igraph_t* g;
            igraph_integer_t eid;
            const std::vector<double>* weights;

        public:
            using value_type = std::pair<Edge, int>;

            const_iterator(const igraph_t* g, igraph_integer_t eid, const std::vector<double>* weights)
                : g(g), eid(eid), weights(weights) {}

            value_type operator*() const {
                igraph_integer_t from, to;
                igraph_edge(g, eid, &from, &to);
                double w = (*weights)[eid];
                int sign = (w == 0.0 && std::signbit(w)) ? -1 : (w >= 0.0 ? 1 : -1);
                return {{static_cast<int>(from), static_cast<int>(to)}, sign};
            }

            const_iterator& operator++() {
                ++eid;
                return *this;
            }

            bool operator!=(const const_iterator& other) const {
                return eid != other.eid;
            }
        };

        const_iterator begin() const { return const_iterator(g, 0, weights); }
        const_iterator end() const { return const_iterator(g, igraph_ecount(g), weights); }

        SignedEdge operator[](igraph_integer_t eid) const;
        SignedEdge operator[](const Edge& e) const;
        int size() const { return igraph_ecount(g); }
    };

    class WeightsView {
    private:
        const igraph_t* g;
        const std::vector<double>* weights;

    public:
        explicit WeightsView(const igraph_t* graph, const std::vector<double>* weights)
            : g(graph), weights(weights) {}

        class const_iterator {
        private:
            const igraph_t* g;
            igraph_integer_t eid;
            const std::vector<double>* weights;

        public:
            using value_type = std::pair<Edge, double>;

            const_iterator(const igraph_t* g, igraph_integer_t eid, const std::vector<double>* weights)
                : g(g), eid(eid), weights(weights) {}

            value_type operator*() const {
                igraph_integer_t from, to;
                igraph_edge(g, eid, &from, &to);
                return {{static_cast<int>(from), static_cast<int>(to)}, (*weights)[eid]};
            }

            const_iterator& operator++() {
                ++eid;
                return *this;
            }

            bool operator!=(const const_iterator& other) const {
                return eid != other.eid;
            }
        };

        const_iterator begin() const { return const_iterator(g, 0, weights); }
        const_iterator end() const { return const_iterator(g, igraph_ecount(g), weights); }

        double operator[](const Edge& e) const;
        int size() const { return igraph_ecount(g); }
    };

    class EdgeIndexesView {
    private:
        const igraph_t* g;
        const MapType& map;

    public:
        explicit EdgeIndexesView(const MapType& m, const igraph_t* graph);

        class const_iterator {
        private:
            const igraph_t* g;
            igraph_integer_t eid;

        public:
            using value_type = std::pair<Edge, int>;

            const_iterator(const igraph_t* g, igraph_integer_t eid)
                : g(g), eid(eid) {}

            value_type operator*() const;
            const_iterator& operator++();
            bool operator!=(const const_iterator& other) const;
        };

        const_iterator begin() const;
        const_iterator end() const;

        int operator[](const Edge& e) const;
    };

    SignedGraph(const std::string& file_path);
    std::unique_ptr<SignedGraph> clone() const;
    virtual ~SignedGraph();

    const SignedEdgesView signs_view() const { return SignedEdgesView(&g, &weights); }
    const WeightsView weights_view() const { return WeightsView(&g, &weights); }
    const EdgeIndexesView edge_index() const;
    // Original & switched edge weights (eid-indexed)
    const std::vector<double>& get_weights() const { return weights; }
    const std::vector<double>& get_switched_weight() const { return switched_weights; }
    // Convenience: key/value iterator view over the current switched signature
    const WeightsView switched_weights_view() const { return WeightsView(&g, &switched_weights); }
    const std::vector<double>& plus_degrees_view() const;
    const std::vector<double>& minus_degrees_view() const;
    const long double max_degree_vertex() const;
    const std::vector<int> neighbors(int u) const;
    int frustrated_edges(const std::vector<int>& partition) const;
    std::vector<Edge> frustrated_edges_keys(const std::vector<int>& partition) const;
    void switching_from_partition(const std::vector<int>& s);
    void single_switching(int u);
    const std::vector<int> greedy_switching();
    // Back-compat helper: uses default (no-kick) options
    std::vector<int> greedy_switching_base(const std::function<int(int,int)>& cmp_fn) {
        return greedy_switching_base(cmp_fn, GreedyKickOptions{});
    }
    void restore_switched_sign();
    bool are_cycles_edge_disjoint(const std::vector<NegativeCycle>& cycles) const;
    bool are_cycles_sign_correct(const std::vector<NegativeCycle>& cycles, bool expect_negative = true) const;

    int vertex_count() const;
    int edge_count() const;
    void print_info() const;
    void save_partition_svg(const std::vector<int>& partition, const std::string& filename, bool custom_layout) const;
};
