// signed_graph.h
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
#include <algorithm>
#include <edge.h>

struct IGraphDeleter {
  void operator()(igraph_t* p) const noexcept { if (p) { igraph_destroy(p); delete p; } }
};

class GraphCore {
protected:
    const std::shared_ptr<igraph_t> g_;   // shared topology
    // Immutable edge sign pattern (±1 ints)
    std::vector<int>                sigma_;
    // Aggregations (by sign counts)
    std::vector<double>             d_pos;
    std::vector<double>             d_neg;

    using MapType = std::unordered_map<Edge, int, EdgeHash>;
    MapType _edge_index;

    // share topology, copy sigma_/aggregations
    GraphCore(const GraphCore&)            = default;
    GraphCore& operator=(const GraphCore&) = default;
    GraphCore(GraphCore&&)                 = default;
    GraphCore& operator=(GraphCore&&)      = default;
    GraphCore(const GraphCore* const other);
    GraphCore(const std::string& file_path);
    GraphCore(const GraphCore* const other, std::vector<double> new_weights); // validates {±1} → sigma_
    GraphCore(const std::shared_ptr<igraph_t> base, std::vector<double> new_weights);     // validates {±1} → sigma_

    // Sign tests (tie-break +0/-0 preserved via double cast)
    inline bool is_pos_(const double w) const {
        return (w > 0.0) || (w == 0.0 && !std::signbit(w));
    }
    inline bool is_neg_(const double w) const {
        return (w < 0.0) || (w == 0.0 &&  std::signbit(w));
    }

    inline void compute_degrees() {
        const int n = vertex_count();
        d_pos.assign(n, 0.0);
        d_neg.assign(n, 0.0);
        for (int eid = 0; eid < edge_count(); ++eid) {
            igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
            if (is_pos_edge(eid)) {
                d_pos[(int)u] += 1.0; d_pos[(int)v] += 1.0;
            } else {
                d_neg[(int)u] += 1.0; d_neg[(int)v] += 1.0;
            }
        }
    }

    static inline void validate_and_fill_sigma(std::vector<int>& dst, const std::vector<double>& src) {
        dst.clear(); dst.reserve(src.size());
        for (double w : src) {
            if (w != 1.0 && w != -1.0) throw std::invalid_argument("Edge signs must be ±1.");
            dst.push_back(w > 0 ? 1 : -1);
        }
    }

public:
    // Hook for salience used by views; base returns 0.0.
    virtual double edge_salience(int u, int v) const noexcept { return 0.0; }

    virtual ~GraphCore() = default;

    class EdgeIndexesView {
    private:
        const igraph_t* g; const MapType& map;
    public:
        explicit EdgeIndexesView(const MapType& m, const igraph_t* graph) : g(graph), map(m) {}
        class const_iterator {
        private: const igraph_t* g; igraph_integer_t eid; 
        public:
            using value_type = std::pair<Edge, int>;
            const_iterator(const igraph_t* g, igraph_integer_t eid) : g(g), eid(eid) {}
            value_type operator*() const { igraph_integer_t u,v; igraph_edge(g, eid, &u, &v); return {{(int)u,(int)v}, (int)eid}; }
            const_iterator& operator++() { ++eid; return *this; }
            bool operator!=(const const_iterator& other) const { return eid != other.eid; }
        };
        const_iterator begin() const { return const_iterator{g, 0}; }
        const_iterator end()   const { return const_iterator{g, igraph_ecount(g)}; }
        int operator[](const Edge& e) const { auto it = map.find(e); if (it == map.end()) throw std::out_of_range("Edge not found"); return it->second; }
    };

    class SignedEdgesView {
    private:
        const igraph_t* g; const GraphCore* core;
    public:
        explicit SignedEdgesView(const GraphCore* core, const igraph_t* graph) : g(graph), core(core) {}
        class const_iterator {
        private: const igraph_t* g; igraph_integer_t eid; const GraphCore* core; 
        public:
            using value_type = SignedEdge;
            const_iterator(const igraph_t* g, igraph_integer_t eid, const GraphCore* core) : g(g), eid(eid), core(core) {}
            value_type operator*() const {
                igraph_integer_t u,v; igraph_edge(g, eid, &u, &v);
                const int sign = core->is_neg_edge((int)eid) ? -1 : +1;
                const double sal = core->edge_salience((int)u, (int)v);
                return { Edge{ (int)u, (int)v }, sign, sal };
            }
            const_iterator& operator++() { ++eid; return *this; }
            bool operator!=(const const_iterator& other) const { return eid != other.eid; }
        };
        const_iterator begin() const { return const_iterator(g, 0, core); }
        const_iterator end()   const { return const_iterator(g, igraph_ecount(g), core); }
        inline SignedEdge operator[](igraph_integer_t e) const {
            igraph_integer_t u,v; igraph_edge(g, e, &u, &v);
            const int sign = core->is_neg_edge((int)e) ? -1 : +1;
            const double sal = core->edge_salience((int)u, (int)v);
            return { Edge{ (int)u, (int)v }, sign, sal };
        }
        inline SignedEdge operator[](const Edge& e) const {
            igraph_integer_t eid; if (igraph_get_eid(g, &eid, e.first, e.second, 0, 0) != IGRAPH_SUCCESS)
                return { e, 0, 0.0 };
            return (*this)[eid];
        }
        int size() const { return igraph_ecount(g); }
    };

    class NegativeEdgesView {
    private: const igraph_t* g; const GraphCore* core; 
    public:
        explicit NegativeEdgesView(const GraphCore* core, const igraph_t* g_) : g(g_), core(core) {}
        struct const_iterator {
            const igraph_t* g; const GraphCore* core; igraph_integer_t eid; igraph_integer_t m;
            void skip_to_next_neg() { while (eid < m) { if (core->is_neg_edge((int)eid)) break; ++eid; } }
            const_iterator(const igraph_t* g_, const GraphCore* core, igraph_integer_t start) : g(g_), core(core), eid(start), m(igraph_ecount(g_)) { skip_to_next_neg(); }
            bool operator!=(const const_iterator& other) const { return eid != other.eid; }
            const_iterator& operator++() { ++eid; skip_to_next_neg(); return *this; }
            std::pair<Edge,int> operator*() const { igraph_integer_t u,v; igraph_edge(g, eid, &u, &v); return { Edge{ (int)u,(int)v }, (int)eid }; }
        };
        const_iterator begin() const { return const_iterator(g, core, 0); }
        const_iterator end()   const { return const_iterator(g, core, igraph_ecount(g)); }
    };

    inline int vertex_count() const { return igraph_vcount(g_.get()); }
    inline int edge_count()   const { return igraph_ecount(g_.get()); }

    inline const EdgeIndexesView edge_index() const { return EdgeIndexesView{_edge_index, g_.get()}; }
    inline SignedEdgesView signs_view() const { return SignedEdgesView(this, g_.get()); }
    inline NegativeEdgesView negative_edges_view() const { return NegativeEdgesView(this, g_.get()); }

    // Degrees
    inline const std::vector<double>& get_pos_degrees() const { return d_pos; }
    inline const std::vector<double>& get_neg_degrees() const { return d_neg; }
    inline double net_degree(int u) const noexcept { return d_pos[u] - d_neg[u]; }
    inline std::vector<double> net_degrees() const { std::vector<double> d(vertex_count()); for (int u = 0; u < vertex_count(); ++u) d[u] = net_degree(u); return d; }
    inline long double max_pos_degree_vertex() const { int n = vertex_count(); long int mx=0, vid=0; for (int i=0;i<n;++i){ if (d_pos[i] > mx) { mx=d_pos[i]; vid=i; } } return vid; }
    inline long double max_degree_vertex() const {
        int n = vertex_count(); igraph_vector_int_t degrees; igraph_real_t mx=0; long int vid=0;
        igraph_vector_int_init(&degrees, n);
        igraph_degree(g_.get(), &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
        for (igraph_integer_t i = 0; i < n; ++i) { igraph_real_t deg = VECTOR(degrees)[i]; if (deg > mx) { mx = deg; vid = i; } }
        igraph_vector_int_destroy(&degrees); return vid;
    }
    inline std::vector<int> neighbors(int u) const { igraph_vector_int_t vec; igraph_vector_int_init(&vec, 0); igraph_neighbors(g_.get(), &vec, u, IGRAPH_ALL); std::vector<int> res(VECTOR(vec), VECTOR(vec)+igraph_vector_int_size(&vec)); igraph_vector_int_destroy(&vec); return res; }
    inline bool are_disjoint(const std::vector<std::vector<Edge>>& edge_sets) const { std::unordered_set<Edge,EdgeHash> seen; for (const auto& es: edge_sets) for (const auto& e: es) if (!seen.insert(e).second) return false; return true; }

    // Cycles must be edge-disjoint (moved from SignedGraph)
    inline bool are_cycles_edge_disjoint(const std::vector<NegativeCycle>& cycles) const {
        std::unordered_set<Edge, EdgeHash> seen_edges;
    
        for (const auto& cycle : cycles) {
            const auto& neg = cycle.neg_edge();
            if (!seen_edges.insert(neg).second) {
                return false;
            }
            for (const auto& e : cycle.pos_edges()) {
                if (!seen_edges.insert(e).second) {
                    return false;
                }
            }
        }
        return true;
    }

    // Sign tests (tie-break +0/-0 preserved via double cast)
    inline virtual bool is_pos_edge(int eid) const {
        const double w = static_cast<double>(sigma_[(int)eid]);
        return is_pos_(w);
    }
    inline virtual bool is_neg_edge(int eid) const {
        const double w = static_cast<double>(sigma_[(int)eid]);
        return is_neg_(w);
    }
    inline virtual bool is_pos_edge(int u, int v) const { igraph_integer_t eid; if (igraph_get_eid(g_.get(), &eid, u, v, 0, 0) != IGRAPH_SUCCESS) return false; return is_pos_edge((int)eid); }
    inline virtual bool is_neg_edge(int u, int v) const { igraph_integer_t eid; if (igraph_get_eid(g_.get(), &eid, u, v, 0, 0) != IGRAPH_SUCCESS) return false; return is_neg_edge((int)eid); }

    inline int negative_edge_count() const { int neg=0, m=edge_count(); for (int eid=0; eid<m; ++eid) if (is_neg_edge(eid)) ++neg; return neg; }
    inline std::vector<Edge> get_negative_edges() const { std::vector<Edge> out; out.reserve(edge_count()/3); for (int eid=0; eid<edge_count(); ++eid) if (is_neg_edge(eid)) { igraph_integer_t u,v; igraph_edge(g_.get(), eid, &u, &v); out.emplace_back((int)u,(int)v);} return out; }
    inline std::vector<int> get_negative_eids() const { std::vector<int> ids; ids.reserve(edge_count()/3); for (int eid=0; eid<edge_count(); ++eid) if (is_neg_edge(eid)) ids.push_back(eid); return ids; }
    inline int positive_edge_count() const { int pos=0,m=edge_count(); for (int eid=0; eid<m; ++eid) if (is_pos_edge(eid)) ++pos; return pos; }
    inline std::vector<Edge> get_positive_edges() const { std::vector<Edge> out; out.reserve(edge_count()/3); for (int eid=0; eid<edge_count(); ++eid) if (is_pos_edge(eid)) { igraph_integer_t u,v; igraph_edge(g_.get(), eid, &u, &v); out.emplace_back((int)u,(int)v);} return out; }
    inline std::vector<int> get_positive_eids() const { std::vector<int> ids; ids.reserve(edge_count()/3); for (int eid=0; eid<edge_count(); ++eid) if (is_pos_edge(eid)) ids.push_back(eid); return ids; }

    inline void print_info() const { std::cout << "Graph has " << vertex_count() << " vertices and " << edge_count() << std::endl; }
};

class SignedGraph : public GraphCore {

protected:
    std::vector<double> s_;        // size |V|, entries in [-1,1]

    // share topology, replace sigma (passed as ±1 doubles), and s
    SignedGraph(const std::shared_ptr<igraph_t> base, std::vector<double> new_sigma, std::vector<double> new_s);
    SignedGraph(const std::string& file_path);

    SignedGraph(const SignedGraph* const other);
    SignedGraph(const SignedGraph* const other, std::vector<double> new_sigma);
    SignedGraph(const SignedGraph* const other, std::vector<double> new_sigma, std::vector<double> new_s);

    void flip_and_refresh_net_degree(int u, std::vector<double>& dtilde);
    std::vector<int> maximal_salience_clique_strict_pos(const std::vector<int>& Z0) const;
    bool has_fractional_switching(double eps = 1e-12) const;
    void single_switching(int u);
    void single_switching(int u, igraph_vector_int_t* incident);

    // compose switching (acts on s_ only; clamp to [-1,1])
    template<class Vec>
    void apply_compose_switching(const Vec& s) {
        for (int u = 0; u < vertex_count(); ++u) {
            s_[u] *= static_cast<double>(s[u]);
            if (s_[u] < -1.0) s_[u] = -1.0; else if (s_[u] > 1.0) s_[u] = 1.0;
        }
        compute_degrees();
    }

private:
    // sign accessors using s[u]*sigma*s[v] with existing tie-break
    inline bool is_pos_edge(int eid, int u, int v) const {
        const double p = s_[(int)u] * static_cast<double>(sigma_[(int)eid]) * s_[(int)v];
        return is_pos_(p);
    }
    inline bool is_neg_edge(int eid, int u, int v) const {
        const double p = s_[(int)u] * static_cast<double>(sigma_[(int)eid]) * s_[(int)v];
        return is_neg_(p);
    }

public:
    // enable copy/move semantics for efficient transfers
    SignedGraph(const SignedGraph&) = default;
    SignedGraph& operator=(const SignedGraph&) = default;
    SignedGraph(SignedGraph&&) = default;
    SignedGraph& operator=(SignedGraph&&) = default;

    using GraphCore::is_pos_edge;
    using GraphCore::is_neg_edge;

    std::unique_ptr<SignedGraph> clone() const;
    virtual ~SignedGraph();

    // Expose raw sigma
    inline const std::vector<int>& get_sigma() const { return sigma_; }
    inline const std::vector<double>& switching_vector() const { return s_; }

    int frustrated_edges() const;
	const std::vector<Edge> frustrated_edges_keys() const;

    inline double vertex_salience(int u) const noexcept { return 1.0 - std::fabs(s_[u]); }
    inline double edge_salience(int u, int v) const noexcept override {
        const double su = s_[u], sv = s_[v];
        double sal = is_neg_edge(u,v) ? 1.0 - std::max(0.0,0.5*(su+sv)) : 0.5 * (1.0 + std::min(su, sv));
        if (sal < 0.0) sal = 0.0; if (sal > 1.0) sal = 1.0; return sal;
    }
    inline std::vector<double> edge_saliences() const {
        const int m = edge_count(); std::vector<double> es(m, 0.0);
        for (int eid = 0; eid < m; ++eid) { igraph_integer_t u,v; igraph_edge(g_.get(), eid, &u, &v); es[eid] = edge_salience((int)u,(int)v); }
        return es;
    }
    bool are_cycles_sign_correct(const std::vector<NegativeCycle>& cycles, bool expect_negative = true) const;
    void save_partition_svg(const std::vector<int>& partition, const std::string& filename, bool custom_layout) const;

    // Override sign accessors using s[u]*sigma*s[v] with existing tie-break
    inline bool is_pos_edge(int eid) const override {
        igraph_integer_t u,v; igraph_edge(g_.get(), eid, &u, &v);
        return is_pos_edge(eid, u, v);
    }
    inline bool is_neg_edge(int eid) const override {
        igraph_integer_t u,v; igraph_edge(g_.get(), eid, &u, &v);
        return is_neg_edge(eid, u, v);
    }

    // Override sign accessors using s[u]*sigma*s[v] with existing tie-break
    inline bool is_pos_edge(int u, int v) const { igraph_integer_t eid; if (igraph_get_eid(g_.get(), &eid, u, v, 0, 0) != IGRAPH_SUCCESS) return false; return is_pos_edge((int)eid, u, v); }
    inline bool is_neg_edge(int u, int v) const { igraph_integer_t eid; if (igraph_get_eid(g_.get(), &eid, u, v, 0, 0) != IGRAPH_SUCCESS) return false; return is_neg_edge((int)eid, u, v); }
};