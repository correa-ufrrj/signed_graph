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
#include <edge.h>


struct IGraphDeleter {
  void operator()(igraph_t* p) const noexcept { if (p) { igraph_destroy(p); delete p; } }
};

class GraphCore {
protected:
	const std::shared_ptr<igraph_t> g_;   // shared topology
	std::vector<double>        		weights_; // this instance’s edge weights (±1 here)
    std::vector<double> 			d_plus;
    std::vector<double> 			d_minus;

    using MapType = std::unordered_map<Edge, int, EdgeHash>;
    MapType _edge_index;
	
	// share topology, copy weights_
	GraphCore(const GraphCore&)            = default;
	GraphCore& operator=(const GraphCore&) = default;
	GraphCore(const GraphCore* const other);
    GraphCore(const std::string& file_path);
	GraphCore(const GraphCore* const other, std::vector<double> new_weights);
	
	// share topology, replace weights
	GraphCore(const std::shared_ptr<igraph_t> base, std::vector<double> new_weights);

	inline void compute_degrees() {
	    int n = vertex_count();
	    d_plus.assign(n, 0.0);
	    d_minus.assign(n, 0.0);
	
	    for (int eid = 0; eid < edge_count(); ++eid) {
	        igraph_integer_t from, to;
	        igraph_edge(g_.get(), eid, &from, &to);
	        double w = weights_[eid];
	        if (is_pos_edge(eid)) {
	            d_plus[from] += w;
	            d_plus[to] += w;
	        } else {
	            d_minus[from] += -w;
	            d_minus[to] += -w;
	        }
	    }
	};

public:
	virtual ~GraphCore() = default;

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

	class SignedEdgesView {
	private:
	    const igraph_t* g;
	    const GraphCore* core; // to call is_pos_edge/is_neg_edge
	
	public:
	    explicit SignedEdgesView(const GraphCore* core, const igraph_t* graph)
	        : g(graph), core(core) {}
	
	    class const_iterator {
	    private:
	        const igraph_t* g;
	        igraph_integer_t eid;
	        const GraphCore* core;
	
	    public:
	        using value_type = std::pair<Edge, int>; // (edge, sign)
	
	        const_iterator(const igraph_t* g, igraph_integer_t eid, const GraphCore* core)
	            : g(g), eid(eid), core(core) {}
	
	        value_type operator*() const {
	            igraph_integer_t from, to;
	            igraph_edge(g, eid, &from, &to);
	            int sign = core->is_neg_edge((int)eid) ? -1 : +1; // zero handled by is_neg_edge
	            return { { static_cast<int>(from), static_cast<int>(to) }, sign };
	        }
	
	        const_iterator& operator++() { ++eid; return *this; }
	        bool operator!=(const const_iterator& other) const { return eid != other.eid; }
	    };
	
	    const_iterator begin() const { return const_iterator(g, 0, core); }
	    const_iterator end()   const { return const_iterator(g, igraph_ecount(g), core); }
	
	    inline SignedEdge operator[](igraph_integer_t e) const {
	        igraph_integer_t u, v; igraph_edge(g, e, &u, &v);
	        int sign = core->is_neg_edge((int)e) ? -1 : +1;
	        return { Edge{ (int)u, (int)v }, sign };
	    }
	    inline SignedEdge operator[](const Edge& e) const {
	        igraph_integer_t eid; if (igraph_get_eid(g, &eid, e.first, e.second, 0, 0) != IGRAPH_SUCCESS)
	            return { e, 0 };
	        return (*this)[eid];
	    }
	
	    int size() const { return igraph_ecount(g); }
	};

	class NegativeEdgesView {
	private:
	    const igraph_t* g;
	    const GraphCore* core; // to call is_pos_edge/is_neg_edge
	
	public:
	    explicit NegativeEdgesView(const GraphCore* core, const igraph_t* g_) : g(g_), core(core) {}
	
	    struct const_iterator {
	        const igraph_t* g;
	        const GraphCore* core;
	        igraph_integer_t eid;
	        igraph_integer_t m;
	
	        void skip_to_next_neg() {
	            while (eid < m) {
	                if (core->is_neg_edge((int)eid)) break;
	                ++eid;
	            }
	        }
	
	        const_iterator(const igraph_t* g_, const GraphCore* core, igraph_integer_t start)
	            : g(g_), core(core), eid(start), m(igraph_ecount(g_)) { skip_to_next_neg(); }
	
	        bool operator!=(const const_iterator& other) const { return eid != other.eid; }
	        const_iterator& operator++() { ++eid; skip_to_next_neg(); return *this; }
	
	        std::pair<Edge,int> operator*() const {
	            igraph_integer_t u, v; igraph_edge(g, eid, &u, &v);
	            return { Edge{ (int)u, (int)v }, (int)eid };
	        }
	    };
	
	    const_iterator begin() const { return const_iterator(g, core, 0); }
	    const_iterator end()   const { return const_iterator(g, core, igraph_ecount(g)); }
	};

	inline int vertex_count() const { return igraph_vcount(g_.get()); };
	inline int edge_count() const { return igraph_ecount(g_.get()); };
    inline const WeightsView weights_view() const { return WeightsView(g_.get(), &weights_); };
	inline const EdgeIndexesView edge_index() const {
	    return GraphCore::EdgeIndexesView{_edge_index, g_.get()};
	};
	inline SignedEdgesView signs_view() const { return SignedEdgesView(this, g_.get()); }
	inline NegativeEdgesView negative_edges_view() const { return NegativeEdgesView(this, g_.get()); }
    inline const std::vector<double>& get_weights() const { return weights_; };
	inline const std::vector<double>& get_plus_degrees() const { return d_plus; };
	inline const std::vector<double>& get_minus_degrees() const { return d_minus; };
	inline const double net_degree(int u) const noexcept { return d_plus[u] - d_minus[u]; };
	inline const std::vector<double> net_degrees() const {
		std::vector<double> d(vertex_count());
		for (int u = 0; u < vertex_count(); u++)
			d[u] = net_degree(u);
		return d;
	};
	inline const long double max_plus_degree_vertex() const {
	    int n = vertex_count();
	
	    long int max_degree = 0;
	    long int max_vertex = 0;
	
	    for (int i = 0; i < n; ++i) {
	        if (d_plus[i] > max_degree) {
	            max_degree = d_plus[i];
	            max_vertex = i;
	        }
	    }
	    return max_vertex;
	};
	inline const long double max_degree_vertex() const {
	    int n = vertex_count();
	
	    igraph_vector_int_t degrees;
	    igraph_real_t max_degree = 0;
	    long int max_vertex = 0;
	
	    igraph_vector_int_init(&degrees, n);
	    igraph_degree(g_.get(), &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
	
	    for (igraph_integer_t i = 0; i < n; ++i) {
	        igraph_real_t deg = VECTOR(degrees)[i];
	        if (deg > max_degree) {
	            max_degree = deg;
	            max_vertex = i;
	        }
	    }
	
	    igraph_vector_int_destroy(&degrees);
	    return max_vertex;
	};
	inline const std::vector<int> neighbors(int u) const {
	    igraph_vector_int_t neighbors;
	    igraph_vector_int_init(&neighbors, 0);
	    igraph_neighbors(g_.get(), &neighbors, u, IGRAPH_ALL);
	
	    std::vector<int> result(VECTOR(neighbors), VECTOR(neighbors) + igraph_vector_int_size(&neighbors));
	    igraph_vector_int_destroy(&neighbors);
	    return result;
	};
    inline bool are_disjoint(const std::vector<std::vector<Edge>>& edge_sets) const {
	    std::unordered_set<Edge, EdgeHash> seen_edges;
	
	    for (const auto& eset : edge_sets) {
	        for (const auto& e : eset) {
	            if (!seen_edges.insert(e).second) {
	                return false;
	            }
	        }
	    }
	    return true;
	};
	inline bool is_pos_edge(int eid) const {
	    const double w = weights_[(int)eid];
	    return (w > 0.0) || (w == 0.0 && !std::signbit(w)); // +0 → positive
	};
	inline bool is_neg_edge(int eid) const {
	    const double w = weights_[(int)eid];
	    return (w < 0.0) || (w == 0.0 && std::signbit(w)); // −0 → negative
	};
	inline bool is_pos_edge(int u, int v) const {
	    igraph_integer_t eid; if (igraph_get_eid(g_.get(), &eid, u, v, 0, 0) != IGRAPH_SUCCESS) return false;
	    return is_pos_edge(eid);
	};
	inline bool is_neg_edge(int u, int v) const {
	    igraph_integer_t eid; if (igraph_get_eid(g_.get(), &eid, u, v, 0, 0) != IGRAPH_SUCCESS) return false;
	    return is_neg_edge(eid);
	};
	inline int negative_edge_count() const {
	    int neg = 0; int m = edge_count();
	    for (int eid = 0; eid < m; ++eid) {
	        if (is_neg_edge(eid)) ++neg;
	    }
	    return neg;
	};
	inline std::vector<Edge> get_negative_edges() const {
	    std::vector<Edge> out;
	    out.reserve(edge_count()/3);
	    for (int eid = 0; eid < edge_count(); ++eid) {
	        if (is_neg_edge(eid)) {
	            igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
	            out.emplace_back((int)u, (int)v);
	        }
	    }
	    return out;
	}	
	inline std::vector<int> get_negative_eids() const {
	    std::vector<int> ids; ids.reserve(edge_count()/3);
	    for (int eid = 0; eid < edge_count(); ++eid) {
	        if (is_neg_edge(eid)) ids.push_back(eid);
	    }
	    return ids;
	}
	inline int positive_edge_count() const {
	    int pos = 0; int m = edge_count();
	    for (int eid = 0; eid < m; ++eid) {
	        if (is_pos_edge(eid)) ++pos;
	    }
	    return pos;
	};
	inline std::vector<Edge> get_positive_edges() const {
	    std::vector<Edge> out;
	    out.reserve(edge_count()/3);
	    for (int eid = 0; eid < edge_count(); ++eid) {
	        if (is_pos_edge(eid)) {
	            igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
	            out.emplace_back((int)u, (int)v);
	        }
	    }
	    return out;
	}	
	inline std::vector<int> get_positive_eids() const {
	    std::vector<int> ids; ids.reserve(edge_count()/3);
	    for (int eid = 0; eid < edge_count(); ++eid) {
	        if (is_pos_edge(eid)) ids.push_back(eid);
	    }
	    return ids;
	}

	inline void print_info() const {
	    std::cout << "Graph has " << vertex_count() << " vertices and " << edge_count() << " edges." << std::endl;
	};
};

class SignedGraph : public GraphCore {
protected:
	std::vector<double> s_;        // size |V|, entries ±1
	bool             	is_switching_{false};
	
	// share topology, replace weights
	SignedGraph(const std::shared_ptr<igraph_t> base, std::vector<double> new_weights, std::vector<double> new_s);
	SignedGraph(const std::string& file_path);

    SignedGraph(const SignedGraph* const other);
    SignedGraph(const SignedGraph* const other, std::vector<double> new_weights);
	SignedGraph(const SignedGraph* const other, std::vector<double> new_weights, std::vector<double> new_s);
            
private:

    void single_switching(int u);
    void single_switching(int u, igraph_vector_int_t* incident);

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
        bool   relax_to_all_pos_if_Z0_empty = true;
        int    delta_m_minus_cap = 512;
        double delta_m_minus_penalty = 0.0;

        // New heuristic knobs
        int R_max = 1;     // round cap
        int K_max = 1;     // max consecutive zero-clique flips per round
        int L_max = 2;     // max integer zero-clique flips in the detour
        int Delta = 0;     // gate on current working signature (deg(u) <= Delta)
    };

    std::unique_ptr<SignedGraph> clone() const;
    virtual ~SignedGraph();

    inline const std::vector<double>& switching_vector() const { return s_; }

    std::vector<int> negative_triangle_count_per_vertex() const;
    int negative_triangle_count_of_vertex(int u) const;

    int frustrated_edges() const;
    std::vector<Edge> frustrated_edges_keys(const std::vector<int>& partition) const;

	// --- Switching interfaces (always apply to *current* state) ---
	// 1) compose a given mask d (±1 or fractional in [-1,1])
	template<class Vec>
	void apply_compose_switching(const Vec& s) {
	  // update node mask
	  for (int u = 0; u < vertex_count(); ++u) s_[u] *= static_cast<double>(s[u]);
	  // update edge weights in place
	  for (int eid = 0; eid < edge_count(); ++eid) {
	    igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
	    weights_[eid] *= static_cast<double>(s[(int)u]) * static_cast<double>(s[(int)v]);
	    // (optional) clamp to [-1,1] if you want hard bounding
	  }
	  compute_degrees();
	}
    template<class Vec>
	SignedGraph compose_switching(const Vec& s) {
		SignedGraph sg(g_, weights_, s_);
		sg.is_switching_ = true;
	   	sg.apply_compose_switching(s);
	   	return sg;
	}

    void apply_greedy_switching(const GreedyKickOptions& opts);
    SignedGraph greedy_switching(const GreedyKickOptions& opts);
    SignedGraph greedy_switching();

	inline double vertex_salience(int u) const noexcept { return 1.0 - std::fabs(s_[u]); } // 1 - abs(2x_u - 1)
	inline double edge_salience(int u, int v) const noexcept {
	    const double su = s_[u];   // in [-1, 1], su = 2*x_u - 1
	    const double sv = s_[v];
	
	    double sal = is_neg_edge(u,v) ? 1.0 - std::max(0.0,0.5*(su+sv)) : 0.5 * (1.0 + std::min(su, sv));
	
	    // numerical safety (should already be in [0,1])
	    if (sal < 0.0) sal = 0.0;
	    if (sal > 1.0) sal = 1.0;
	    return sal;
	}

	inline std::vector<double> edge_saliences() const {
	    const int m = edge_count();
	
	    std::vector<double> es(m, 0.0);
	    const auto& sv = signs_view();  // full-eid order
	
	    for (int eid = 0; eid < m; ++eid) {
	        const int u = (int)sv[eid].points.first;
	        const int v = (int)sv[eid].points.second;
	        es[eid] = edge_salience(u, v);  // keep your current salience formula
	        // If you prefer salience in [0,1], use:
	        // es[eid] = std::min(1.0, 0.25 * edge_salience(u, v, yuv));
	    }
	    return es;
	}
    SignedGraph integer_projection() const;

    // === Promoted helpers from apply_greedy_switching ===
    std::vector<int> maximal_salience_clique_strict_pos(const std::vector<int>& Z0) const;
    bool has_fractional_switching(double eps = 1e-12) const;
    void flip_and_refresh_net_degree(int u, std::vector<double>& dtilde);

    bool are_cycles_edge_disjoint(const std::vector<NegativeCycle>& cycles) const;
    bool are_cycles_sign_correct(const std::vector<NegativeCycle>& cycles, bool expect_negative = true) const;

    void save_partition_svg(const std::vector<int>& partition, const std::string& filename, bool custom_layout) const;
};
