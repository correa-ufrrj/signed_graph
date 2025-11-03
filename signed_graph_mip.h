// signed_graph_mip.h
#pragma once

#include "signed_graph.h"
#include <vector>
#include <optional>

// Forward-declare the stream class
class NegativeCycleBatch;

// Mutable MIP layer over SignedGraph.
// Spec:
// (i)  SignedGraph is immutable for external users.
// (ii) SignedGraphForMIP derives from SignedGraph and is mutable.
// (iii) Parameterized by xhat ∈ [0,1]^V.
// (iv) s[u] = 2*xhat[u] − 1.
// (v)  Default xhat is unit (all ones) ⇒ s ≡ +1.
class SignedGraphForMIP : public SignedGraph {
private:
	// Helper: compute x = (s + 1)/2
    inline double x_from_s_(double s) const { return (s + 1.0) / 2.0; }
	// Helper: round a switching vector to ±1 in-place (0 → +1)
	inline double round_pm1_(double s) const {
	    return is_pos_(s) ? +1.0 : -1.0;
	}
	inline void round_pm1_inplace(std::vector<double>& s) const {
	    for (double& x : s) x = round_pm1_(x);
	}

public:
    using SignedGraph::apply_compose_switching;

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
        int R_max = 1; int Delta = 0;
    };

    // Read-only live view of fractional edge polarities p_{uv} = s[u] * sigma[uv] * s[v] ∈ [-1,1]
    class EdgePolarityView {
    private:
        const igraph_t* g;
        const SignedGraphForMIP* sg;
    public:
        struct const_iterator {
            const igraph_t* g; igraph_integer_t eid; const SignedGraphForMIP* sg;
            using value_type = std::pair<Edge, double>;
            value_type operator*() const {
                igraph_integer_t u,v; igraph_edge(g, eid, &u, &v);
                const double p = sg->switching_vector()[(int)u]
                                 * static_cast<double>(sg->get_sigma()[(int)eid])
                                 * sg->switching_vector()[(int)v];
                return { Edge{ (int)u,(int)v }, p };
            }
            const_iterator& operator++() { ++eid; return *this; }
            bool operator!=(const const_iterator& o) const { return eid != o.eid; }
        };
        explicit EdgePolarityView(const SignedGraphForMIP* sg_, const igraph_t* graph) : g(graph), sg(sg_) {}
        const_iterator begin() const { return const_iterator{ g, 0, sg }; }
        const_iterator end()   const { return const_iterator{ g, igraph_ecount(g), sg }; }
        inline double operator[](igraph_integer_t e) const {
            igraph_integer_t u,v; igraph_edge(g, e, &u, &v);
            return sg->switching_vector()[(int)u]
                 * static_cast<double>(sg->get_sigma()[(int)e])
                 * sg->switching_vector()[(int)v];
        }
        inline double operator[](const Edge& e) const {
            igraph_integer_t eid; if (igraph_get_eid(g, &eid, e.first, e.second, 0, 0) != IGRAPH_SUCCESS) return 0.0;
            return (*this)[eid];
        }
        int size() const { return igraph_ecount(g); }
    };

    // ---- Constructors (declarations only; no inline definitions here) ----
    // Load from file; default xhat ≡ 1 → s ≡ +1
    SignedGraphForMIP(const std::string& file_path);
    // Construct from an existing SignedGraph (immutable base)
    SignedGraphForMIP(SignedGraph&& base);
    // Copy constructor
    SignedGraphForMIP(const SignedGraphForMIP* const other);
    // Replace sigma (±1), keep current xhat → s recomputed from stored xhat
    SignedGraphForMIP(const SignedGraphForMIP* const other, std::vector<double> new_sigma);
    // Replace sigma (±1) and xhat ∈ [0,1]^V; s = 2*xhat − 1
    SignedGraphForMIP(const SignedGraphForMIP* const other, std::vector<double> new_sigma,
                      const std::vector<double>& xhat);
    ~SignedGraphForMIP();

    // Integer projection on s (round to ±1). Returns a new MIP graph.
    void apply_integer_projection();
    SignedGraphForMIP integer_projection() const;

    // Locally update s[u] from xu ∈ [0,1] via s[u] = 2*xu − 1 (affects only incident edges).
    void align_switching(int u, double xu);

    // Greedy switching interfaces (return-by-value); this instance remains mutable too.
    void apply_greedy_switching(const GreedyKickOptions& opts);
    SignedGraphForMIP greedy_switching(const GreedyKickOptions& opts) const;
    SignedGraphForMIP greedy_switching() const;

    // Rebuild current switching from a new x̂ (overwrites current s_; sigma unchanged)
    void reseed_switching(const std::vector<double>& xhat);

    // Compose switching (returns a NEW MIP graph)
    template<class Vec>
    SignedGraphForMIP compose_switching(const Vec& s) {
        SignedGraphForMIP out(this);   // copy MIP directly (reuses g_, sigma_, s_)
        out.apply_compose_switching(s); // updates s_ and recomputes degrees
        return out;
    }

    // Convenience accessors in x/y space
    inline double get_x(int u) const { return x_from_s_(s_[u]); }
    inline double get_y(int u, int v) const { return get_x(u) * get_x(v); }
    inline double get_y(Edge uv) const { return get_y(uv.first, uv.second); }
    inline double get_rounded_x(int u) const { return x_from_s_(round_pm1_(s_[u])); }
    inline double get_rounded_y(int u, int v) const { return get_rounded_x(u) * get_rounded_x(v); }
    inline double get_rounded_y(Edge uv) const { return get_rounded_y(uv.first, uv.second); }

    // Edge polarity helpers (read-only live view and eager snapshot)
    inline EdgePolarityView edge_polarities_view() const { return EdgePolarityView(this, g_.get()); }
    const std::vector<double> edge_polarities() const;

    // Stream factory
    NegativeCycleBatch open_negative_cycle_stream(bool cover, bool use_triangle_order = false) const;
};