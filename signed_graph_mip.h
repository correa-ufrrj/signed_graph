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
// Fractional polarity layer:
//  - polarity_of(u,v,eid) = s[u]*sigma[eid]*s[v] ∈ [−1, +1] (s ∈ [−1,1]^V).
//  - dtilde_pos_[u] = Σ_v max(0, polarity_of(u,v)), dtilde_neg_[u] = Σ_v max(0, −polarity_of(u,v)).
//  - on_switch_sign_changed_(u, from, to, ...):
//       * first calls base SignedGraph::on_switch_sign_changed_ for (+/−) structure.
//       * then updates (dtilde_pos_, dtilde_neg_) by per-edge deltas using σ via GraphCore::sign_of.
//       * does not touch sigma_ directly; no division by 'from' or 'to'.
class SignedGraphForMIP : public SignedGraph {
private:
    // Cached polarity-degree decompositions:
    // d̃⁺[u] = Σ_v max(0, p_uv), d̃⁻[u] = Σ_v max(0, -p_uv), p_uv = s[u]*σ_uv*s[v]
    std::vector<double> dtilde_pos_;
    std::vector<double> dtilde_neg_;

	// Helper: compute x = (s + 1)/2
    inline double x_from_s_(double s) const { return (s + 1.0) / 2.0; }

	// polarity_if(u,v,su): counterfactual polarity for edge (u,v)
	// if s[u] were su, keeping s[v] fixed. Uses raw σ via GraphCore::sign_of.
    inline double polarity_if(int u, int v, double su) const {
        return su * GraphCore::sign_of(u, v) * s_[v];
    }
    inline double polarity_if(int u, int v, int eid, double su) const {
        return su * GraphCore::sign_of(u, v, eid) * s_[v];
    }
    inline double polarity_of(int u, int v, int eid) const {
        return s_[u] * GraphCore::sign_of(u, v, eid) * s_[v];
    }
    inline double polarity_of(const Edge& e) const {
        return polarity_if(e.first, e.second, s_[e.first]);
    }
	void on_switch_sign_changed_(int u, double from, double to, igraph_vector_int_t* incident) override;

public:
    using SignedGraph::compose_switching_inplace;

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
        const SignedGraphForMIP* core;
    public:
        struct const_iterator {
            const igraph_t* g; igraph_integer_t eid; const SignedGraphForMIP* core;
            using value_type = EdgePolarity;
            value_type operator*() const {
                igraph_integer_t u,v; igraph_edge(g, eid, &u, &v);
                return { Edge{ (int)u, (int)v }, core->polarity_of(u, v, eid) };
            }
            const_iterator& operator++() { ++eid; return *this; }
            bool operator!=(const const_iterator& o) const { return eid != o.eid; }
        };
        explicit EdgePolarityView(const SignedGraphForMIP* sg_, const igraph_t* graph) : g(graph), core(sg_) {}
        const_iterator begin() const { return const_iterator{ g, 0, core }; }
        const_iterator end()   const { return const_iterator{ g, igraph_ecount(g), core }; }
        inline EdgePolarity operator[](igraph_integer_t eid) const {
            igraph_integer_t u,v; igraph_edge(g, eid, &u, &v);
            return { Edge{ (int)u, (int)v }, core->polarity_of(u, v, eid) };
        }
        inline EdgePolarity operator[](const Edge& e) const {
            return { e, core->polarity_of(e) };
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
        out.compose_switching_inplace(s); // updates s_ and recomputes degrees
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
    NegativeCycleBatch open_negative_cycle_stream() const;
};
