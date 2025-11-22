// signed_graph_mip.h
#pragma once

#include "signed_graph.h"
#include <vector>
#include <optional>
#include <cassert>
#include <utility>

// Forward-declare the stream class
class ShortestPathGraph;
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
    mutable std::unique_ptr<ShortestPathGraph> sp_graph_;

	explicit SignedGraphForMIP(int n, double p, double q, uint64_t seed);
	explicit SignedGraphForMIP(int n, double p, double q);

	// Helper: compute x = (s + 1)/2
    inline double x_from_s_(double s) const { return (s + 1.0) / 2.0; }
    
	// Rebuild d̃⁺/d̃⁻ from scratch (linear in |E|).
	// d̃⁺[u] = Σ_v max(0, p_uv), d̃⁻[u] = Σ_v max(0, −p_uv), p_uv = s[u]·σ_uv·s[v] ∈ [−1,1].
	inline void recompute_polarity_degrees_() {
	    const int n = vertex_count();
	    const int m = edge_count();
	    dtilde_pos_.assign(n, 0.0);
	    dtilde_neg_.assign(n, 0.0);
	
	    for (int eid = 0; eid < m; ++eid) {
	        igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
	        const double p = polarity_of((int)u, (int)v, eid);  // uses s_ and GraphCore::sign_of
	        const double pp = (p >= 0.0 ?  p : 0.0);
	        const double pn = (p <= 0.0 ? -p : 0.0);
	        dtilde_pos_[(int)u] += pp; dtilde_pos_[(int)v] += pp;
	        dtilde_neg_[(int)u] += pn; dtilde_neg_[(int)v] += pn;
	    }
	}

	// polarity_if(u,v,su): counterfactual polarity for edge (u,v)
	// if s[u] were su, keeping s[v] fixed. Uses raw σ via GraphCore::sign_of.
    inline double polarity_if(const int u, const int v, double su) const {
//		igraph_integer_t eid; igraph_get_eid(g_.get(), &eid, u, v, 0, 0);
		int eid = _edge_index.at(Edge{u,v});
        return su * GraphCore::sign_of(u, v, eid) * s_[v];
    }
    inline double polarity_if(const int u, const int v, int eid, double su) const {
        return su * GraphCore::sign_of(u, v, eid) * s_[v];
    }
    inline double polarity_of(const int u, const int v, int eid) const {
        return polarity_if(u, v, eid, s_[u]);
    }
    inline double polarity_of(const int u, const int v) const {
        return polarity_if(u, v, s_[u]);
    }
    inline double polarity_of(const Edge& e) const {
        return polarity_if(e.first, e.second, s_[e.first]);
    }
    size_t count_crossing_neg_edges_(const GraphCore::Bitmap& anchor) const;
	void on_vertex_flip_(int u, double from, double to) override;
	void on_neg_flip_batch_(const GraphCore::Bitmap& anchor) override;

public:
    using SignedGraph::compose_switching_inplace;

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
                      
    SignedGraphForMIP(const SignedGraphForMIP&) = delete;
    SignedGraphForMIP& operator=(const SignedGraphForMIP&) = delete;

    // Movable: do NOT move the SP; recreate lazily on first use in the new owner
    SignedGraphForMIP(SignedGraphForMIP&&) noexcept;
    SignedGraphForMIP& operator=(SignedGraphForMIP&&) noexcept = delete;

    // Convenience factory
	static std::unique_ptr<SignedGraphForMIP>
	make_random(int n, double p, double q, uint64_t seed);
	static std::unique_ptr<SignedGraphForMIP>
	make_random(int n, double p, double q);

    ~SignedGraphForMIP();

	inline double pos_cum_polarity(int u) const { return dtilde_pos_[u]; }
	inline double neg_cum_polarity(int u) const { return dtilde_neg_[u]; }
	inline double net_polarity(int u) const { return dtilde_pos_[u] - dtilde_neg_[u]; }

    // Integer projection on s (round to ±1). Returns a new MIP graph.
    void apply_integer_projection();
    SignedGraphForMIP integer_projection() const;

    // Locally update s[u] from xu ∈ [0,1] via s[u] = 2*xu − 1 (affects only incident edges).
    void align_switching(int u, double xu);

    // Greedy switching interfaces (return-by-value); this instance remains mutable too.
    void apply_greedy_switching();
    SignedGraphForMIP greedy_switching() const;

    // Rebuild current switching from a new x̂ (overwrites current s_; sigma unchanged)
    void reseed_switching(const std::vector<double>& xhat);

    // Compose switching (returns a NEW MIP graph)
    template<class Vec>
    SignedGraphForMIP compose_switching(const Vec& s) {
        SignedGraphForMIP out(this);      // copy MIP directly (reuses g_, sigma_, s_)
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
    inline const double edge_polarity(int u, int v) const { return polarity_of(u,v); };
    
    ShortestPathGraph& shortest_path_graph() const;
};

class ShortestPathGraph final {
public:
    struct Path {
    private:
        std::vector<int>   nodes_;      // endpoints are nodes_.front(), nodes_.back()
        std::vector<int>   pos_peids_;  // edge ids in the POS graph
        double             cost_ = 0.0; // sum of working weights
        bool               reachable_ = false;

        friend class ShortestPathGraph;
    public:
        // light, read-only view
        inline const std::vector<int>&  nodes()     const { return nodes_; }
        inline std::vector<Edge> 		edges()		const {
			std::vector<Edge> edges_; edges_.reserve(std::max(0, length() - 1)); 
	        const auto& ns = nodes_;
	        for (int i = 1; i < (int)ns.size(); ++i){
	            edges_.emplace_back(ns[i-1], ns[i]);
	        }
	        return edges_;
	    }
        inline double                   cost()      const { return cost_; }
        inline bool                     reachable() const { return reachable_; }
        inline int                      length()    const { return (int)nodes_.size(); }
    };

    // Creation — by value from the owner (only stores a reference to owner).
    explicit ShortestPathGraph(const SignedGraphForMIP& owner,
                            double eps = 1e-12, double max_cap = 0.0)
        : G_(owner), w_eps(eps), w_cap(max_cap) {
        // Silence igraph "Couldn't reach some vertices" warnings just for this batch.
        struct IgraphWarningSilencer {
            igraph_warning_handler_t* prev = nullptr;
            static void noop(const char*, const char*, int) {}
            IgraphWarningSilencer() : prev(igraph_set_warning_handler(noop)) {}
            ~IgraphWarningSilencer() { igraph_set_warning_handler(prev); }
        } _silence_igraph_warnings_guard;
		igraph_vector_int_init(&sp_nodes_, 0);
		sp_nodes_init_ = true;
    }

    ~ShortestPathGraph() {
        if (built_) igraph_destroy(&g_pos_);
        if (w_pos_init_) igraph_vector_destroy(&w_pos_);
        if (sp_nodes_init_) igraph_vector_int_destroy(&sp_nodes_);
    }
    ShortestPathGraph(const ShortestPathGraph&) = delete;
    ShortestPathGraph& operator=(const ShortestPathGraph&) = delete;
    // Custom move to avoid double-free of igraph objects
    ShortestPathGraph(ShortestPathGraph&& other) noexcept
        : G_(other.G_)
        , g_pos_(other.g_pos_)
        , built_(other.built_)
        , w_pos_(other.w_pos_)
        , w_pos_init_(other.w_pos_init_)
        , sp_nodes_(other.sp_nodes_)
        , sp_nodes_init_(other.sp_nodes_init_)
        , w_eps(other.w_eps)
        , w_cap(other.w_cap)
        , full2pos_eid_(std::move(other.full2pos_eid_))
        , pos2full_eid_(std::move(other.pos2full_eid_)) {
        other.built_ = false;
        other.w_pos_init_ = false;
        other.sp_nodes_init_ = false;
        // leave other's PODs in a safe-to-destroy state
    }
    ShortestPathGraph& operator=(ShortestPathGraph&& other) noexcept {
        if (this == &other) return *this;
        if (built_)          igraph_destroy(&g_pos_);
        if (w_pos_init_)     igraph_vector_destroy(&w_pos_);
        if (sp_nodes_init_)  igraph_vector_int_destroy(&sp_nodes_);
        // move-transfer
        g_pos_ = other.g_pos_; built_ = other.built_;
        w_pos_ = other.w_pos_; w_pos_init_ = other.w_pos_init_;
        sp_nodes_ = other.sp_nodes_; sp_nodes_init_ = other.sp_nodes_init_;
        w_eps = other.w_eps; w_cap = other.w_cap;
        full2pos_eid_ = std::move(other.full2pos_eid_);
        pos2full_eid_ = std::move(other.pos2full_eid_);
        // null out source
        other.built_ = false;
        other.w_pos_init_ = false;
        other.sp_nodes_init_ = false;
        return *this;
    }

    // ===== Weights (endpoint-based) =====
    // NOTE: Per request, no guards: callers must ensure (u,v) is a POS edge.
    inline double weight(int u, int v) const {
        // Map (u,v) → full_eid via owner's edge index
        const int full_eid = G_.edge_index(u, v);             // -1 if absent

        // Pull the current working weight (Dijkstra weights vector)
        return weight(full_eid);
    }
    inline double weight(const int full_eid) const {
        // Map full_eid → pos_eid; negative edges have -1
        const igraph_integer_t peid = full2pos_eid_[(size_t)full_eid];

        // Pull the current working weight (Dijkstra weights vector)
        return weight_(peid);
    }
    inline void set_weight(int u, int v, double w) {
        const int fe = G_.edge_index(u, v);
        const int pe = full2pos_eid_[(size_t)fe];
        set_weight(pe, w);
    }
    // Bump with clamps (eps, max_cap=0 → no cap)
    inline void bump_weight(int fe, double delta) {
        const int pe = full2pos_eid_[(size_t)fe];
        double w = weight(pe) + delta;
        if (w < w_eps) w = w_eps;
        if (w_cap > 0.0 && w > w_cap) w = w_cap;
        set_weight(pe, w);
    }
    inline void bump_weight(int u, int v, double delta) {
        const int fe = G_.edge_index(u, v);
		bump_weight(fe, delta);
    }
    // Bump with clamps (eps, max_cap=0 → no cap)
    inline void bump_weights(const Path& p, double delta) {
        for (int pe : p.pos_peids_) {
	        double w = weight(pe) + delta;
	        if (w < w_eps) w = w_eps;
	        if (w_cap > 0.0 && w > w_cap) w = w_cap;
	        set_weight(pe, w);
        }
    }

    // Reseed from a FULL-edge “base” vector (size == |E_full|); only POS edges are read.
    void reseed_weights(const std::vector<double>& base) {
        const long mpos = (long)pos2full_eid_.size();
        assert((long)base.size() >= (long)full2pos_eid_.size());
        for (long pe = 0; pe < mpos; ++pe) {
            const int fe = pos2full_eid_[pe];
            double w = base[(size_t)fe];
            set_weight(pe, w);
        }
    }

    // ===== Shortest paths =====
    // Dijkstra(s,t) using current working weights; returns Path with nodes, pos_peids, cost, reachable.
    Path dijkstra(int s, int t) const;

    // Visit edges on a Path by endpoints (delegates topology; caller never sees maps).
    template <class F>
    void for_each_edge(const Path& p, F&& f) const {
        const auto& ns = p.nodes_;
        for (int i = 1; i < (int)ns.size(); ++i) f(ns[i-1], ns[i]);
    }
	template<class F>
	void for_each_edge(F&& f) const {
	    for (const auto& e : G_.get_edges(+1)) f(e.first, e.second);
	}

    // visit full-eids on a Path (mapping hidden; caller gets IDs only).
    template <class F>
    void for_each_eid(const Path& p, F&& f) const {
        for (int pe : p.pos_peids_) f(pos2full_eid_[(size_t)pe]);
    }
    template <class F>
    void for_each_eid(F&& f) const {
        for (int pe : pos2full_eid_) f(pe);
    }

    // visit full-eids on a Path (mapping hidden; caller gets IDs only).
    template <class F>
    void for_each_weight(const Path& p, F&& f) const {
        for (int pe : p.pos_peids_) f(weight(pe));
    }

private:
    friend class SignedGraphForMIP; // only the owner may trigger rebuilds
    const SignedGraphForMIP& G_;

    // POS-only graph + working weights
    igraph_t        g_pos_{};
    bool            built_{false};
    igraph_vector_t w_pos_{};
    bool            w_pos_init_{false};
    mutable bool    sp_nodes_init_{false};
    
    double w_eps{1e-12};
    double w_cap{0.0};

    // FULL↔POS maps (internal)
    std::vector<int> full2pos_eid_; // size |E_full| ; -1 for non-POS
    std::vector<int> pos2full_eid_; // size |E_pos|
    
    // to reuse in dijkstra (mutable: written in const dijkstra())
    mutable igraph_vector_int_t sp_nodes_;

    inline double weight_(const int peid) const {
        // Pull the current working weight (Dijkstra weights vector)
        return (double)VECTOR(w_pos_)[peid];
    }
    inline void set_weight(const int peid, double w) {
        VECTOR(w_pos_)[peid] = w;
    }
    void build_graph_();
    // Rebuild internal POS graph + maps after a switching in the owner.
    inline void notify_signs_changed() { build_graph_(); }

};
