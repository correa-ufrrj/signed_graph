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

// GraphCore.hpp
#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include <iterator>

class GraphCore {
public:
    // ===== Immutable 64-bit bitmap (externally) =====
    class Bitmap {
    public:
        using word_t = uint64_t;
        static constexpr size_t kWordBits = 64;

        // ---- queries ----
        inline size_t size()  const { return nbits_; }
        inline bool   empty() const { return nbits_ == 0; }
        inline bool contains(size_t i) const {
            return (w_[i >> 6] >> (i & 63)) & 1ULL;
        }

        // ---- popcount ----
        inline size_t count() const {
            size_t s = 0;
            for (word_t x : w_) s += __builtin_popcountll(x);
            return s;
        }

        // ===== Iterator over set indices =====
        class const_iterator {
        public:
            using value_type = size_t;
            using difference_type = std::ptrdiff_t;
            using iterator_category = std::input_iterator_tag;

            const_iterator() = default;
            value_type operator*() const { return cur_; }
            const_iterator& operator++() { advance_(); return *this; }
            const_iterator operator++(int) { auto t = *this; advance_(); return t; }
            bool operator==(const const_iterator& o) const {
                return b_ == o.b_ && wi_ == o.wi_ && rem_ == o.rem_ && cur_ == o.cur_;
            }
            bool operator!=(const const_iterator& o) const { return !(*this == o); }

        private:
            friend class Bitmap;
            const Bitmap* b_ = nullptr;
            size_t wi_ = 0;
            word_t rem_ = 0;
            size_t cur_ = 0;

            void init_(const Bitmap* b, bool end) {
                b_ = b;
                if (!b_ || end || b_->w_.empty()) {
                    wi_ = b_ ? b_->w_.size() : 0; rem_ = 0; cur_ = b ? b->nbits_ : 0; return;
                }
                wi_ = 0; rem_ = b_->w_[0];
                skip_zero_words_();
                cur_ = (wi_ < b_->w_.size()) ? ((wi_ << 6) + __builtin_ctzll(rem_)) : b_->nbits_;
            }
            inline void advance_() {
                rem_ &= (rem_ - 1);
                if (rem_) { cur_ = (wi_ << 6) + __builtin_ctzll(rem_); return; }
                ++wi_; skip_zero_words_();
                cur_ = (wi_ < b_->w_.size()) ? ((wi_ << 6) + __builtin_ctzll(rem_)) : b_->nbits_;
            }
            inline void skip_zero_words_() {
                const size_t n = b_->w_.size();
                while (wi_ < n) {
                    rem_ = b_->w_[wi_];
                    if (rem_) break;
                    ++wi_;
                }
            }
        };

        const_iterator begin() const { const_iterator it; it.init_(this, false); return it; }
        const_iterator end()   const { const_iterator it; it.init_(this, true ); return it; }

        // ---- materialize indices ----
        std::vector<uint32_t> to_vector() const {
            std::vector<uint32_t> out;
            out.reserve(count());
            for (size_t i : *this) out.push_back(static_cast<uint32_t>(i));
            return out;
        }

    private:
        friend class GraphCore; // GraphCore can construct and mutate internally
        size_t nbits_ = 0;
        std::vector<word_t> w_;

        // Construction restricted to GraphCore
        Bitmap(size_t nbits, size_t nwords) : nbits_(nbits), w_(nwords, 0) {}
        explicit Bitmap(size_t nbits)
        : nbits_(nbits), w_((nbits + kWordBits - 1) / kWordBits, 0) {}

        static Bitmap from_indices(size_t nbits, const uint32_t* idx, size_t len) {
            Bitmap bm(nbits);
            for (size_t k = 0; k < len; ++k) {
                const size_t i = idx[k];
                bm.w_[i >> 6] |= (word_t(1) << (i & 63));
            }
            return bm;
        }
        static Bitmap from_words(size_t nbits, std::vector<word_t> words) {
            Bitmap bm(nbits, words.size());
            bm.w_ = std::move(words);                // tail assumed zeroed by caller
            return bm;
        }

        // ===== In-place set algebra (private) =====
        // All assume same shape (true by construction).
        inline void ior (const Bitmap& b)  { for (size_t i = 0, n = w_.size(); i < n; ++i) w_[i] |=  b.w_[i]; }
        inline void iand(const Bitmap& b)  { for (size_t i = 0, n = w_.size(); i < n; ++i) w_[i] &=  b.w_[i]; }
        inline void ixor(const Bitmap& b)  { for (size_t i = 0, n = w_.size(); i < n; ++i) w_[i] ^=  b.w_[i]; }
        inline void isub(const Bitmap& b)  { for (size_t i = 0, n = w_.size(); i < n; ++i) w_[i] &= ~b.w_[i]; }
    };

private:
    std::vector<Bitmap> N_pos_;   // + neighbors
    std::vector<Bitmap> N_neg_;   // - neighbors

    static Bitmap make_bitmap_from_indices(size_t nbits, const uint32_t* idx, size_t len) {
        return Bitmap::from_indices(nbits, idx, len);
    }
    static Bitmap make_bitmap_from_words(size_t nbits, std::vector<Bitmap::word_t> words) {
        return Bitmap::from_words(nbits, std::move(words));
    }
    static Bitmap make_empty_bitmap(size_t nbits) { return Bitmap(nbits); }

    // helpers to build N_pos_/N_neg_
    void init_empty_bitmaps_(int n) {
        N_pos_.clear(); N_neg_.clear();
        N_pos_.reserve(n); N_neg_.reserve(n);
        for (int i = 0; i < n; ++i) { N_pos_.emplace_back(make_empty_bitmap((size_t)n)); N_neg_.emplace_back(make_empty_bitmap((size_t)n)); }
    }
    void rebuild_bitmaps_from_sigma_() {
        const int n = vertex_count(), m = edge_count();
        init_empty_bitmaps_(n);
        for (int eid = 0; eid < m; ++eid) {
            igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
            auto& A = (sigma_[eid] > 0) ? N_pos_ : N_neg_;
            // direct bit set (GraphCore is friend of Bitmap)
            A[(int)u].w_[(int)v >> 6] |= (Bitmap::word_t(1) << ((int)v & 63));
            A[(int)v].w_[(int)u >> 6] |= (Bitmap::word_t(1) << ((int)u & 63));
        }
    }

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
        d_pos.resize(n); d_neg.resize(n);
        for (int u = 0; u < n; ++u) {
            d_pos[u] = static_cast<double>(neighbors_bm((size_t)u, +1).count());
            d_neg[u] = static_cast<double>(neighbors_bm((size_t)u, -1).count());
        }
    }

    static inline void validate_and_fill_sigma(std::vector<int>& dst, const std::vector<double>& src) {
        dst.clear(); dst.reserve(src.size());
        for (double w : src) {
            if (w != 1.0 && w != -1.0) throw std::invalid_argument("Edge signs must be ±1.");
            dst.push_back(w > 0 ? 1 : -1);
        }
    }

    // === Bitmap update hooks used by SignedGraph ===
    void flip_incident_edges_bitmaps_only_(int u);
    void rebuild_bitmaps_from_sigma_and_switching_(const std::vector<double>& s);

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
    // sign>0 → positive neighbors; sign<0 → negative neighbors
    inline const Bitmap& neighbors_bm(size_t u, int sign) const {
        return (sign > 0) ? N_pos_[u] : N_neg_[u];
    }
    // anti = neighbors of the opposite sign (useful for P2 with sign flip)
    inline const Bitmap& antineighbors_bm(size_t u, int sign) const {
        return (sign > 0) ? N_neg_[u] : N_pos_[u];
    }
    // union (by value) for callers that want sign-agnostic neighborhood
    inline Bitmap neighbors_bm(size_t u) const {
        Bitmap out = N_pos_[u];
        out.ior(N_neg_[u]);
        return out;
    }

    // materialized neighbors (sign-agnostic) via bitmap fast path
    inline std::vector<int> neighbors(int u) const {
        Bitmap both = neighbors_bm((size_t)u);               // union (+) | (-)
        std::vector<uint32_t> vv = both.to_vector();
        return std::vector<int>(vv.begin(), vv.end());
    }

    // common neighbors (sign-agnostic): (N⁺(u)∪N⁻(u)) ∩ (N⁺(v)∪N⁻(v))
    inline Bitmap common_neighbors(size_t u, size_t v) const {
        Bitmap a = N_pos_[u]; a.ior(N_neg_[u]);
        Bitmap b = N_pos_[v]; b.ior(N_neg_[v]);
        a.iand(b);
        return a;
    }
    // symmetric difference of sign-agnostic neighborhoods
    inline Bitmap common_antineighbors(size_t u, size_t v) const {
        Bitmap a = N_pos_[u]; a.ior(N_neg_[u]);
        Bitmap b = N_pos_[v]; b.ior(N_neg_[v]);
        a.ixor(b);
        return a;
    }
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
    inline virtual bool is_pos_edge(int u, int v) const {
        return neighbors_bm((size_t)u, +1).contains((size_t)v);
    }
    inline virtual bool is_neg_edge(int u, int v) const {
        return neighbors_bm((size_t)u, -1).contains((size_t)v);
    }
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
    // Guarded switching vector: all writes update N_pos_/N_neg_ bitmaps
    class SwitchingVector {
    public:
        struct Ref {
            SwitchingVector* owner; int u;
            operator double() const noexcept { return owner->vals_[u]; }
            Ref& operator=(double x) { owner->assign_(u, x); return *this; }
            Ref& operator*=(double a) { owner->assign_(u, owner->vals_[u] * a); return *this; }
            Ref& operator+=(double a) { owner->assign_(u, owner->vals_[u] + a); return *this; }
            Ref& operator-=(double a) { owner->assign_(u, owner->vals_[u] - a); return *this; }
        };

        SwitchingVector() = default;
        SwitchingVector(SignedGraph* g, int n, double init_val = 1.0)
            : g_(g), vals_(n, init_val) {
            // Build effective bitmaps from sigma and the (uniform) switching
            if (g_ && g_->vertex_count() == n) {
                if (init_val == 1.0) { /* nothing to flip */ }
                else if (g_->is_neg_(init_val)) {
                    for (int u = 0; u < n; ++u) g_->flip_incident_edges_bitmaps_only_(u);
                } else {
                    // +0 handled as positive by is_pos_ tie-break; no flips
                }
            }
        }
        SwitchingVector(SignedGraph* g, const std::vector<double>& init_vals, bool rebuild_bitmaps)
            : g_(g), vals_(init_vals) {
            if (g_ && rebuild_bitmaps) g_->rebuild_bitmaps_from_sigma_and_switching_(vals_);
        }

        Ref operator[](int u) { return Ref{this,u}; }
        double operator[](int u) const { return vals_[u]; }
        int size() const noexcept { return (int)vals_.size(); }
        const std::vector<double>& readonly() const noexcept { return vals_; }

        void begin_batch() { batching_ = true; dirty_.clear(); }
        void end_batch() {
            if (!batching_) return;
            batching_ = false;
            if (dirty_.empty() || !g_) { dirty_.clear(); return; }
            std::sort(dirty_.begin(), dirty_.end());
            dirty_.erase(std::unique(dirty_.begin(), dirty_.end()), dirty_.end());
            for (int u : dirty_) g_->flip_incident_edges_bitmaps_only_(u);
            dirty_.clear();
        }

        // Replace whole vector and rebuild effective bitmaps in O(m)
        void set_all(const std::vector<double>& s_new) {
            vals_ = s_new;
            if (g_) g_->rebuild_bitmaps_from_sigma_and_switching_(vals_);
        }

    private:
        friend class SignedGraph;
        SignedGraph* g_ = nullptr;
        std::vector<double> vals_;
        bool batching_ = false;
        std::vector<int> dirty_;

        inline static double clamp_(double x) {
            if (x < -1.0) return -1.0; if (x > 1.0) return 1.0; return x;
        }
        void assign_(int u, double x) {
            x = clamp_(x);
            const bool was_neg = g_->is_neg_(vals_[u]);
            const bool now_neg = g_->is_neg_(x);
            vals_[u] = x;
            if (was_neg != now_neg) {
                if (batching_) dirty_.push_back(u);
                else if (g_) g_->flip_incident_edges_bitmaps_only_(u);
            }
        }
    };

    struct BatchGuard { SwitchingVector& sv; explicit BatchGuard(SwitchingVector& x) : sv(x) { sv.begin_batch(); } ~BatchGuard() { sv.end_batch(); } };

    SwitchingVector s_;        // size |V|, entries in [-1,1]

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
        BatchGuard guard(s_);                 // coalesce bitmap flips
        const int n = vertex_count();
        for (int u = 0; u < n; ++u) {
            s_[u] *= static_cast<double>(s[u]); // guarded writes
        }
        // guard dtor applies flips; recompute degrees once
        compute_degrees();
    }

private:
    // sign accessors using s[u]*sigma*s[v] with existing tie-break
    inline bool is_pos_edge(int eid, int u, int v) const {
        const double p = s_.readonly()[(int)u] * static_cast<double>(sigma_[(int)eid]) * s_.readonly()[(int)v];
        return is_pos_(p);
    }
    inline bool is_neg_edge(int eid, int u, int v) const {
        const double p = s_.readonly()[(int)u] * static_cast<double>(sigma_[(int)eid]) * s_.readonly()[(int)v];
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

    // === Override eid-based sign queries to reflect switching via bitmaps ===
    inline bool is_pos_edge(int eid) const override {
        igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
        return neighbors_bm((size_t)u, +1).contains((size_t)v);
    }
    inline bool is_neg_edge(int eid) const override {
        igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
        return neighbors_bm((size_t)u, -1).contains((size_t)v);
    }

    // === Batch helpers for efficient guarded updates ===
    inline void flip_vertices_batch(const std::vector<int>& U) {
        BatchGuard guard(s_);               // coalesce bitmap flips
        for (int u : U) s_[u] = -s_[u];     // each write is guarded
        compute_degrees();                  // recompute degrees once
    }
    inline void set_switching_vector(const std::vector<double>& s_new) {
        s_.set_all(s_new);                  // rebuilds N_pos_/N_neg_ from scratch
        compute_degrees();                  // keep degree views in sync
    }

    std::unique_ptr<SignedGraph> clone() const;
    virtual ~SignedGraph();

    // Expose raw sigma
    inline const std::vector<int>& get_sigma() const { return sigma_; }
    inline const std::vector<double>& switching_vector() const { return s_.readonly(); }

    int frustrated_edges() const;
    const std::vector<Edge> frustrated_edges_keys() const;

    inline double vertex_salience(int u) const noexcept { return 1.0 - std::fabs(s_.readonly()[u]); }
    inline double edge_salience(int u, int v) const noexcept override {
        const double su = s_.readonly()[u], sv = s_.readonly()[v];
        return 0.5 * (vertex_salience(u) + vertex_salience(v));
    }
    inline std::vector<double> edge_saliences() const {
        const int m = edge_count(); std::vector<double> es(m, 0.0);
        for (int eid = 0; eid < m; ++eid) { igraph_integer_t u,v; igraph_edge(g_.get(), eid, &u, &v); es[eid] = edge_salience((int)u,(int)v); }
        return es;
    }
    bool are_cycles_sign_correct(const std::vector<NegativeCycle>& cycles, bool expect_negative = true) const;
    void save_partition_svg(const std::vector<int>& partition, const std::string& filename, bool custom_layout) const;

    
};
