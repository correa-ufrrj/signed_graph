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
#include <cstdint>
#include <cstddef>
#include <vector>
#include <iterator>


struct IGraphDeleter {
  void operator()(igraph_t* p) const noexcept { if (p) { igraph_destroy(p); delete p; } }
};


// Sign convention (immutable base):
//  - sigma_[eid] ∈ {−1, +1} is the raw, unswitched sign of edge eid.
//  - GraphCore::sign_of(u,v,eid) returns sigma_[eid] (no switching).
//  - Bitmap rows N_pos_/N_neg_ are maintained to reflect the CURRENT effective sign;
//    in GraphCore they are built from sigma only (no switching).
class GraphCore {
public:
    // ===== Immutable 64-bit bitmap (externally) =====
    class Bitmap {
    public:
        using word_t = uint64_t;
        static constexpr size_t kWordBits = 64;


        // ---- queries ----
        inline size_t extent()  const { return nbits_; }
        inline bool   empty() const {
            for (word_t x : w_) if (x) return false;
            return true;
        }
        inline bool contains(size_t i) const {
            return (w_[i >> 6] >> (i & 63)) & 1ULL;
        }
		// Removes and returns the smallest element (index) in the set.
		// Precondition: !empty(); returns -1 only if called on an empty bitmap.
		inline int pop_front() {
		    for (size_t wi = 0, n = w_.size(); wi < n; ++wi) {
		        word_t word = w_[wi];
		        if (word) {
		            const unsigned tz = (unsigned)__builtin_ctzll(word); // index of lowest 1-bit
		            w_[wi] = word & (word - 1);                          // clear that bit
		            return int((wi << 6) + tz);
		        }
		    }
		    return -1; // unreachable if precondition holds
		}

        // ---- popcount ----
        inline size_t count() const {
            size_t s = 0;
            for (word_t x : w_) s += __builtin_popcountll(x);
            return s;
        }
        
        static inline void swap(Bitmap& one, Bitmap& other) {
        	std::swap(one.w_, other.w_);
		};


        // ===== In-place set algebra (private) =====
        // All assume same shape (true by construction).
        inline void ineg ()  {
			for (size_t i = 0, n = w_.size(); i < n; ++i) w_[i] = ~w_[i];
		    if (!w_.empty()) {
		        const size_t rem = nbits_ & 63;
		        if (rem) {
		            const Bitmap::word_t mask = (Bitmap::word_t(1) << rem) - 1;
		            w_.back() &= mask;
		        }
		    }
		 }
        inline void ior (const Bitmap& b)  { for (size_t i = 0, n = w_.size(); i < n; ++i) w_[i] |=  b.w_[i]; }
        inline void iand(const Bitmap& b)  { for (size_t i = 0, n = w_.size(); i < n; ++i) w_[i] &=  b.w_[i]; }
        inline void ixor(const Bitmap& b)  { for (size_t i = 0, n = w_.size(); i < n; ++i) w_[i] ^=  b.w_[i]; }
        inline void isub(const Bitmap& b)  { for (size_t i = 0, n = w_.size(); i < n; ++i) w_[i] &= ~b.w_[i]; }
        inline void iclear(const size_t u) { const size_t wu = u >> 6; const Bitmap::word_t umask = Bitmap::word_t(1) << ((size_t)u & 63); w_[wu] &= ~umask; }
        inline void iset(const size_t u) { const size_t wu = u >> 6; const Bitmap::word_t umask = Bitmap::word_t(1) << ((size_t)u & 63); w_[wu] |= umask; }
		inline Bitmap complement() const { Bitmap out = *this; out.ineg(); return out; }

        // ===== Iterator over set indices, increasing order =====
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
    };


private:
    // Immutable edge sign pattern (±1 ints)
    std::vector<int>    sigma_;


    static Bitmap make_bitmap_from_indices(size_t nbits, const uint32_t* idx, size_t len) {
        return Bitmap::from_indices(nbits, idx, len);
    }
    static Bitmap make_bitmap_from_words(size_t nbits, std::vector<Bitmap::word_t> words) {
        return Bitmap::from_words(nbits, std::move(words));
    }
    static Bitmap make_empty_bitmap(size_t nbits) { return Bitmap(nbits); }


    // helpers to build N_pos_/N_neg_
    inline void init_empty_bitmaps_(int n) {
        N_pos_.clear(); N_neg_.clear();
        N_pos_.reserve(n); N_neg_.reserve(n);
        for (int i = 0; i < n; ++i) { N_pos_.emplace_back(make_empty_bitmap((size_t)n)); N_neg_.emplace_back(make_empty_bitmap((size_t)n)); }
    }
	inline void validate_and_fill_sigma(std::vector<int>& dst,
	                                           const std::vector<double>& src) {
	    dst.clear();
	    dst.reserve(src.size());
	    for (double w : src) {
	        if (w != 1.0 && w != -1.0)
	            throw std::invalid_argument("Edge signs must be ±1.");
	        dst.push_back(w > 0.0 ? 1 : -1);
	    }
	}
    inline void rebuild_bitmaps_from_sign_() {
        const int n = vertex_count(), m = edge_count();
        init_empty_bitmaps_(n);
        for (int eid = 0; eid < m; ++eid) {
            igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
            auto& A = is_pos_edge(eid) ? N_pos_ : N_neg_;
            // direct bit set (GraphCore is friend of Bitmap)
            A[(int)u].w_[(int)v >> 6] |= (Bitmap::word_t(1) << ((int)v & 63));
            A[(int)v].w_[(int)u >> 6] |= (Bitmap::word_t(1) << ((int)u & 63));
        }
    }


protected:
    const std::shared_ptr<igraph_t> g_;   // shared topology
    std::vector<Bitmap> N_pos_;   // + neighbors
    std::vector<Bitmap> N_neg_;   // - neighbors
    // Aggregations (by sign counts)
    std::vector<int>    d_pos;
    std::vector<int>    d_neg;
    int    				m_pos;
    int    				m_neg;


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
    inline virtual int sign_of(const int u, const int v, int eid) const {
        return sigma_[(int)eid];
    }

    inline void compute_degrees() {
        const int n = vertex_count();
        d_pos.resize(n); d_neg.resize(n);
        m_pos = 0; m_neg = 0;
        for (int u = 0; u < n; ++u) {
            d_pos[u] = static_cast<double>(neighbors_bm((size_t)u, +1).count());
            d_neg[u] = static_cast<double>(neighbors_bm((size_t)u, -1).count());
            m_pos += d_pos[u];
            m_neg += d_neg[u];
        }
        m_pos >>= 1;
        m_neg >>= 1;
    }

	void on_vertex_flip_(int u);
	virtual void on_neg_flip_batch_(const GraphCore::Bitmap& anchor);

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
                return { Edge{ (int)u, (int)v }, core->sign_of(u, v, eid) };
            }
            const_iterator& operator++() { ++eid; return *this; }
            bool operator!=(const const_iterator& other) const { return eid != other.eid; }
        };
        const_iterator begin() const { return const_iterator(g, 0, core); }
        const_iterator end()   const { return const_iterator(g, igraph_ecount(g), core); }
        inline SignedEdge operator[](igraph_integer_t eid) const {
            igraph_integer_t u,v; igraph_edge(g, eid, &u, &v);
            return { Edge{ (int)u, (int)v }, core->sign_of(eid) };
        }
        inline SignedEdge operator[](const Edge& e) const {
            return { e, core->sign_of(e) };
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
    inline const int edge_index(int u, int v) const { return _edge_index.at(Edge{u,v}); }
    inline SignedEdgesView signs_view() const { return SignedEdgesView(this, g_.get()); }
    inline NegativeEdgesView negative_edges_view() const { return NegativeEdgesView(this, g_.get()); }

    // Degrees
    inline int degree(int u, int sign) const { return (sign > 0) ? d_pos[u] : d_neg[u]; }
    inline int net_degree(int u) const noexcept { return d_pos[u] - d_neg[u]; }
    inline std::vector<double> net_degrees() const { std::vector<double> d(vertex_count()); for (int u = 0; u < vertex_count(); ++u) d[u] = net_degree(u); return d; }
    inline long double max_pos_degree_vertex() const { int n = vertex_count(); long int mx=0, vid=0; for (int i=0;i<n;++i){ if (d_pos[i] > mx) { mx=d_pos[i]; vid=i; } } return vid; }

    // antineighbors: not neighbors
    inline const Bitmap antineighbors_bm(size_t u) const {
        Bitmap out = neighbors_bm(u);
        out.ineg();
        out.iclear(u);
        return out;
    }
    // common antineighborhoods
	inline Bitmap common_antineighbors(size_t u, size_t v) const {
	    Bitmap a = neighbors_bm(u);  // N(u) = N⁺(u) ∪ N⁻(u)
	    a.ior(neighbors_bm(v));      // N(u) ∪ N(v)
	    a.ineg();                    // complement: V \ (N(u) ∪ N(v))
	    a.iclear(u); a.iclear(v);    // never include endpoints
	    return a;
	}
    
    // sign>0 → positive neighbors; sign<0 → negative neighbors
    // enumeration in increasing order of neighbors
    inline const Bitmap& neighbors_bm(size_t u, int sign) const {
        return (sign > 0) ? N_pos_[u] : N_neg_[u];
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

    // common positive neighbors: N⁺(u) ∩ N⁺(v)
    inline Bitmap common_neighbors(size_t u, size_t v, int sign) const {
        if (sign > 0) {
			Bitmap a = N_pos_[u];
	        a.iand(N_pos_[v]);
	        return a;
		}
		Bitmap a = N_neg_[u];
    	a.iand(N_neg_[v]);			
        return a;
    }

    // common positive neighbors: N⁺(u) ∩ N⁺(v)
    inline Bitmap common_neighbors(Edge e, int sign) const {
		return common_neighbors((size_t)e.first, (size_t)e.second, sign);
    }

    // common neighbors (sign-agnostic): (N⁺(u)∪N⁻(u)) ∩ (N⁺(v)∪N⁻(v))
    inline Bitmap common_neighbors(size_t u, size_t v) const {
        Bitmap a = neighbors_bm(u);
        a.iand(neighbors_bm(v));
        return a;
    }
    
    std::vector<GraphCore::Bitmap> connected_components(int sign) const;
    
	inline std::vector<int> maximal_clique_strict_pos(const std::vector<int>& order) const {
	    if (order.empty()) return {};
	
	    std::vector<int> Q;
	    Q.reserve(order.size());
	
	    const int seed = order[0];
	    Q.push_back(seed);
	
	    // Running intersection C = ⋂_{w∈Q} N⁺(w)
	    Bitmap C = neighbors_bm((size_t)seed, +1);  // copy; safe to mutate
	
	    for (size_t i = 1; i < order.size(); ++i) {
	        const int u = order[i];
	        if (!C.contains((size_t)u)) continue;   // u is not adjacent to all in Q
	        Q.push_back(u);
	        C.iand(N_pos_[(size_t)u]);              // GraphCore can call Bitmap::iand
	        if (C.empty()) break;              // no further vertex can join
	    }
	    return Q;
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

    // Sign value
    inline virtual int sign_of(int eid) const {
        return sigma_[(int)eid];
    }
    inline int sign_of(const int u, const int v) const {
//        igraph_integer_t eid; igraph_get_eid(g_.get(), &eid, u, v, 0, 0);
		int eid = _edge_index.at(Edge{u,v});
        return sign_of(u, v, eid);
    }
    inline int sign_of(const Edge& e) const {
        return sign_of(e.first, e.second);
    }


    // Sign tests (tie-break +0/-0 preserved via double cast)
    inline bool is_pos_edge(int eid) const {
        return is_pos_(sign_of((int)eid));
    }
    inline bool is_neg_edge(int eid) const {
        return is_neg_(sign_of((int)eid));
    }
    inline bool is_pos_edge(int u, int v) const {
        return neighbors_bm((size_t)u, +1).contains((size_t)v);
    }
    inline bool is_neg_edge(int u, int v) const {
        return neighbors_bm((size_t)u, -1).contains((size_t)v);
    }
	inline int edge_count(int sign) const { return (sign > 0) ? m_pos : m_neg; }
	// --- enumerate only negative edges via bitmaps (v > u to avoid duplicates) ---
	inline std::vector<Edge> get_edges(int sign) const { std::vector<Edge> out; out.reserve(std::max(1, edge_count(sign))); const int n = vertex_count(); for (int u = 0; u < n; ++u) { const auto& bm = neighbors_bm(static_cast<size_t>(u), sign); for (size_t vs : bm) { int v = static_cast<int>(vs); if (v > u) out.emplace_back(u, v); }} return out; }

	inline igraph_vector_int_t get_edge_eids(int sign) const {
		igraph_vector_int_t ids; igraph_vector_int_init(&ids, static_cast<int>(std::max(1, edge_count(sign))));
		const int n = vertex_count();
		int v;
		int i = 0;
		for (int u = 0; u < n; ++u) {
			const auto& bm = neighbors_bm(static_cast<size_t>(u), sign);
			for (size_t vs : bm) {
				v = static_cast<int>(vs);
				if (v < u) {
					VECTOR(ids)[i++] = _edge_index.at({u,v});
				}
				else break;
			}
		}
		return ids; // caller owns and must destroy
	}

    inline void print_info() const { std::cout << "Graph has " << vertex_count() << " vertices and " << edge_count() << " edges." << std::endl; }
};


// Effective sign under (integer) switching s ∈ {−1,+1}^V:
//  - sign_of(u,v,eid) overrides GraphCore and returns round_pm1( s[u]*sigma[eid]*s[v] ) ∈ {−1,+1}.
//  - is_pos_edge/is_neg_edge (eid) use virtual sign_of → they reflect CURRENT switching.
//  - on_switch_sign_changed_(u, from, to, ...):
//       * called BEFORE s[u] is assigned the new value.
//       * updates degree buckets d_pos/d_neg and bitmap rows N_pos_/N_neg_ to the POST-switch state.
//       * complexity O(deg(u)). No direct access to sigma_.
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
        };

		SwitchingVector(const SwitchingVector& src, SignedGraph* new_owner)
		: g_(new_owner), vals_(src.vals_) {}
		
		SwitchingVector(SwitchingVector&& src, SignedGraph* new_owner) noexcept
		: g_(new_owner), vals_(std::move(src.vals_)) {}
		
        SwitchingVector(SignedGraph* g, double init_val = 1.0)
            : g_(g), vals_(g->vertex_count(), init_val) {
        }
        
		// forbid plain copy/move to prevent accidental stale owner
		SwitchingVector(const SwitchingVector&) = delete;
		SwitchingVector& operator=(const SwitchingVector&) = delete;
		SwitchingVector(SwitchingVector&&) = delete;
		SwitchingVector& operator=(SwitchingVector&&) = delete;

		inline Ref operator[](int u) { return Ref{ this, u }; }
        inline double operator[](int u) const { return vals_[u]; }
        inline int size() const noexcept { return (int)vals_.size(); }
        inline const std::vector<double>& readonly() const noexcept { return vals_; }

		inline void flip_all_no_bitmaps() noexcept {
		    for (double& x : vals_) x = -x;
		}

    private:
        friend class SignedGraph;
        SignedGraph* g_ = nullptr;
        std::vector<double> vals_;

        inline static double clamp_(double x) {
            if (x < -1.0) return -1.0; if (x > 1.0) return 1.0; return x;
        }
		inline void assign_(int u, double x) {
		    if (vals_[u] == clamp_(x)) return;
	    	g_->on_vertex_flip_(u, vals_[u], clamp_(x));
		    vals_[u] = clamp_(x);
		}
		inline void neg_flip_batch_(const GraphCore::Bitmap& anchor) {
	    	g_->on_neg_flip_batch_(anchor);
	    	for (int u : anchor)
		    	vals_[u] = -vals_[u];
		}
    };

    SwitchingVector s_;        // size |V|, entries in [-1,1]

    // share topology, replace sigma (passed as ±1 doubles), and s
    SignedGraph(const std::shared_ptr<igraph_t> base, std::vector<double> new_sigma, std::vector<double> new_s);
    SignedGraph(const std::string& file_path);

    SignedGraph(const SignedGraph* const other);
    SignedGraph(const SignedGraph* const other, std::vector<double> new_sigma);
    SignedGraph(const SignedGraph* const other, std::vector<double> new_sigma, std::vector<double> new_s);

	// Helper: round a switching vector to ±1 in-place (0 → +1)
	inline double round_pm1_(double s) const {
	    return is_pos_(s) ? +1.0 : -1.0;
	}
	inline void round_pm1_inplace(std::vector<double>& s) const {
	    for (double& x : s) x = round_pm1_(x);
	}

    inline int sign_of(const int u, const int v, int eid) const override {
        return round_pm1_(s_[u] * GraphCore::sign_of(u, v, eid) * s_[v]);
    }

    std::vector<int> maximal_salience_clique_strict_pos(const std::vector<int>& Z0) const;
    bool has_fractional_switching(double eps = 1e-12) const;
    void single_switching(int u);
	// Called when s[u] changes sign (from and to have opposite signs).
	// snapshot before changing s[u]
	// Updates degrees and (+/−) bitmaps; O(deg(u)).
	virtual void on_vertex_flip_(int u, double from, double to);
	inline void neg_flip_batch_(const GraphCore::Bitmap& anchor) {
		s_.neg_flip_batch_(anchor); // C.3 helper
	}

    // compose switching (acts on s_ only; clamp to [-1,1])
    template<class Vec>
    void compose_switching_inplace(const Vec& s) {
        // Minimal & efficient: coalesce bitmap flips; guarded writes only flip when sign changes
        {
            const int n = vertex_count();
            for (int u = 0; u < n; ++u) {
                s_[u] *= static_cast<double>(s[u]);
            }
        }
	}

public:
    // enable copy/move semantics for efficient transfers
	SignedGraph(const SignedGraph& o)
	: GraphCore(o), s_(o.s_, this) {}
	
	SignedGraph& operator=(const SignedGraph& o) = delete;
	
	SignedGraph(SignedGraph&& o) noexcept
	: GraphCore(std::move(o)), s_(std::move(o.s_), this) {}
	
	SignedGraph& operator=(SignedGraph&& o) = delete;


    using GraphCore::is_pos_edge;
    using GraphCore::is_neg_edge;
    using GraphCore::sign_of;


    // === Override eid-based sign queries to reflect switching via bitmaps ===
    inline int sign_of(int eid) const override {
        igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
        return sign_of(u, v, eid);
    }


    std::unique_ptr<SignedGraph> clone() const;
    virtual ~SignedGraph();


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


