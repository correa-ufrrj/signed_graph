// signed_graph.cpp
#include "signed_graph.h"
#include <igraph/igraph_strvector.h>
#include <iostream>
#include <limits>
#include <boost/heap/pairing_heap.hpp>
#include <ostream>
#include <algorithm>
#include <cstdio>

#ifndef SG_DEBUG
#define SG_DEBUG 1   // flip to 0 to disable all probes
#endif

#if SG_DEBUG
  #include <chrono>
  #define SGLOG(...) do { std::fprintf(stderr, __VA_ARGS__); std::fprintf(stderr, "\n"); } while(0)
  struct SGTick { std::chrono::high_resolution_clock::time_point t0; const char* label;
    SGTick(const char* L): t0(std::chrono::high_resolution_clock::now()), label(L) {}
    ~SGTick(){ auto t1=std::chrono::high_resolution_clock::now();
      auto ms=std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
      std::fprintf(stderr,"[SGT] %s: %lld ms\n", label, (long long)ms);
    }
  };
  // --- forward declaration of the debug helper used earlier in this file ---
  static int cc_on_strict_pos(const igraph_t* g, const std::vector<double>& switched_weights,
                              std::vector<int>& comp_id, std::vector<int>& comp_size);
#else
  #define SGLOG(...)
  struct SGTick { SGTick(const char*){} };
#endif

// Single helper: numeric sign as +1, 0, -1 (used only where we truly want number sign)
static inline int sign_of_num(double w) noexcept {
    return (w > 0.0) - (w < 0.0);
}

std::ostream& operator<<(std::ostream& os, const Edge& e) {
    os << "(" << e.first << ", " << e.second << ")";
    return os;
}

double GraphCore::WeightsView::operator[](const Edge& e) const {
    igraph_integer_t eid;
    if (igraph_get_eid(g, &eid, e.first, e.second, 0, 0) != IGRAPH_SUCCESS) {
        throw std::out_of_range("Edge not found");
    }
    double w = (*weights)[eid];
    return w;
}

GraphCore::EdgeIndexesView::EdgeIndexesView(const MapType& m, const igraph_t* graph) : map(m), g(graph) {}

auto GraphCore::EdgeIndexesView::begin() const -> const_iterator {
    return const_iterator{g, 0};
}

auto GraphCore::EdgeIndexesView::end() const -> const_iterator {
    return const_iterator{g, igraph_ecount(g)};
}

auto GraphCore::EdgeIndexesView::const_iterator::operator*() const -> value_type {
    igraph_integer_t from, to;
    igraph_edge(g, eid, &from, &to);
    return {{static_cast<int>(from), static_cast<int>(to)}, static_cast<int>(eid)};
}

auto GraphCore::EdgeIndexesView::const_iterator::operator++() -> const_iterator& {
    ++eid;
    return *this;
}

bool GraphCore::EdgeIndexesView::const_iterator::operator!=(const const_iterator& other) const {
    return eid != other.eid;
}

int GraphCore::EdgeIndexesView::operator[](const Edge& e) const {
    auto it = map.find(e);
    if (it == map.end()) throw std::out_of_range("Edge not found");
    return it->second;
}

NegativeCycle::NegativeCycle(const Edge& neg_edge, std::vector<Edge>&& pos_edges)
    : neg_edge_(neg_edge), pos_edges_(std::move(pos_edges)) {
    if (pos_edges_.empty()) {
        throw std::runtime_error("Invalid negative cycle: no path found.");
    }
    if (!neg_edge_.is_adjacent_to(pos_edges_.front()) || !neg_edge_.is_adjacent_to(pos_edges_.back())) {
        std::cout << "Neg edge: " << neg_edge << "\nPos path: ";
        for (auto e : pos_edges_) std::cout << e << " ";
        std::cout << "\n";
        throw std::runtime_error("Invalid cycle: path endpoints do not match negative edge.");
    }
}

const Edge& NegativeCycle::neg_edge() const { return neg_edge_; }
const std::vector<Edge>& NegativeCycle::pos_edges() const { return pos_edges_; }

std::ostream& operator<<(std::ostream& os, const SignedEdge& e) {
    return os << e.points << ": " << e.sign;
}

std::ostream& operator<<(std::ostream& os, const NegativeCycle& nc) {
    os << "Negative Edge: " << nc.neg_edge() << "\n";
    os << "Positive Edges: ";
    for (const auto& e : nc.pos_edges()) {
        os << e;
    }
    return os;
}

// Share topology; copy trivial fields
GraphCore::GraphCore(const GraphCore* const other)
    : g_(other->g_),                    // share the igraph_t
      weights_(other->weights_),
      d_plus(other->d_plus),
      d_minus(other->d_minus),
      _edge_index(other->_edge_index) {}

// Share topology; replace weights; recompute degrees for new weights
GraphCore::GraphCore(const GraphCore* const other, std::vector<double> new_weights)
    : g_(other->g_),
      weights_(std::move(new_weights)),
      _edge_index(other->_edge_index)
{
    const auto m = igraph_ecount(g_.get());
    if ((int)weights_.size() != m)
        throw std::invalid_argument("Weight vector size does not match edge count.");
    compute_degrees();
}

// Share topology supplied by caller; replace weights
GraphCore::GraphCore(const std::shared_ptr<igraph_t> base,
                     std::vector<double> new_weights)
: g_(base),
  weights_(std::move(new_weights))
{
    const int m = igraph_ecount(g_.get());
    if ((int)weights_.size() != m)
        throw std::invalid_argument("Weight vector size does not match edge count.");

    // Build eid map (same topology ⇒ deterministic)
    _edge_index.clear();
    _edge_index.reserve(m);
    for (int eid = 0; eid < m; ++eid) {
        igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
        _edge_index.emplace(Edge{(int)u, (int)v}, eid);
    }

    compute_degrees();
}

// Build a new topology from file; allocate g_ in initializer list
GraphCore::GraphCore(const std::string& file_path)
    : g_(std::shared_ptr<igraph_t>(new igraph_t, IGraphDeleter{}))  // allocate BEFORE igraph_* calls
{
    std::ifstream infile(file_path);
    if (!infile.is_open()) {
        throw std::runtime_error("Error opening file: " + file_path);
    }

    igraph_set_attribute_table(&igraph_cattribute_table);

    std::vector<igraph_integer_t> edges_flat;
    std::string line;
    int max_vertex = 0;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string token;
        int u, v, sign;
        std::getline(ss, token, ','); u = std::stoi(token);
        std::getline(ss, token, ','); v = std::stoi(token);
        std::getline(ss, token, ','); sign = std::stoi(token);

        if (sign != 1 && sign != -1)
            throw std::runtime_error("Invalid sign in line: " + line);

        edges_flat.push_back(u);
        edges_flat.push_back(v);
        weights_.push_back(sign);
        max_vertex = std::max(max_vertex, std::max(u, v));
    }

    // Build topology in g_
    igraph_empty(g_.get(), max_vertex + 1, IGRAPH_UNDIRECTED);
    igraph_vector_int_t edges_vec;
    igraph_vector_int_view(&edges_vec, edges_flat.data(), (igraph_integer_t)edges_flat.size());
    igraph_add_edges(g_.get(), &edges_vec, nullptr);

    compute_degrees();

    // Build edge index
    _edge_index.clear();
    int id = 0;
    for (const auto& [e, _] : weights_view()) {
        _edge_index[e] = id++;
    }

#if SG_DEBUG
    {
        int m = edge_count();
        int pos=0, neg=0, zpos=0, zneg=0, other=0;
        for (int eid=0; eid<m; ++eid) {
            double w = weights_[eid];
            if (w == 0.0) { if (std::signbit(w)) ++zneg; else ++zpos; }
            else if (w > 0.0) ++pos;
            else if (w < 0.0) ++neg;
            else ++other;
        }
        SGLOG("[SG-PROBE] load: m=%d pos=%d neg=%d +0=%d -0=%d other=%d", m,pos,neg,zpos,zneg,other);
        SGLOG("[SG-PROBE] edge_index size = %zu (expect m)", _edge_index.size());
    }
#endif
}

SignedGraph::SignedGraph(const std::shared_ptr<igraph_t> base,
                         std::vector<double> new_weights,
                         std::vector<double> new_s)
: GraphCore(base, std::move(new_weights)),
  s_(std::move(new_s))
{
    const int n = igraph_vcount(base.get());
    if ((int)s_.size() != n)
        throw std::invalid_argument("Switching vector size does not match vertex count.");
    is_switching_ = true;
}

SignedGraph::SignedGraph(const SignedGraph* const other)
    : GraphCore(static_cast<const GraphCore* const>(other)),
      s_(other->s_),
      is_switching_(other->is_switching_) {}

SignedGraph::SignedGraph(const SignedGraph* const other, std::vector<double> new_weights)
    : GraphCore(static_cast<const GraphCore* const>(other), std::move(new_weights)),
      s_(other->s_),
      is_switching_(other->is_switching_) {}

SignedGraph::SignedGraph(const std::string& file_path)
  : GraphCore(file_path), s_(std::vector<double>(vertex_count(), 1.0)) {}
  
// in SignedGraph (and similarly already present in SignedGraphForMIP)
SignedGraph::SignedGraph(const SignedGraph* const other,
                         std::vector<double> new_weights,
                         std::vector<double> new_s)
: GraphCore(other, std::move(new_weights)),
  s_(std::move(new_s))
{
    if ((int)s_.size() != igraph_vcount(g_.get()))
        throw std::invalid_argument("Switching vector size does not match vertex count.");
    is_switching_ = true;
}

std::unique_ptr<SignedGraph> SignedGraph::clone() const {
    return std::make_unique<SignedGraph>(*this);
}

SignedGraph::~SignedGraph() {
}

int SignedGraph::frustrated_edges() const {
    int count = 0;
    auto sign = [&](int u){ return s_[u] >= 0 ? 1 : -1; };

    for (int eid = 0; eid < igraph_ecount(g_.get()); ++eid) {
        igraph_integer_t from, to;
        igraph_edge(g_.get(), eid, &from, &to);
        double w = weights_[eid];
        bool same_partition = sign(from) == sign(to);
        bool is_positive = (w >= 0.0 && !(w == 0.0 && std::signbit(w)));
        if ((is_positive && !same_partition) || (!is_positive && same_partition)) {
            count++;
        }
    }
    return count;
}

std::vector<Edge> SignedGraph::frustrated_edges_keys(const std::vector<int>& partition) const {
    std::vector<Edge> frustrated;
    for (int eid = 0; eid < igraph_ecount(g_.get()); ++eid) {
        igraph_integer_t from, to;
        igraph_edge(g_.get(), eid, &from, &to);
        double w = weights_[eid];
        bool same_partition = partition[from] == partition[to];
        bool is_positive = (w >= 0.0 && !(w == 0.0 && std::signbit(w)));
        if ((is_positive && !same_partition) || (!is_positive && same_partition)) {
            frustrated.emplace_back(from, to);
        }
    }
    return frustrated;
}

void SignedGraph::single_switching(int u, igraph_vector_int_t* incident) {
    const igraph_t* G = g_.get();
    // Ensure incident list is available
    if (incident == nullptr) {
        igraph_vector_int_t tmp; igraph_vector_int_init(&tmp, 0);
        single_switching(u, &tmp);
        igraph_vector_int_destroy(&tmp);
        return;
    }
    if (igraph_vector_int_size(incident) == 0) {
        igraph_incident(G, incident, u, IGRAPH_ALL);
    }

    // Flip vertex mask
    s_[u] = -s_[u];

    // All incident edges invert sign => swap degree buckets for u
    std::swap(d_plus[u], d_minus[u]);

    // Flip incident edges and update neighbor degree buckets
    for (int k = 0; k < (int)igraph_vector_int_size(incident); ++k) {
        const int eid = VECTOR(*incident)[k];
        const int v   = IGRAPH_OTHER(G, eid, u);
        const double w_old = weights_[eid];

        // Move contribution on neighbor between plus/minus buckets
        if (is_pos_edge(eid)) {
            d_plus[v]  -= w_old;
            d_minus[v] += w_old;
        } else {
            const double a = -w_old; // |w_old|
            d_minus[v] -= a;
            d_plus[v]  += a;
        }

        // Flip edge weight in working signature
        weights_[eid] = -w_old;
    }
}

void SignedGraph::single_switching(int u) {
    igraph_vector_int_t incident; igraph_vector_int_init(&incident, 0);
    single_switching(u, &incident);
    igraph_vector_int_destroy(&incident);
}

void SignedGraph::apply_greedy_switching(const GreedyKickOptions& opts)
{
    const igraph_t* G = g_.get();
    const int n = vertex_count();
    const double EPS = 1e-12;

    // initialize d~ from precomputed degree buffers
    std::vector<double> dtilde = net_degrees();

    // incumbent
    int best_mminus = negative_edge_count();
	int cur_mminus = best_mminus + 1;
    std::vector<double> best_s = s_;
    std::vector<double> best_w = weights_;

    // --- chooser for step B (weighted by salience). If frac_y==nullptr, fall back to vertex_salience ---
    auto pick_u_star_in_Q = [&](const std::vector<int>& Q) -> int {
        if (Q.empty()) return -1;
        igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
        double best_score = -std::numeric_limits<double>::infinity();
        int    best_u     = -1;
        for (int u : Q) {
            double score = 0.0;
            igraph_incident(G, &inc, u, IGRAPH_ALL);
            for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
                const int eid = VECTOR(inc)[ii];
                const int v   = IGRAPH_OTHER(G, eid, u);
                if (!is_pos_edge(u, v)) continue; // strictly positive only
                score += edge_salience(u, v);
            }
            igraph_vector_int_clear(&inc);
            if (score > best_score || (score == best_score && u < best_u)) {
                best_score = score; best_u = u;
            }
        }
        igraph_vector_int_destroy(&inc);
        return best_u;
    };

    struct Key { double d; double sal; int u; int ver; };
    auto cmp = [](const Key& a, const Key& b){
        if (a.d   != b.d  ) return a.d   > b.d;   // min-heap
        if (a.sal != b.sal) return a.sal < b.sal; // prefer higher salience
        return a.u > b.u;
    };

    const int R = std::max(1, opts.R_max);
    bool advanced = true;
    for (int r = 0; advanced && r < R; ++r) {
        advanced = false;

        // A) Fractional greedy flips (min-heap on d~, tie-break by salience)
        {
            std::priority_queue<Key, std::vector<Key>, decltype(cmp)> pq(cmp);
            std::vector<int> ver(n, 0);
            for (int u = 0; u < n; ++u) if (dtilde[u] < -EPS) pq.push({dtilde[u], vertex_salience(u), u, ver[u]});

            igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
            while (!pq.empty() && pq.top().d < -EPS) {
//				std::cout << "[GREEDY-PQ] top=" << pq.top().u
//						  << " top.d=" << pq.top().d
//						  << "\n";
                const int u = pq.top().u;
                const int veru = pq.top().ver;
				pq.pop();
                if (ver[u] != veru) continue; // stale key
//				std::cout << "[GREEDY-PQ] u=" << u
//						  << " d=" << dtilde[u]
//						  << "\n";

                flip_and_refresh_net_degree(u, dtilde);
                advanced = true;

                // refresh keys for u and its neighbors
//                igraph_vector_int_init(&inc, 0);
                igraph_incident(G, &inc, u, IGRAPH_ALL);
                ++ver[u];
//                pq.push({dtilde[u], vertex_salience(u), u, ver[u]});
                for (int k = 0; k < (int)igraph_vector_int_size(&inc); ++k) {
                    const int v = IGRAPH_OTHER(G, VECTOR(inc)[k], u);
                    ++ver[v];
                    if (dtilde[v] < -EPS) {
	                    pq.push({dtilde[v], vertex_salience(v), v, ver[v]});
                    }
                }
            }
            igraph_vector_int_destroy(&inc);
            if (advanced) {
	            cur_mminus = negative_edge_count();
	            if (cur_mminus <= best_mminus) { best_mminus = cur_mminus; best_s = s_; best_w = weights_; }
            }
        }
//        if (advanced) continue;

        // B) Zero-clique step (strictly positive edges; weighted by salience; cap K_max)
        if (opts.K_max > 0) {
            std::vector<int> Z0; Z0.reserve(n);
            for (int u = 0; u < n; ++u) if (std::fabs(dtilde[u]) <= EPS) Z0.push_back(u);
            
            std::cout << "[GREEDY-B] |Z0|=" << Z0.size()
            		  << " cur_minus=" << cur_mminus
            		  << " best_minus=" << best_mminus
            		  << " project=" << integer_projection().negative_edge_count()
            		  << "\n";

            if (!Z0.empty()) {
                auto Q = maximal_salience_clique_strict_pos(Z0);
            
	            std::cout << "[GREEDY-B] |Q|=" << Q.size()
	            		  << "\n";

                if ((int)Q.size() >= 2) {
		            for (int v: Q) {
						flip_and_refresh_net_degree(v, dtilde);
					}
                    advanced = true;
//                    const int ustar = pick_u_star_in_Q(Q);
//                    if (ustar >= 0) {
//                        flip_and_refresh_net_degree(ustar, dtilde);
//                        advanced = true;
//                    }
                }
            }
        }
        if (advanced) continue;

        // C) Integer detour (only if s is nonintegral; length ≤ L_max; gate d~(u*) ≤ Δ)
        if (opts.L_max > 0 && has_fractional_switching()) {
            // 1) integer projection and a short integer greedy pass there
            SignedGraph SGint = integer_projection();
            GreedyKickOptions int_opts = opts; int_opts.R_max = 1; int_opts.K_max = 0; int_opts.L_max = 0;
            SGint.apply_greedy_switching(int_opts); // SGint has integer s_, so this naturally skips its own Step C

            // 2) compute d_int on SGint from its degree buffers
            std::vector<double> d_int(SGint.vertex_count(), 0.0);
            for (int u = 0; u < SGint.vertex_count(); ++u) d_int[u] = SGint.d_plus[u] - SGint.d_minus[u];

            // 3) short replay sequence S from a weighted Z0-int clique (strictly positive edges in SGint)
            std::vector<int> S; S.reserve(opts.L_max);
            for (int t = 0; t < opts.L_max; ++t) {
                std::vector<int> Z0i; Z0i.reserve(n);
                for (int u = 0; u < n; ++u) if (std::fabs(d_int[u]) <= EPS) Z0i.push_back(u);
                if (Z0i.empty()) break;
                auto Qi = SGint.maximal_salience_clique_strict_pos(Z0i);
                if ((int)Qi.size() < 2) break;
                const int ustar = Qi.front();

                // gate on CURRENT fractional net degree
                if (dtilde[ustar] > (double)opts.Delta) break;

                // flip u* in SGint and refresh d_int locally
                igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
                S.push_back(ustar);
                SGint.single_switching(ustar, &inc);
                d_int[ustar] = SGint.d_plus[ustar] - SGint.d_minus[ustar];
                for (int k = 0; k < (int)igraph_vector_int_size(&inc); ++k) {
                    const int v  = IGRAPH_OTHER(G, VECTOR(inc)[k], ustar);
                    d_int[v] = SGint.d_plus[v] - SGint.d_minus[v];
                }
                igraph_vector_int_destroy(&inc);
            }

            // 4) Replay S on the working state; accept iff any d~ becomes negative
            bool created_neg = false; std::vector<int> replayed; replayed.reserve(S.size());
            for (int u : S) {
                flip_and_refresh_net_degree(u, dtilde);
                replayed.push_back(u);
                if (dtilde[u] < -EPS) { created_neg = true; break; }
                igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
                igraph_incident(G, &inc, u, IGRAPH_ALL);
                for (int k = 0; k < (int)igraph_vector_int_size(&inc); ++k) {
                    const int v = IGRAPH_OTHER(G, VECTOR(inc)[k], u);
                    if (dtilde[v] < -EPS) { created_neg = true; break; }
                }
                igraph_vector_int_destroy(&inc);
                if (created_neg) break;
            }

            if (created_neg) {
                advanced = true;
                const int cur_mminus = negative_edge_count();
                if (cur_mminus < best_mminus) { best_mminus = cur_mminus; best_s = s_; best_w = weights_; }
            } else {
                for (int i = (int)replayed.size() - 1; i >= 0; --i) flip_and_refresh_net_degree(replayed[i], dtilde);
            }
        }

        if (!advanced) break;
    }

    // commit incumbent
    s_ = std::move(best_s);
    weights_ = std::move(best_w);
    compute_degrees();
    is_switching_ = true;
}

SignedGraph SignedGraph::greedy_switching(const GreedyKickOptions& opts) {
    SignedGraph sg(g_, weights_, s_);
    sg.is_switching_ = true;
    sg.apply_greedy_switching(opts);
    return sg;
}

SignedGraph SignedGraph::greedy_switching() {
    GreedyKickOptions opts; return greedy_switching(opts);
}

//SignedGraph SignedGraph::integer_projection() const {
//    const int n = vertex_count();
//    const int m = edge_count();
//
//    std::vector<double> s_int(n);
//    for (int u = 0; u < n; ++u) {
//        s_int[u] = (s_[u] >= 0.0) ? +1.0 : -1.0; // round to {±1}
//    }
//
//	const auto s = signs_view();
//    std::vector<double> sigma_int(m);
//    for (int eid = 0; eid < m; ++eid) {
//		sigma_int[eid] = s[eid].sign; // σ_{s_int} = s_int ⊙ σ_s ⊙ s_int
//    }
//
//    SignedGraph sg(this, sigma_int, s_int);
//    sg.is_switching_ = true;
//    return sg;
//}
//// ---------- SignedGraph ----------
SignedGraph SignedGraph::integer_projection() const {
    const int n = vertex_count();
    const int m = edge_count();

    // 1) round node mask to ±1
    std::vector<double> s_int(n);
    for (int u = 0; u < n; ++u)
        s_int[u] = (s_[u] >= 0.0) ? +1.0 : -1.0;

    // 2) recompute edge signs from CURRENT weights_ by sign-only switching
    std::vector<double> w_proj(weights_);
    for (int eid = 0; eid < m; ++eid) {
        igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
        w_proj[eid] = w_proj[eid] * s_int[(int)u] * s_int[(int)v];
        // w_proj stays in {±1} if weights_ were in {±1}
    }

    // 3) return a graph with s_=s_int and σ = s_int ⊙ σ_current ⊙ s_int
    SignedGraph out(this, std::move(w_proj), std::move(s_int));
    out.is_switching_ = true;
    return out;
}


bool SignedGraph::are_cycles_edge_disjoint(const std::vector<NegativeCycle>& cycles) const {
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

bool SignedGraph::are_cycles_sign_correct(const std::vector<NegativeCycle>& cycles, bool expect_negative) const {
    auto signs_v = signs_view();
    for (const auto& cycle : cycles) {
        int product = signs_v[cycle.neg_edge()].sign;
        for (const auto& e : cycle.pos_edges()) {
            product *= signs_v[e].sign;
        }
        bool is_negative = (product < 0);
        if (is_negative != expect_negative) {
            std::cerr << "[Error] Unexpected sign in cycle: "
                      << "neg_edge = (" << cycle.neg_edge().first << ", " << cycle.neg_edge().second << "), "
                      << "product = " << product << "\n";
            return false;
        }
    }
    return true;
}

std::vector<int> SignedGraph::maximal_salience_clique_strict_pos(const std::vector<int>& Z0) const {
    if (Z0.empty()) return {};
    std::vector<int> cand = Z0;
    std::sort(cand.begin(), cand.end(), [&](int a, int b){
        const double sa = const_cast<SignedGraph*>(this)->vertex_salience(a);
        const double sb = const_cast<SignedGraph*>(this)->vertex_salience(b);
        if (sa != sb) return sa > sb;
        return a < b;
    });
    std::vector<int> Q; Q.reserve(cand.size());
    Q.push_back(cand.front());
    for (size_t i = 1; i < cand.size(); ++i) {
        int u = cand[i]; bool ok = true;
        for (int w : Q) { if (!is_pos_edge(u, w)) { ok = false; break; } }
        if (ok) Q.push_back(u);
    }
    return Q;
}

bool SignedGraph::has_fractional_switching(double eps) const {
    for (double su : s_) if (std::fabs(su) < 1.0 - eps) return true;
    return false;
}
void SignedGraph::flip_and_refresh_net_degree(int u, std::vector<double>& dtilde) {
    const igraph_t* G = g_.get();
    igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
    igraph_incident(G, &inc, u, IGRAPH_ALL);
    std::vector<int> nbrs; nbrs.reserve((int)igraph_vector_int_size(&inc));
    for (int k = 0; k < (int)igraph_vector_int_size(&inc); ++k) {
        const int eid = VECTOR(inc)[k];
        const int v   = IGRAPH_OTHER(G, eid, u);
        nbrs.push_back(v);
    }
    single_switching(u, &inc); // updates s_, weights_, d_plus/d_minus

    dtilde[u] = net_degree(u);
    for (int v : nbrs) dtilde[v] = net_degree(v);

    igraph_vector_int_destroy(&inc);
}

void SignedGraph::save_partition_svg(const std::vector<int>& partition, const std::string& filename, bool custom_layout) const {
    igraph_vector_ptr_t v_attrs;
    igraph_vector_ptr_init(&v_attrs, 0);
    igraph_vector_ptr_t e_attrs;
    igraph_vector_ptr_init(&e_attrs, 0);

    igraph_strvector_t colors;
    igraph_strvector_init(&colors, 2);
    igraph_strvector_set(&colors, 0, "red");
    igraph_strvector_set(&colors, 1, "blue");

    igraph_vector_t color_attr;
    igraph_vector_init(&color_attr, igraph_vcount(g_.get()));
    for (int i = 0; i < igraph_vcount(g_.get()); ++i) {
        VECTOR(color_attr)[i] = partition[i];
    }
    igraph_cattribute_VAN_setv(const_cast<igraph_t*>(g_.get()), "color", &color_attr);
    igraph_vector_destroy(&color_attr);

    igraph_strvector_t ecolors;
    igraph_strvector_init(&ecolors, igraph_ecount(g_.get()));

    int i = 0;
    for (const auto& [edge, sign] : signs_view()) {
        bool is_frustrated = (sign == 1 && partition[edge.first] != partition[edge.second]) ||
                             (sign == -1 && partition[edge.first] == partition[edge.second]);
        igraph_strvector_set(&ecolors, i++, is_frustrated ? "red" : "blue");
    }
    igraph_cattribute_EAS_setv(const_cast<igraph_t*>(g_.get()), "color", &ecolors);

    igraph_t gcopy;
    igraph_copy(&gcopy, g_.get());

    igraph_matrix_t layout;
    if (custom_layout) {
        igraph_matrix_init(&layout, igraph_vcount(g_.get()), 2);
        int above = 0, below = 0;
        for (int i2 = 0; i2 < igraph_vcount(g_.get()); ++i2) {
            double x = (partition[i2] == 0 ? above++ : below++);
            double y = (partition[i2] == 0 ? 1.0 : 0.0);
            MATRIX(layout, i2, 0) = x;
            MATRIX(layout, i2, 1) = y;
        }
    }
    std::string dot_filename = filename;
    if (dot_filename.size() > 4 && dot_filename.substr(dot_filename.size() - 4) == ".svg") {
        dot_filename = dot_filename.substr(0, dot_filename.size() - 4) + ".dot";
    }
    std::cout << "[INFO] Writing graph to: " << dot_filename << std::endl;
    std::cout << "[INFO] To convert to SVG: dot -Tsvg " << dot_filename << " -o " << dot_filename.substr(0, dot_filename.size() - 4) + ".svg" << std::endl;
    FILE* fout = fopen(filename.c_str(), "w");
    if (custom_layout) {
        igraph_write_graph_dot(&gcopy, fout);
        igraph_matrix_destroy(&layout);
    } else {
        igraph_write_graph_dot(&gcopy, fout);
    }
    fclose(fout);
    igraph_destroy(&gcopy);

    igraph_strvector_destroy(&colors);
    igraph_strvector_destroy(&ecolors);
    igraph_vector_ptr_destroy_all(&v_attrs);
    igraph_vector_ptr_destroy_all(&e_attrs);
}

namespace { inline bool __is_pos_weight(double w) { return (w > 0.0) || (w == 0.0 && !std::signbit(w)); } }

std::vector<int> SignedGraph::negative_triangle_count_per_vertex() const {
    const int n = igraph_vcount(g_.get());
    const int m = igraph_ecount(g_.get());
    std::vector<int> tri(n, 0);
    if (n == 0 || m == 0) return tri;

    std::vector<std::vector<int>> adj(n);
    for (igraph_integer_t eid = 0; eid < m; ++eid) {
        igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
        int a = static_cast<int>(u), b = static_cast<int>(v);
        adj[a].push_back(b);
        adj[b].push_back(a);
    }
    for (int u = 0; u < n; ++u) std::sort(adj[u].begin(), adj[u].end());

    std::vector<char> mark(n, 0);

    for (int u = 0; u < n; ++u) {
        for (int v : adj[u]) mark[v] = 1;
        for (int v : adj[u]) if (v > u) {
            for (int w : adj[v]) if (w > v && mark[w]) {
                Edge uv(u,v), vw(v,w), uw(u,w);
                auto it_uv = _edge_index.find(uv);
                auto it_vw = _edge_index.find(vw);
                auto it_uw = _edge_index.find(uw);
                if (it_uv == _edge_index.end() || it_vw == _edge_index.end() || it_uw == _edge_index.end()) continue;
                const double w_uv = weights_[it_uv->second];
                const double w_vw = weights_[it_vw->second];
                const double w_uw = weights_[it_uw->second];
                const int s_uv = __is_pos_weight(w_uv) ? 1 : -1;
                const int s_vw = __is_pos_weight(w_vw) ? 1 : -1;
                const int s_uw = __is_pos_weight(w_uw) ? 1 : -1;
                if (s_uv * s_vw * s_uw < 0) { ++tri[u]; ++tri[v]; ++tri[w]; }
            }
        }
        for (int v : adj[u]) mark[v] = 0;
    }
    return tri;
}

int SignedGraph::negative_triangle_count_of_vertex(int u) const {
    auto vec = negative_triangle_count_per_vertex();
    if (u < 0 || u >= (int)vec.size()) return 0;
    return vec[u];
}

#if SG_DEBUG
// --- definition for the forward-declared helper ---
static int cc_on_strict_pos(const igraph_t* g, const std::vector<double>& switched_weights,
                            std::vector<int>& comp_id, std::vector<int>& comp_size) {
    const int n = igraph_vcount(g);
    const int m = igraph_ecount(g);
    comp_id.assign(n, -1);
    comp_size.clear();
    std::vector<std::vector<int>> adj(n);
    for (int eid=0; eid<m; ++eid) {
        if (switched_weights[eid] > 0.0) {
            igraph_integer_t u,v; igraph_edge(g, eid, &u, &v);
            adj[(int)u].push_back((int)v);
            adj[(int)v].push_back((int)u);
        }
    }
    int cid=0;
    std::vector<int> st; st.reserve(n);
    for (int s=0; s<n; ++s) if (comp_id[s] < 0) {
        comp_size.push_back(0);
        st.clear(); st.push_back(s);
        comp_id[s]=cid;
        while(!st.empty()){
            int u=st.back(); st.pop_back();
            ++comp_size.back();
            for (int w: adj[u]) if (comp_id[w] < 0) { comp_id[w]=cid; st.push_back(w); }
        }
        ++cid;
    }
    return cid;
}
#endif
