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

SignedEdge SignedGraph::SignedEdgesView::operator[](igraph_integer_t eid) const {
    double w = (*weights)[eid];
    int sign = (w == 0.0 && std::signbit(w)) ? -1 : (w >= 0.0 ? 1 : -1);
    igraph_integer_t from, to;
    igraph_edge(g, eid, &from, &to);
    return SignedEdge{Edge{static_cast<int>(from), static_cast<int>(to)}, sign};
}

SignedEdge SignedGraph::SignedEdgesView::operator[](const Edge& e) const {
    igraph_integer_t eid;
    if (igraph_get_eid(g, &eid, e.first, e.second, 0, 0) != IGRAPH_SUCCESS) {
        throw std::out_of_range("Edge not found");
    }
    return (*this)[eid];
}

double SignedGraph::WeightsView::operator[](const Edge& e) const {
    igraph_integer_t eid;
    if (igraph_get_eid(g, &eid, e.first, e.second, 0, 0) != IGRAPH_SUCCESS) {
        throw std::out_of_range("Edge not found");
    }
    double w = (*weights)[eid];
    return w;
}

SignedGraph::EdgeIndexesView::EdgeIndexesView(const MapType& m, const igraph_t* graph) : map(m), g(graph) {}

auto SignedGraph::EdgeIndexesView::begin() const -> const_iterator {
    return const_iterator{g, 0};
}

auto SignedGraph::EdgeIndexesView::end() const -> const_iterator {
    return const_iterator{g, igraph_ecount(g)};
}

auto SignedGraph::EdgeIndexesView::const_iterator::operator*() const -> value_type {
    igraph_integer_t from, to;
    igraph_edge(g, eid, &from, &to);
    return {{static_cast<int>(from), static_cast<int>(to)}, static_cast<int>(eid)};
}

auto SignedGraph::EdgeIndexesView::const_iterator::operator++() -> const_iterator& {
    ++eid;
    return *this;
}

bool SignedGraph::EdgeIndexesView::const_iterator::operator!=(const const_iterator& other) const {
    return eid != other.eid;
}

int SignedGraph::EdgeIndexesView::operator[](const Edge& e) const {
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

SignedGraph::SignedGraph(const SignedGraph* const other) {
    igraph_copy(&g, &other->g);
    weights = other->weights;
    switched_weights = other->switched_weights;
    d_plus = other->d_plus;
    d_minus = other->d_minus;
    switched_d_plus = other->switched_d_plus;
    switched_d_minus = other->switched_d_minus;
    _edge_index = other->_edge_index;
};

SignedGraph::SignedGraph(const SignedGraph* const other, std::vector<double> new_weights) {
    igraph_copy(&g, &other->g);
    if (new_weights.size() != static_cast<size_t>(igraph_ecount(&g))) {
        throw std::invalid_argument("Weight vector size does not match edge count in graph.");
    }
    weights = std::move(new_weights);
    switched_weights = other->switched_weights;
    d_plus = other->d_plus;
    d_minus = other->d_minus;
    switched_d_plus = other->switched_d_plus;
    switched_d_minus = other->switched_d_minus;
    _edge_index = other->_edge_index;
}

SignedGraph::SignedGraph(const std::string& file_path) {
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
        weights.push_back(sign);
        max_vertex = std::max(max_vertex, std::max(u, v));
    }

    igraph_empty(&g, max_vertex + 1, IGRAPH_UNDIRECTED);
    igraph_vector_int_t edges_vec;
    igraph_vector_int_view(&edges_vec, edges_flat.data(), edges_flat.size());
    igraph_add_edges(&g, &edges_vec, nullptr);

    switched_weights = weights;

    compute_degrees();
    switched_d_plus = d_plus;
    switched_d_minus = d_minus;

    _edge_index.clear();
    int id = 0;
    for (const auto& [e, _] : signs_view()) {
        _edge_index[e] = id++;
    }

#if SG_DEBUG
    {
        int m = igraph_ecount(&g);
        int pos=0, neg=0, zpos=0, zneg=0, other=0;
        for (int eid=0; eid<m; ++eid) {
            double w = weights[eid];
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

std::unique_ptr<SignedGraph> SignedGraph::clone() const {
    return std::make_unique<SignedGraph>(*this);
}

SignedGraph::~SignedGraph() {
    igraph_destroy(&g);
}

const std::vector<int> SignedGraph::neighbors(int u) const {
    igraph_vector_int_t neighbors;
    igraph_vector_int_init(&neighbors, 0);
    igraph_neighbors(&g, &neighbors, u, IGRAPH_ALL);

    std::vector<int> result(VECTOR(neighbors), VECTOR(neighbors) + igraph_vector_int_size(&neighbors));
    igraph_vector_int_destroy(&neighbors);
    return result;
}

void SignedGraph::compute_degrees() {
    int n = igraph_vcount(&g);
    d_plus.assign(n, 0.0);
    d_minus.assign(n, 0.0);

    for (int eid = 0; eid < igraph_ecount(&g); ++eid) {
        igraph_integer_t from, to;
        igraph_edge(&g, eid, &from, &to);
        double w = weights[eid];
        if (w >= 0.0 && !(w == 0.0 && std::signbit(w))) {
            d_plus[from] += w;
            d_plus[to] += w;
        } else {
            d_minus[from] += -w;
            d_minus[to] += -w;
        }
    }
}

inline bool SignedGraph::strictly_positive(double w) { return w > 0.0; }
inline double SignedGraph::vertex_salience(double x) {
    double a = std::fabs(x - 0.5) * 2.0;
    if (a > 1.0) a = 1.0;
    return 1.0 - a;
}
inline double SignedGraph::edge_salience_from_y(double y) {
    double a = std::fabs(y - 0.5) * 2.0;
    if (a > 1.0) a = 1.0;
    return 1.0 - a;
}
inline double SignedGraph::tau_from_x(double xu, double xv) {
    return 4.0 * xu * xv - 2.0 * xu - 2.0 * xv + 1.0;
}

// Maximal clique inside Z0 using strictly positive edges from switched_frac_edge_sign
std::vector<int> SignedGraph::maximal_salience_clique_strict_pos(
    const std::vector<int>& z0_vertices,
    const std::vector<char>& /*inZ0*/,
    const std::vector<double>& switched_frac_edge_sign,
    const std::vector<double>* frac_x
) const {
    std::vector<int> cand = z0_vertices;
    std::sort(cand.begin(), cand.end(), [&](int a, int b){
        double sa = frac_x ? vertex_salience((*frac_x)[a]) : 0.0;
        double sb = frac_x ? vertex_salience((*frac_x)[b]) : 0.0;
        if (sa != sb) return sa > sb;
        return a < b;
    });

    std::vector<int> Q;
    Q.reserve(cand.size());
    for (int u : cand) {
        bool ok = true;
        for (int v : Q) {
            igraph_integer_t eid;
            if (igraph_get_eid(const_cast<igraph_t*>(&g), &eid, u, v, /*directed=*/0, /*error=*/0) != IGRAPH_SUCCESS) {
                ok = false; break;
            }
            double w = switched_frac_edge_sign[(int)eid];
            if (!strictly_positive(w)) { ok = false; break; }
        }
        if (ok) Q.push_back(u);
    }
    if ((int)Q.size() < 2) Q.clear();
    return Q;
}

int SignedGraph::pick_u_star_in_Q(
    const std::vector<int>& Q,
    const std::vector<double>& switched_frac_edge_sign,
    const std::vector<double>* frac_y
) const {
    if (Q.empty()) return -1;

    igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
    double best_score = -1.0;
    int best_u = -1;

    for (int u : Q) {
        igraph_incident(const_cast<igraph_t*>(&g), &inc, u, IGRAPH_ALL);
        double score = 0.0;
        for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
            int eid = VECTOR(inc)[ii];
            double w = switched_frac_edge_sign[eid];
            if (!strictly_positive(w)) continue;
            if (frac_y) {
                score += edge_salience_from_y((*frac_y)[eid]);
            }
        }
        if (score > best_score) { best_score = score; best_u = u; }
    }
    igraph_vector_int_destroy(&inc);
    return best_u;
}

// REFAC: project x to s_int ∈ {±1} and compute sigma_int_edge_sign[eid] = s_int[u]*base_edge_sign[eid]*s_int[v]
void SignedGraph::integer_projection_from_x(
    const std::vector<double>& x,
    const std::vector<double>& base_edge_sign,
    std::vector<int>& s_int,
    std::vector<double>& sigma_int_edge_sign) const
{
    const int n = igraph_vcount(&g);
    const int m = igraph_ecount(&g);

    s_int.resize(n);
    for (int u = 0; u < n; ++u) {
        // s = 1 - 2*round(x)  (x>=.5 -> s=-1; else +1)
        s_int[u] = (x[u] >= 0.5) ? -1 : +1;
    }

    sigma_int_edge_sign.resize(m);
    for (int eid = 0; eid < m; ++eid) {
        igraph_integer_t uu, vv;
        igraph_edge(&g, eid, &uu, &vv);
        int u = static_cast<int>(uu), v = static_cast<int>(vv);
        // base_edge_sign is the current working signature at eid (±1)
        sigma_int_edge_sign[eid] =
            static_cast<double>(s_int[u]) * base_edge_sign[eid] * static_cast<double>(s_int[v]);
    }
}

void SignedGraph::integer_greedy_pass_minheap(
    std::vector<double>& sigma_int_edge_sign,
    std::vector<int>& s_int,
    std::vector<int>& d_int) const
{
    const int n = (int)s_int.size();
    d_int.assign(n, 0);

    igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
    for (int u = 0; u < n; ++u) {
        igraph_incident(const_cast<igraph_t*>(&g), &inc, u, IGRAPH_ALL);
        int sum = 0;
        for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
            int eid = VECTOR(inc)[ii];
            sum += (sigma_int_edge_sign[eid] > 0.0) ? 1 : -1;
        }
        d_int[u] = sum;
    }

    struct Node { int u; int key; };
    auto cmp = [](const Node& a, const Node& b){ return a.key > b.key; };
    std::priority_queue<Node,std::vector<Node>,decltype(cmp)> pq(cmp);
    for (int u = 0; u < n; ++u) pq.push({u, d_int[u]});

    while (!pq.empty() && pq.top().key < 0) {
        int u = pq.top().u; pq.pop();
        if (d_int[u] >= 0) continue;

        s_int[u] = -s_int[u];

        igraph_incident(const_cast<igraph_t*>(&g), &inc, u, IGRAPH_ALL);
        d_int[u] = -d_int[u];

        for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
            int eid = VECTOR(inc)[ii];
            int v = IGRAPH_OTHER(const_cast<igraph_t*>(&g), eid, u);

            sigma_int_edge_sign[eid] = -sigma_int_edge_sign[eid];
            if (sigma_int_edge_sign[eid] > 0.0) {
                d_int[v] += 2;
            } else {
                d_int[v] -= 2;
            }
            pq.push({v, d_int[v]});
        }
        pq.push({u, d_int[u]});
    }
    igraph_vector_int_destroy(&inc);
}

const SignedGraph::EdgeIndexesView SignedGraph::edge_index() const {
    return SignedGraph::EdgeIndexesView{_edge_index, &g};
}

const std::vector<double>& SignedGraph::plus_degrees_view() const { return d_plus; }
const std::vector<double>& SignedGraph::minus_degrees_view() const { return d_minus; }

const long double SignedGraph::max_degree_vertex() const {
    int n = vertex_count();

    igraph_vector_int_t degrees;
    igraph_real_t max_degree = 0;
    long int max_vertex = 0;

    igraph_vector_int_init(&degrees, n);
    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

    for (igraph_integer_t i = 0; i < n; ++i) {
        igraph_real_t deg = VECTOR(degrees)[i];
        if (deg > max_degree) {
            max_degree = deg;
            max_vertex = i;
        }
    }

    igraph_vector_int_destroy(&degrees);
    return max_vertex;
}

void SignedGraph::single_switching(int u, igraph_vector_int_t* incident) {
    igraph_incident(&g, incident, u, IGRAPH_ALL);

    for (int i = 0; i < igraph_vector_int_size(incident); ++i) {
        int eid = VECTOR(*incident)[i];
        double& w = switched_weights[eid];
        w = -w; // toggle sign (keeps ±0 convention)
    }

    std::swap(switched_d_plus[u], switched_d_minus[u]);

    for (int i = 0; i < igraph_vector_int_size(incident); ++i) {
        int eid = VECTOR(*incident)[i];
        int v = IGRAPH_OTHER(&g, eid, u);
        double w = switched_weights[eid];
        if (w >= 0.0 && !(w == 0.0 && std::signbit(w))) {
            switched_d_plus[v] += w;
            switched_d_minus[v] -= w;
        } else {
            switched_d_plus[v] -= -w;
            switched_d_minus[v] += -w;
        }
    }
}

int SignedGraph::frustrated_edges(const std::vector<int>& partition) const {
    int count = 0;
    for (int eid = 0; eid < igraph_ecount(&g); ++eid) {
        igraph_integer_t from, to;
        igraph_edge(&g, eid, &from, &to);
        double w = weights[eid];
        bool same_partition = partition[from] == partition[to];
        bool is_positive = (w >= 0.0 && !(w == 0.0 && std::signbit(w)));
        if ((is_positive && !same_partition) || (!is_positive && same_partition)) {
            count++;
        }
    }
    return count;
}

std::vector<Edge> SignedGraph::frustrated_edges_keys(const std::vector<int>& partition) const {
    std::vector<Edge> frustrated;
    for (int eid = 0; eid < igraph_ecount(&g); ++eid) {
        igraph_integer_t from, to;
        igraph_edge(&g, eid, &from, &to);
        double w = weights[eid];
        bool same_partition = partition[from] == partition[to];
        bool is_positive = (w >= 0.0 && !(w == 0.0 && std::signbit(w)));
        if ((is_positive && !same_partition) || (!is_positive && same_partition)) {
            frustrated.emplace_back(from, to);
        }
    }
    return frustrated;
}

void SignedGraph::switching_from_partition(const std::vector<int>& s) {
    igraph_integer_t edge_count = igraph_ecount(&g);

    switched_weights = weights;
    switched_d_plus.assign(d_plus.begin(), d_plus.end());
    switched_d_minus.assign(d_minus.begin(), d_minus.end());

    for (igraph_integer_t eid = 0; eid < edge_count; ++eid) {
        igraph_integer_t from, to;
        igraph_edge(&g, eid, &from, &to);

        if (s[from] != s[to]) {
            double& w = switched_weights[eid];
            w = -w;

            if (w >= 0.0 && !(w == 0.0 && std::signbit(w))) {
                switched_d_plus[from] += w;
                switched_d_plus[to] += w;
                switched_d_minus[from] -= w;
                switched_d_minus[to] -= w;
            } else {
                switched_d_plus[from] -= -w;
                switched_d_plus[to] -= -w;
                switched_d_minus[from] += -w;
                switched_d_minus[to] += -w;
            }
        }
    }
}

void SignedGraph::single_switching(int u) {
    igraph_vector_int_t incident;
    igraph_vector_int_init(&incident, 0);
    single_switching(u, &incident);
    igraph_vector_int_destroy(&incident);
}

void SignedGraph::restore_switched_sign() {
    switched_weights = weights;
    switched_d_plus = d_plus;
    switched_d_minus = d_minus;
}

std::vector<int> SignedGraph::greedy_switching_base(
    const std::function<int(int,int)>& /*cmp_fn*/,
    const GreedyKickOptions& opts)
{
    const int n = vertex_count();
    const int m = igraph_ecount(&g);

    // Working switching s (start at +1)
    std::vector<int> s(n, 1);

    // Edge signs under current switching: initialize from switched_weights sign
    std::vector<double> sigma_s_edge_sign(m, 1.0);
    for (int eid = 0; eid < m; ++eid) {
        double w = switched_weights[eid];
        sigma_s_edge_sign[eid] = static_cast<double>( sign_of_num(w) );
    }

#if SG_DEBUG
    {
        std::vector<int> comp_id, comp_sz;
        int comps = cc_on_strict_pos(&g, switched_weights, comp_id, comp_sz);
        int nloc = vertex_count();
        long long sum=0; int biggest=0;
        for (int szv: comp_sz){ sum+=szv; biggest=std::max(biggest,szv); }
        SGLOG("[SG-PROBE] +edge CCs=%d, n=%d, largest=%d", comps, nloc, biggest);
        if (comps>1) {
            int small=0; for (int szv: comp_sz) if (szv<10) ++small;
            SGLOG("[SG-PROBE] +edge CC small comps(<10): %d", small);
        }
    }
#endif

#if SG_DEBUG
    {
        // verify sigma_s_edge_sign matches switched_weights’ sign
        int ml = igraph_ecount(&g);
        int mismatch=0; int neg_edges=0;
        for (int eid=0; eid<ml; ++eid) {
            int a = (switched_weights[eid] > 0.0) ? +1 : (switched_weights[eid] < 0.0 ? -1 : 0);
            int b = (sigma_s_edge_sign[eid] > 0.0) ? +1 : (sigma_s_edge_sign[eid] < 0.0 ? -1 : 0);
            if (a != b) ++mismatch;
            if (b < 0) ++neg_edges;
        }
        SGLOG("[SG-PROBE] init sigma_s check: mismatch=%d, negE=%d", mismatch, neg_edges);
    }
#endif

    const bool has_frac = (opts.frac_x != nullptr);
    const auto* X = opts.frac_x;
    const auto* Y = opts.frac_y;

    // tilde_sigma[eid] = sigma_s(uv) * tau(x_u, x_v)
    std::vector<double> tilde_sigma(m, 0.0);

    auto rebuild_tilde_sigma = [&](){
        if (!has_frac) {
            for (int eid = 0; eid < m; ++eid) tilde_sigma[eid] = sigma_s_edge_sign[eid];
            return;
        }
        for (int eid = 0; eid < m; ++eid) {
            igraph_integer_t uu, vv; igraph_edge(&g, eid, &uu, &vv);
            int u = static_cast<int>(uu), v = static_cast<int>(vv);
            double tau = tau_from_x((*X)[u], (*X)[v]);
            tilde_sigma[eid] = sigma_s_edge_sign[eid] * tau;
        }
    };

    // fractional net degrees
    std::vector<double> dtilde(n, 0.0);
    auto rebuild_dtilde = [&](){
        std::fill(dtilde.begin(), dtilde.end(), 0.0);
        igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
        for (int u = 0; u < n; ++u) {
            igraph_incident(&g, &inc, u, IGRAPH_ALL);
            double sum = 0.0;
            for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
                int eid = VECTOR(inc)[ii];
                sum += tilde_sigma[eid];
            }
            dtilde[u] = sum;
        }
        igraph_vector_int_destroy(&inc);
    };

    auto count_mminus = [&](){
        int neg = 0;
        for (int eid = 0; eid < m; ++eid) if (sigma_s_edge_sign[eid] < 0.0) ++neg;
        return neg;
    };

    auto flip_working = [&](int u){
        s[u] = -s[u];
        igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
        igraph_incident(&g, &inc, u, IGRAPH_ALL);
        for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
            int eid = VECTOR(inc)[ii];
            sigma_s_edge_sign[eid] = -sigma_s_edge_sign[eid];
            tilde_sigma[eid]       = -tilde_sigma[eid];
        }
        igraph_vector_int_destroy(&inc);

#if SG_DEBUG
        int cur_mminus = count_mminus();
        SGLOG("[SG-PROBE] flip A: u=%d, m_minus=%d, min(dtilde)≈%.6f", u, cur_mminus,
              *std::min_element(dtilde.begin(), dtilde.end()));
#endif
    };

    rebuild_tilde_sigma();
    rebuild_dtilde();

#if SG_DEBUG
    if (opts.frac_x) {
        int ml = igraph_ecount(&g);
        int neg_tau=0, zero_tau=0, pos_tau=0;
        for (int eid=0; eid<ml; ++eid) {
            igraph_integer_t u,v; igraph_edge(&g, eid, &u, &v);
            double tau = tau_from_x((*opts.frac_x)[(int)u], (*opts.frac_x)[(int)v]);
            if (tau < 0) ++neg_tau; else if (tau > 0) ++pos_tau; else ++zero_tau;
        }
        int neg_dtilde=0, z_dtilde=0;
        for (double d : dtilde) { if (d < 0) ++neg_dtilde; else if (d == 0.0) ++z_dtilde; }
        SGLOG("[SG-PROBE] tau dist: neg=%d zero=%d pos=%d; vertices: neg_dtilde=%d zero_dtilde=%d",
              neg_tau, zero_tau, pos_tau, neg_dtilde, z_dtilde);
    }
#endif

    int best_mminus = count_mminus();
    std::vector<int> s_best = s;

#if SG_DEBUG
    SGLOG("[SG-PROBE] start greedy: R_max=%d K_max=%d L_max=%d Delta=%d has_frac=%d",
          opts.R_max, opts.K_max, opts.L_max, opts.Delta, (opts.frac_x?1:0));
    SGLOG("[SG-PROBE] start m_minus=%d", best_mminus);
#endif

    for (int r = 0; r < opts.R_max; ++r) {
        int z = 0; // consecutive zero-clique flips

        // A) Fractional greedy flips (min-heap on dtilde)
        {
            struct Node { int u; double key; double sal; };
            auto cmp = [](const Node& a, const Node& b){
                if (a.key != b.key) return a.key > b.key;
                return a.sal < b.sal;
            };
            std::priority_queue<Node,std::vector<Node>,decltype(cmp)> pq(cmp);

            for (int u = 0; u < n; ++u) {
                double sal_u = has_frac ? vertex_salience((*X)[u]) : 0.0;
                pq.push({u, dtilde[u], sal_u});
            }

            while (!pq.empty() && pq.top().key < 0.0) {
                int u = pq.top().u; pq.pop();
                if (dtilde[u] >= 0.0) continue;

                flip_working(u);

                igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
                igraph_incident(&g, &inc, u, IGRAPH_ALL);

                dtilde[u] = -dtilde[u];
                for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
                    int eid = VECTOR(inc)[ii];
                    int v = IGRAPH_OTHER(&g, eid, u);

                    double sumv = 0.0;
                    igraph_vector_int_t incv; igraph_vector_int_init(&incv, 0);
                    igraph_incident(&g, &incv, v, IGRAPH_ALL);
                    for (int jj = 0; jj < (int)igraph_vector_int_size(&incv); ++jj) {
                        int e2 = VECTOR(incv)[jj];
                        sumv += tilde_sigma[e2];
                    }
                    dtilde[v] = sumv;
                    igraph_vector_int_destroy(&incv);

                    double sal_v = has_frac ? vertex_salience((*X)[v]) : 0.0;
                    pq.push({v, dtilde[v], sal_v});
                }
                igraph_vector_int_destroy(&inc);

                int cur_mminus = count_mminus();
                if (cur_mminus < best_mminus) {
                    best_mminus = cur_mminus;
                    s_best = s;
                }
                double sal_u = has_frac ? vertex_salience((*X)[u]) : 0.0;
                pq.push({u, dtilde[u], sal_u});
            }
        }

#if SG_DEBUG
        int z0_count = 0; for (double d: dtilde) if (d==0.0) ++z0_count;
        SGLOG("[SG-PROBE] B: Z0 size=%d (cap z<K_max? %d<%d)", z0_count, z, opts.K_max);
#endif

        // B) Fractional zero-clique (strictly positive edges); limited by K_max
        bool advanced = false;
        if (z < opts.K_max) {
            std::vector<int> z0_vertices; z0_vertices.reserve(n);
            for (int u = 0; u < n; ++u) if (dtilde[u] == 0.0) z0_vertices.push_back(u);

            if (!z0_vertices.empty()) {
                auto Q = maximal_salience_clique_strict_pos(z0_vertices, /*inZ0*/{}, tilde_sigma, X);

#if SG_DEBUG
                SGLOG("[SG-PROBE] B: Q.size=%zu", Q.size());
#endif

                if (Q.size() >= 2) {
                    int ustar = pick_u_star_in_Q(Q, tilde_sigma, Y);

#if SG_DEBUG
                    SGLOG("[SG-PROBE] B: u*=%d", ustar);
#endif

                    if (ustar != -1) {
                        flip_working(ustar);

                        igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
                        igraph_incident(&g, &inc, ustar, IGRAPH_ALL);
                        dtilde[ustar] = -dtilde[ustar];
                        for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
                            int eid = VECTOR(inc)[ii];
                            int v = IGRAPH_OTHER(&g, eid, ustar);
                            double sumv = 0.0;
                            igraph_vector_int_t incv; igraph_vector_int_init(&incv, 0);
                            igraph_incident(&g, &incv, v, IGRAPH_ALL);
                            for (int jj = 0; jj < (int)igraph_vector_int_size(&incv); ++jj) sumv += tilde_sigma[VECTOR(incv)[jj]];
                            dtilde[v] = sumv;
                            igraph_vector_int_destroy(&incv);
                        }
                        igraph_vector_int_destroy(&inc);

                        int cur_mminus = count_mminus();
                        if (cur_mminus < best_mminus) { best_mminus = cur_mminus; s_best = s; }

                        ++z;
                        advanced = true;
                    }
                }
            }
        }
        if (advanced) continue;

        // C) Integer detour
        {
            std::vector<int> s_int;
            std::vector<double> sigma_int_edge_sign;
            if (has_frac) integer_projection_from_x(*X, /*base*/sigma_s_edge_sign, s_int, sigma_int_edge_sign);
            else {
                std::vector<double> X1(n, 1.0);
                integer_projection_from_x(X1, /*base*/sigma_s_edge_sign, s_int, sigma_int_edge_sign);
            }

            std::vector<int> d_int;
            integer_greedy_pass_minheap(sigma_int_edge_sign, s_int, d_int);

            std::vector<int> S_replay; S_replay.reserve(opts.L_max);
            for (int t = 0; t < opts.L_max; ++t) {
#if SG_DEBUG
                SGLOG("[SG-PROBE] C: entering detour L_max=%d", opts.L_max);
                int z0i=0; for (int uu=0; uu<n; ++uu) if (d_int[uu]==0) ++z0i;
                SGLOG("[SG-PROBE] C: Z0_int=%d", z0i);
#endif
                std::vector<int> z0_vertices; z0_vertices.reserve(n);
                for (int u = 0; u < n; ++u) if (d_int[u] == 0) z0_vertices.push_back(u);
                if (z0_vertices.empty()) break;

                std::vector<double> sigma_int_as_tilde = sigma_int_edge_sign;
                auto Q = maximal_salience_clique_strict_pos(z0_vertices, /*inZ0*/{}, sigma_int_as_tilde, X);
                if (Q.size() < 2) break;

                int ustar = pick_u_star_in_Q(Q, sigma_int_as_tilde, Y);

#if SG_DEBUG
                SGLOG("[SG-PROBE] C: Q.size=%zu, u*=%d, gateDelta=%d (deg_u computed next)", Q.size(), ustar, opts.Delta);
#endif

                if (ustar == -1) break;

                // Gate: deg(u*) on current working signature (sigma_s) <= Delta
                int deg_u = 0;
                igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
                igraph_incident(&g, &inc, ustar, IGRAPH_ALL);
                for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
                    int eid = VECTOR(inc)[ii];
                    if (sigma_s_edge_sign[eid] > 0.0) ++deg_u; else --deg_u;
                }
                igraph_vector_int_destroy(&inc);

#if SG_DEBUG
                SGLOG("[SG-PROBE] C: deg(u*) on working sig = %d (<=Delta? %s)", deg_u, (deg_u<=opts.Delta?"yes":"no"));
#endif

                if (deg_u > opts.Delta) break;

                S_replay.push_back(ustar);

                igraph_vector_int_t inc2; igraph_vector_int_init(&inc2, 0);
                igraph_incident(&g, &inc2, ustar, IGRAPH_ALL);
                s_int[ustar] = -s_int[ustar];
                d_int[ustar] = -d_int[ustar];
                for (int ii = 0; ii < (int)igraph_vector_int_size(&inc2); ++ii) {
                    int eid = VECTOR(inc2)[ii];
                    int v = IGRAPH_OTHER(&g, eid, ustar);
                    double &sig = sigma_int_edge_sign[eid];
                    sig = -sig;
                    if (sig > 0.0) d_int[v] += 2; else d_int[v] -= 2;
                }
                igraph_vector_int_destroy(&inc2);
            }

            bool created_neg = false;
            std::vector<int> replayed;
            for (int u : S_replay) {
#if SG_DEBUG
                SGLOG("[SG-PROBE] C: replay |S|=%zu", S_replay.size());
#endif
                flip_working(u);

                igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
                igraph_incident(&g, &inc, u, IGRAPH_ALL);
                dtilde[u] = -dtilde[u];
                for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
                    int eid = VECTOR(inc)[ii];
                    int v = IGRAPH_OTHER(&g, eid, u);
                    double sumv = 0.0;
                    igraph_vector_int_t incv; igraph_vector_int_init(&incv, 0);
                    igraph_incident(&g, &incv, v, IGRAPH_ALL);
                    for (int jj = 0; jj < (int)igraph_vector_int_size(&incv); ++jj) sumv += tilde_sigma[VECTOR(incv)[jj]];
                    dtilde[v] = sumv;
                    igraph_vector_int_destroy(&incv);
                    if (dtilde[v] < 0.0) created_neg = true;
                }
                igraph_vector_int_destroy(&inc);

                replayed.push_back(u);
                if (dtilde[u] < 0.0) created_neg = true;

                if (created_neg) break;
            }

            if (created_neg) {
                int cur_mminus = count_mminus();
                if (cur_mminus < best_mminus) { best_mminus = cur_mminus; s_best = s; }
                continue;
            } else {
                for (int i = (int)replayed.size()-1; i >= 0; --i) {
                    flip_working(replayed[i]); // revert
                }
            }
        }

        break; // no advance in this round
    }

#if SG_DEBUG
    {
        std::vector<int> comp_id, comp_sz;
        int comps = cc_on_strict_pos(&g, switched_weights, comp_id, comp_sz);
        int biggest = 0; for (int szv: comp_sz) biggest = std::max(biggest, szv);
        SGLOG("[SG-PROBE] end: +edge CCs=%d, largest=%d", comps, biggest);
    }
#endif

    return s_best;
}

const std::vector<int> SignedGraph::greedy_switching() {
    auto trivial_cmp = [](int, int) { return false; };
    return greedy_switching_base(trivial_cmp);
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

int SignedGraph::vertex_count() const { return igraph_vcount(&g); }
int SignedGraph::edge_count() const { return igraph_ecount(&g); }

void SignedGraph::print_info() const {
    std::cout << "Graph has " << vertex_count() << " vertices and " << edge_count() << " edges." << std::endl;
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
    igraph_vector_init(&color_attr, igraph_vcount(&g));
    for (int i = 0; i < igraph_vcount(&g); ++i) {
        VECTOR(color_attr)[i] = partition[i];
    }
    igraph_cattribute_VAN_setv(const_cast<igraph_t*>(&g), "color", &color_attr);
    igraph_vector_destroy(&color_attr);

    igraph_strvector_t ecolors;
    igraph_strvector_init(&ecolors, igraph_ecount(&g));

    int i = 0;
    for (const auto& [edge, sign] : signs_view()) {
        bool is_frustrated = (sign == 1 && partition[edge.first] != partition[edge.second]) ||
                             (sign == -1 && partition[edge.first] == partition[edge.second]);
        igraph_strvector_set(&ecolors, i++, is_frustrated ? "red" : "blue");
    }
    igraph_cattribute_EAS_setv(const_cast<igraph_t*>(&g), "color", &ecolors);

    igraph_t gcopy;
    igraph_copy(&gcopy, &g);

    igraph_matrix_t layout;
    if (custom_layout) {
        igraph_matrix_init(&layout, igraph_vcount(&g), 2);
        int above = 0, below = 0;
        for (int i2 = 0; i2 < igraph_vcount(&g); ++i2) {
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
    const int n = igraph_vcount(&g);
    const int m = igraph_ecount(&g);
    std::vector<int> tri(n, 0);
    if (n == 0 || m == 0) return tri;

    std::vector<std::vector<int>> adj(n);
    for (igraph_integer_t eid = 0; eid < m; ++eid) {
        igraph_integer_t u, v; igraph_edge(&g, eid, &u, &v);
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
                const double w_uv = switched_weights[it_uv->second];
                const double w_vw = switched_weights[it_vw->second];
                const double w_uw = switched_weights[it_uw->second];
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
