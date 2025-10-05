// File: signed_graph.cpp
#include "signed_graph.h"
#include <igraph/igraph_strvector.h>
#include <iostream>
#include <limits>
#include <boost/heap/pairing_heap.hpp>
#include <ostream>

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
    _edge_index = other->_edge_index;
};

// SignedGraph weight-based clone constructor
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

// --- Signed zeros on edges (already used style: rely on +0/-0 via signbit) ---
inline bool SignedGraph::strictly_positive(double w) { return w > 0.0; }

// Vertex salience from x (peaks at 0.5)
inline double SignedGraph::vertex_salience(double x) {
    double a = std::fabs(x - 0.5) * 2.0;
    if (a > 1.0) a = 1.0;
    return 1.0 - a; // in [0,1]
}

// Edge salience from y (if provided); otherwise 0.
inline double SignedGraph::edge_salience_from_y(double y) {
    double a = std::fabs(y - 0.5) * 2.0;
    if (a > 1.0) a = 1.0;
    return 1.0 - a; // in [0,1]
}

// Compute product-surrogate tau = 4 x_u x_v - 2 x_u - 2 x_v + 1
inline double SignedGraph::tau_from_x(double xu, double xv) {
    return 4.0 * xu * xv - 2.0 * xu - 2.0 * xv + 1.0;
}

// Build a maximal clique Q inside Z0 on strictly-positive edges,
// greedily by vertex salience. Returns list of vertex ids (may be empty).
std::vector<int> SignedGraph::maximal_salience_clique_strict_pos(
    const std::vector<int>& z0_vertices,
    const std::vector<char>& inZ0,
    const std::vector<double>& switched_frac_edge_sign, // \widetilde\sigma on edges (signed zeros allowed)
    const std::vector<double>* frac_x // for vertex salience
) const {
    const int n = (int)inZ0.size();
    std::vector<int> cand = z0_vertices;
    // Order candidates by descending vertex salience
    std::sort(cand.begin(), cand.end(), [&](int a, int b){
        double sa = frac_x ? vertex_salience((*frac_x)[a]) : 0.0;
        double sb = frac_x ? vertex_salience((*frac_x)[b]) : 0.0;
        if (sa != sb) return sa > sb;
        return a < b;
    });

    std::vector<int> Q;
    Q.reserve(cand.size());

    // (Simple maximal clique by greedy grow)
    for (int u : cand) {
        bool ok = true;
        for (int v : Q) {
            igraph_integer_t eid;
            if (igraph_get_eid(const_cast<igraph_t*>(&g_), &eid, u, v, /*directed=*/0, /*error=*/0) != IGRAPH_SUCCESS) {
                ok = false; break;
            }
            double w = switched_frac_edge_sign[(int)eid];
            if (!strictly_positive(w)) { ok = false; break; } // STRICT > 0, no +0 allowed
        }
        if (ok) Q.push_back(u);
    }

    // Ensure |Q| >= 2, otherwise empty
    if ((int)Q.size() < 2) Q.clear();
    return Q;
}

int SignedGraph::pick_u_star_in_Q(
    const std::vector<int>& Q,
    const std::vector<double>& switched_frac_edge_sign,
    const std::vector<double>* frac_y // edge salience source
) const {
    if (Q.empty()) return -1;

    igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
    double best_score = -1.0;
    int best_u = -1;

    for (int u : Q) {
        igraph_incident(const_cast<igraph_t*>(&g_), &inc, u, IGRAPH_ALL);
        double score = 0.0;
        for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
            int eid = VECTOR(inc)[ii];
            double w = switched_frac_edge_sign[eid];
            if (!strictly_positive(w)) continue; // STRICT
            if (frac_y) {
                score += edge_salience_from_y((*frac_y)[eid]);
            }
        }
        if (score > best_score) { best_score = score; best_u = u; }
    }
    igraph_vector_int_destroy(&inc);
    return best_u;
}

void SignedGraph::integer_projection_from_x(
    const std::vector<double>& x,
    std::vector<int>& s_int,           // out: s_int[u] = 1 - 2*round(x_u>=0.5)
    std::vector<double>& sigma_int_edge_sign // out: +1 / -1 on edges after switching by s_int
) const {
    const int n = (int)x.size();
    s_int.assign(n, 1);
    for (int u = 0; u < n; ++u) {
        int xu = (x[u] >= 0.5) ? 1 : 0;
        s_int[u] = 1 - 2 * xu; // in {-1,+1}
    }
    // sigma_int(uv) = s_int[u] * sigma_s(uv) * s_int[v]
    sigma_int_edge_sign.resize(m_); // m_ is edge count
    for (int eid = 0; eid < m_; ++eid) {
        int u = edges_[eid].first, v = edges_[eid].second;
        double sig = (double)signature_[eid]; // \sigma_s, assume {-1,+1}
        sigma_int_edge_sign[eid] = (double)s_int[u] * sig * (double)s_int[v];
    }
}

void SignedGraph::integer_greedy_pass_minheap(
    std::vector<double>& sigma_int_edge_sign, // +1/-1
    std::vector<int>& s_int,                  // flips update this
    std::vector<int>& d_int                   // out: integer net degree per vertex
) const {
    const int n = (int)s_int.size();
    d_int.assign(n, 0);

    // Compute integer net degree d_u = sum_{v} sigma_int(uv)
    igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
    for (int u = 0; u < n; ++u) {
        igraph_incident(const_cast<igraph_t*>(&g_), &inc, u, IGRAPH_ALL);
        int sum = 0;
        for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
            int eid = VECTOR(inc)[ii];
            sum += (sigma_int_edge_sign[eid] > 0.0) ? 1 : -1;
        }
        d_int[u] = sum;
    }

    // min-heap by d_int[u]
    struct Node { int u; int key; };
    auto cmp = [](const Node& a, const Node& b){ return a.key > b.key; };
    std::priority_queue<Node,std::vector<Node>,decltype(cmp)> pq(cmp);
    for (int u = 0; u < n; ++u) pq.push({u, d_int[u]});

    // Greedy flips while min key < 0
    while (!pq.empty() && pq.top().key < 0) {
        int u = pq.top().u; pq.pop();

        // Revalidate (lazy heap)
        if (d_int[u] >= 0) continue;

        // Flip u in integer view
        s_int[u] = -s_int[u];

        // update incident edges & neighbor degrees
        igraph_incident(const_cast<igraph_t*>(&g_), &inc, u, IGRAPH_ALL);
        d_int[u] = -d_int[u]; // flipping u flips the sign of all incident terms

        for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
            int eid = VECTOR(inc)[ii];
            int v = IGRAPH_OTHER(const_cast<igraph_t*>(&g_), eid, u);

            // sigma_int(uv) flips sign
            sigma_int_edge_sign[eid] = -sigma_int_edge_sign[eid];

            // Update neighbor degree: one incident term flipped; adjust +/-2 accordingly
            if (sigma_int_edge_sign[eid] > 0.0) {
                // became +1, neighbor gains +2 (from -1 to +1)
                d_int[v] += 2;
            } else {
                // became -1, neighbor loses 2 (from +1 to -1)
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

const std::vector<double>& SignedGraph::plus_degrees_view() const {
    return d_plus;
}

const std::vector<double>& SignedGraph::minus_degrees_view() const {
    return d_minus;
}

const long double SignedGraph::max_degree_vertex() const {
    int n = vertex_count();

    // Initialize and compute degrees
    igraph_vector_int_t degrees;
    igraph_real_t max_degree = 0;
    long int max_vertex = 0;

    igraph_vector_int_init(&degrees, n);
    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

    // Find vertex with maximum degree
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
        w = -w; // switch sign of weight, including +0.0 to -0.0 and vice versa
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

// signed_graph.cpp
std::vector<int> SignedGraph::greedy_switching_base(
    const std::function<int(int,int)>& cmp_fn,
    const GreedyKickOptions& opts)
{
    const int n = n_;
    const int m = m_;

    // --- Working switching s (start at +1) ---
    std::vector<int> s(n, 1);
    // If you maintain an external current switching, seed from it.

    // --- Edge signs under current switching: sigma_s(uv) in {-1,+1} ---
    std::vector<double> sigma_s_edge_sign(m, 1.0);
    for (int eid = 0; eid < m; ++eid) {
        int u = edges_[eid].first, v = edges_[eid].second;
        sigma_s_edge_sign[eid] = (double)s[u] * (double)signature_[eid] * (double)s[v];
    }

    // --- Fractional surrogate per-edge: \tilde\sigma = sigma_s * tau(x) ---
    const bool has_frac = (opts.frac_x != nullptr);
    const auto* X = opts.frac_x;
    const auto* Y = opts.frac_y; // optional, only for salience
    std::vector<double> tilde_sigma(m, 0.0);

    auto rebuild_tilde_sigma = [&](){
        if (!has_frac) {
            // integer-only: tilde_sigma coincides with sigma_s
            for (int eid = 0; eid < m; ++eid) tilde_sigma[eid] = sigma_s_edge_sign[eid];
            return;
        }
        for (int eid = 0; eid < m; ++eid) {
            int u = edges_[eid].first, v = edges_[eid].second;
            double tau = tau_from_x((*X)[u], (*X)[v]); // product surrogate
            // signed zeros honored naturally via IEEE if tau==0
            tilde_sigma[eid] = sigma_s_edge_sign[eid] * tau;
        }
    };

    // --- Fractional net degrees \tilde d_u (sum of tilde_sigma over incident edges) ---
    std::vector<double> dtilde(n, 0.0);
    auto rebuild_dtilde = [&](){
        std::fill(dtilde.begin(), dtilde.end(), 0.0);
        igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
        for (int u = 0; u < n; ++u) {
            igraph_incident(&g_, &inc, u, IGRAPH_ALL);
            double sum = 0.0;
            for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
                int eid = VECTOR(inc)[ii];
                sum += tilde_sigma[eid]; // numerically computed; signed zero tag only on edges
            }
            dtilde[u] = sum;
        }
        igraph_vector_int_destroy(&inc);
    };

    // --- Count negatives on current signature (objective tracker) ---
    auto count_mminus = [&](){
        int neg = 0;
        for (int eid = 0; eid < m; ++eid) if (sigma_s_edge_sign[eid] < 0.0) ++neg;
        return neg;
    };

    // --- Apply one flip in working switching s and update sigma_s_edge_sign + tilde ---
    auto flip_working = [&](int u){
        s[u] = -s[u];
        igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
        igraph_incident(&g_, &inc, u, IGRAPH_ALL);
        for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
            int eid = VECTOR(inc)[ii];
            sigma_s_edge_sign[eid] = -sigma_s_edge_sign[eid];
            // tilde flips sign because sigma_s flips
            tilde_sigma[eid] = -tilde_sigma[eid];
        }
        igraph_vector_int_destroy(&inc);
    };

    // --- Initialize state ---
    rebuild_tilde_sigma();
    rebuild_dtilde();
    int best_mminus = count_mminus();
    std::vector<int> s_best = s;

    // --- Round loop ---
    for (int r = 0; r < opts.R_max; ++r) {
        int z = 0; // consecutive zero-clique flips in this round

        // A) Fractional greedy flips: min-heap by dtilde[u]
        {
            struct Node { int u; double key; double sal; };
            auto cmp = [](const Node& a, const Node& b){
                if (a.key != b.key) return a.key > b.key;      // min-heap
                return a.sal < b.sal;                          // tie: larger sal first
            };
            std::priority_queue<Node,std::vector<Node>,decltype(cmp)> pq(cmp);

            for (int u = 0; u < n; ++u) {
                double sal_u = has_frac ? vertex_salience((*X)[u]) : 0.0;
                pq.push({u, dtilde[u], sal_u});
            }

            while (!pq.empty() && pq.top().key < 0.0) {
                int u = pq.top().u; double k = pq.top().key; pq.pop();
                if (dtilde[u] >= 0.0) continue; // lazy recheck

                // Flip u in working state
                flip_working(u);

                // update dtilde locally
                igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
                igraph_incident(&g_, &inc, u, IGRAPH_ALL);

                // flipping u flips the sign of all incident edge contributions
                dtilde[u] = -dtilde[u];
                for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
                    int eid = VECTOR(inc)[ii];
                    int v = IGRAPH_OTHER(&g_, eid, u);
                    // For neighbor v, its incident sum changes by +/-2*tilde contribution of that eid before flip.
                    // But since we don't store the old per-edge, a safe route: recompute v's dtilde incrementally by:
                    // dtilde[v] += (-tilde_old_eid - tilde_old_eid) = -2 * tilde_old_eid.
                    // We can't get tilde_old_eid now (already negated). So rebuild neighbor's dtilde by scanning incident edges.
                    // To stay simple & robust, rebuild only neighbors:
                    double sumv = 0.0;
                    igraph_vector_int_t incv; igraph_vector_int_init(&incv, 0);
                    igraph_incident(&g_, &incv, v, IGRAPH_ALL);
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

                // Track best signature
                int cur_mminus = count_mminus();
                if (cur_mminus < best_mminus) {
                    best_mminus = cur_mminus;
                    s_best = s;
                }
                // Push u back
                double sal_u = has_frac ? vertex_salience((*X)[u]) : 0.0;
                pq.push({u, dtilde[u], sal_u});
            }
        }

        // B) Fractional zero-clique step (strictly-positive edges, |Q|>=2), capped per-round
        bool advanced = false;
        if (z < opts.K_max) {
            // Collect Z0 by numeric zero
            std::vector<char> inZ0(n, 0);
            std::vector<int> z0_vertices; z0_vertices.reserve(n);
            for (int u = 0; u < n; ++u) if (dtilde[u] == 0.0) { inZ0[u] = 1; z0_vertices.push_back(u); }

            if (!z0_vertices.empty()) {
                // Find a maximal clique Q with high vertex salience
                auto Q = maximal_salience_clique_strict_pos(z0_vertices, inZ0, tilde_sigma, X);
                if (Q.size() >= 2) {
                    int ustar = pick_u_star_in_Q(Q, tilde_sigma, Y);
                    if (ustar != -1) {
                        flip_working(ustar);

                        // recompute dtilde only at ustar neighbors (simple & safe: rebuild neighbors)
                        igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
                        igraph_incident(&g_, &inc, ustar, IGRAPH_ALL);
                        dtilde[ustar] = -dtilde[ustar];
                        for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
                            int eid = VECTOR(inc)[ii];
                            int v = IGRAPH_OTHER(&g_, eid, ustar);
                            double sumv = 0.0;
                            igraph_vector_int_t incv; igraph_vector_int_init(&incv, 0);
                            igraph_incident(&g_, &incv, v, IGRAPH_ALL);
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
        if (advanced) {
            // go back to next round iteration (which restarts A)
            continue;
        }

        // C) Integer detour (only if B cannot act)
        {
            // 1) projection
            std::vector<int> s_int;
            std::vector<double> sigma_int_edge_sign;
            if (has_frac) integer_projection_from_x(*X, s_int, sigma_int_edge_sign);
            else {
                // if no fractional input, project from all-ones X (degenerate case)
                std::vector<double> X1(n, 1.0);
                integer_projection_from_x(X1, s_int, sigma_int_edge_sign);
            }

            // 2) integer greedy pass
            std::vector<int> d_int;
            integer_greedy_pass_minheap(sigma_int_edge_sign, s_int, d_int);

            // 3) integer zero-clique (bounded; mirrors fractional case)
            std::vector<int> S_replay; S_replay.reserve(opts.L_max);
            for (int t = 0; t < opts.L_max; ++t) {
                // Build Z0^int: d_int[u] == 0
                std::vector<char> inZ0(n, 0);
                std::vector<int> z0_vertices; z0_vertices.reserve(n);
                for (int u = 0; u < n; ++u) if (d_int[u] == 0) { inZ0[u] = 1; z0_vertices.push_back(u); }

                if (z0_vertices.empty()) break;

                // Maximal clique Q on strictly positive edges in sigma_int
                // Reuse tilde_sigma slot temporarily to avoid extra array: build from sigma_int
                std::vector<double> sigma_int_as_tilde = sigma_int_edge_sign; // +1/-1
                auto Q = maximal_salience_clique_strict_pos(z0_vertices, inZ0, sigma_int_as_tilde, X);
                if (Q.size() < 2) break;

                int ustar = pick_u_star_in_Q(Q, sigma_int_as_tilde, Y);
                if (ustar == -1) break;

                // Gate: deg(u*) on current working signature (sigma_s) <= Delta
                // deg(u) = (#pos incident - #neg incident) under sigma_s
                int deg_u = 0;
                igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
                igraph_incident(&g_, &inc, ustar, IGRAPH_ALL);
                for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
                    int eid = VECTOR(inc)[ii];
                    if (sigma_s_edge_sign[eid] > 0.0) ++deg_u; else --deg_u;
                }
                igraph_vector_int_destroy(&inc);
                if (deg_u > opts.Delta) {
                    // try next best in Q
                    // (simple: remove ustar and retry once)
                    // You can implement a full ranking if needed.
                    // For now, stop if gate fails at top choice.
                    break;
                }

                // Accept in integer view
                S_replay.push_back(ustar);

                // Flip ustar in integer view
                // This flips sigma_int_edge_sign and updates d_int as in greedy pass:
                // reuse the logic from integer_greedy_pass_minheap for a single vertex.
                {
                    igraph_vector_int_t inc2; igraph_vector_int_init(&inc2, 0);
                    igraph_incident(&g_, &inc2, ustar, IGRAPH_ALL);
                    // flip s_int
                    s_int[ustar] = -s_int[ustar];
                    // degree update at ustar
                    d_int[ustar] = -d_int[ustar];
                    for (int ii = 0; ii < (int)igraph_vector_int_size(&inc2); ++ii) {
                        int eid = VECTOR(inc2)[ii];
                        int v = IGRAPH_OTHER(&g_, eid, ustar);
                        double &sig = sigma_int_edge_sign[eid];
                        double old = sig;
                        sig = -sig; // flip
                        if (sig > 0.0) d_int[v] += 2; else d_int[v] -= 2;
                    }
                    igraph_vector_int_destroy(&inc2);
                }
            }

            // 4) Replay S_replay on fractional state; stop early if we create a fractional negative
            bool created_neg = false;
            std::vector<int> replayed;
            for (int u : S_replay) {
                flip_working(u);

                // Refresh dtilde for u and neighbors (local rebuild)
                igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
                igraph_incident(&g_, &inc, u, IGRAPH_ALL);
                dtilde[u] = -dtilde[u];
                for (int ii = 0; ii < (int)igraph_vector_int_size(&inc); ++ii) {
                    int eid = VECTOR(inc)[ii];
                    int v = IGRAPH_OTHER(&g_, eid, u);
                    double sumv = 0.0;
                    igraph_vector_int_t incv; igraph_vector_int_init(&incv, 0);
                    igraph_incident(&g_, &incv, v, IGRAPH_ALL);
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
                // Continue next round (go back to A in the for-loop)
                continue;
            } else {
                // rollback flips from replay (no effect)
                for (int i = (int)replayed.size()-1; i >= 0; --i) {
                    flip_working(replayed[i]); // revert
                }
                // nothing advanced; round ends
            }
        } // end detour

        // If we reach here with no advance, break early
        break;
    } // rounds

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
            return false; // Duplicate edge found
        }

        for (const auto& e : cycle.pos_edges()) {
            if (!seen_edges.insert(e).second) {
                return false; // Duplicate edge found
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

int SignedGraph::vertex_count() const {
    return igraph_vcount(&g);
}

int SignedGraph::edge_count() const {
    return igraph_ecount(&g);
}

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
        for (int i = 0; i < igraph_vcount(&g); ++i) {
            double x = (partition[i] == 0 ? above++ : below++);
            double y = (partition[i] == 0 ? 1.0 : 0.0);
            MATRIX(layout, i, 0) = x;
            MATRIX(layout, i, 1) = y;
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
        igraph_write_graph_dot(&gcopy, fout);  // replaced SVG layout with DOT export
        igraph_matrix_destroy(&layout);
    } else {
        igraph_write_graph_dot(&gcopy, fout);  // replaced SVG output with DOT export
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

    // Build adjacency lists
    std::vector<std::vector<int>> adj(n);
    for (igraph_integer_t eid = 0; eid < m; ++eid) {
        igraph_integer_t u, v; igraph_edge(&g, eid, &u, &v);
        int a = static_cast<int>(u), b = static_cast<int>(v);
        adj[a].push_back(b);
        adj[b].push_back(a);
    }
    for (int u = 0; u < n; ++u) std::sort(adj[u].begin(), adj[u].end());

    // Mark array for fast intersection
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
