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

std::vector<int> SignedGraph::greedy_switching_base(const std::function<int(int, int)>& cmp_fn,
                                                     const GreedyKickOptions& opts) {
    int n = igraph_vcount(&g);
    std::vector<double> d(n);
    for (int i = 0; i < n; ++i)
        d[i] = switched_d_plus[i] - switched_d_minus[i];

    std::vector<int> s(n, 0);

    using Heap = boost::heap::pairing_heap<int, boost::heap::compare<std::function<bool(int, int)>> >;
    std::unordered_map<int, Heap::handle_type> handles;
    std::unordered_set<int> in_heap;

    auto heap_cmp = [&](int a, int b) {
        int cmp = cmp_fn(a, b);
        if (d[a] < 0.0 && d[b] < 0.0 && cmp != 0) return cmp > 0;
        if (d[a] == d[b]) return std::signbit(d[b]) && !std::signbit(d[a]);
        return d[a] > d[b];
    };

    Heap pq(heap_cmp);
    for (int i = 0; i < n; ++i) {
        handles[i] = pq.push(i);
        in_heap.insert(i);
    }

    igraph_vector_int_t incident;
    igraph_vector_int_init(&incident, 0);

    auto find_max_clique = [&](std::unordered_set<int>& zero_nodes,
                               const std::unordered_map<int, std::unordered_set<int>>& adj) -> std::vector<int> {
        std::vector<int> best;
        for (auto it = zero_nodes.begin(); it != zero_nodes.end(); ) {
            int u = *it;
            auto it_u = adj.find(u);
            if (it_u == adj.end()) {
                it = zero_nodes.erase(it);
                continue;
            }

            std::vector<int> current_clique = {u};
            for (int v : it_u->second) {
                bool connected = true;
                auto it_v = adj.find(v);
                if (it_v == adj.end()) continue;
                for (int w : current_clique) {
                    if (!it_v->second.count(w)) {
                        connected = false;
                        break;
                    }
                }
                if (connected) current_clique.push_back(v);
            }
            if (current_clique.size() > best.size()) best = std::move(current_clique);
            it++;
        }
        return best;
    };

    // --- KICK summary acc (per-run & global) ---
    static int RUNS = 0;
    static long long SUM_KICKS = 0;
    static long long SUM_CREATED_NEG = 0;
    static double SUM_AW = 0.0;
    static bool KICK_HDR = false;

    int kicks_used = 0;
    int created_neg_this_run = 0;
    double aw_sum_this_run = 0.0;

    while (!pq.empty()) {
        int u = pq.top(); pq.pop();
        if (d[u] > 0.0 || (d[u] == 0.0 && !std::signbit(d[u]))) {
	        in_heap.erase(u);
			std::unordered_set<int> zero_nodes;
			for (int i = 0; i < n; ++i)
			    if (d[i] == 0.0 && !std::signbit(d[i])) zero_nodes.insert(i);
			
			// Fast exit: if Z0 is empty, skip building the adjacency-over-Z0 (O(m)).
			std::unordered_map<int, std::unordered_set<int>> adj;
			if (!zero_nodes.empty()) {
				for (int eid = 0; eid < igraph_ecount(&g); ++eid) {
		            double& w = switched_weights[eid];
		            if (w < 0.0 || (w == 0.0 && std::signbit(w))) continue;
	                igraph_integer_t from, to;
	                igraph_edge(&g, eid, &from, &to);
	                int u = static_cast<int>(from);
	                int v = static_cast<int>(to);
	                if (zero_nodes.count(u) && zero_nodes.count(v)) {
	                    adj[u].insert(v);
	                    adj[v].insert(u);
	                }
	            }
            }

			if (!zero_nodes.empty() && !adj.empty()) {
                std::vector<int> clique = find_max_clique(zero_nodes, adj);
                if (!clique.empty()) {
                    u = clique.front();
                    single_switching(u, &incident);
                    s[u] = 1 - s[u];
                    d[u] = switched_d_plus[u] - switched_d_minus[u];
                    for (int v : clique) {
                        if (v == u) continue;
                        d[v] = switched_d_plus[v] - switched_d_minus[v];
                        if (in_heap.count(v)) pq.update(handles[v], v);
                        else {
                            handles[v] = pq.push(v);
                            in_heap.insert(v);
                        }
                    }
                    continue;
                }
            }

            // ---- No negative-d vertex and no zero-clique switch: try guarded KICK ----
            auto is_neg = [&](double w){ return (w < 0.0) || (w == 0.0 && std::signbit(w)); };
            int m = igraph_ecount(&g);
            int mminus = 0;
            for (int eid = 0; eid < m; ++eid) if (is_neg(switched_weights[eid])) ++mminus;
            bool gate = false;
            if (opts.neg_edge_threshold_abs >= 0 && mminus > opts.neg_edge_threshold_abs) gate = true;
            if (opts.neg_edge_threshold_frac >= 0 && (double)mminus / (double)m > opts.neg_edge_threshold_frac) gate = true;

			if (kicks_used < opts.max_kicks && gate) {
			    // --- Build Z0 (zero net-degree vertices) ---
			    std::vector<char> inZ0(n, 0);
			    int Z0count = 0;
			    for (int i = 0; i < n; ++i)
			        if (d[i] == 0.0 && !std::signbit(d[i])) { inZ0[i] = 1; ++Z0count; }
			
			    auto positive_weight = [&](int eid){
			        double w = switched_weights[eid];
			        return (w > 0.0) || (w == 0.0 && !std::signbit(w));
			    };
			
			    // Fast path: if Z0 is empty and the policy says “don’t relax”, don’t try KICK.
			    if (Z0count == 0 && !opts.relax_to_all_pos_if_Z0_empty) {
			        break; // exit the greedy loop as before
			    }
			    const bool restrict_to_Z0 = (Z0count > 0);  // if Z0 exists, keep the classic restriction
			
			    int best_u = -1; double best_score = 0.0; double best_aw = 0.0; int best_created = 0;
			
			    igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
			    for (int cand = 0; cand < n; ++cand) {
			        igraph_incident(&g, &inc, cand, IGRAPH_ALL);
			        const int deg  = (int)igraph_vector_int_size(&inc);
			        const int scan = std::min(deg, opts.neighbor_cap);
			
			        // --- O(scan) estimate of Δm⁻ for flipping 'cand' (pos→neg minus neg→pos) ---
			        // Sample the first 'scan' incident edges without any Z0 restriction.
			        int sample_pos = 0, sample_neg = 0;
			        for (int ii = 0; ii < scan; ++ii) {
			            int eid = VECTOR(inc)[ii];
			            double w = switched_weights[eid];
			            if ((w > 0.0) || (w == 0.0 && !std::signbit(w))) ++sample_pos;
			            else                                            ++sample_neg;
			        }
			        // Scale the sample to a conservative integer estimate of the full degree effect.
			        // Using ceil() to be cautious; if scan==deg this equals the exact value.
			        const double scale = (scan > 0 ? (double)deg / (double)scan : 1.0);
			        const int delta_m_minus_est = (scan > 0)
			                                    ? (int)std::ceil((sample_pos - sample_neg) * scale)
			                                    : 0;
			
			        // Hard gate: skip vertices whose estimated Δm⁻ exceeds the cap.
			        if (delta_m_minus_est > opts.delta_m_minus_cap) {
			            continue;
			        }
			
			        // --- First pass: score using allowed neighbors (Z0 or all-positive when relaxed) ---
			        double Aw = 0.0; int created = 0;
			        std::vector<int> posN; posN.reserve(scan);
			
			        for (int ii = 0; ii < scan; ++ii) {
			            int eid = VECTOR(inc)[ii];
			            if (!positive_weight(eid)) continue; // only positive neighbors pre-flip
			            int v = IGRAPH_OTHER(&g, eid, cand);
			            if (restrict_to_Z0 && !inZ0[v]) continue; // restrict to Z0 if requested
			
			            double w = switched_weights[eid];
			            double sal = 0.0;
			            if (opts.edge_salience && eid < (int)opts.edge_salience->size()) {
			                sal = (*(opts.edge_salience))[eid];
			                if      (sal < 0.0) sal = 0.0;
			                else if (sal > 1.0) sal = 1.0;
			            }
			
			            // salience-nudged contribution (nonnegative, preserves weighting mode)
			            const double contrib = (opts.use_weighted_degree ? w : 1.0) * (1.0 + opts.kick_salience_bias * sal);
			            Aw += contrib;
			            ++created;
			            posN.push_back(v);
			        }
			        if (Aw <= 0.0) continue;
			
			        double score = Aw;
			
			        // Optional soft penalty (kept for parity; you can set penalty=0 to skip)
			        if (opts.delta_m_minus_penalty > 0.0 && delta_m_minus_est > 0)
			            score -= opts.delta_m_minus_penalty * (double)delta_m_minus_est;
			
			        // Triangle tiebreak on the (already bounded) posN
			        if (opts.use_triangle_tiebreak && !posN.empty()) {
			            int cap = opts.triangle_cap_per_u;
			            int counted = 0, neg_pairs = 0;
			            for (size_t i1 = 0; i1 < posN.size() && counted < cap; ++i1) {
			                for (size_t i2 = i1 + 1; i2 < posN.size() && counted < cap; ++i2) {
			                    int v = posN[i1], wv = posN[i2];
			                    igraph_integer_t eid_vw;
			                    if (igraph_get_eid(&g, &eid_vw, v, wv, 0, 0) == IGRAPH_SUCCESS) {
			                        double w = switched_weights[eid_vw];
			                        if ((w < 0.0) || (w == 0.0 && std::signbit(w))) ++neg_pairs;
			                    }
			                    ++counted;
			                }
			            }
			            score += opts.triangle_beta * (double)neg_pairs;
			        }
			
			        if (score > best_score) {
			            best_score = score; best_u = cand; best_aw = Aw; best_created = created;
			        }
			    }
			    igraph_vector_int_destroy(&inc);
			
			    if (best_u != -1 && best_score > 0.0 /*kick applies*/) {
			        single_switching(best_u, &incident);
			        s[best_u] = 1 - s[best_u];
			        d[best_u] = switched_d_plus[best_u] - switched_d_minus[best_u];
			        for (int i = 0; i < igraph_vector_int_size(&incident); ++i) {
			            int v = IGRAPH_OTHER(&g, VECTOR(incident)[i], best_u);
			            d[v] = switched_d_plus[v] - switched_d_minus[v];
			            if (in_heap.count(v)) pq.update(handles[v], v);
			            else if (d[v] < 0.0 || (d[v] == 0.0 && std::signbit(d[v]))) {
			                handles[v] = pq.push(v);
			                in_heap.insert(v);
			            }
			        }
			        ++kicks_used;
			        created_neg_this_run += best_created;
			        aw_sum_this_run += best_aw;
			        std::cout << "[KICK] M-=" << mminus << "/" << m
			                  << ", |Z0|=" << Z0count
			                  << ", u*=" << best_u
			                  << ", Aw=" << best_aw
			                  << ", created=" << best_created
			                  << std::endl;
			        continue; // resume greedy; fresh negatives likely exist
			    }
			}
            break;
        }

        in_heap.erase(u);
        single_switching(u, &incident);
        s[u] = 1 - s[u];
        d[u] = switched_d_plus[u] - switched_d_minus[u];

        for (int i = 0; i < igraph_vector_int_size(&incident); ++i) {
            int v = IGRAPH_OTHER(&g, VECTOR(incident)[i], u);
            d[v] = switched_d_plus[v] - switched_d_minus[v];

            if (in_heap.count(v))
                pq.update(handles[v], v);
            else if (d[v] < 0.0 || (d[v] == 0.0 && std::signbit(d[v]))) {
                handles[v] = pq.push(v);
                in_heap.insert(v);
            }
        }
    }

    igraph_vector_int_destroy(&incident);

    // --- Per-run KICK summary aggregation & periodic print ---
    ++RUNS;
    SUM_KICKS += kicks_used;
    SUM_CREATED_NEG += created_neg_this_run;
    SUM_AW += aw_sum_this_run;
    if (RUNS % 50 == 0) {
        if (!KICK_HDR) {
            std::cout << std::setw(10) << "Runs"
                      << std::setw(12) << "Kicks"
                      << std::setw(16) << "AvgAw/Run"
                      << std::setw(16) << "AvgCreated"
                      << std::endl;
            KICK_HDR = true;
        }
        double avgAw = (RUNS ? (SUM_AW / RUNS) : 0.0);
        double avgCr = (RUNS ? ((double)SUM_CREATED_NEG / RUNS) : 0.0);
        std::cout << std::setw(10) << RUNS
                  << std::setw(12) << SUM_KICKS
                  << std::setw(16) << avgAw
                  << std::setw(16) << avgCr
                  << std::endl;
    }

    return s;
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
