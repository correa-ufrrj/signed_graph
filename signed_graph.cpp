// signed_graph.cpp
#include "signed_graph.h"
#include <igraph/igraph_strvector.h>
#include <limits>
#include <cstdio>
#include <sstream>

#include <iomanip>

std::ostream& operator<<(std::ostream& os, const Edge& e) {
    os << "(" << e.first << ", " << e.second << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const SignedEdge& e) {
    std::ostringstream tmp;
    tmp.setf(std::ios::fixed); tmp<<std::setprecision(3)<<e.salience;
    os << e.points << ": " << e.sign << " [salience=" << tmp.str() << "]";
    return os;
}

// ===== GraphCore ctors =====
GraphCore::GraphCore(const GraphCore* const other)
    : g_(other->g_), sigma_(other->sigma_), d_pos(other->d_pos), d_neg(other->d_neg), _edge_index(other->_edge_index) {}

GraphCore::GraphCore(const GraphCore* const other, std::vector<double> new_weights)
    : g_(other->g_), _edge_index(other->_edge_index) {
    validate_and_fill_sigma(sigma_, new_weights);
    compute_degrees();
}

GraphCore::GraphCore(const std::shared_ptr<igraph_t> base, std::vector<double> new_sigma)
    : g_(base) {
    validate_and_fill_sigma(sigma_, new_sigma);
    const int m = igraph_ecount(g_.get());
    _edge_index.clear();
    _edge_index.reserve(m);
    for (int eid = 0; eid < m; ++eid) {
        igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
        _edge_index.emplace(Edge{(int)u, (int)v}, eid);
    }
    compute_degrees();
}

GraphCore::GraphCore(const std::string& file_path)
    : g_(std::shared_ptr<igraph_t>(new igraph_t, IGraphDeleter{})) {
    std::ifstream infile(file_path);
    if (!infile.is_open()) throw std::runtime_error("Error opening file: " + file_path);
    igraph_set_attribute_table(&igraph_cattribute_table);

    std::vector<igraph_integer_t> edges_flat;
    std::string line;
    int max_v = 0;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string t;
        int u, v, sign;
        std::getline(ss, t, ','); u = std::stoi(t);
        std::getline(ss, t, ','); v = std::stoi(t);
        std::getline(ss, t, ','); sign = std::stoi(t);
        if (sign != 1 && sign != -1) throw std::runtime_error("Invalid sign in line: " + line);
        edges_flat.push_back(u);
        edges_flat.push_back(v);
        sigma_.push_back(sign);
        max_v = std::max(max_v, std::max(u, v));
    }

    igraph_empty(g_.get(), max_v + 1, IGRAPH_UNDIRECTED);
    igraph_vector_int_t edges_vec;
    igraph_vector_int_view(&edges_vec, edges_flat.data(), (igraph_integer_t)edges_flat.size());
    igraph_add_edges(g_.get(), &edges_vec, nullptr);

    _edge_index.clear();
    const int m = edge_count();
    _edge_index.reserve(m);
    for (int eid = 0; eid < m; ++eid) {
        igraph_integer_t u, v; igraph_edge(g_.get(), eid, &u, &v);
        _edge_index.emplace(Edge{(int)u, (int)v}, eid);
    }

    compute_degrees();
}

// ===== SignedGraph ctors =====
SignedGraph::SignedGraph(const std::shared_ptr<igraph_t> base,
                         std::vector<double> new_sigma,
                         std::vector<double> new_s)
    : GraphCore(base, std::move(new_sigma)),
      s_(std::move(new_s)) {
    const int n = igraph_vcount(base.get());
    if ((int)s_.size() != n) throw std::invalid_argument("Switching vector size does not match vertex count.");
}

SignedGraph::SignedGraph(const SignedGraph* const other)
    : GraphCore(static_cast<const GraphCore* const>(other)),
      s_(other->s_) {}

SignedGraph::SignedGraph(const SignedGraph* const other, std::vector<double> new_sigma)
    : GraphCore(static_cast<const GraphCore* const>(other), std::move(new_sigma)),
      s_(other->s_) {}

SignedGraph::SignedGraph(const std::string& file_path)
    : GraphCore(file_path), s_(std::vector<double>(vertex_count(), 1.0)) {
}

SignedGraph::SignedGraph(const SignedGraph* const other,
                         std::vector<double> new_sigma,
                         std::vector<double> new_s)
    : GraphCore(other, std::move(new_sigma)),
      s_(std::move(new_s)) {
    if ((int)s_.size() != igraph_vcount(g_.get()))
        throw std::invalid_argument("Switching vector size does not match vertex count.");
}

std::unique_ptr<SignedGraph> SignedGraph::clone() const {
    return std::make_unique<SignedGraph>(*this);
}

SignedGraph::~SignedGraph() {}

int SignedGraph::frustrated_edges() const {
    return negative_edge_count();
}

const std::vector<Edge> SignedGraph::frustrated_edges_keys() const {
    return get_negative_edges();
}

bool SignedGraph::has_fractional_switching(double eps) const {
    for (double su : s_) if (std::fabs(su) < 1.0 - eps) return true;
    return false;
}

// signed_graph.cpp
void SignedGraph::single_switching(int u, igraph_vector_int_t* incident) {
    const igraph_t* G = g_.get();

    igraph_vector_int_t local;
    bool owned = false;
    if (incident == nullptr) {
        igraph_vector_int_init(&local, 0);
        incident = &local;
        owned = true;
    }
    if (igraph_vector_int_size(incident) == 0) {
        igraph_incident(G, incident, u, IGRAPH_ALL);
    }

    const int deg = (int)igraph_vector_int_size(incident);
    std::vector<char> was_pos(deg, 0);

    // Snapshot sign & salience BEFORE the flip
    for (int k = 0; k < deg; ++k) {
        const int eid = VECTOR(*incident)[k];
        const int v   = IGRAPH_OTHER(G, eid, u);
        was_pos[k] = !is_neg_edge(eid);
    }

    // Flip and clamp s[u]
    s_[u] = -s_[u];

    // All incident edges invert sign ⇒ swap degree buckets at u
    std::swap(d_pos[u], d_neg[u]);

    // Update neighbor degree buckets and incident salience only for edges (u,v)
    for (int k = 0; k < deg; ++k) {
        const int eid = VECTOR(*incident)[k];
        const int v   = IGRAPH_OTHER(G, eid, u);

        if (was_pos[k]) { // became negative after flip
            d_pos[v] -= 1.0;
            d_neg[v] += 1.0;
        } else {           // became positive after flip
            d_neg[v] -= 1.0;
            d_pos[v] += 1.0;
        }
    }

    if (owned) {
        igraph_vector_int_destroy(&local);
    }
}

void SignedGraph::single_switching(int u) {
    igraph_vector_int_t incident; igraph_vector_int_init(&incident, 0);
    single_switching(u, &incident);
    igraph_vector_int_destroy(&incident);
}

void SignedGraph::flip_and_refresh_net_degree(int u, std::vector<double>& dtilde) {
    const igraph_t* G = g_.get();

    igraph_vector_int_t inc; igraph_vector_int_init(&inc, 0);
    igraph_incident(G, &inc, u, IGRAPH_ALL);

    const int deg = (int)igraph_vector_int_size(&inc);
    std::vector<int> nbrs; nbrs.reserve(deg);
    std::vector<double> p0; p0.reserve(deg);

    // Record pre-flip contributions p0 = s[u]*sigma_e*s[v]
    for (int k = 0; k < deg; ++k) {
        const int eid = VECTOR(inc)[k];
        const int v   = IGRAPH_OTHER(G, eid, u);
        nbrs.push_back(v);
        p0.push_back(s_[u] * static_cast<double>(sigma_[eid]) * s_[v]);
    }

    // Local update on structure (degrees + salience)
    single_switching(u, &inc);

    // Update fractional net degrees locally
    dtilde[u] = -dtilde[u];
    for (int i = 0; i < deg; ++i) {
        const int v = nbrs[i];
        dtilde[v] -= 2.0 * p0[i];
    }

    igraph_vector_int_destroy(&inc);
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
        for (int w : Q) { 
			if (!is_pos_edge(u, w)) { ok = false; break; } }
        if (ok) Q.push_back(u);
    }
    return Q;
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
    for (const auto& [edge, sign, _] : signs_view()) {
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
