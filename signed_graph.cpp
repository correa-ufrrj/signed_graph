// signed_graph.cpp
#include "signed_graph.h"
#include <igraph/igraph_strvector.h>
#include <limits>
#include <cstdio>
#include <sstream>
#include <iomanip>

namespace {

// --- Fast 64-bit hash for sign masks (no external deps) ---------------------
// SplitMix64 mix (public-domain quality mix; good enough for probes)
static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x = x ^ (x >> 31);
    return x;
}

} // namespace

std::ostream& operator<<(std::ostream& os, const Edge& e) {
    os << "(" << e.first << ", " << e.second << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const SignedEdge& e) {
    std::ostringstream tmp;
    os << e.points << ": " << e.sign << "]";
    return os;
}

std::shared_ptr<igraph_t> GraphCore::make_empty_graph_ptr_(int n) {
    igraph_t* raw = new igraph_t;
    if (igraph_empty(raw, n, /*directed=*/0) != IGRAPH_SUCCESS) {
        delete raw;
        throw std::runtime_error("igraph_empty failed");
    }
    return std::shared_ptr<igraph_t>(raw, [](igraph_t* p){
        if (p) { igraph_destroy(p); delete p; }
    });
}

// ===== GraphCore ctors =====
GraphCore::GraphCore(const GraphCore* const other)
    : g_(other->g_),
      sigma_(other->sigma_),
      d_pos(other->d_pos),
      d_neg(other->d_neg),
      m_pos(other->m_pos),m_neg(other->m_neg),
      _edge_index(other->_edge_index),
      N_pos_(other->N_pos_),
      N_neg_(other->N_neg_) {}

GraphCore::GraphCore(const GraphCore* const other, std::vector<double> new_weights)
    : g_(other->g_), _edge_index(other->_edge_index) {
    validate_and_fill_sigma(sigma_, new_weights);
    rebuild_bitmaps_from_sign_();
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
    rebuild_bitmaps_from_sign_();
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

    rebuild_bitmaps_from_sign_();
    compute_degrees();
}

GraphCore::GraphCore(int n, double p, uint64_t seed, bool /*fast_tag*/)
    : g_(make_empty_graph_ptr_(n))
{
    if (n < 3) throw std::invalid_argument("n must be >= 3");
    if (!(p >= 0.0 && p <= 1.0)) throw std::invalid_argument("p must be in [0,1]");

    // (1) Sample activated triples via geometric gaps
    const uint64_t M = choose3_(static_cast<uint64_t>(n));
    std::mt19937_64 rng(seed);
    std::geometric_distribution<uint64_t> geo(p/n);
    std::bernoulli_distribution coin(0.5);

    std::vector<uint64_t> idx;
    idx.reserve(static_cast<size_t>(std::min<uint64_t>(M, 8 + static_cast<uint64_t>(M * p * 1.1))));
    for (uint64_t i = 0; i < M; ) {
        if (p == 0.0) break;
        uint64_t gap = geo(rng);
        if (i + gap >= M) break;
        i += gap;
        idx.push_back(i);
        ++i;
    }

    // (2) Build sign map by completing only missing edges to keep odd negatives in each activated triple
    std::unordered_map<Edge, int, EdgeHash> sign_map;
    sign_map.reserve(static_cast<size_t>(n) * 6u);

    auto get = [&](int a, int b)->int {
        auto it = sign_map.find(Edge{a,b});
        return it == sign_map.end() ? 0 : it->second; // 0=undef
    };
    auto set_if_undef = [&](int a, int b, int s){
        Edge e{a,b};
        if (sign_map.find(e) == sign_map.end())
            sign_map.emplace(std::move(e), s); // never overwrite
    };

    int u, v, w;
    for (uint64_t t : idx) {
        index_to_triple_(static_cast<uint32_t>(n), t, u, v, w);

        const int suv = get(u,v);
        const int svw = get(v,w);
        const int suw = get(u,w);

        const int def = (suv != 0) + (svw != 0) + (suw != 0);
        const int neg = (suv < 0) + (svw < 0) + (suw < 0);

        if (def == 0) {
            const int pick = std::uniform_int_distribution<int>(0,2)(rng);
            set_if_undef(u,v, pick==0 ? -1 : +1);
            set_if_undef(v,w, pick==1 ? -1 : +1);
            set_if_undef(u,w, pick==2 ? -1 : +1);
        } else if (def == 1) {
            if (suv != 0) {
                if (suv > 0) {
                    if (coin(rng)) { set_if_undef(v,w, -1); set_if_undef(u,w, +1); }
                    else           { set_if_undef(v,w, +1); set_if_undef(u,w, -1); }
                } else {
                    const int s = +1; // coin(rng) ? +1 : -1;
                    set_if_undef(v,w, s);
                    set_if_undef(u,w, s);
                }
            } else if (svw != 0) {
                if (svw > 0) {
                    if (coin(rng)) { set_if_undef(u,v, -1); set_if_undef(u,w, +1); }
                    else           { set_if_undef(u,v, +1); set_if_undef(u,w, -1); }
                } else {
                    const int s = +1; // coin(rng) ? +1 : -1;
                    set_if_undef(u,v, s);
                    set_if_undef(u,w, s);
                }
            } else { // suw defined
                if (suw > 0) {
                    if (coin(rng)) { set_if_undef(u,v, -1); set_if_undef(v,w, +1); }
                    else           { set_if_undef(u,v, +1); set_if_undef(v,w, -1); }
                } else {
                    const int s = +1; // coin(rng) ? +1 : -1;
                    set_if_undef(u,v, s);
                    set_if_undef(v,w, s);
                }
            }
        } else { // def == 2
            if (suv == 0)      set_if_undef(u,v, (neg % 2 == 0) ? -1 : +1);
            else if (svw == 0) set_if_undef(v,w, (neg % 2 == 0) ? -1 : +1);
            else               set_if_undef(u,w, (neg % 2 == 0) ? -1 : +1);
        }
        // def==3 cannot occur by construction here.
    }

    // (3) Materialize to igraph + sigma
    std::vector<Edge> edges; edges.reserve(sign_map.size());
    std::vector<double> sigma_d; sigma_d.reserve(sign_map.size());
    for (const auto& kv : sign_map) { edges.push_back(kv.first); sigma_d.push_back(static_cast<double>(kv.second)); }

    igraph_vector_int_t evec;
    igraph_vector_int_init(&evec, static_cast<long>(2 * edges.size()));
    for (size_t i = 0; i < edges.size(); ++i) {
        VECTOR(evec)[2*i + 0] = static_cast<igraph_real_t>(edges[i].first);
        VECTOR(evec)[2*i + 1] = static_cast<igraph_real_t>(edges[i].second);
    }
	if (igraph_add_edges(g_.get(), &evec, /*attr=*/nullptr) != IGRAPH_SUCCESS) {
	    igraph_vector_int_destroy(&evec);
	    throw std::runtime_error("igraph_add_edges failed");
	}
    igraph_vector_int_destroy(&evec);

    // (4) Fill sigma_ (validated to ±1), index map, bitmaps, degrees
    validate_and_fill_sigma(sigma_, sigma_d);

    _edge_index.clear(); _edge_index.reserve(edges.size());
    const igraph_integer_t m = igraph_ecount(g_.get());
    for (igraph_integer_t eid = 0; eid < m; ++eid) {
        igraph_integer_t a,b; igraph_edge(g_.get(), eid, &a, &b);
        _edge_index.emplace(Edge{(int)a,(int)b}, (int)eid);
    }

    rebuild_bitmaps_from_sign_();
    compute_degrees();
}

std::vector<GraphCore::Bitmap> GraphCore::connected_components(int sign) const {
    const int n = vertex_count();

    // visited bitmap over all vertices
    Bitmap visited = make_empty_bitmap((size_t)n);

    std::vector<Bitmap> components;
    components.reserve(4); // heuristic reserve

    // Start a new component at s
    Bitmap q = make_empty_bitmap((size_t)n);

    for (int s = 0; s < n; ++s) {
        if (visited.contains((size_t)s)) continue;

        // Start a new component at s
        Bitmap comp = make_empty_bitmap((size_t)n);

        visited.iset((size_t)s);
        comp.iset((size_t)s);
        q.iset((size_t)s);

        while (!q.empty()) {
            const int u = q.pop_front();
            
			// inside the BFS loop
			const Bitmap& Nb = neighbors_bm((size_t)u, sign);
			Bitmap add = Nb;         // copy
			add.isub(visited);       // keep only *new* vertices
			visited.ior(add);        // mark them visited
			comp.ior(add);           // add to current component
			q.ior(add);              // enqueue them
        }

        components.emplace_back(std::move(comp));
    }

    return components;
}

void GraphCore::on_neg_flip_batch_(const GraphCore::Bitmap& anchor) {
    // C̄ = V \ A
    GraphCore::Bitmap comp = anchor.complement();

    size_t moved_total = 0; // number of edges moved from N⁻ to N⁺ (counted once)

    // For each vertex v in the complement, convert all N⁻(v) ∩ A to positive
    for (size_t vs : comp) {
        const int v = static_cast<int>(vs);

        // av = set of anchor endpoints u with (u,v) currently negative
        GraphCore::Bitmap av = neighbors_bm((size_t)v, -1); // copy
        av.iand(anchor);

        const size_t k = av.count();
        if (k == 0) continue;
        moved_total += k;

        // Move those edges for v in one shot: N⁻(v) -= av, N⁺(v) |= av
        N_neg_[(size_t)v].isub(av);
        N_pos_[(size_t)v].ior(av);

        // Update degree buckets for v
        d_neg[v] -= static_cast<int>(k);
        d_pos[v] += static_cast<int>(k);

        for (size_t us : av) {
            const int u = static_cast<int>(us);

            // Move v within u's adjacency bitmaps (single-bit ops)
            N_neg_[(size_t)u].iclear((size_t)v);
            N_pos_[(size_t)u].iset  ((size_t)v);

            // Update degrees for u
            d_neg[u] -= 1;
            d_pos[u] += 1;
        }
    }

    // Update global edge counters once (each crossing edge counted exactly once)
    m_neg -= static_cast<int>(moved_total);
    m_pos += static_cast<int>(moved_total);
}

// snapshot before changing s[u]
void GraphCore::on_vertex_flip_(int u)
{
	for (int v : neighbors_bm(u)) {
		int old_sign = sign_of(u,v);
	    d_pos[v] -= old_sign; d_neg[v] += old_sign;
	    if (is_pos_(old_sign)) {
	        N_pos_[(size_t)v].iclear((size_t)u);  // clear (v,u) from positive row
	        m_pos--;
	        N_neg_[(size_t)v].iset((size_t)u);  // set   (v,u) in  negative row
	        m_neg++;
	    } else {
	        N_neg_[(size_t)v].iclear((size_t)u);
	        m_neg--;
	        N_pos_[(size_t)v].iset((size_t)u);
	        m_pos++;
	    }
    }

    // Swap degree buckets at u and swap its bitmap rows (+ ↔ −)
    std::swap(d_pos[u], d_neg[u]);
    Bitmap::swap(N_pos_[(size_t)u], N_neg_[(size_t)u]);
}

// ===== SignedGraph ctors =====

// Build from an existing igraph + new σ and a switching vector new_s.
// GraphCore(base, new_sigma) builds bitmaps for raw σ → we must apply switching.
SignedGraph::SignedGraph(const std::shared_ptr<igraph_t> base,
                         std::vector<double> new_sigma,
                         std::vector<double> new_s)
    : GraphCore(base, std::move(new_sigma)),
      s_(this)
{
    compose_switching_inplace(new_s);   // flips bitmaps/degree buckets to match new_s
}

// Pointer-copy from another SignedGraph.
// GraphCore(other) copies already-switched bitmaps/degree buckets → do NOT flip again.
SignedGraph::SignedGraph(const SignedGraph* const other)
    : GraphCore(static_cast<const GraphCore* const>(other)),
      s_(this)
{
    s_.vals_ = other->s_.readonly();    // silent copy; no structural updates
}

// Replace σ, keep other's switching.
// GraphCore(other, new_sigma) is raw σ → we must apply other's s to rebuild effective signs.
SignedGraph::SignedGraph(const SignedGraph* const other, std::vector<double> new_sigma)
    : GraphCore(static_cast<const GraphCore* const>(other), std::move(new_sigma)),
      s_(this)
{
    compose_switching_inplace(other->s_.readonly());  // bring structure in sync with s_copy
}

// Load from file: σ-only; s_ defaults to +1, which already matches the structure.
SignedGraph::SignedGraph(const std::string& file_path)
    : GraphCore(file_path),
      s_(this)
{
    // nothing else to do
}

// Replace σ and s explicitly.
// GraphCore(other, new_sigma) is raw σ → apply new_s to structure.
SignedGraph::SignedGraph(const SignedGraph* const other,
                         std::vector<double> new_sigma,
                         std::vector<double> new_s)
    : GraphCore(static_cast<const GraphCore* const>(other), std::move(new_sigma)),
      s_(this)
{
    compose_switching_inplace(new_s);
}

SignedGraph::SignedGraph(int n, double p, uint64_t seed)
    : GraphCore(n, p, seed, /*fast_tag*/true)
    , s_(this, /*init_val=*/1.0)  // s ≡ +1; no bitmap flips necessary
{
    // WHY: base ctor already built sigma_, bitmaps, degrees.
    // s_ is initialized to +1 so effective signs equal base sigma_.
}

std::unique_ptr<SignedGraph>
SignedGraph::make_random(int n, double p, uint64_t seed) {
    return std::unique_ptr<SignedGraph>(new SignedGraph(n, p, seed));
}

std::unique_ptr<SignedGraph> SignedGraph::clone() const {
    return std::make_unique<SignedGraph>(*this);
}

SignedGraph::~SignedGraph() {}

// SignedGraph::hash_sign_mask — sign-only hash of s_ (SwitchingVector)
uint64_t SignedGraph::hash_sign_mask() const {
    const std::vector<double>& sv = s_.readonly();
    const size_t n = sv.size();

    // Seed with length to make empty/non-empty distinct.
    uint64_t h = 0x9e3779b97f4a7c15ULL ^ static_cast<uint64_t>(n);

    uint64_t word = 0;
    int bit = 0;

    // Pack 64 sign-bits per word: bit=1 iff round_pm1_(s[u]) < 0.
    for (size_t i = 0; i < n; ++i) {
        const bool is_neg = is_neg_(round_pm1_(sv[i]));
        word |= (static_cast<uint64_t>(is_neg) << bit);
        if (++bit == 64) {
            h ^= splitmix64(word ^ static_cast<uint64_t>(i));
            h = (h << 27) | (h >> (64 - 27));
            h = h * 0x165667919E3779F9ULL + 0x9E3779B97F4A7C15ULL;
            word = 0;
            bit  = 0;
        }
    }

    // Tail (if n % 64 != 0)
    if (bit != 0) {
        h ^= splitmix64(word ^ static_cast<uint64_t>(n));
        h = (h << 27) | (h >> (64 - 27));
        h = h * 0x165667919E3779F9ULL + 0x9E3779B97F4A7C15ULL;
    }

    // Murmur3-style finalizer
    h ^= (h >> 33);
    h *= 0xff51afd7ed558ccdULL;
    h ^= (h >> 33);
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= (h >> 33);
    return h;
}


// snapshot before changing s[u]
void SignedGraph::on_vertex_flip_(int u, double from, double to)
{
    // If signs did not change, nothing to do.
    if (is_neg_(from) == is_neg_(to)) return;
	GraphCore::on_vertex_flip_(u);
}

int SignedGraph::frustrated_edges() const {
    return edge_count(-1);
}

const std::vector<Edge> SignedGraph::frustrated_edges_keys() const {
    return get_edges(-1);
}

bool SignedGraph::has_fractional_switching(double eps) const {
    for (double su : s_.readonly()) if (std::fabs(su) < 1.0 - eps) return true;
    return false;
}

void SignedGraph::single_switching(int u) {
//	single_switching(u, nullptr);
	s_[u] = -s_[u];
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
    return maximal_clique_strict_pos(cand);
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
