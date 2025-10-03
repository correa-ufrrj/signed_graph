// File: signed_graph_mip.cpp
#include "signed_graph_mip.h"
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <set>
#include <boost/heap/pairing_heap.hpp>
#include <numeric>

// Correct signature for igraph ≥ 0.10.x
void silent_warning_handler(const char* reason, const char* file, int line) {
    // Suppress all warnings
}

// SignedGraphForMIP default clone constructor
SignedGraphForMIP::SignedGraphForMIP(const SignedGraph* const other)
  : SignedGraph(other), frac_weights(other->edge_count(), 0.0),
    mask_weights(other->edge_count(), 0.0) {
    const auto& edge_map = this->edge_index();
    for (const auto& [edge, eid] : edge_map) {
        edge_to_eid[edge] = eid;
    }
}

// SignedGraphForMIP weight-based clone constructor
SignedGraphForMIP::SignedGraphForMIP(const SignedGraph* const other, std::vector<double> new_weights)
  : SignedGraph(other, std::move(new_weights)), frac_weights(other->edge_count(), 0.0),
    mask_weights(other->edge_count(), 0.0) {
    const auto& edge_map = this->edge_index();
    for (const auto& [edge, eid] : edge_map) {
        edge_to_eid[edge] = eid;
    }
}

SignedGraphForMIP::SignedGraphForMIP(const std::string& file_path)
    : SignedGraph(file_path)
{
    for (int eid = 0; eid < edge_count(); ++eid) {
        igraph_integer_t u, v;
        igraph_edge(&g, eid, &u, &v);
        edge_to_eid[{static_cast<int>(u), static_cast<int>(v)}] = eid;
    }

    frac_weights.resize(weights.size());
	mask_weights.resize(weights.size());
	std::transform(weights.begin(), weights.end(), frac_weights.begin(),
	               [](double v){ return static_cast<double>(v); });
	std::fill(mask_weights.begin(), mask_weights.end(), 0.0);
}

SignedGraphForMIP::~SignedGraphForMIP() {
}

bool SignedGraphForMIP::weighting_from_fractional(const std::vector<double>& x,
                                                  const std::vector<double>& y) {
    igraph_integer_t ecount = edge_count();
    bool has_changed = false;
    auto signs = signs_view();

    if (mask_weights.size() != static_cast<size_t>(ecount))
        mask_weights.assign(ecount, 0.0);

    for (igraph_integer_t eid = 0; eid < ecount; ++eid) {
        igraph_integer_t from, to; igraph_edge(&g, eid, &from, &to);
        const int old_w = signs[eid].sign; // ±1 (or ±0.0 treated consistently)
        // frac_weights in [0,2] since tau∈[-1,1]; map to [0,1] via 0.5×
        const double tau = 4.0 * y[eid] - 2.0 * x[from] - 2.0 * x[to] + 1.0; // ∈[-1,1]
        frac_weights[eid] = 1.0 - tau * old_w;            // ∈[0,2]
        const double t01 = 0.5 * frac_weights[eid];       // ∈[0,1]
        // salience: 1 at 0.5, fades to 0 at 0 or 1
        mask_weights[eid] = 1.0 - std::min(1.0, 2.0 * std::abs(t01 - 0.5));

        has_changed |= std::abs(y[eid] - x[from] * x[to]) > 1e-5;
    }
    return has_changed;
}

const std::vector<int> SignedGraphForMIP::greedy_switching() {
    auto trivial_cmp = [](int, int) { return 0; };
    auto s = greedy_switching_base(trivial_cmp);

    // Copy "switched_weights" → "frac_weights"
    std::transform(switched_weights.begin(), switched_weights.end(), frac_weights.begin(),
                   [](double val) { return static_cast<double>(val); });
	std::cout << "[SignedGraph] greedy_switching finished." << std::endl;

    return s;
}

// ---- fractional_greedy_switching overloads (SignedGraphForMIP) ----
// Note: if these definitions already exist elsewhere, keep only one definition.
// They route the fractional weighting priority into the greedy with kick options.

std::optional<std::shared_ptr<const std::vector<int>>>
SignedGraphForMIP::fractional_greedy_switching(const SignedGraph::GreedyKickOptions& opts_in) {
    int n = igraph_vcount(&g);
    std::vector<double> w_plus(n, 0.0), w_minus(n, 0.0);

    for (igraph_integer_t eid = 0; eid < igraph_ecount(&g); ++eid) {
        igraph_integer_t u, v;
        igraph_edge(&g, eid, &u, &v);
        double weight = 1 - frac_weights[eid];
        if (weight >= 0) {
            w_plus[u] += weight;
            w_plus[v] += weight;
        } else {
            w_minus[u] += -weight;
            w_minus[v] += -weight;
        }
    }

    std::vector<double> w(n);
    std::transform(w_plus.begin(), w_plus.end(), w_minus.begin(), w.begin(), std::minus<double>());

    auto cmp_fn = [&](int a, int b) -> int {
        if (w[a] > w[b]) return 1;
        if (w[a] < w[b]) return -1;
        return 0;
    };

    auto opts = opts_in; // copy, then augment
	opts.edge_salience = &edge_salience_view();   // uses salience_full_ if set, else mask_weights
    opts.kick_salience_bias = 0.5;       // or expose as parameter / heuristic

    auto s = greedy_switching_base(cmp_fn, opts);
    if (std::all_of(s.begin(), s.end(), [](int v){ return v == 0; }))
        return std::nullopt;
    return std::make_shared<const std::vector<int>>(std::move(s));
}

// Backwards-compatible overload
std::optional<std::shared_ptr<const std::vector<int>>>
SignedGraphForMIP::fractional_greedy_switching() {
    return fractional_greedy_switching(SignedGraph::GreedyKickOptions{});
}

// --- SignedGraphForMIP wrappers over the stream ---
std::vector<NegativeCycle>
SignedGraphForMIP::find_switched_lower_bound(bool cover) {
    std::vector<NegativeCycle> flat;
    auto stream = open_negative_cycle_stream(cover);
    std::vector<NegativeCycle> batch;
    while (stream.next(batch)) {
        flat.insert(flat.end(),
                    std::make_move_iterator(batch.begin()),
                    std::make_move_iterator(batch.end()));
    }
    return flat;
}

std::vector<std::vector<NegativeCycle>>
SignedGraphForMIP::find_switched_lower_bound_grouped(bool cover) const {
    std::vector<std::vector<NegativeCycle>> groups;
    auto stream = open_negative_cycle_stream(cover);
    std::vector<NegativeCycle> batch;
    while (stream.next(batch)) {
        groups.emplace_back();
        groups.back().swap(batch);
    }
    return groups;
}

NegativeCycleBatchStream SignedGraphForMIP::open_negative_cycle_stream(bool cover, bool use_triangle_order) const {
    return NegativeCycleBatchStream(*this, cover, use_triangle_order);
}

// --- NegativeCycleBatchStream definitions ---
NegativeCycleBatchStream::NegativeCycleBatchStream(const SignedGraphForMIP& G, bool cover, bool use_triangle_order)
    : G_(G),
      vcount_(0),
      ecount_(0),
      med_base_pos_(1.0),
      saved_weights_init_(false),
      cover_(cover),
      finished_(false),
      total_found_(0),
      batches_emitted_(0),
      g_pos_built_(false),
      saved_weights_pos_init_(false)
       {
    use_tri_order_ = use_triangle_order;
    if (use_tri_order_) {
        try { neg_tri_vert_ = G_.negative_triangle_count_per_vertex(); }
        catch(...) { neg_tri_vert_.assign(igraph_vcount(&G_.g), 0); }
    } else {
        neg_tri_vert_.clear();
    }
    build_initial_state_();
}

int  NegativeCycleBatchStream::total_cycles_emitted() const { return static_cast<int>(total_found_); }
int  NegativeCycleBatchStream::batches_emitted()      const { return batches_emitted_; }

void NegativeCycleBatchStream::build_mask_for_batch_() {
    // Reset within-batch masks to base weights and clear reuse accumulators
    if (!saved_weights_pos_init_) return; // build_initial_state_ not called yet
    for (long pe = 0; pe < (long)pos2full_eid_.size(); ++pe) {
        igraph_integer_t fe = pos2full_eid_[pe];
        if (use_lp_weights_ && fe < (igraph_integer_t)lp_score_full_.size()) {
            const double base = base_pos_[(size_t)fe];
            const double hint = lp_score_full_[(size_t)fe];
            VECTOR(saved_weights_pos_)[pe] = std::max(1e-12, lp_alpha_ * base / (1.0 + lp_beta_ * hint));
        } else {
            VECTOR(saved_weights_pos_)[pe] = base_pos_[(size_t)fe];
        }
    }
    if (saved_weights_init_) {
        for (long i = 0; i < (long)ecount_; ++i) VECTOR(saved_weights_)[i] = 0.0;
    }
    std::fill(reuse_accum_.begin(), reuse_accum_.end(), 0.0);
}

void NegativeCycleBatchStream::set_lp_scores_full_edges(const std::vector<double>& s, double alpha, double beta) {
    // WHY: steer SPs toward violated cycles while keeping Dijkstra non-negative
    lp_score_full_ = s; use_lp_weights_ = true; lp_alpha_ = alpha; lp_beta_ = beta;
    if (saved_weights_pos_init_) {
        // apply immediately for current batch too
        for (long pe = 0; pe < (long)pos2full_eid_.size(); ++pe) {
            igraph_integer_t fe = pos2full_eid_[pe];
            const double base = base_pos_[(size_t)fe];
            const double hint = (fe < (igraph_integer_t)s.size() ? s[(size_t)fe] : 0.0);
            VECTOR(saved_weights_pos_)[pe] = std::max(1e-12, lp_alpha_ * base / (1.0 + lp_beta_ * hint));
        }
    }
}

inline bool NegativeCycleBatchStream::edge_is_pos(igraph_integer_t eid) const {
    // Patch A: decide by switched signs (robust), not raw weights or neg_hard_
    // Assumes G_.switched_signs[eid] ∈ {+1, -1}
    return G_.switched_weights[eid] > 0.0;
}

void NegativeCycleBatchStream::build_initial_state_() {
    vcount_ = igraph_vcount(&G_.g);
    ecount_ = igraph_ecount(&G_.g);

    base_pos_.assign((size_t)ecount_, 1.0);
    reuse_accum_.assign((size_t)ecount_, 0.0);
    neg_edges_.clear(); neg_edges_.reserve((size_t)ecount_);
    neg_deg_.assign((size_t)vcount_, 0);
    pos_deg_.assign((size_t)vcount_, 0);

    // Set base_pos_ from absolute weights on positive edges and collect negative edges.
    std::vector<double> pos_bases; pos_bases.reserve((size_t)ecount_);
    for (igraph_integer_t eid = 0; eid < ecount_; ++eid) {
        igraph_integer_t u, v; igraph_edge(&G_.g, eid, &u, &v);
        const double w = G_.switched_weights[eid];
        if (w > 0.0) {
            base_pos_[(size_t)eid] = std::max(1e-12, std::fabs(w));
            pos_bases.push_back(base_pos_[(size_t)eid]);
            ++pos_deg_[(size_t)u]; ++pos_deg_[(size_t)v];
        } else if (w < 0.0) {
            neg_edges_.emplace_back((int)u, (int)v);
            ++neg_deg_[(size_t)u]; ++neg_deg_[(size_t)v];
        } else { // zero weight: treat as tiny positive to allow paths
            base_pos_[(size_t)eid] = 1e-12;
            pos_bases.push_back(1e-12);
            ++pos_deg_[(size_t)u]; ++pos_deg_[(size_t)v];
        }
    }
    
    // Adaptive per-vertex cap (store to member)
    tri_cap_per_v_.assign((size_t)vcount_, tri_cap_per_vertex_);
    for (int v = 0; v < (int)vcount_; ++v) {
	    // WHY: high-degree vertices can safely host more triangles
	    // cap = base + ceil(deg⁺ / 64), clamped to [tri_cap, 64]
	    int bump = (pos_deg_[(size_t)v] + 63) / 64;
	    tri_cap_per_v_[(size_t)v] = std::min(64, tri_cap_per_vertex_ + bump);
	}

    // Median of positive base weights
    if (!pos_bases.empty()) {
        const size_t mid = pos_bases.size()/2;
        std::nth_element(pos_bases.begin(), pos_bases.begin()+mid, pos_bases.end());
        med_base_pos_ = std::max(1e-12, pos_bases[mid]);
    } else {
        med_base_pos_ = 1.0;
    }

    // Build positive-only graph and mappings (pos2full / full2pos).
    std::vector<igraph_integer_t> pos_edges; pos_edges.reserve((size_t)ecount_ * 2);
    full2pos_eid_.assign((size_t)ecount_, -1);
    pos2full_eid_.clear(); pos2full_eid_.reserve((size_t)ecount_);

    for (igraph_integer_t eid = 0; eid < ecount_; ++eid) {
        if (!edge_is_pos(eid)) continue;
        igraph_integer_t u, v; igraph_edge(&G_.g, eid, &u, &v);
        full2pos_eid_[(size_t)eid] = (igraph_integer_t)pos2full_eid_.size();
        pos2full_eid_.push_back(eid);
        pos_edges.push_back(u); pos_edges.push_back(v);
    }

    if (g_pos_built_) igraph_destroy(&g_pos_);
    igraph_vector_int_t edges_vec; igraph_vector_int_init(&edges_vec, (long)pos_edges.size());
    for (long i = 0; i < (long)pos_edges.size(); ++i) VECTOR(edges_vec)[i] = (igraph_integer_t)pos_edges[i];
    igraph_create(&g_pos_, &edges_vec, vcount_, /*directed=*/IGRAPH_UNDIRECTED);
    igraph_vector_int_destroy(&edges_vec);
    g_pos_built_ = true;

    // Initialize saved weights for both graphs
    if (saved_weights_pos_init_) igraph_vector_destroy(&saved_weights_pos_);
    igraph_vector_init(&saved_weights_pos_, (long)pos2full_eid_.size());
    for (long pe = 0; pe < (long)pos2full_eid_.size(); ++pe) {
        igraph_integer_t fe = pos2full_eid_[pe];
        VECTOR(saved_weights_pos_)[pe] = base_pos_[(size_t)fe];
    }
    saved_weights_pos_init_ = true;

    if (!saved_weights_init_) {
        igraph_vector_init(&saved_weights_, (long)ecount_);
        for (long i = 0; i < (long)ecount_; ++i) VECTOR(saved_weights_)[i] = 0.0; // only used as within-batch mask
        saved_weights_init_ = true;
    }

    // Triangle usage caps and coverage flags
    tri_used_per_vertex_.assign((size_t)vcount_, 0);
    neg_edge_covered_.assign((size_t)ecount_, 0);

    // Light stats line (optional)
    int negE = (int)neg_edges_.size();
    double avgPosDeg = 0.0; int posV = 0; int hist[8] = {0};
    for (int v = 0; v < (int)vcount_; ++v) {
        const int d = pos_deg_[(size_t)v];
        if (d > 0) { avgPosDeg += d; ++posV; }
        int b = (d<=2?0:d<=4?1:d<=8?2:d<=16?3:d<=32?4:d<=64?5:d<=128?6:7);
        hist[b]++;
    }
    if (posV>0) avgPosDeg /= (double)posV;
    K_tri_per_neg_ = std::max(3, std::min(8, (int)std::sqrt(std::max(1, (int)avgPosDeg))));
    std::cout << "[TRI-DEG] negE=" << negE
              << ", V=" << vcount_ << ", E=" << ecount_
              << ", E[minDeg+]=" << std::fixed << std::setprecision(2) << avgPosDeg
              << ", hist=[0-2]=" << hist[0]
              << " [3-4]=" << hist[1]
              << " [5-8]=" << hist[2]
              << " [9-16]=" << hist[3]
              << " [17-32]=" << hist[4]
              << " [33-64]=" << hist[5]
              << " [65-128]=" << hist[6]
              << " [129+]=" << hist[7]
              << std::endl;
}

bool NegativeCycleBatchStream::next(std::vector<NegativeCycle>& out) {
    out.clear();
    if (finished_) return false;
    
    // --- PROBE counters/timers (per call) ---
    using clock = std::chrono::steady_clock;
    long long ms_dijkstra = 0, ms_check = 0, ms_emit = 0, ms_misc = 0;
    size_t neg_edges_scanned = 0, cycles_emitted_now = 0, triangles_emitted_now = 0;
    long long path_nodes_scanned = 0, pos_edges_on_paths = 0;

    build_mask_for_batch_();

    const size_t before_total = total_found_;

	// === Fast triangle-only first batch when requested ===
    if (use_tri_order_ && batches_emitted_ == 0) {
        struct TriCand { int u,v,w; igraph_integer_t e_aw, e_bw; igraph_integer_t neg_eid; double base_score; int bucket; };
        std::vector<std::vector<TriCand>> per_neg_cands; per_neg_cands.resize(neg_edges_.size());
        std::vector<int> per_neg_best_idx; per_neg_best_idx.resize(neg_edges_.size(), -1);

        // For overlap-aware scoring
        std::vector<int> pos_use_count(pos2full_eid_.size(), 0);
        auto penalized_score = [&](const TriCand& c)->double{
            igraph_integer_t paw = full2pos_eid_[c.e_aw];
            igraph_integer_t pbw = full2pos_eid_[c.e_bw];
            int reuse = 0;
            if (paw >= 0) reuse += pos_use_count[(size_t)paw];
            if (pbw >= 0) reuse += pos_use_count[(size_t)pbw];
            return c.base_score - overlap_penalty_gamma_ * (double)reuse;
        };

        // Helper: bucket by min positive degree of endpoints
        auto bucket_of = [&](int u, int v){
            int md = std::min(pos_deg_[(size_t)u], pos_deg_[(size_t)v]);
            if (md <= 2) return 0; if (md <= 4) return 1; if (md <= 8) return 2; if (md <= 16) return 3;
            if (md <= 32) return 4; if (md <= 64) return 5; if (md <= 128) return 6; return 7;
        };

        // === 1) Build up to K candidates per negative edge ===
        const double EPS = 1e-12; const double ALPHA = 1.0; const double BETA = 0.25; const double lambda_negdeg = 0.02;
        int Tmax = 0; for (int t : neg_tri_vert_) Tmax = std::max(Tmax, t);

        std::vector<char> mark(vcount_, 0); std::vector<int> touched; touched.reserve(128);
        for (size_t idx = 0; idx < neg_edges_.size(); ++idx) {
            const auto& e_neg = neg_edges_[idx];
            const int u = e_neg.first, v = e_neg.second;
            igraph_integer_t neg_eid = -1; igraph_get_eid(&G_.g, &neg_eid, u, v, 0, 1);
            double neg_mag = std::fabs(G_.switched_weights[neg_eid]);

            igraph_vector_int_t nu, nv; igraph_vector_int_init(&nu, 0); igraph_vector_int_init(&nv, 0);
            igraph_neighbors(&G_.g, &nu, u, IGRAPH_ALL);
            igraph_neighbors(&G_.g, &nv, v, IGRAPH_ALL);
            const bool mark_u = (igraph_vector_int_size(&nu) <= igraph_vector_int_size(&nv));
            igraph_vector_int_t &nmark = mark_u ? nu : nv;
            igraph_vector_int_t &nprobe = mark_u ? nv : nu;
            const int a = mark_u ? u : v; const int b = mark_u ? v : u;

            // mark positives from a
            for (int i = 0; i < igraph_vector_int_size(&nmark); ++i) {
                int w = VECTOR(nmark)[i]; if (w == b) continue;
                igraph_integer_t eaw; if (igraph_get_eid(&G_.g, &eaw, a, w, 0, 0) == IGRAPH_SUCCESS && edge_is_pos(eaw)) { mark[w] = 1; touched.push_back(w); }
            }

            std::vector<TriCand> local; local.reserve(8);
            for (int i = 0; i < igraph_vector_int_size(&nprobe); ++i) {
                const int w = VECTOR(nprobe)[i]; if (!mark[w] || w == a) continue;
                igraph_integer_t ebw; if (igraph_get_eid(&G_.g, &ebw, b, w, 0, 0) != IGRAPH_SUCCESS || !edge_is_pos(ebw)) continue;
                igraph_integer_t eaw; if (igraph_get_eid(&G_.g, &eaw, a, w, 0, 0) != IGRAPH_SUCCESS || !edge_is_pos(eaw)) continue;
                
				double lp_aw = 0.0, lp_bw = 0.0, lp_neg = 0.0;
				if (use_lp_weights_) {
				    lp_aw = (eaw < (igraph_integer_t)lp_score_full_.size() ? lp_score_full_[(size_t)eaw] : 0.0);
				    lp_bw = (ebw < (igraph_integer_t)lp_score_full_.size() ? lp_score_full_[(size_t)ebw] : 0.0);
				    igraph_integer_t ne; igraph_get_eid(&G_.g, &ne, u, v, 0, 1);
				    if (ne >= 0 && ne < (igraph_integer_t)lp_score_full_.size()) lp_neg = lp_score_full_[(size_t)ne];
				}

                // base score (no overlap penalty yet)
                double cheap = 1.0/(EPS + base_pos_[eaw]) + 1.0/(EPS + base_pos_[ebw]);
                double central = (Tmax > 0 ? (double)(neg_tri_vert_[u] + neg_tri_vert_[v] + neg_tri_vert_[w])/(double)Tmax : 0.0);
                double negdeg_bonus = lambda_negdeg * ((double)neg_deg_[u] + (double)neg_deg_[v] + (double)neg_deg_[w]);
                double cover_bonus = (!neg_edge_covered_[(size_t)neg_eid] ? 0.25 : 0.0);
				double lp_bonus = 0.4 * (lp_neg + lp_aw + lp_bw);
				double base_score = neg_mag + ALPHA*cheap + BETA*central + negdeg_bonus + cover_bonus + lp_bonus;

                TriCand c{u,v,w,eaw,ebw,neg_eid,base_score,bucket_of(u,v)};
                local.push_back(c);
            }

            // keep top-K by base_score
            std::nth_element(local.begin(), local.begin()+std::min((int)local.size(), K_tri_per_neg_), local.end(),
                             [](const TriCand& A, const TriCand& B){ return A.base_score > B.base_score; });
            if ((int)local.size() > K_tri_per_neg_) local.resize(K_tri_per_neg_);
            per_neg_cands[idx] = std::move(local);
            if (!per_neg_cands[idx].empty()) per_neg_best_idx[idx] = 0;

            for (int w : touched) mark[w] = 0; touched.clear();
            igraph_vector_int_destroy(&nu); igraph_vector_int_destroy(&nv);
        }

        // Build candidate pool with indices
        struct PoolIdx { size_t neg_idx; int cand_idx; double base_score; int bucket; };
        std::vector<PoolIdx> pool; pool.reserve(neg_edges_.size()* (size_t)K_tri_per_neg_);
        for (size_t i = 0; i < per_neg_cands.size(); ++i) for (int j = 0; j < (int)per_neg_cands[i].size(); ++j)
            pool.push_back({i,j,per_neg_cands[i][j].base_score, per_neg_cands[i][j].bucket});
        if (pool.empty()) { ++batches_emitted_; finished_ = true; return false; }

        // Global soft budget B (looser than hard top-B): clamp(ceil(0.12*N), v/2, 4500)
        const int N = (int)pool.size();
        int B = std::clamp((int)std::ceil(0.12 * (double)N), (int)(vcount_/2), 4500);

        // Threshold θ from max/median base_score
        std::vector<double> sc; sc.reserve(pool.size()); for (auto &p: pool) sc.push_back(p.base_score);
        std::nth_element(sc.begin(), sc.begin() + (sc.size()*7)/10, sc.end()); // 70th percentile
		const double theta = sc[(sc.size()*7)/10]; // keep top ~30% by (penalized) score

        // Per-bucket caps (Bucketed top-B)
        int nonempty_b = 0; std::array<int,8> bucket_total{}; bucket_total.fill(0);
        for (auto &p: pool) { bucket_total[(size_t)p.bucket]++; }
        for (int b = 0; b < 8; ++b) if (bucket_total[b] > 0) ++nonempty_b;
        const int per_bucket_cap = std::max(1, (B + nonempty_b - 1) / std::max(1, nonempty_b));
        std::array<int,8> bucket_used{}; bucket_used.fill(0);

		auto violates_vertex_cap = [&](const TriCand& c){
		    if (K_tri_per_neg_ <= 1) return false; // relax when strict per-neg-edge regime
		    return tri_used_per_vertex_[c.u] >= tri_cap_per_v_[c.u] ||
                   tri_used_per_vertex_[c.v] >= tri_cap_per_v_[c.v] ||
                   tri_used_per_vertex_[c.w] >= tri_cap_per_v_[c.w];
		};

        // Round 1: per-neg-edge floor – try to take best per negative edge
        int accepted = 0;
        std::vector<char> pos_edge_hard_used(ecount_, 0); // only used when K=1
        for (size_t i = 0; i < per_neg_cands.size() && accepted < B; ++i) {
            if (per_neg_cands[i].empty()) continue;
            const TriCand &c = per_neg_cands[i][0];

            // Enforce per-vertex cap
			if (violates_vertex_cap(c)) continue;

            // Disjointness handling
            bool overlaps = false;
            if (K_tri_per_neg_ == 1) {
                if (pos_edge_hard_used[(size_t)c.e_aw] || pos_edge_hard_used[(size_t)c.e_bw]) overlaps = true;
            }
            if (overlaps) continue;

            // Accept
            std::vector<Edge> cyc_path; cyc_path.emplace_back(c.u, c.w); cyc_path.emplace_back(c.w, c.v);
            out.emplace_back(neg_edges_[i], std::move(cyc_path));
            ++total_found_; ++accepted; ++triangles_emitted_now;
            ++tri_used_per_vertex_[c.u]; ++tri_used_per_vertex_[c.v]; ++tri_used_per_vertex_[c.w];
            neg_edge_covered_[(size_t)c.neg_eid] = 1;

            // Update usage/penalties
            igraph_integer_t paw = full2pos_eid_[c.e_aw]; igraph_integer_t pbw = full2pos_eid_[c.e_bw];
            if (paw >= 0) ++pos_use_count[(size_t)paw];
            if (pbw >= 0) ++pos_use_count[(size_t)pbw];
            if (K_tri_per_neg_ == 1) { pos_edge_hard_used[(size_t)c.e_aw] = 1; pos_edge_hard_used[(size_t)c.e_bw] = 1; }
            // cross-batch bump
            bump_cross_batch_(c.e_aw, 3); bump_cross_batch_(c.e_bw, 3);
            ++bucket_used[(size_t)c.bucket];
        }

        // Round 2: Thresholded greedy over remaining pool with bucket caps
        std::sort(pool.begin(), pool.end(), [&](const PoolIdx& A, const PoolIdx& Bx){ return A.base_score > Bx.base_score; });
        for (auto &p : pool) {
            if (accepted >= B) break;
            const TriCand &c = per_neg_cands[p.neg_idx][p.cand_idx];
            double s = penalized_score(c);
            if (s < theta) continue;
            if (bucket_used[(size_t)c.bucket] >= per_bucket_cap) continue;

			if (violates_vertex_cap(c)) continue;

            bool overlaps = false;
            if (K_tri_per_neg_ == 1) {
                if (pos_edge_hard_used[(size_t)c.e_aw] || pos_edge_hard_used[(size_t)c.e_bw]) overlaps = true;
            }
            if (overlaps) continue;

            // Accept
            std::vector<Edge> cyc_path; cyc_path.emplace_back(c.u, c.w); cyc_path.emplace_back(c.w, c.v);
            out.emplace_back(neg_edges_[p.neg_idx], std::move(cyc_path));
            ++total_found_; ++accepted; ++triangles_emitted_now;
            ++tri_used_per_vertex_[c.u]; ++tri_used_per_vertex_[c.v]; ++tri_used_per_vertex_[c.w];
            neg_edge_covered_[(size_t)c.neg_eid] = 1;
            igraph_integer_t paw = full2pos_eid_[c.e_aw]; igraph_integer_t pbw = full2pos_eid_[c.e_bw];
            if (paw >= 0) ++pos_use_count[(size_t)paw];
            if (pbw >= 0) ++pos_use_count[(size_t)pbw];
            if (K_tri_per_neg_ == 1) { pos_edge_hard_used[(size_t)c.e_aw] = 1; pos_edge_hard_used[(size_t)c.e_bw] = 1; }
            bump_cross_batch_(c.e_aw, 3); bump_cross_batch_(c.e_bw, 3);
            ++bucket_used[(size_t)c.bucket];
        }

        ++batches_emitted_;
        const bool made_progress = (total_found_ > before_total);
        if (!made_progress) {
            finished_ = true;
        } else if (cover_) {
            // Keep only negative edges not yet covered by a triangle (first round)
            std::vector<Edge> residual; residual.reserve(neg_edges_.size());
            for (size_t i = 0; i < neg_edges_.size(); ++i) {
                igraph_integer_t ne; igraph_get_eid(&G_.g, &ne, neg_edges_[i].first, neg_edges_[i].second, 0, 1);
                if (!neg_edge_covered_[(size_t)ne]) residual.push_back(neg_edges_[i]);
            }
            neg_edges_.swap(residual);
            if (neg_edges_.empty()) finished_ = true; // all covered by triangles
        } else {
            finished_ = true; // only one batch when not in cover mode
        }
        return !out.empty();
    }

    std::vector<Edge> new_disconnected; new_disconnected.reserve(neg_edges_.size());

	// allow k alternate shortest paths per neg edge (soft mask)
	int K_sp_per_neg = 3;                  // knob
	double alt_path_bump = med_base_pos_;  // small positive bump to encourage diversity
    for (const auto& e_neg : neg_edges_) {
        const int uu = e_neg.first, vv = e_neg.second;

		int accepted_here = 0;
		// Track unique positive-paths between (uu,vv) within this batch
        std::unordered_set<uint64_t> seen_paths;  // signature per path
	    for (int rep = 0; rep < K_sp_per_neg; ++rep) {
	        igraph_vector_int_t path; igraph_vector_int_init(&path, 0);
	        const auto t_dij0 = clock::now();
			igraph_get_shortest_path_dijkstra(&g_pos_, &path, nullptr, uu, vv, &saved_weights_pos_, IGRAPH_ALL);
			if (igraph_vector_int_size(&path) < 2) { 
			    igraph_vector_int_destroy(&path); 
			    // Truly disconnected under current weights; stop trying more reps
			    new_disconnected.push_back(e_neg);
			    break; 
			}
	        ms_dijkstra += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t_dij0).count();
	        ++neg_edges_scanned;
	        
			// --- Build edge list once (needed both for duplication check and emit) ---
            std::vector<Edge>                cycle_path_tmp; 
            std::vector<igraph_integer_t>    path_pos_eids; 
            std::vector<igraph_integer_t>    path_full_eids;
            const int nodes_on_path = igraph_vector_int_size(&path);
            cycle_path_tmp.reserve(nodes_on_path);
            path_pos_eids.reserve(nodes_on_path);
            path_full_eids.reserve(nodes_on_path);
            for (int i = 1; i < nodes_on_path; ++i) {
                const int a = VECTOR(path)[i - 1], b = VECTOR(path)[i];
                cycle_path_tmp.emplace_back(a, b);
                igraph_integer_t peid;
                if (igraph_get_eid(&g_pos_, &peid, a, b, 0, 0) == IGRAPH_SUCCESS) {
                    path_pos_eids.push_back(peid);
                    igraph_integer_t feid = pos2full_eid_[peid];
                    path_full_eids.push_back(feid);
                }
            }
            // --- Duplicate detection (signature of pos-edge ids) ---
            const auto t_chk0 = clock::now();
            bool valid = true;
            uint64_t sig = 1469598103934665603ull; // FNV-1a
            for (auto peid : path_pos_eids) { sig ^= (uint64_t)peid; sig *= 1099511628211ull; }
            auto ins = seen_paths.insert(sig);
            if (!ins.second) valid = false; // duplicate path, reject
            ms_check += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t_chk0).count();
 
            if (valid) {        
            	const auto t_emit0 = clock::now();

            	path_nodes_scanned += nodes_on_path;

            	for (size_t i = 0; i < path_pos_eids.size(); ++i) ++pos_edges_on_paths;

            	// Triangle test in pos graph: 2 edges on path
            	bool is_triangle = (nodes_on_path == 3);
            	bool accepted_triangle = false; // track acceptance
            	if (is_triangle) {
            		// Recover third vertex w of triangle (uu --pe1-- w --pe2-- vv)
            		igraph_integer_t pe1 = path_pos_eids[0], pe2 = path_pos_eids[1];
            		igraph_integer_t e1  = pos2full_eid_[pe1], e2 = pos2full_eid_[pe2];

            		igraph_integer_t u1, v1, u2, v2;
            		igraph_edge(&G_.g, e1, &u1, &v1);
            		igraph_edge(&G_.g, e2, &u2, &v2);
            		// w is the common endpoint of e1 and e2 different from uu and vv
            		igraph_integer_t w = (u1 == u2 || u1 == v2) ? u1 : v1;
            		if (w == uu || w == vv) w = (u1 == uu || u1 == vv) ? v1 : u1;

            		// Per-vertex cap vs USAGE
            		if (tri_used_per_vertex_[(int)uu] < tri_cap_per_v_[(int)uu] &&
            				tri_used_per_vertex_[(int)vv] < tri_cap_per_v_[(int)vv] &&
            				tri_used_per_vertex_[(int)w ] < tri_cap_per_v_[(int)w ]) {
            			++tri_used_per_vertex_[(int)uu];
            			++tri_used_per_vertex_[(int)vv];
            			++tri_used_per_vertex_[(int)w];
            			accepted_triangle = true;
            		}
            		// Mark coverage of this neg edge (used by bonuses in future selections)
            		igraph_integer_t neg_eid; igraph_get_eid(&G_.g, &neg_eid, uu, vv, /*directed=*/0, /*error=*/1);
            		neg_edge_covered_[(size_t)neg_eid] = 1;
            	}

            	// Cross-batch penalty: bump base (full) and mirror on pos
            	const int L = nodes_on_path;
            	for (size_t i = 0; i < path_pos_eids.size(); ++i) {
            		igraph_integer_t feid = path_full_eids[i];
            		igraph_integer_t peid = path_pos_eids[i];
            		bump_cross_batch_(feid, L);
            	}

            	out.emplace_back(e_neg, std::move(cycle_path_tmp));
            	++total_found_; ++cycles_emitted_now; ++accepted_here;
            	if (accepted_triangle) ++triangles_emitted_now;

            	// mark neg edge covered for all accepted cycles
            	igraph_integer_t neg_eid; 
            	igraph_get_eid(&G_.g, &neg_eid, uu, vv, /*directed=*/0, /*error=*/1);
            	neg_edge_covered_[(size_t)neg_eid] = 1;

            	for (auto peid : path_pos_eids) {
            		const igraph_integer_t eid = pos2full_eid_[peid]; // map pos→full
            		reuse_accum_[eid] += base_pos_[eid];              // these are all positive by construction
            	}
            	const double RESET_TH = 0.5 * static_cast<double>(ecount_) * med_base_pos_;
            	bool do_reset = false;
            	for (auto peid : path_pos_eids) {
            		const igraph_integer_t eid = pos2full_eid_[peid];
            		if (reuse_accum_[eid] >= RESET_TH) { do_reset = true; break; }
            	}
            	if (do_reset) std::fill(reuse_accum_.begin(), reuse_accum_.end(), 0.0);

            	// Soft mask: encourage a different path next rep
            	for (auto peid : path_pos_eids) VECTOR(saved_weights_pos_)[peid] += alt_path_bump;
            	ms_emit += std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t_emit0).count();
            } else {
            	// Duplicate path -> escalate soft bump to escape tie
            	for (auto peid : path_pos_eids) VECTOR(saved_weights_pos_)[peid] += 2.0 * alt_path_bump;
            }

	        igraph_vector_int_destroy(&path);
	        // Stop early if nothing new could be accepted this edge
	        if (accepted_here == 0) break;
	    }
    }

    ++batches_emitted_;

    const bool made_progress = (total_found_ > before_total);
    if (!made_progress) {
        finished_ = true;
    } else if (cover_) {
        neg_edges_.swap(new_disconnected);
        if (neg_edges_.empty()) finished_ = true;
    } else {
        finished_ = true;
    }
    
    // --- PROBE: per-call profile ---
    // Note: "intersections" are not performed in this enumerator (shortest-path based),
    // so we proxy work via path length sums and number of SP calls.
    std::cout << "[TRI-PROFILE] negE_scanned=" << neg_edges_scanned
              << ", cycles_out=" << cycles_emitted_now
              << ", tri_out=" << triangles_emitted_now
              << ", sum_path_nodes=" << path_nodes_scanned
              << ", pos_edges_on_paths=" << pos_edges_on_paths
              << ", ms_dijkstra=" << ms_dijkstra
              << ", ms_check=" << ms_check
              << ", ms_emit=" << ms_emit
              << ", ms_misc=" << ms_misc
              << ", batches=" << batches_emitted_
              << ", batch_size=" << out.size()
              << std::endl;
	size_t rej_cap=0, rej_bucket=0, rej_overlap=0, rej_theta=0; // TODO: increment in branches
	std::cout << "[TRI-REJECTS] cap=" << rej_cap 
	          << " bucket=" << rej_bucket 
	          << " overlap=" << rej_overlap 
	          << " theta=" << rej_theta << std::endl;

    return !out.empty();
}
