// separation_pipeline.cpp
#include "separation_pipeline.h"
#include <unordered_map>
#include <limits>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <sstream>

// Forward-declare the C hooks with correct linkage so RAII scope can call them.
extern "C" {
void TBB_set_active(void* ctx,
                    void (*emit)(void*, int, double),
                    void (*accept)(void*, int, double),
                    int  (*budget)(void*, int));
void TBB_clear_active();
// NOTE: Do NOT define TBB_budget_override here (provided elsewhere) to avoid ODR conflicts.
}

// TriangleBucketBatch → generic hook bridge (single strong defs, thread-local)
namespace {
struct TBB_VTable {
    void* ctx{nullptr};
    void (*emit)(void*, int, double){nullptr};
    void (*accept)(void*, int, double){nullptr};
    int  (*budget)(void*, int){nullptr};
};
static thread_local TBB_VTable g_tbb_vt{};

// Minimal RAII scope to wire/unwire hooks safely
struct TBBHookScope {
    TBBHookScope(void* ctx,
                 void (*emit)(void*, int, double),
                 void (*accept)(void*, int, double),
                 int  (*budget)(void*, int)) {
        TBB_set_active(ctx, emit, accept, budget);
    }
    ~TBBHookScope() {
        TBB_clear_active();
    }
};

// Small helpers
static inline std::tuple<double,double,double> min_mean_max(const std::vector<double>& v) {
    if (v.empty()) return {0.0,0.0,0.0};
    auto [mn_it, mx_it] = std::minmax_element(v.begin(), v.end());
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / static_cast<double>(v.size());
    return {*mn_it, mean, *mx_it};
}

static inline double median_of(std::vector<double> v) {
    if (v.empty()) return 0.0;
    size_t n = v.size();
    size_t mid = n / 2;
    std::nth_element(v.begin(), v.begin()+mid, v.end());
    double med = v[mid];
    if (n % 2 == 0) {
        auto max_lo = *std::max_element(v.begin(), v.begin()+mid);
        med = 0.5 * (med + max_lo);
    }
    return med;
}

// Put near other utils in this TU (static internal helpers)
static double percentile_of(std::vector<double> v, double q01) {
    if (v.empty()) return 0.0;
    q01 = std::clamp(q01, 0.0, 1.0);
    const size_t k = (size_t)std::floor(q01 * (v.size()-1));
    std::nth_element(v.begin(), v.begin()+k, v.end());
    return v[k];
}

} // namespace

// SP stage: stash only the count here
thread_local int g_sp_cycles_accepted = 0;

extern "C" {
// Called by TriangleBucketBatch/NCB during enumeration/selection
void TBB_on_emit(int edge_id, double used_density) {
    if (g_tbb_vt.emit) g_tbb_vt.emit(g_tbb_vt.ctx, edge_id, used_density);
}
void TBB_on_accept(int edge_id, double density) {
    if (g_tbb_vt.accept) g_tbb_vt.accept(g_tbb_vt.ctx, edge_id, density);
}
// Registration API for callers (pipeline, NCB, etc.)
void TBB_set_active(void* ctx,
                    void (*emit)(void*, int, double),
                    void (*accept)(void*, int, double),
                    int  (*budget)(void*, int)) {
    g_tbb_vt.ctx    = ctx;
    g_tbb_vt.emit   = emit;
    g_tbb_vt.accept = accept;
    g_tbb_vt.budget = budget;
}
void TBB_clear_active() {
    g_tbb_vt = TBB_VTable{};
}
} // extern "C"

// ─────────────────────────────────────────────────────────────
// TriangleCyclePipeline — driver
// ─────────────────────────────────────────────────────────────

// One-shot guard to reuse the *current* switching once (skip re-switching in next round).
// Thread-local so callbacks/threads don't interfere.
static thread_local bool g_skip_switching_once = false;

TriangleCyclePipeline::TriangleCyclePipeline(SeparationPersistent& persistent,
                                             SeparationConfig cfg)
    : S_(persistent), C_(cfg)
{
    // Ensure persistent edge-aligned arrays are sized
//    S_.init_sizes_if_needed(G_);
}

TriangleCyclePipeline::Result
TriangleCyclePipeline::run_round(const SignedGraphForMIP& G_)
{
    // Keep persistent arrays in sync with |E| (defensive)
    S_.init_sizes_if_needed(G_);

    // Reset per-round scratch.
    round_.clear();
	tri_selected_.clear();
	// Initialize uncovered set = all negative full-eids
	uncovered_neg_eids_ = G_.get_edge_eids_bm(-1);

    // Read modes from config
    const bool triangles_enabled = C_.modes.use_triangles;
    const bool sp_enabled        = C_.modes.use_negcycles;
	DBG_PIPE(std::cout << "[PIPELINE] modes: triangles=" << (triangles_enabled?"on":"off")
	                   << " negcycles=" << (sp_enabled?"on":"off") << "\n";);

    const int m = G_.edge_count();
    round_.selected_count_full.assign(m, 0.0);
	round_.used_in_batch_full.assign(m, 0.0);

    g_sp_cycles_accepted = 0;

    // Clear per-round exported keys vector
    accepted_keys_.clear();
	dbg_sample_feids_.clear();

    // === Round-setup & switching probes ===
	// SWITCH probe (including top-5 degrees & sign hash)
	DBG_SWITCH({
        const int pos_edges = G_.edge_count(+1);
        const int neg_edges = G_.edge_count(-1);
        const int m_guard   = m;
        const int delta     = (pos_edges + neg_edges) - m_guard;
	    const uint64_t sign_hash = G_.hash_sign_mask();
        std::cout << "[SWITCH-PROBE] |E+|=" << pos_edges
                  << " |E-|=" << neg_edges
                  << " m="   << m_guard
                  << " guard_delta=" << delta
                  << " sign_hash=0x"  << std::hex << sign_hash << std::dec << "\n";

        // Per-vertex positive/negative degree top-5
        const int n = G_.vertex_count();
        auto top5 = [&](const int sign){
            std::vector<std::pair<int,int>> vv; vv.reserve(n);
            for (int i=0;i<n;++i) vv.emplace_back(G_.degree(i, sign), i);
            std::partial_sort(vv.begin(), vv.begin()+std::min<int>(5,(int)vv.size()), vv.end(),
                              [](auto& a, auto& b){ return (a.first>b.first) || (a.first==b.first && a.second<b.second); });
            std::ostringstream oss;
            int k = std::min<int>(5, (int)vv.size());
            for (int i=0;i<k;++i) oss << " (" << vv[i].second << ":" << vv[i].first << ")";
            return oss.str();
        };
        std::cout << "[SWITCH-PROBE] top5 deg_pos:" << top5(+1) << "\n";
        std::cout << "[SWITCH-PROBE] top5 deg_neg:" << top5(-1) << "\n";
	});

    // === Phase-difference probe (persistent arrays & LP knob) ===
	DBG_PHASE(std::cout << "[PHASE-PROBE] lambda_LP=" << C_.ranking.lambda_LP
                    << " m_pos=" << G_.edge_count(+1) << " m=" << m << "\n";);

    // === 3) Build edge salience (FULL-edge array) ===
    round_.sal_full = G_.edge_saliences();

	// Salience stats (median, zeros%) – compute only if probing
	DBG_PHASE({
	    double sal_min, sal_mean, sal_max;
	    auto [mn,mean,mx] = min_mean_max(round_.sal_full);
	    size_t nz = 0, no = 0; for (double v : round_.sal_full){ if (v>0.0) ++nz; if (v>0.0 && v<1.0) ++no; }
	    std::vector<double> tmp = round_.sal_full;
	    const double med = median_of(tmp);
	    const double p0 = round_.sal_full.empty()?0.0 : 100.0 * (double)(round_.sal_full.size()-nz) / (double)round_.sal_full.size();
	    const double f0 = round_.sal_full.empty()?0.0 : 100.0 * (double)(no) / (double)round_.sal_full.size();
	    std::cout << "[PHASE-PROBE] salience stats min=" << mn
	              << " mean=" << mean << " median=" << med
	              << " max=" << mx << " zeros%=" << p0 << "% fracs%=" << f0 << "%\n";
	});

    // === 4) Build working weights ω′ on E⁺ (includes λ_LP blend) ===
    seed_sp_weights_(G_);

    // === Run stages ===
    if (triangles_enabled && C_.budget.B_tri > 0) {
        triangle_stage_(G_);
    }
    if (sp_enabled) {
        sp_stage_(G_);
    }

    // === Commit persistent updates (ω/H) ===
    if ((triangles_enabled && C_.budget.B_tri > 0) || sp_enabled) {
        DBG_SP({ double sel_sum=0.0; int sel_nz=0; for (double v: round_.selected_count_full){ sel_sum+=v; if (v!=0.0) ++sel_nz; }
          double used_sum=0.0; int used_nz=0; for (double v: round_.used_in_batch_full){ used_sum+=v; if (v!=0.0) ++used_nz; }
          std::cout << "[SP-PROBE] sel_nz="<<sel_nz<<" sel_sum="<<sel_sum
                    << " used_nz_full="<<used_nz<<" used_sum_full="<<used_sum << "\n"; });
    	commit_stage_(G_);
    }

    Result R{};
    R.triangles_accepted = (int)tri_selected_.size();
    R.cycles_accepted    = g_sp_cycles_accepted;
    R.accepted_keys      = std::move(accepted_keys_);
    return R;
}

void TriangleCyclePipeline::seed_sp_weights_(const SignedGraphForMIP& G_) {
    auto& SP = G_.shortest_path_graph();

    const double eps   = std::max(1e-12, C_.weights.omega_eps);
    const double w_cap = (C_.weights.omega_max > 0.0 ? C_.weights.omega_max
                                                     : std::numeric_limits<double>::infinity());

    const auto& omega = S_.omega;        // |E|
    const auto& H     = S_.H;            // |E|
    const auto& sal   = round_.sal_full; // |E|

    // Compile-out all telemetry in experiment mode (why: avoid runtime cost/IO noise).
    DBG_PHASE_DECL(
        double sum_withLP = 0.0, sum_noLP = 0.0, lp_term_sum = 0.0;
        double wb_sum = 0.0, wb_min = +INFINITY, wb_max = -INFINITY; // ω base
        double ht_sum = 0.0, ht_min = +INFINITY, ht_max = -INFINITY; // H[log1p]*λ_hist
        double lp_sum = 0.0, lp_min = +INFINITY, lp_max = -INFINITY; // λ_LP·sal
        double wp_sum = 0.0, wp_min = +INFINITY, wp_max = -INFINITY; // ω′
        int m_pos = 0;
        int samples = 0;
    );

    SP.for_each_edge([&](int u, int v){
        const int fe = G_.edge_index(u, v);

        const double base = (fe >= 0 && fe < (int)omega.size()) ? omega[fe] : 1.0;
        const double rep0 = (fe >= 0 && fe < (int)H.size())     ? std::log1p(std::max(0.0, H[fe])) : 0.0;
        const double rep  = C_.ranking.lambda_hist * rep0;
        const double lp   = (fe >= 0 && fe < (int)sal.size())   ? C_.ranking.lambda_LP * sal[fe]   : 0.0;

        double w_noLP = base + rep;
        double w = std::max(eps, std::min(w_cap, w_noLP - lp));

        SP.set_weight(u, v, w);

        DBG_PHASE({
            ++m_pos;
            sum_withLP += w;
            sum_noLP   += std::max(eps, std::min(w_cap, w_noLP));
            lp_term_sum += lp;

            wb_sum += base; wb_min = std::min(wb_min, base); wb_max = std::max(wb_max, base);
            ht_sum += rep;  ht_min = std::min(ht_min,  rep ); ht_max = std::max(ht_max,  rep );
            lp_sum += lp;   lp_min = std::min(lp_min,  lp  ); lp_max = std::max(lp_max,  lp  );
            wp_sum += w;    wp_min = std::min(wp_min,  w   ); wp_max = std::max(wp_max,  w   );

            if (samples < 3) {
                std::cout << "[OMEGA'-SAMPLE]"
                          << " (u="<<u<<", v="<<v<<", feid="<<fe
                          << ", ω="<<base<<", H="<< (fe<(int)H.size()?H[fe]:0.0)
                          << ", sal="<<(fe<(int)sal.size()?sal[fe]:0.0)
                          << ", ω′="<<w<<")\n";
                ++samples;
            }
        });
    });

    DBG_PHASE({
        const double mean_withLP = (m_pos ? sum_withLP / m_pos : 0.0);
        const double mean_noLP   = (m_pos ? sum_noLP   / m_pos : 0.0);
        const double wb_mean = (m_pos ? wb_sum / m_pos : 0.0);
        const double ht_mean = (m_pos ? ht_sum / m_pos : 0.0);
        const double lp_mean = (m_pos ? lp_sum / m_pos : 0.0);
        const double wp_mean = (m_pos ? wp_sum / m_pos : 0.0);

        std::cout << "[PHASE-PROBE] omega' mean_withLP=" << mean_withLP
                  << " mean_noLP=" << mean_noLP
                  << " mean_delta(withLP-noLP)=" << (mean_withLP - mean_noLP)
                  << " lp_term_sum=" << lp_term_sum
                  << " m_pos=" << m_pos << "\n";

        std::cout << "[OMEGA'-STATS] ω base    min=" << wb_min << " mean=" << wb_mean << " max=" << wb_max << "\n";
        std::cout << "[OMEGA'-STATS] H[log1p]*λ_hist min=" << ht_min << " mean=" << ht_mean << " max=" << ht_max << "\n";
        std::cout << "[OMEGA'-STATS] λ_LP·sal  min=" << lp_min << " mean=" << lp_mean << " max=" << lp_max << "\n";
        std::cout << "[OMEGA'-STATS] ω'        min=" << wp_min << " mean=" << wp_mean << " max=" << wp_max << "\n";
    });
}

void TriangleCyclePipeline::triangle_stage_(const SignedGraphForMIP& G_) {
    auto& SP = G_.shortest_path_graph();

    const int n = G_.vertex_count();
    const int m = G_.edge_count();

    // --------- build anchors & per-anchor intersections as bitmaps ----------
    const auto neg_edges_raw = G_.get_edges(-1);
    const int neg_raw = (int)neg_edges_raw.size();

    // Gather saliences per full-eid to gate anchors
    std::vector<double> neg_sals;
    neg_sals.reserve(neg_raw);
    for (const auto& e : neg_edges_raw) {
        const int feid = G_.edge_index(e);
        neg_sals.push_back(round_.sal_full[feid]);
    }

    double mean = 0.0;
    for (double s : neg_sals) mean += s;
    mean = (neg_raw ? mean / neg_raw : 0.0);

    const double p50 = percentile_of(neg_sals, 0.50);
    const double sal_min_fallback = std::max(p50, mean);

    std::vector<Edge> neg_edges;
    neg_edges.reserve(neg_raw);
    TriangleBucketBatch::CommonAdj pos_adj;                     // vector<GraphCore::Bitmap>
    pos_adj.reserve(neg_raw);

    int nsel = 0;
    for (int i = 0; i < neg_raw; i++) {
        if (neg_sals[i] < sal_min_fallback) continue;

		const auto e = neg_edges_raw[i];
        neg_edges.emplace_back(e);
        // Precompute W = N⁺(u) ∩ N⁺(v) as a bitmap (fast later enumeration)
		auto tri_e = G_.common_neighbors(e, +1);
		tri_e.ior(G_.common_neighbors(e, -1).greater_than(e.second));
        pos_adj.emplace_back(tri_e);
        ++nsel;
    }

    DBG_TBB(std::cout << "[TBB-FILTER] negE_raw=" << neg_raw
                      << " negE_keep=" << nsel
                      << " sal_fallback=" << sal_min_fallback
                      << " (p50=" << p50 << ", mean=" << mean << ")\n";);

    // Params from config (respect defaults)
    TriangleBucketBatch::Params P = C_.to_tbb_params();

    // Annealed triangle budget (effective), unchanged
    if (P.B_tri <= 0) P.B_tri = std::max(1, (int)neg_edges.size());
    const int base_B = C_.budget.B_tri;
    const int B_eff_this_round = (S_.B_tri_cur > 0 ? std::min(base_B, S_.B_tri_cur) : base_B);
    P.B_tri = std::min(P.B_tri, B_eff_this_round);

    TriangleBucketBatch tbb(neg_edges, pos_adj, G_.edge_index(), P);

    auto scorer = [&](TriangleBucketBatch::Candidate& c){
		c.all_neg = G_.is_neg_edge(c.pos_eid_uw);
        const double p = std::max(1.0, C_.ranking.inv_power);
        double inv_hmean = 0.0;
		double w1, w2;
		if (c.all_neg) {
	        w1 = round_.sal_full[c.pos_eid_uw];
	        w2 = round_.sal_full[c.pos_eid_wv];
		} else {
	        w1 = SP.weight(c.pos_eid_uw);
	        w2 = SP.weight(c.pos_eid_wv);
		}
        if (p == 1.0) {
            inv_hmean = 0.5 * (1.0 / w1 + 1.0 / w2);
        } else {
            inv_hmean = 0.5 * (std::pow(w1, -p) + std::pow(w2, -p));
        }
        c.score_primary   = C_.ranking.alpha * inv_hmean + C_.ranking.theta * round_.sal_full[c.neg_eid];
        c.score_secondary = 0.0;
    };

    // Build buckets and probe TBB stats before selection
    tbb.build_buckets(scorer);
    DBG_TBB({
        const auto& ST = tbb.stats();
        std::cout << "[TBB-STATS] negE=" << ST.neg_edges
                  << " buckets_nonempty=" << ST.buckets_nonempty
                  << " tri_candidates=" << ST.candidates
                  << " B_tri=" << ST.B_tri
                  << " B_eff(pre)~" << P.B_tri
                  << "\n";
    });

    // Selection with hook wiring (register → select → clear)
    auto emit_cb = [](void* ctx, int eid, double d) {
        static_cast<TriangleCyclePipeline*>(ctx)->on_emit_(eid, d);
    };
    auto accept_cb = [](void* ctx, int eid, double d) {
        static_cast<TriangleCyclePipeline*>(ctx)->on_accept_(eid, d);
    };
    auto budget_cb = [](void* ctx, int base) -> int {
        return static_cast<TriangleCyclePipeline*>(ctx)->override_budget_(base);
    };

    int selected_count = 0;
    int covered_count  = 0;
    {
        TBBHookScope scope(static_cast<void*>(this), emit_cb, accept_cb, budget_cb);
        #if FM_PROBE_TBB
            // Build only when probes are ON.
            std::vector<int> anchors_with_candidates;            // telemetry-only
        	tri_selected_ = tbb.select(anchors_with_candidates);
        #else
            // Probes OFF: avoid building the anchors vector.
            // We still call select with a dummy sink, but select() will compile-out the pushes (see patch below).
            static std::vector<TriangleBucketBatch::EdgeId> __unused_sink; // static → no per-round allocs
            tri_selected_ = tbb.select(__unused_sink);
        #endif
        
        // Persist SELECTION and compute the set of SELECTED anchors for gating
        DBG_TBB({
        	for (const auto& c : tri_selected_) {
                auto keep_dbg = [&](int fe){
                    if (fe >= 0 && (int)dbg_sample_feids_.size() < 12) dbg_sample_feids_.push_back(fe);
                };
                keep_dbg(c.pos_eid_uw);
                keep_dbg(c.pos_eid_wv);
	        }
        });
        selected_count = (int)tri_selected_.size();
        
        // covered_by_tri := anchors of triangles that were actually SELECTED this round
		for (const auto& c : tri_selected_) {
		    // Mark anchor as covered → clear bit from the "uncovered" set
		    uncovered_neg_eids_.iclear(static_cast<size_t>(c.neg_eid));
		}
        
		DBG_TBB_DECL(
            const int candidates_count = (int)anchors_with_candidates.size();
            const int selected_anchor_count = neg_raw-(int)uncovered_neg_eids_.count();
		)
        DBG_TBB({
            std::cout << "[TBB-GATE-PROBE] candidates=" << candidates_count
                      << " selected_anchors=" << selected_anchor_count
                      << " (delta=" << (candidates_count - selected_anchor_count) << ")\n";
        });

        // After selection, log final TBB stats and acceptance rate
        DBG_TBB({
            const auto& ST = tbb.stats();
            const int B_eff = (ST.B_eff > 0 ? ST.B_eff : P.B_tri);
            const double acc_rate_vs_budget = (B_eff > 0) ? (100.0 * (double)selected_count / (double)B_eff) : 0.0;
            const double acc_rate_vs_buckets = (ST.buckets_nonempty > 0) ? (100.0 * (double)selected_count / (double)ST.buckets_nonempty) : 0.0;
            std::cout << "[TBB-SELECT] selected=" << selected_count
                      << " B_eff=" << B_eff
                      << " acc_rate_vs_budget=" << acc_rate_vs_budget << "% "
                      << "acc_rate_vs_buckets=" << acc_rate_vs_buckets << "% "
                      << "anchors_with_candidates=" << candidates_count
                      << " anchors_selected=" << selected_anchor_count
                      << "\n";
        });
    }

    round_.emitted_triangles = selected_count;

    // ── Update annealing state for next round (EMA on a proxy 'violation') ──
    {
        // Use mean phi as a lightweight proxy (bounded in [0,1]) when true viol not available.
        double v_round = 0.0;
        if (!tri_selected_.empty()) {
            double sum_phi = 0.0;
            for (const auto& c : tri_selected_) sum_phi += round_.sal_full[c.neg_eid];
            v_round = sum_phi / (double)tri_selected_.size();
            if (v_round < 0.0) v_round = 0.0;
            if (v_round > 1.0) v_round = 1.0;
        }
        
        // EMA on proxy violation
        const auto& A   = C_.anneal_tri;
        const double tau = std::clamp(A.tau, 0.0, 1.0);
        S_.bar_v_triangle = (1.0 - tau) * S_.bar_v_triangle + tau * v_round;
        
        // Previous effective budget (boot from base if unset)
        const int base_B = C_.budget.B_tri;
        const int prev_B = (S_.B_tri_cur > 0 ? S_.B_tri_cur : base_B);
        
        // Map to γ in (gamma_min, gamma_max) with higher v → smaller γ
        const double ratio = std::clamp(S_.bar_v_triangle / std::max(1e-12, A.v0), 0.0, 1.0);
        double gamma = A.gamma_max - (A.gamma_max - A.gamma_min) * ratio;
        gamma = std::min(gamma, 0.9995); // never ≥ 1
        
        // Extra shrink if the batch was underutilized
        double util = (P.B_tri > 0) ? (double)selected_count / (double)P.B_tri : 0.0;
        if (util < 0.30) gamma *= 0.97;   // +3% extra shrink
        if (util < 0.15) gamma *= 0.94;   // +6% extra shrink
        
        // Gentle multiplicative step
        int b_next = (int)std::floor(prev_B * gamma + 1e-9);
        
        // --- soft floor + patience-based hard-off ---
        static int zero_streak = 0;
        if (selected_count == 0) ++zero_streak; else zero_streak = 0;
        
        // Keep a small soft floor for a while
        const int soft_floor = std::max(A.B_min, std::max(1, base_B / 64)); // e.g., 512→16
        b_next = std::max(soft_floor, b_next);
        
        // If we keep selecting nothing at/near the floor, turn it off
        const int patience = 1; // rounds at/near floor with zero selections
        if (b_next <= soft_floor && zero_streak >= patience) {
            S_.B_tri_cur = 0; // hard-off
            DBG_TBB(std::cout << "[ANNEAL] hard-off: floor=" << soft_floor
                              << " zero_streak=" << zero_streak << "\n";);
        } else {
            S_.B_tri_cur = b_next;
        }
        
        DBG_TBB(std::cout << "[ANNEAL] prev_B=" << prev_B
                          << " selected=" << selected_count
                          << " v_round~phi=" << v_round
                          << " bar_v_triangle=" << S_.bar_v_triangle
                          << " gamma=" << gamma
                          << " B_tri_next=" << S_.B_tri_cur << "\n";);
    }
  
    // ── Export accepted triangles as fmkey::CycleKey and deduplicate ──
    {
        accepted_keys_.reserve(accepted_keys_.size() + tri_selected_.size());
        for (const auto& r : tri_selected_) {
            fmkey::CycleKey k = fmkey::make_from_triangle(r.u, r.v, r.w);
            // Dedup against persistent sets
            if (S_.in_model_keys.find(k) != S_.in_model_keys.end()) continue;
            if (S_.recent_keys.find(k)   != S_.recent_keys.end())   continue;
            accepted_keys_.push_back(std::move(k));
        }
    }
}

void TriangleCyclePipeline::sp_stage_(const SignedGraphForMIP& G_) {
    g_sp_cycles_accepted = 0;

    // Build uncovered negatives = graph negatives \ covered_neg_eids_
    const auto edge_idx = G_.edge_index();

    auto neg_eids_all = G_.get_edge_eids(-1);
    std::vector<double> neg_sals;
    const int neg_raw = G_.edge_count(-1);
    DBG_SP(std::cout << "--------------------------------------------- neg_raw=" << neg_raw << "\n";);
    neg_sals.reserve(neg_raw);
    double mean = 0.0;
    double max = 0.0, min = 1.0;
    for (int eid = 0; eid < neg_raw; eid++) {
        const int feid = neg_eids_all[eid];
        double s = round_.sal_full[feid];
        neg_sals.push_back(s);
        mean += s;
        if (s > max) max = s;
        if (s < min) min = s;
    }
    mean = (neg_raw ? mean / neg_raw : 0.0);
    
    const double p50 = percentile_of(neg_sals, 0.5);
    const double sal_min_fallback = std::max(p50, mean);
    
    // Build uncovered list:
    const auto& sv = G_.signs_view();
    // ---- order anchors by descending salience ----
    struct Anchor { Edge e; int fe; double w; double sal; };
	std::vector<Anchor> anchors;
	anchors.reserve(static_cast<size_t>(uncovered_neg_eids_.count())); // only uncovered
	for (size_t sfe : uncovered_neg_eids_) {
	    const int fe = static_cast<int>(sfe);
	    const double sal = round_.sal_full[fe];
	    if (sal + 1e-7 < sal_min_fallback) continue;  // numeric guard preserved
	
	    const Edge e = sv[(size_t)fe].points;
	    const double w = (fe >= 0 && fe < (int)S_.omega.size()) ? S_.omega[fe] : 1.0;
	    anchors.push_back(Anchor{e, fe, w, sal});
	}
    std::sort(anchors.begin(), anchors.end(),
              [](const Anchor& A, const Anchor& B){
                  if (A.w != B.w) return A.w < B.w;   // ascending ω
                  if (A.sal != B.sal) return A.sal > B.sal;   // descending salience
                  return A.fe < B.fe;                 // deterministic tie-break
              });
    // Build sorted uncovered list in order
    int kept = 0;
    std::vector<Edge> neg_uncov; neg_uncov.reserve(neg_raw);
    for (size_t i = 0; i < anchors.size(); ++i) {
        neg_uncov.push_back(anchors[i].e);
        ++kept;
    }
    
    DBG_SP({
		const int uncovered = static_cast<int>(uncovered_neg_eids_.count());
		const int covered   = neg_raw - uncovered;
		const double denom  = std::max(1.0, static_cast<double>(uncovered));
        std::cout << "[SP-GATE] neg_raw=" << neg_raw
                  << " covered_by_tri=" << covered
                  << " kept_for_sp=" << kept
                  << " sal_min=" << sal_min_fallback
                  << " keep_rate=" << std::fixed << std::setprecision(2) << (100.0 * (kept / denom))
                  << "% (min=" << min << ", max=" << max << ", p50=" << p50 << ", mean=" << mean << ")\n";
    });
    if (kept == 0) { DBG_SP(std::cout << "[SP-GATE] kept_for_sp=0 → skipping SP stage.\n";); return; }
    
    // --- debug: accepted/used edges sample (FULL eids) for this SP stage ---
    DBG_SP_DECL(std::vector<int> dbg_sp_samples;);
    DBG_SP(dbg_sp_samples.reserve(16););
	DBG_SP_DECL(
        auto dbg_keep = [&](int fe){
            if (fe >= 0 && (int)dbg_sp_samples.size() < 12) dbg_sp_samples.push_back(fe);
        };
    );

    NegativeCycleBatch ncb(G_, neg_uncov, C_.to_ncb_params());

    std::vector<NegativeCycle> batch;
    while (ncb.next(batch)) {
        DBG_SP_DECL(
			std::unordered_map<int,int> lnodes_hist;
        	int accepted_in_batch = 0;
        );
        const int emitted_batch = (int)batch.size();
        g_sp_cycles_accepted += (int)batch.size();

        accepted_keys_.reserve(accepted_keys_.size() + batch.size());
        // --- Accumulate commit signals for SP acceptances ---
        for (const auto& cyc : batch) {
            fmkey::CycleKey k = fmkey::make_from_cycle(cyc);
            if (S_.in_model_keys.find(k) != S_.in_model_keys.end()) continue;
            if (S_.recent_keys.find(k)   != S_.recent_keys.end())   continue;
            accepted_keys_.push_back(std::move(k));

            // Density ~ 1/L where L = number of nodes on the pos path
            const int L_nodes = (int)cyc.pos_edges().size() + 1;
            double sal_nodes = G_.vertex_salience(cyc.neg_edge().first) + G_.vertex_salience(cyc.neg_edge().second);
            for (const auto& pe : cyc.pos_edges()) {
                sal_nodes += G_.vertex_salience(pe.first) + G_.vertex_salience(pe.second);
            }
            DBG_SP({
				++lnodes_hist[L_nodes];
            	++accepted_in_batch;
            });
            const double dens = sal_nodes / (2.0 * std::max(1, L_nodes));

            // Bump along positive path edges (FULL-ids) and mirror to POS-domain used mask
            for (const auto& pe : cyc.pos_edges()) {
                auto feid = edge_idx[pe];
                round_.selected_count_full[(size_t)feid] += dens;
            }

            // Also give density credit to the negative anchor (u,v) reconstructed from path endpoints
            auto neg_eid = edge_idx[cyc.neg_edge()];
            round_.selected_count_full[(size_t)neg_eid] += dens;
        }
        ncb.flush_emitted_to([&](int full_eid, double d){
            if (full_eid >= 0 && full_eid < (int)round_.used_in_batch_full.size())
                round_.used_in_batch_full[(size_t)full_eid] += d;
            // debug sample
                DBG_SP(dbg_keep(full_eid););
        });
        
        DBG_SP({
            std::vector<std::pair<int,int>> hist_pairs;
            hist_pairs.reserve(lnodes_hist.size());
            for (auto &kv : lnodes_hist) hist_pairs.emplace_back(kv.first, kv.second);
            std::sort(hist_pairs.begin(), hist_pairs.end(), [](const auto& a, const auto& b){ return a.first < b.first; });
            std::ostringstream oss;
            oss << "[SP-BATCH] accepted=" << accepted_in_batch << " emitted=" << emitted_batch << " L_nodes_hist=";
            for (const auto& kv : hist_pairs) oss << " " << kv.first << ":" << kv.second;
            std::cout << oss.str() << "\n";
        });
        batch.clear();
    }
    
    // --- debug: print a few accepted/used pos edges with rich ω′ terms ---
    DBG_SP({
        if (!dbg_sp_samples.empty()) {
            auto& SPref       = G_.shortest_path_graph();
            const auto& sv    = G_.signs_view();       // FULL-aligned
            const auto& omega = S_.omega;              // FULL-aligned
            const auto& H     = S_.H;                  // FULL-aligned
            const auto& sal   = round_.sal_full;       // FULL-aligned
            const double lh   = C_.ranking.lambda_hist;
            const double llp  = C_.ranking.lambda_LP;
        
            int shown = 0;
            for (int fe : dbg_sp_samples) {
                if (shown++ >= 6) break;                       // keep it short
                if (fe < 0 || fe >= (int)sv.size()) continue;
        
                const int u = (int)sv[(size_t)fe].points.first;
                const int v = (int)sv[(size_t)fe].points.second;
        
                // best-effort POS id; fall back to -1 if not exposed
                int pos_eid = -1;
                // if your SP exposes a mapping, uncomment one of these:
                // pos_eid = SPref.full2pos(fe);
                // pos_eid = SPref.eid(u, v);
                // otherwise we leave -1 (still printing feid and weights)
        
                const double base     = (fe < (int)omega.size()) ? omega[fe] : 1.0;
                const double H_raw    = (fe < (int)H.size())     ? H[fe]     : 0.0;
                const double H_scaled = lh * std::log1p(std::max(0.0, H_raw));
                const double s        = (fe < (int)sal.size())   ? sal[fe]   : 0.0;
                const double lp_term  = llp * s;
        
                double wprime = 0.0;
                wprime = SPref.weight(u, v); // ω′ as currently seeded/bumped
        
                std::cout << "[OMEGA'-SAMPLE/ACCEPTED]"
                          << " (pos_eid=" << pos_eid << ", feid=" << fe << ")"
                          << " H_raw="    << H_raw
                          << " H_scaled=" << H_scaled
                          << " ω="        << base
                          << " LP="       << lp_term
                          << " ω′="       << wprime
                          << "\n";
            }
        }
    });

    round_.emitted_cycles = g_sp_cycles_accepted;
}

void TriangleCyclePipeline::commit_stage_(const SignedGraphForMIP& G_) {
	auto& SP = G_.shortest_path_graph();

    // ── Persistent updates (ω, H via EMA, pool_count) ───────────
    const int m = G_.edge_count();
    S_.init_sizes_if_needed(G_);

    // Ensure round arrays have size m
    if ((int)round_.selected_count_full.size() != m)
        round_.selected_count_full.assign(m, 0.0);
    if ((int)S_.omega.size() != m)      S_.omega.resize(m, 1.0);
    if ((int)S_.H.size() != m)          S_.H.resize(m, 0.0);

    const double beta_sel = C_.weights.beta_sel;
    const double eps      = std::max(1e-12, C_.weights.omega_eps);
    const double w_cap    = (C_.weights.omega_max > 0.0 ? C_.weights.omega_max
                                                        : std::numeric_limits<double>::infinity());

    const double delta = C_.ema.delta;
    const double kappa = C_.ema.kappa;

    // (A) ω drift across rounds using accepted density
    DBG_PHASE_DECL(
        size_t omega_inc_nz = 0;
        double omega_inc_sum = 0.0;
        size_t omega_clamp_min = 0, omega_clamp_max = 0;
    );

    for (int eid = 0; eid < m; ++eid) {
        double inc = round_.selected_count_full[eid];
        if (inc != 0.0) {
            double w_prev = S_.omega[eid];
            double w = w_prev + beta_sel * inc;
            if (w < eps) { w = eps; DBG_PHASE(++omega_clamp_min;); }
            if (w > w_cap) { w = w_cap; DBG_PHASE(++omega_clamp_max;); }
            S_.omega[eid] = w;
            DBG_PHASE(++omega_inc_nz; omega_inc_sum += (S_.omega[eid] - w_prev););
        } else {
            // keep ω within clamps
            if (S_.omega[eid] < eps) { S_.omega[eid] = eps; DBG_PHASE(++omega_clamp_min;); }
            if (S_.omega[eid] > w_cap) { S_.omega[eid] = w_cap; DBG_PHASE(++omega_clamp_max;); }
        }
	}

    // (B) H EMA: H ← δ·H + ρ·accepted + κ·(emitted−accepted)_pos
    // Start with decay
    std::vector<double> Hnew = S_.H;
    for (double& v : Hnew) v *= delta;
    
    // add acceptance on FULL ids
	const double rho = C_.ema.rho; // new knob, try 0.05–0.1
	for (int eid = 0; eid < m; ++eid)
	    Hnew[eid] += rho * round_.selected_count_full[eid];
    
    // add rejection on FULL ids, only for current positive edges
    const auto& sv = G_.signs_view();
    DBG_PHASE_DECL(size_t H_rej_nz = 0; double H_rej_sum = 0.0;);
    SP.for_each_eid([&](int eid){
        double emitted  = round_.used_in_batch_full[eid];
        double accepted = round_.selected_count_full[eid];
        double rejected = emitted - accepted;
        if (rejected > 0.0) {
            Hnew[eid] += kappa * rejected;
            DBG_PHASE(++H_rej_nz; H_rej_sum += kappa * rejected;);
        }
    });
    S_.H.swap(Hnew);

    // ── PROBE: summarize update magnitudes
    DBG_PHASE(std::cout << "[COMMIT] omega_inc_nz=" << omega_inc_nz
                        << " omega_inc_sum=" << omega_inc_sum
                        << " clamp_min=" << omega_clamp_min
                        << " clamp_max=" << omega_clamp_max
                        << " H_rej_nz=" << H_rej_nz
                        << " H_rej_sum=" << H_rej_sum
                        << "\n";);
}

// ─────────────────────────────────────────────────────────────
// TBB hook handlers → update per-round scratch only
// ─────────────────────────────────────────────────────────────
void TriangleCyclePipeline::on_emit_(int full_eid, double used_density) {
    if (full_eid < 0 || full_eid >= (int)round_.used_in_batch_full.size()) return;
    // Within-batch “used” counter is FULL-aligned now.
    round_.used_in_batch_full[(size_t)full_eid] += used_density;
    // No ω′ bump here — SP working weights are managed by ShortestPathGraph.
}

void TriangleCyclePipeline::on_accept_(int full_eid, double density) {
    if (full_eid < 0 || full_eid >= (int)round_.selected_count_full.size()) return;
    // Per-round density counter (∑ 1/|C| for edges in accepted cuts)
    round_.selected_count_full[full_eid] += density;
}

int TriangleCyclePipeline::override_budget_(int base) const {
    // Use annealed B_tri if provided; otherwise fall back to base
    if (S_.B_tri_cur > 0) return std::min(S_.B_tri_cur, base);
    return base;
}

void TriangleCyclePipeline::resize_round_full_arrays_(int m_full) {
    round_.selected_count_full.assign(m_full, 0.0);
}
