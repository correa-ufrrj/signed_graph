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
#include "triangle_bucket_batch.h"
#include "negative_cycle_batch.h"   // for SP stage
#include "separation_pipeline_tls.h"

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

static inline void log_vec_stats(const char* tag, const std::vector<double>& v) {
    auto [mn,mean,mx] = min_mean_max(v);
    std::cout << tag << " min=" << mn << " mean=" << mean << " max=" << mx << " size=" << v.size() << "\n";
}

template <typename T>
static inline size_t count_nz(const std::vector<T>& a) {
    size_t c=0; for (const auto& x : a) if (x!=T{}) ++c; return c;
}

} // namespace

// SP stage: stash only the count here
thread_local int g_sp_cycles_accepted = 0;

// Reheat pool keys (non-violated cuts to retry on the next fractional round)
thread_local ReheatPool g_reheat_pool;
thread_local std::unordered_set<fmkey::CycleKey, fmkey::CycleKeyHash, fmkey::CycleKeyEq> g_reheat_inflight;

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

TriangleCyclePipeline::TriangleCyclePipeline(SignedGraphForMIP& G,
                                             SeparationPersistent& persistent,
                                             SeparationConfig cfg)
    : G_(G), S_(persistent), C_(cfg)
{
    // Ensure persistent edge-aligned arrays are sized
    S_.init_sizes_if_needed(G_);
    rebuild_pos_maps();
}

// Public: request that the next run_round() *reuses* the graph's current switching.
// Implementation: set thread-local guard; run_round will skip restore+greedy.
void TriangleCyclePipeline::reuse_current_switching_once() {
    g_skip_switching_once = true;
}

void TriangleCyclePipeline::rebuild_pos_maps() {
    const int m = G_.edge_count();

    // Reset view
    Gt_.clear();
    Gt_.full2pos.assign(m, -1);
    Gt_.pos_mask.assign(m, 0);
    Gt_.neg_mask.assign(m, 0);

    // Build masks & maps from the CURRENT switched signature.
    Gt_.pos2full.reserve(m);
    const auto& swv = G_.get_switched_weight();
    for (int eid = 0; eid < m; ++eid) {
        const double sw = (eid < (int)swv.size() ? swv[eid] : 0.0);
        if (sw > 0.0) {
            Gt_.full2pos[eid] = (int)Gt_.pos2full.size();
            Gt_.pos2full.push_back(eid);
            Gt_.pos_mask[eid] = 1;
        } else if (sw < 0.0) {
            Gt_.neg_mask[eid] = 1;
        }
    }
    Gt_.m_pos = (int)Gt_.pos2full.size();

    // ── Map-construction probes
    const int pos_nz = (int)count_nz(Gt_.pos_mask);
    const int neg_nz = (int)count_nz(Gt_.neg_mask);
    std::cout << "[MAP-PROBE] m=" << m
              << " full2pos.size=" << Gt_.full2pos.size()
              << " pos2full.size=" << Gt_.pos2full.size()
              << " pos_mask.nz=" << pos_nz
              << " neg_mask.nz=" << neg_nz
              << " m_pos=" << Gt_.m_pos << "\n";

    // Samples: first/last 3 full→pos where defined
    {
        std::vector<std::pair<int,int>> fp;
        fp.reserve(Gt_.m_pos);
        for (int eid=0; eid<m; ++eid) if (Gt_.full2pos[eid] >= 0) fp.emplace_back(eid, Gt_.full2pos[eid]);
        auto print_pairs = [&](const char* tag){
            std::ostringstream oss;
            oss << tag << " ";
            int take_front = std::min<int>(3, (int)fp.size());
            int take_back  = std::min<int>(3, (int)fp.size() - take_front);
            oss << "front:";
            for (int i=0;i<take_front;++i) oss << " (" << fp[i].first << "→" << fp[i].second << ")";
            oss << " back:";
            for (int i=(int)fp.size()-take_back;i<(int)fp.size();++i) if (i>=0) oss << " (" << fp[i].first << "→" << fp[i].second << ")";
            std::cout << oss.str() << "\n";
        };
        print_pairs("[MAP-PROBE] full→pos");
    }
    {
        std::ostringstream oss;
        oss << "[MAP-PROBE] pos→full front:";
        int take_front = std::min<int>(3, (int)Gt_.pos2full.size());
        for (int i=0;i<take_front;++i) oss << " (" << i << "→" << Gt_.pos2full[i] << ")";
        oss << " back:";
        int take_back = std::min<int>(3, (int)Gt_.pos2full.size()-take_front);
        for (int i=(int)Gt_.pos2full.size()-take_back;i<(int)Gt_.pos2full.size();++i) if (i>=0) oss << " (" << i << "→" << Gt_.pos2full[i] << ")";
        std::cout << oss.str() << "\n";
    }

    // Size per-round arrays accordingly
    resize_round_pos_arrays_(Gt_.m_pos);
    resize_round_full_arrays_(m);
}

TriangleCyclePipeline::Result
TriangleCyclePipeline::run_round(const std::vector<double>* x_hat,
                                 const std::vector<double>* y_hat,
                                 Phase phase)
{
    // Keep persistent arrays in sync with |E| (defensive)
    S_.init_sizes_if_needed(G_);

    // Reset per-round scratch.
    round_.clear();
	tri_selected_.clear();
	covered_neg_eids_.clear();

    const int m = G_.edge_count();
    round_.selected_count_full.assign(m, 0.0);
    g_sp_cycles_accepted = 0;

    // Clear per-round exported keys vector
    accepted_keys_.clear();
	
	// --- Reheat pre-stage: move pool keys to the head of this round ---
	// Copy keys from the reheat pool so the callback evaluates them first at this LP.
	std::size_t reheat_staged = 0;
    g_reheat_inflight.clear();
	if (!g_reheat_pool.empty()) {
	    // Optionally: put a simple ordering (e.g., higher EMA first). For now, copy as-is.
	    for (const auto& kv : g_reheat_pool) {
	        accepted_keys_.push_back(kv.first);
		    g_reheat_inflight.insert(kv.first);
	        ++reheat_staged;
	    }
	    // Decrement TTL for all items we just selected; erase those that hit zero.
	    auto [dec_cnt, erase_cnt] = g_reheat_pool.prune_by_ttl(); // TTL-- and erase zeros
	    std::cout << "[REHEAT] staged=" << reheat_staged
	              << " pool_after_stage=" << g_reheat_pool.size()
	              << " ttl_dec=" << dec_cnt
	              << " ttl_erased=" << erase_cnt << "\n";
	}

    // === 1) Choose switching s and apply it ===
	bool wff_called = false; bool wff_changed = false;
	bool skip_reswitch = g_skip_switching_once;
	if (skip_reswitch) {
	    // One-shot reuse of the *current* switching (set by caller).
	    // Do NOT restore/sign or recompute; just consume and clear the flag.
	    std::cout << "[SWITCH-GUARD] reuse_current_switching_once=1 (skipping re-switch)\n";
	    g_skip_switching_once = false;
	} else {
	    G_.restore_switched_sign();
	    if (phase == Phase::Fractional && x_hat && y_hat) {
	        wff_called = true;
	        wff_changed = G_.weighting_from_fractional(*x_hat, *y_hat);
	    }
	    auto s = G_.greedy_switching();
	    G_.switching_from_partition(s);
	}

    // === 2) Rebuild positive-subgraph maps under the new switching ===
    rebuild_pos_maps();

    // === Round-setup & switching probes ===
    {
        const int pos_edges = Gt_.m_pos;
        const int neg_edges = (int)count_nz(Gt_.neg_mask);
        const int m_guard   = m;
        const int delta     = (pos_edges + neg_edges) - m_guard;
        std::cout << "[SWITCH-PROBE] |E+|=" << pos_edges
                  << " |E-|=" << neg_edges
                  << " m="   << m_guard
                  << " guard_delta=" << delta << "\n";

        // Per-vertex positive/negative degree top-5
        const int n = G_.vertex_count();
        std::vector<int> deg_pos(n,0), deg_neg(n,0);
        for (const auto& kv : G_.edge_index()) {
            const Edge& e = kv.first; int eid = kv.second;
            if (eid < 0 || eid >= m) continue;
            if (Gt_.pos_mask[eid]) { ++deg_pos[e.first]; ++deg_pos[e.second]; }
            else if (Gt_.neg_mask[eid]) { ++deg_neg[e.first]; ++deg_neg[e.second]; }
        }
        auto top5 = [&](const std::vector<int>& d){
            std::vector<std::pair<int,int>> vv; vv.reserve(d.size());
            for (int i=0;i<(int)d.size();++i) vv.emplace_back(d[i], i);
            std::partial_sort(vv.begin(), vv.begin()+std::min<int>(5,(int)vv.size()), vv.end(),
                              [](auto& a, auto& b){ return (a.first>b.first) || (a.first==b.first && a.second<b.second); });
            std::ostringstream oss;
            int k = std::min<int>(5, (int)vv.size());
            for (int i=0;i<k;++i) oss << " (" << vv[i].second << ":" << vv[i].first << ")";
            return oss.str();
        };
        std::cout << "[SWITCH-PROBE] top5 deg_pos:" << top5(deg_pos) << "\n";
        std::cout << "[SWITCH-PROBE] top5 deg_neg:" << top5(deg_neg) << "\n";
    }

    // === Phase-difference probe (persistent arrays & LP knob) ===
    std::cout << "[PHASE-PROBE] phase=" << (phase==Phase::Build?"Build":"Fractional")
              << " lambda_LP=" << C_.ranking.lambda_LP
              << " Gt.m_pos=" << Gt_.m_pos << " m=" << m
              << (wff_called?" wff_called=1":" wff_called=0")
              << (wff_called?(wff_changed?" wff_changed=1":" wff_changed=0"):"")
              << "\n";
    if (phase == Phase::Build) {
        log_vec_stats("[PHASE-PROBE] omega", S_.omega);
        log_vec_stats("[PHASE-PROBE] H", S_.H);
        log_vec_stats("[PHASE-PROBE] pool_count", S_.pool_count);
    }

    // === 3) Build salience (FULL-edge array) ===
    build_salience_(x_hat, y_hat);

    // Compute salience stats (with median & %zero)
    double sal_min, sal_mean, sal_max;
    {
        auto [mn,mean,mx] = min_mean_max(round_.sal_full);
        sal_min = mn; sal_mean = mean; sal_max = mx;
        size_t nz = 0; for (double v : round_.sal_full) if (v>0.0) ++nz;
        const double p0 = round_.sal_full.empty()?0.0 : 100.0 * (double)(round_.sal_full.size()-nz) / (double)round_.sal_full.size();
        std::vector<double> tmp = round_.sal_full;
        double med = median_of(tmp);
        std::cout << "[PHASE-PROBE] salience stats min=" << sal_min
                  << " mean=" << sal_mean
                  << " median=" << med
                  << " max=" << sal_max
                  << " zeros%=" << p0 << "%\n";
    }

    // === 4) Build working weights ω′ on E⁺ (includes λ_LP blend) ===
    build_omega_prime_pos_();

    // Build-phase gating confirmation: λ_LP should have no effect (λ=0 or sal all zeros)
    if (phase == Phase::Build) {
        bool lambda_zero = (std::abs(C_.ranking.lambda_LP) <= 0.0);
        bool sal_all_zero = std::all_of(round_.sal_full.begin(), round_.sal_full.end(),
                                        [](double v){ return v==0.0; });
        std::cout << "[BUILD-GUARD] LP_effect_off=" << ((lambda_zero||sal_all_zero)?1:0)
                  << " reason=" << (lambda_zero ? "lambda_LP=0" : (sal_all_zero ? "zero-salience" : "LP-active-WARNING"))
                  << "\n";
    }

    // === Run stages ===
    pre_enumeration_reheat_();
    triangle_stage_();
    sp_stage_();

    // === Commit persistent updates (ω/H/pool_count) ===
    commit_stage_();

    Result R{};
    R.triangles_accepted = (int)tri_selected_.size();
    R.cycles_accepted    = g_sp_cycles_accepted;
    R.accepted_keys      = std::move(accepted_keys_);
    return R;
}

void TriangleCyclePipeline::build_salience_(const std::vector<double>* x_hat,
                                            const std::vector<double>* y_hat)
{
    const int m = G_.edge_count();
    round_.sal_full.assign(m, 0.0);

    // If both x̂ (per-vertex) and ŷ (per-edge) are available, compute:
    // sal(uv) = min{1, 4 * | x[u]*x[v] - y[uv] | }.
    if (x_hat && y_hat) {
        const auto& x = *x_hat;
        const auto& y = *y_hat;

        const int n = G_.vertex_count();
        if ((int)x.size() >= n && (int)y.size() >= m) {
            const auto& sv = G_.signs_view();
            for (int eid = 0; eid < m; ++eid) {
                const int u = (int)sv[eid].points.first;
                const int v = (int)sv[eid].points.second;

                // safe fetch + light clamp to [0,1] (helps with tiny numerical drift)
                const double xu  = std::min(1.0, std::max(0.0, x[u]));
                const double xv  = std::min(1.0, std::max(0.0, x[v]));
                const double yuv = std::min(1.0, std::max(0.0, y[eid]));

                double s = 4.0 * std::fabs(xu * xv - yuv);
                if (s > 1.0) s = 1.0;            // clamp to [0,1]
                // s >= 0 by construction
                round_.sal_full[eid] = s;
            }
            std::cout << "[SALIENCE] source=spec (4*|x_u x_v - y_uv|, clamped)\n";
            return;
        } else {
            std::cout << "[SALIENCE-WARN] insufficient sizes for x̂/ŷ "
                      << "(x=" << (int)x.size() << " vs n=" << G_.vertex_count()
                      << ", y=" << (int)y.size() << " vs m=" << m
                      << "). Falling back.\n";
        }
    }

    // Fallback: use graph-prepared salience if present (e.g., from weighting_from_fractional()).
    const std::vector<double>& gsal = G_.edge_salience_view();
    if ((int)gsal.size() == m) {
        round_.sal_full = gsal; // copy
        std::cout << "[SALIENCE] source=edge_salience_view() (graph-prepared)\n";
        return;
    }

    // Last resort: keep zeros.
    std::cout << "[SALIENCE] source=none (zeros)\n";
}

void TriangleCyclePipeline::build_omega_prime_pos_()
{
    // ω′(uv) = max{ ε, ω(uv) + λ_hist·log(1+H(uv)) − λ_LP·sal(uv) } on E⁺.
    const int m_pos = Gt_.m_pos;
    if (m_pos <= 0) {
        round_.omega_prime_pos.clear();
        round_.used_in_batch_pos.clear();
        return;
    }

    round_.omega_prime_pos.assign(m_pos, 1.0);
    round_.used_in_batch_pos.assign(m_pos, 0.0);

    const double eps   = std::max(1e-12, C_.weights.omega_eps);
    const double w_cap = (C_.weights.omega_max > 0.0 ? C_.weights.omega_max
                                                     : std::numeric_limits<double>::infinity());

    const auto& omega = S_.omega;        // |E|
    const auto& H     = S_.H;            // |E|
    const auto& sal   = round_.sal_full; // |E|


	// --- Probe: clamp & blend parameters
	std::cout << std::setprecision(6) << std::fixed;
	std::cout << "[OMEGA'-CLAMP] eps=" << eps
	          << " omega_max=" << (std::isfinite(w_cap) ? w_cap : -1.0) << (std::isfinite(w_cap) ? "" : " (inf)")
	          << " lambda_hist=" << C_.ranking.lambda_hist
	          << " lambda_LP="   << C_.ranking.lambda_LP
	          << " m_pos="       << m_pos << "\n";

	// Probing accumulators & per-term collections (for min/mean/max)
	double sum_withLP = 0.0, sum_noLP = 0.0, lp_term_sum = 0.0;
	
	std::vector<double> w_base;    w_base.reserve(m_pos);
	std::vector<double> h_term;    h_term.reserve(m_pos);  // λ_hist·log1p(H)
	std::vector<double> lp_term;   lp_term.reserve(m_pos); // λ_LP·sal
	std::vector<double> w_prime;   w_prime.reserve(m_pos);
	
	for (int pid = 0; pid < m_pos; ++pid) {
	    const int eid = Gt_.pos2full[pid];
	
	    const double base = (eid >= 0 && eid < (int)omega.size()) ? omega[eid] : 1.0;
	    const double rep0 = (eid >= 0 && eid < (int)H.size())     ? std::log1p(std::max(0.0, H[eid])) : 0.0;
	    const double rep  = C_.ranking.lambda_hist * rep0;
	    const double lp   = (eid >= 0 && eid < (int)sal.size())   ? C_.ranking.lambda_LP * sal[eid] : 0.0;
	
	    double w_noLP = base + rep;
	    double w = w_noLP - lp;
	    if (w < eps) w = eps;
	    if (w > w_cap) w = w_cap;
	
	    round_.omega_prime_pos[pid] = w;
	
	    // For stats
	    w_base.push_back(base);
	    h_term.push_back(rep);
	    lp_term.push_back(lp);
	    w_prime.push_back(w);
	
	    sum_withLP += w;
	    sum_noLP   += std::max(eps, std::min(w_cap, w_noLP));
	    lp_term_sum += lp;
	}

    const double mean_withLP = (m_pos>0) ? sum_withLP / (double)m_pos : 0.0;
    const double mean_noLP   = (m_pos>0) ? sum_noLP   / (double)m_pos : 0.0;

    std::cout << "[PHASE-PROBE] omega' mean_withLP=" << mean_withLP
              << " mean_noLP=" << mean_noLP
              << " mean_delta(withLP-noLP)=" << (mean_withLP - mean_noLP)
              << " lp_term_sum=" << lp_term_sum
              << " m_pos=" << m_pos << "\n";

	// --- Probe: term-wise ranges
	auto [wb_min, wb_mean, wb_max] = min_mean_max(w_base);
	auto [ht_min, ht_mean, ht_max] = min_mean_max(h_term);
	auto [lp_min, lp_mean, lp_max] = min_mean_max(lp_term);
	auto [wp_min, wp_mean, wp_max] = min_mean_max(w_prime);
	
	std::cout << "[OMEGA'-STATS] ω base    min=" << wb_min << " mean=" << wb_mean << " max=" << wb_max << "\n";
	std::cout << "[OMEGA'-STATS] H[log1p]*λ min=" << ht_min << " mean=" << ht_mean << " max=" << ht_max << "\n";
	std::cout << "[OMEGA'-STATS] λ_LP·sal  min=" << lp_min << " mean=" << lp_mean << " max=" << lp_max << "\n";
	std::cout << "[OMEGA'-STATS] ω'        min=" << wp_min << " mean=" << wp_mean << " max=" << wp_max << "\n";
	
	// --- Probe: blend deltas
	std::cout << "[OMEGA'-BLEND] mean_withLP=" << mean_withLP
	          << " mean_noLP=" << mean_noLP
	          << " mean_delta(withLP-noLP)=" << (mean_withLP - mean_noLP)
	          << " lp_term_sum=" << lp_term_sum
	          << " m_pos=" << m_pos << "\n";

	// --- Probe: sample tuples (eid, ω, H, sal, ω′)
	{
	    std::ostringstream oss;
	    oss << "[OMEGA'-SAMPLE]";
	    int take = std::min(3, m_pos);
	    for (int i = 0; i < take; ++i) {
	        int eid = Gt_.pos2full[i];
	        double w0 = (eid >=0 && eid < (int)omega.size()? omega[eid] : 0.0);
	        double H0 = (eid >=0 && eid < (int)H.size()? H[eid] : 0.0);
	        double s0 = (eid >=0 && eid < (int)sal.size()? sal[eid] : 0.0);
	        double wp = round_.omega_prime_pos[i];
	        oss << " (eid=" << eid << ", ω=" << w0 << ", H=" << H0 << ", sal=" << s0 << ", ω′=" << wp << ")";
	    }
	    std::cout << oss.str() << "\n";
	}
	
	// --- Domain check: ω′ only over E⁺ and sane values
	size_t nan_cnt=0, inf_cnt=0, below_eps=0, above_cap=0;
	for (double w : round_.omega_prime_pos) {
	    if (!std::isfinite(w)) { if (std::isnan(w)) ++nan_cnt; else ++inf_cnt; }
	    if (w < eps - 1e-12) ++below_eps;
	    if (std::isfinite(w_cap) && w > w_cap + 1e-12) ++above_cap;
	}
	const bool size_ok = ((int)round_.omega_prime_pos.size() == m_pos);
	int map_anomalies = 0;
	// verify every pos edge has a valid mapping
	for (int pid = 0; pid < m_pos; ++pid) {
	    int eid = Gt_.pos2full[pid];
	    if (eid < 0 || eid >= (int)Gt_.full2pos.size() || Gt_.full2pos[eid] != pid) ++map_anomalies;
	}
	std::cout << "[OMEGA'-DOMAIN] size=" << round_.omega_prime_pos.size()
	          << " m_pos=" << m_pos
	          << " size_ok=" << (size_ok?1:0)
	          << " map_anomalies=" << map_anomalies
	          << " nan=" << nan_cnt
	          << " inf=" << inf_cnt
	          << " below_eps=" << below_eps
	          << " above_cap=" << above_cap << "\n";
}

void TriangleCyclePipeline::pre_enumeration_reheat_() {
// (Reserved for later steps.)
}

void TriangleCyclePipeline::triangle_stage_() {
    // === Triangle-first via TriangleBucketBatch ===
    tri_selected_.clear();
    covered_neg_eids_.clear();

    const int n = G_.vertex_count();
    const int m = G_.edge_count();

    // Edge-index map (min(u,v),max(u,v)) -> full eid
    std::unordered_map<long long, int> key2eid;
    key2eid.reserve((size_t)m * 2);
    auto key64 = [](int a, int b) -> long long {
        if (a > b) std::swap(a,b);
        return ( (static_cast<long long>(a) << 32) ^ static_cast<unsigned long long>(b) );
    };
    for (const auto& kv : G_.edge_index()) {
        const Edge& e = kv.first; int eid = kv.second;
        key2eid[key64(e.first, e.second)] = eid;
    }

    // Positive adjacency + list of negative anchors under current switching
    TriangleBucketBatch::PosAdj pos_adj(n);
    std::vector<std::pair<int,int>> neg_edges;
    neg_edges.reserve(m/2);
    for (const auto& kv : G_.edge_index()) {
        const Edge& e = kv.first; int eid = kv.second;
        if (eid < 0 || eid >= (int)Gt_.pos_mask.size()) continue;
        if (Gt_.pos_mask[eid]) {
            pos_adj[e.first].push_back(e.second);
            pos_adj[e.second].push_back(e.first);
        } else if (Gt_.neg_mask[eid]) {
            neg_edges.emplace_back(e.first, e.second);
        }
    }

    // Params from config (respect defaults)
    TriangleBucketBatch::Params P = C_.to_tbb_params();

    // Annealed triangle budget (effective)
    if (P.B_tri <= 0) P.B_tri = std::max(1, (int)neg_edges.size());
    const int P_before = P.B_tri;
    if (S_.B_tri_cur > 0) P.B_tri = std::min(P.B_tri, S_.B_tri_cur);

	const int base_B = C_.budget.B_tri;
	// Use annealed budget if available; otherwise fall back to base on first round(s)
	const int B_eff_this_round = (S_.B_tri_cur > 0 ? std::min(base_B, S_.B_tri_cur) : base_B);
	
	P.B_tri = std::min(P.B_tri, B_eff_this_round); // keep existing min w/ TBB param

    TriangleBucketBatch tbb(neg_edges, pos_adj, key2eid, P);

	// Scorer: primary = α * (1 / Hmean(ω′ on positive edges of tri)) + θ·φ
	auto scorer = [&](TriangleBucketBatch::Candidate& c){
	    const int pid_uw = (c.pos_eid_uw >= 0 && c.pos_eid_uw < (int)Gt_.full2pos.size())
	                       ? Gt_.full2pos[c.pos_eid_uw] : -1;
	    const int pid_wv = (c.pos_eid_wv >= 0 && c.pos_eid_wv < (int)Gt_.full2pos.size())
	                       ? Gt_.full2pos[c.pos_eid_wv] : -1;
	
	    const auto inv = [&](int pid)->double{
	        if (pid < 0 || pid >= (int)round_.omega_prime_pos.size()) return 0.0;
	        return 1.0 / std::max(C_.weights.omega_eps, round_.omega_prime_pos[pid]);
	    };
	
	    // inverse harmonic mean (average of inverses)
	    double inv_sum = 0.0; int k = 0;
	    if (pid_uw >= 0) { inv_sum += inv(pid_uw); ++k; }
	    if (pid_wv >= 0) { inv_sum += inv(pid_wv); ++k; }
	    const double inv_hmean = (k > 0) ? (inv_sum / (double)k) : 0.0;
	
	    auto sal = [&](int eid){
	        return (eid >= 0 && eid < (int)round_.sal_full.size()) ? round_.sal_full[eid] : 0.0;
	    };
	    c.phi = (sal(c.neg_eid) + sal(c.pos_eid_uw) + sal(c.pos_eid_wv)) / 3.0;
	
	    c.score_primary   = C_.ranking.alpha * inv_hmean + C_.ranking.theta * c.phi;
	    c.score_secondary = 0.0;
	    c.viol = 0.0; // reserved for later violation-aware scoring
	};

    // Build buckets and probe TBB stats before selection
    tbb.build_buckets(scorer);
    {
        const auto& ST = tbb.stats();
        std::cout << "[TBB-STATS] negE=" << ST.neg_edges
                  << " buckets_nonempty=" << ST.buckets_nonempty
                  << " tri_candidates=" << ST.candidates
                  << " B_tri=" << ST.B_tri
                  << " B_eff(pre)~" << P.B_tri
                  << "\n";
    }

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
        std::vector<int> covered; // full eids of negative anchors with nonempty buckets
        const auto& sel = tbb.select(covered);

        // Persist outputs in round storage
        covered_neg_eids_.assign(covered.begin(), covered.end());
        covered_count = (int)covered_neg_eids_.size();

        tri_selected_.reserve(sel.size());
        for (const auto& c : sel) {
            TriAcc r;
            r.u = c.u; r.v = c.v; r.w = c.w;
            r.neg_eid = c.neg_eid;
            r.pos_eid_uw = c.pos_eid_uw;
            r.pos_eid_wv = c.pos_eid_wv;
            r.score_primary = c.score_primary;
            r.score_secondary = c.score_secondary;
            r.viol = c.viol; r.phi = c.phi;
            tri_selected_.push_back(std::move(r));
        }
        selected_count = (int)tri_selected_.size();
    }

    // After selection, log final TBB stats and acceptance rate
    {
        const auto& ST = tbb.stats();
        const int B_eff = (ST.B_eff > 0 ? ST.B_eff : P.B_tri);
        const double acc_rate_vs_budget = (B_eff > 0) ? (100.0 * (double)selected_count / (double)B_eff) : 0.0;
        const double acc_rate_vs_buckets = (ST.buckets_nonempty > 0) ? (100.0 * (double)selected_count / (double)ST.buckets_nonempty) : 0.0;
        std::cout << "[TBB-SELECT] selected=" << selected_count
                  << " B_eff=" << B_eff
                  << " acc_rate_vs_budget=" << acc_rate_vs_budget << "% "
                  << "acc_rate_vs_buckets=" << acc_rate_vs_buckets << "% "
                  << "covered_neg=" << covered_count
                  << "\n";
    }

    round_.emitted_triangles = selected_count;

    // ── Update annealing state for next round (EMA on a proxy 'violation') ──
    {
        // Use mean phi as a lightweight proxy (bounded in [0,1]) when true viol not available.
        double v_round = 0.0;
        if (!tri_selected_.empty()) {
            double sum_phi = 0.0;
            for (const auto& r : tri_selected_) sum_phi += r.phi;
            v_round = sum_phi / (double)tri_selected_.size();
            if (v_round < 0.0) v_round = 0.0;
            if (v_round > 1.0) v_round = 1.0;
        }
		
		// EMA on proxy violation
		const auto& A   = C_.anneal_tri;
		const double tau = std::clamp(A.tau, 0.0, 1.0);
		S_.bar_v_triangle = (1.0 - tau) * S_.bar_v_triangle + tau * v_round;
		
		// Previous effective budget (boot from base if unset)
		const int   base_B = C_.budget.B_tri;
		const int prev_B = (S_.B_tri_cur > 0 ? S_.B_tri_cur : base_B);
		
		// Map to γ in (gamma_min, gamma_max) with higher v → smaller γ
		const double ratio = std::clamp(S_.bar_v_triangle / std::max(1e-12, A.v0), 0.0, 1.0);
		double gamma = A.gamma_max - (A.gamma_max - A.gamma_min) * ratio;
		gamma = std::min(gamma, 0.9995); // never ≥ 1
		
		// Gentle multiplicative step
		int b_next = (int)std::floor(prev_B * gamma + 1e-9);
		
		// --- soft floor + patience-based hard-off ---
		static int zero_streak = 0;
		if (selected_count == 0) ++zero_streak; else zero_streak = 0;
		
		// Keep a small soft floor for a while
		const int soft_floor = std::max(A.B_min, std::max(1, base_B / 32)); // e.g., 512→16
		b_next = std::max(soft_floor, b_next);
		
		// If we keep selecting nothing at/near the floor, turn it off
		const int patience = 3; // rounds at/near floor with zero selections
		if (b_next <= soft_floor && zero_streak >= patience) {
		    S_.B_tri_cur = 0; // hard-off
		    std::cout << "[ANNEAL] hard-off: floor=" << soft_floor
		              << " zero_streak=" << zero_streak << "\n";
		} else {
		    S_.B_tri_cur = b_next;
		}
		
		std::cout << "[ANNEAL] prev_B=" << prev_B
		          << " selected=" << selected_count
		          << " v_round~phi=" << v_round
		          << " bar_v_triangle=" << S_.bar_v_triangle
		          << " gamma=" << gamma
		          << " B_tri_next=" << S_.B_tri_cur << "\n";
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

void TriangleCyclePipeline::sp_stage_() {
    g_sp_cycles_accepted = 0;

    // Build uncovered negatives = graph negatives \ covered_neg_eids_
    const auto edge_idx = G_.edge_index();
    std::unordered_set<int> covered_set(covered_neg_eids_.begin(), covered_neg_eids_.end());

    std::vector<Edge> neg_uncov;
    neg_uncov.reserve(covered_set.size() + 16);

    for (const auto& kv : G_.edge_index()) {
        const Edge& e = kv.first; int feid = kv.second;
        if (!Gt_.neg_mask[feid]) continue;      // only negatives under current switching
        if (covered_set.count(feid)) continue;  // skip anchors covered by triangle stage
        neg_uncov.push_back(e);
    }

    NegativeCycleBatch ncb(G_, neg_uncov, /*cover=*/false, C_.to_ncb_params());

    std::vector<NegativeCycle> batch;
    while (ncb.next(batch)) {
        g_sp_cycles_accepted += (int)batch.size();

        accepted_keys_.reserve(accepted_keys_.size() + batch.size());
        for (const auto& r : batch) {
            fmkey::CycleKey k = fmkey::make_from_cycle(r);
            if (S_.in_model_keys.find(k) != S_.in_model_keys.end()) continue;
            if (S_.recent_keys.find(k)   != S_.recent_keys.end())   continue;
            accepted_keys_.push_back(std::move(k));
        }
        batch.clear();
    }

    round_.emitted_cycles = g_sp_cycles_accepted;
}

void TriangleCyclePipeline::commit_stage_() {
    // ── Persistent updates (ω, H via EMA, pool_count) ───────────
    const int m = G_.edge_count();
    S_.init_sizes_if_needed(G_);

    // Ensure round arrays have size m
    if ((int)round_.selected_count_full.size() != m)
        round_.selected_count_full.assign(m, 0.0);
    if ((int)S_.omega.size() != m)      S_.omega.resize(m, 1.0);
    if ((int)S_.H.size() != m)          S_.H.resize(m, 0.0);
    if ((int)S_.pool_count.size() != m) S_.pool_count.resize(m, 0.0);

    const double beta_sel = C_.weights.beta_sel;
    const double eps      = std::max(1e-12, C_.weights.omega_eps);
    const double w_cap    = (C_.weights.omega_max > 0.0 ? C_.weights.omega_max
                                                        : std::numeric_limits<double>::infinity());

    const double delta = S_.ema_delta;
    const double rho   = S_.ema_rho;
    const double kappa = S_.ema_kappa;

    // (A) ω drift across rounds using accepted density
    size_t omega_inc_nz = 0;
    double omega_inc_sum = 0.0;
    size_t omega_clamp_min = 0, omega_clamp_max = 0;

    for (int eid = 0; eid < m; ++eid) {
        double inc = round_.selected_count_full[eid];
        if (inc != 0.0) {
            double w_prev = S_.omega[eid];
            double w = w_prev + beta_sel * inc;
            if (w < eps) { w = eps; ++omega_clamp_min; }
            if (w > w_cap) { w = w_cap; ++omega_clamp_max; }
            S_.omega[eid] = w;
            ++omega_inc_nz;
            omega_inc_sum += (S_.omega[eid] - w_prev);
        } else {
            // keep ω within clamps
            if (S_.omega[eid] < eps) { S_.omega[eid] = eps; ++omega_clamp_min; }
            if (S_.omega[eid] > w_cap) { S_.omega[eid] = w_cap; ++omega_clamp_max; }
        }
    }

    // (B) H EMA: H ← δ·H + ρ·accepted + κ·(emitted−accepted)_pos
    // Start with decay
    std::vector<double> Hnew = S_.H;
    for (double& v : Hnew) v *= delta;

    // + ρ·accepted over ALL edges (neg & pos)
    size_t H_acc_nz = 0;
    for (int eid = 0; eid < m; ++eid) {
        double acc = round_.selected_count_full[eid];
        if (acc != 0.0) { Hnew[eid] += rho * acc; ++H_acc_nz; }
    }

    // Compute (emitted − accepted) on POS edges via pos mapping
    size_t H_rej_nz = 0;
    double H_rej_sum = 0.0;
    if (Gt_.m_pos > 0 && !round_.used_in_batch_pos.empty()) {
        // Map accepted FULL-edge densities to POS indices
        std::vector<double> acc_pos(Gt_.m_pos, 0.0);
        for (int pid = 0; pid < Gt_.m_pos; ++pid) {
            int eid = Gt_.pos2full[pid];
            if (eid >= 0 && eid < m) acc_pos[pid] = round_.selected_count_full[eid];
        }
        for (int pid = 0; pid < Gt_.m_pos; ++pid) {
            const double emitted  = (pid < (int)round_.used_in_batch_pos.size()) ? round_.used_in_batch_pos[pid] : 0.0;
            const double accepted = acc_pos[pid];
            double rejected = emitted - accepted;
            if (rejected < 0.0) rejected = 0.0;
            if (rejected > 0.0) {
                int eid = Gt_.pos2full[pid];
                if (eid >= 0 && eid < m) {
                    Hnew[eid] += kappa * rejected;
                    ++H_rej_nz;
                    H_rej_sum += (kappa * rejected);
                }
            }
        }
    }
    S_.H.swap(Hnew);

    // (C) pool_count: accumulate accepted density on edges of accepted cuts
    size_t pool_inc_nz = 0;
    double pool_inc_sum = 0.0;
    for (int eid = 0; eid < m; ++eid) {
        double inc = round_.selected_count_full[eid];
        if (inc != 0.0) { S_.pool_count[eid] += inc; ++pool_inc_nz; pool_inc_sum += inc; }
    }

    // ── PROBE: summarize update magnitudes
    std::cout << "[COMMIT] omega_inc_nz=" << omega_inc_nz
              << " omega_inc_sum=" << omega_inc_sum
              << " clamp_min=" << omega_clamp_min
              << " clamp_max=" << omega_clamp_max
              << " H_acc_nz=" << H_acc_nz
              << " H_rej_nz=" << H_rej_nz
              << " H_rej_sum=" << H_rej_sum
              << " pool_inc_nz=" << pool_inc_nz
              << " pool_inc_sum=" << pool_inc_sum
              << "\n";
}

// ─────────────────────────────────────────────────────────────
// TBB hook handlers → update per-round scratch only
// ─────────────────────────────────────────────────────────────
void TriangleCyclePipeline::on_emit_(int full_eid, double used_density) {
    if (full_eid < 0 || full_eid >= (int)Gt_.full2pos.size()) return;
    const int pid = Gt_.full2pos[full_eid];
    if (pid < 0 || pid >= (int)round_.omega_prime_pos.size()) return;
    round_.used_in_batch_pos[pid] = used_density;
    // ω′ ← ω′ + β_emit · used_density (clamped) — working (within-round) update
    double w = round_.omega_prime_pos[pid] + C_.weights.beta_emit * used_density;
    if (w < C_.weights.omega_eps) w = C_.weights.omega_eps;
    if (C_.weights.omega_max > 0.0 && w > C_.weights.omega_max) w = C_.weights.omega_max;
    round_.omega_prime_pos[pid] = w;
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

void TriangleCyclePipeline::resize_round_pos_arrays_(int m_pos) {
    round_.omega_prime_pos.assign(m_pos, 1.0);
    round_.used_in_batch_pos.assign(m_pos, 0.0);
}

void TriangleCyclePipeline::resize_round_full_arrays_(int m_full) {
    round_.selected_count_full.assign(m_full, 0.0);
    round_.sal_full.assign(m_full, 0.0);
}
