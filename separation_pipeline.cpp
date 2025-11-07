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

// Put near other utils in this TU (static internal helpers)
static double percentile_of(std::vector<double> v, double q01) {
    if (v.empty()) return 0.0;
    q01 = std::clamp(q01, 0.0, 1.0);
    const size_t k = (size_t)std::floor(q01 * (v.size()-1));
    std::nth_element(v.begin(), v.begin()+k, v.end());
    return v[k];
}

static double safe_div(double num, double den) {
    return (den > 0.0) ? (num / den) : 0.0;
}

static inline void log_vec_stats(const char* tag, const std::vector<double>& v) {
    auto [mn,mean,mx] = min_mean_max(v);
    std::cout << tag << " min=" << mn << " mean=" << mean << " max=" << mx << " size=" << v.size() << "\n";
}

template <typename T>
static inline size_t count_nz(const std::vector<T>& a) {
    size_t c=0; for (const auto& x : a) if (x!=T{}) ++c; return c;
}

// --- Fast 64-bit hash for sign masks (no external deps) ---------------------
// SplitMix64 mix (public-domain quality mix; good enough for probes)
static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x = x ^ (x >> 31);
    return x;
}

// Generic reader: treats any mask/container with size() and operator[] as bytes (0/1)
template <class Mask>
uint64_t hash_sign_mask(const Mask& mask) {
    const size_t n = mask.size();
    const uint8_t* base = nullptr;

    // If mask stores bytes contiguously, use them directly; otherwise copy cheaply.
    std::vector<uint8_t> tmp; tmp.reserve(n);
    if constexpr (std::is_pointer_v<decltype(mask.data())>) {
        // Try to use data() if available and looks like bytes; otherwise fall back below.
        using byte_t = std::remove_pointer_t<decltype(mask.data())>;
        if constexpr (std::is_same_v<byte_t, uint8_t> || std::is_same_v<byte_t, unsigned char> || std::is_same_v<byte_t, char>) {
            base = reinterpret_cast<const uint8_t*>(mask.data());
        } else {
            for (size_t i = 0; i < n; ++i) tmp.push_back(static_cast<uint8_t>(!!mask[i]));
            base = tmp.data();
        }
    } else {
        for (size_t i = 0; i < n; ++i) tmp.push_back(static_cast<uint8_t>(!!mask[i]));
        base = tmp.data();
    }

    // Chunk fold 8 bytes → 64-bit and mix; then finalize like xxHash-style avalanche.
    uint64_t h = 0x9e3779b97f4a7c15ULL ^ static_cast<uint64_t>(n);
    size_t i = 0;
    while (i + 8 <= n) {
        uint64_t w =
            (uint64_t)base[i+0]        |
            (uint64_t)base[i+1] <<  8  |
            (uint64_t)base[i+2] << 16  |
            (uint64_t)base[i+3] << 24  |
            (uint64_t)base[i+4] << 32  |
            (uint64_t)base[i+5] << 40  |
            (uint64_t)base[i+6] << 48  |
            (uint64_t)base[i+7] << 56;
        h ^= splitmix64(w + i);
        // cheap rotate + multiply drift
        h = (h << 27) | (h >> (64 - 27));
        h = h * 0x165667919E3779F9ULL + 0x9E3779B97F4A7C15ULL;
        i += 8;
    }
    // Tail
    uint64_t tail = 0; int sh = 0;
    while (i < n) { tail |= (uint64_t)base[i++] << sh; sh += 8; }
    h ^= splitmix64(tail + n);

    // Final avalanche (Murmur3 finalizer)
    h ^= (h >> 33);
    h *= 0xff51afd7ed558ccdULL;
    h ^= (h >> 33);
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= (h >> 33);
    return h;
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
//    rebuild_pos_maps();
}

void TriangleCyclePipeline::rebuild_pos_maps(const SignedGraphForMIP& G_) {
    const int m = G_.edge_count();

    // Reset view
    Gt_.clear();
    Gt_.full2pos.assign(m, -1);
    Gt_.pos_mask.assign(m, 0);
    Gt_.neg_mask.assign(m, 0);

    // Build masks & maps from the CURRENT switched signature.
    Gt_.pos2full.reserve(m);
    for (int eid = 0; eid < m; ++eid) {
        if (G_.is_pos_edge(eid)) {
            Gt_.full2pos[eid] = (int)Gt_.pos2full.size();
            Gt_.pos2full.push_back(eid);
            Gt_.pos_mask[eid] = 1;
        } else if (G_.is_neg_edge(eid)) {
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
TriangleCyclePipeline::run_round(const SignedGraphForMIP& G_)
{
    // Keep persistent arrays in sync with |E| (defensive)
    S_.init_sizes_if_needed(G_);

    // Reset per-round scratch.
    round_.clear();
	tri_selected_.clear();
	covered_neg_eids_.clear();

    // Read modes from config
    const bool triangles_enabled = C_.modes.use_triangles;
    const bool sp_enabled        = C_.modes.use_negcycles;
    std::cout << "[PIPELINE] modes: triangles=" << (triangles_enabled?"on":"off")
              << " negcycles=" << (sp_enabled?"on":"off") << "\n";

    const int m = G_.edge_count();
    round_.selected_count_full.assign(m, 0.0);
    g_sp_cycles_accepted = 0;

    // Clear per-round exported keys vector
    accepted_keys_.clear();
	
    // === 1) Choose switching s and apply it ===
    // switching has been applied in the frustration model

    // === 2) Rebuild positive-subgraph maps under the new switching ===
    rebuild_pos_maps(G_);

    // === Round-setup & switching probes ===
    {
        const int pos_edges = Gt_.m_pos;
        const int neg_edges = (int)count_nz(Gt_.neg_mask);
        const int m_guard   = m;
        const int delta     = (pos_edges + neg_edges) - m_guard;
	    const uint64_t sign_hash = hash_sign_mask(Gt_.neg_mask);
        std::cout << "[SWITCH-PROBE] |E+|=" << pos_edges
                  << " |E-|=" << neg_edges
                  << " m="   << m_guard
                  << " guard_delta=" << delta
                  << " sign_hash=0x"  << std::hex << sign_hash << std::dec << "\n";

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
    std::cout << "[PHASE-PROBE] lambda_LP=" << C_.ranking.lambda_LP
              << " Gt.m_pos=" << Gt_.m_pos << " m=" << m
              << "\n";

    // === 3) Build edge salience (FULL-edge array) ===
    round_.sal_full = G_.edge_saliences();

    // Compute salience stats (with median & %zero)
    double sal_min, sal_mean, sal_max;
    {
        auto [mn,mean,mx] = min_mean_max(round_.sal_full);
        sal_min = mn; sal_mean = mean; sal_max = mx;
        size_t nz = 0; for (double v : round_.sal_full) if (v>0.0) ++nz;
        size_t no = 0; for (double v : round_.sal_full) if (v>0.0 && v<1.0) ++no;
        const double p0 = round_.sal_full.empty()?0.0 : 100.0 * (double)(round_.sal_full.size()-nz) / (double)round_.sal_full.size();
        const double f0 = round_.sal_full.empty()?0.0 : 100.0 * (double)(no) / (double)round_.sal_full.size();
        std::vector<double> tmp = round_.sal_full;
        double med = median_of(tmp);
        std::cout << "[PHASE-PROBE] salience stats min=" << sal_min
                  << " mean=" << sal_mean
                  << " median=" << med
                  << " max=" << sal_max
                  << " zeros%=" << p0 << "%"
                  << " fracs%=" << f0 << "%\n";
    }

    // === 4) Build working weights ω′ on E⁺ (includes λ_LP blend) ===
    build_omega_prime_pos_();

    // === Run stages ===
    if (triangles_enabled && C_.budget.B_tri > 0) {
        triangle_stage_(G_);
    }
    if (sp_enabled) {
        sp_stage_(G_);
    }

    // === Commit persistent updates (ω/H) ===
    if ((triangles_enabled && C_.budget.B_tri > 0) || sp_enabled) {
        { double sel_sum=0.0; int sel_nz=0; for (double v: round_.selected_count_full){ sel_sum+=v; if (v!=0.0) ++sel_nz; }
          double used_sum=0.0; int used_nz=0; for (double v: round_.used_in_batch_pos){ used_sum+=v; if (v!=0.0) ++used_nz; }
          std::cout << "[SP-PROBE] sel_nz="<<sel_nz<<" sel_sum="<<sel_sum
                    << " used_nz_pos="<<used_nz<<" used_sum_pos="<<used_sum << "\n"; }
    	commit_stage_(G_);
    }

    Result R{};
    R.triangles_accepted = (int)tri_selected_.size();
    R.cycles_accepted    = g_sp_cycles_accepted;
    R.accepted_keys      = std::move(accepted_keys_);
    return R;
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
	std::cout << "[OMEGA'-STATS] H[log1p]*λ_hist min=" << ht_min << " mean=" << ht_mean << " max=" << ht_max << "\n";
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
	std::cout << "[OMEGA'-DOMAIN] size=" << round_.omega_prime_pos.size()
	          << " m_pos=" << m_pos
	          << " size_ok=" << (size_ok?1:0)
	          << " nan=" << nan_cnt
	          << " inf=" << inf_cnt
	          << " below_eps=" << below_eps
	          << " above_cap=" << above_cap << "\n";
}

void TriangleCyclePipeline::triangle_stage_(const SignedGraphForMIP& G_) {
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
    {
        // Build E⁺ adjacency from helper (current switched signature)
        auto pos_edges = G_.get_positive_edges();
        for (const auto& e : pos_edges) {
            pos_adj[e.first].push_back(e.second);
            pos_adj[e.second].push_back(e.first);
        }
    }
    // Negative anchors (full edges) under current switching
	auto neg_eids = G_.get_negative_eids();
	std::vector<double> neg_sals;
	neg_sals.reserve(neg_eids.size());
	double mean = 0.0;
	for (int feid : neg_eids) { double s = round_.sal_full[feid]; neg_sals.push_back(s); mean += s; }
	const int neg_raw = (int)neg_eids.size();
	mean = (neg_raw ? mean / neg_raw : 0.0);
	
	// F: p70 gate + mean fallback
	const double p70 = percentile_of(neg_sals, 0.70);
	const double sal_min_fallback = std::max(p70, 0.5 * mean);
	
	std::vector<Edge> neg_edges;                 // FIX: no pre-size
	neg_edges.reserve(neg_eids.size() / 1.5);          // capacity only
	const auto& sv = G_.signs_view();
	int nsel = 0;
	for (int feid : neg_eids) {
	    if (round_.sal_full[feid] < sal_min_fallback) continue;  // skip anchors not fractional enough
	    const auto& p = sv[feid].points;
	    neg_edges.emplace_back((int)p.first, (int)p.second);
	    ++nsel;
	}
	
	std::cout << "[TBB-FILTER] negE_raw=" << neg_raw
	          << " negE_keep=" << nsel
	          << " sal_fallback=" << sal_min_fallback << " (p70=" << p70 << ", mean=" << mean << ")\n";

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
	
		const double p = std::max(1.0, C_.ranking.inv_power);
		auto w_at = [&](int pid) {
		    return std::max(C_.weights.omega_eps, round_.omega_prime_pos[pid]);
		};
		
		double inv_hmean = 0.0;
		if (pid_uw >= 0 && pid_wv >= 0) {
		    const double w1 = w_at(pid_uw), w2 = w_at(pid_wv);
		    // If you want exactly α / Hmean(w): set C_.ranking.inv_power = 1.0
		    if (p == 1.0) {
		        inv_hmean = 0.5 * (1.0 / w1 + 1.0 / w2);   // = 1 / Hmean(w1,w2)
		    } else {
		        inv_hmean = 0.5 * (std::pow(w1, -p) + std::pow(w2, -p)); // stronger tilt to small ω′
		    }
		}
		c.score_primary = C_.ranking.alpha * inv_hmean + C_.ranking.theta * c.phi;	    
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
		std::vector<int> anchors_with_candidates;            // for telemetry only
		const auto& sel = tbb.select(anchors_with_candidates);
		
		// Persist SELECTION and compute the set of SELECTED anchors for gating
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
		
		// covered_by_tri := anchors of triangles that were actually SELECTED this round
		covered_neg_eids_.clear();
		covered_neg_eids_.reserve(tri_selected_.size());
		{
		    std::unordered_set<int> seen;
		    for (const auto& r : tri_selected_) {
		        if (seen.insert(r.neg_eid).second)
		            covered_neg_eids_.push_back(r.neg_eid);
		    }
		}
		
		// Keep candidates count for probes
		const int candidates_count = (int)anchors_with_candidates.size();
		const int selected_anchor_count = (int)covered_neg_eids_.size();
		std::cout << "[TBB-GATE-PROBE] candidates=" << candidates_count
		          << " selected_anchors=" << selected_anchor_count
		          << " (delta=" << (candidates_count - selected_anchor_count) << ")\n";

	    // After selection, log final TBB stats and acceptance rate
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
		const int base_B = C_.budget.B_tri;
		const int prev_B = (S_.B_tri_cur > 0 ? S_.B_tri_cur : base_B);
		
		// Map to γ in (gamma_min, gamma_max) with higher v → smaller γ
		const double ratio = std::clamp(S_.bar_v_triangle / std::max(1e-12, A.v0), 0.0, 1.0);
		double gamma = A.gamma_max - (A.gamma_max - A.gamma_min) * ratio;
		gamma = std::min(gamma, 0.9995); // never ≥ 1
		
		// Extra shrink if the batch was underutilized
		double util = (P.B_tri > 0) ? (double)selected_count / (double)P.B_tri : 0.0;
		// If very few were selected relative to the budget, shrink harder.
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

void TriangleCyclePipeline::sp_stage_(const SignedGraphForMIP& G_) {
    g_sp_cycles_accepted = 0;

    // Build uncovered negatives = graph negatives \ covered_neg_eids_
    const auto edge_idx = G_.edge_index();
    std::unordered_set<int> covered_set(covered_neg_eids_.begin(), covered_neg_eids_.end());
    
	auto neg_eids_all = G_.get_negative_eids();
	std::vector<double> neg_sals;
	neg_sals.reserve(neg_eids_all.size());
	double mean = 0.0;
	double max = 0.0, min = 1.0;
	for (int feid : neg_eids_all) {
		double s = round_.sal_full[feid];
		neg_sals.push_back(s);
		mean += s;
		if (s > max) max = s;
		if (s < min) min = s;
	}
	const int neg_raw = (int)neg_sals.size();
	mean = (neg_raw ? mean / neg_raw : 0.0);
	
	const double p70 = percentile_of(neg_sals, 0.7);
	const double sal_min_fallback = std::max(p70, 0.5 * mean);
	
	// Build uncovered list:
	std::vector<Edge> neg_uncov; neg_uncov.reserve(neg_eids_all.size());
	const auto& sv = G_.signs_view();
	int kept = 0;
	for (int feid : neg_eids_all) {
	    if (covered_set.count(feid)) continue;
	    if (round_.sal_full[feid] < sal_min_fallback - 1e-07) continue;
	    const auto& p = sv[feid].points;
	    if (fabs(G_.vertex_salience((int)p.first)+G_.vertex_salience((int)p.second)) < 1e-07) continue;
	    neg_uncov.emplace_back((int)p.first, (int)p.second);
	    ++kept;
	}
	
	const int covered = (int)covered_set.size();
	const double denom = std::max(1.0, (double)(neg_raw - covered));
	std::cout << "[SP-GATE] neg_raw=" << neg_raw
	          << " covered_by_tri=" << covered
	          << " kept_for_sp=" << kept
	          << " sal_min=" << sal_min_fallback
	          << " keep_rate=" << std::fixed << std::setprecision(2) << (100.0 * (kept / denom))
	          << "% (min=" << min << ", max=" << max << ", p70=" << p70 << ", mean=" << mean << ")\n";
    if (kept == 0) { std::cout << "[SP-GATE] kept_for_sp=0 → skipping SP stage.\n"; return; }

    NegativeCycleBatch ncb(G_, neg_uncov, /*cover=*/false, C_.to_ncb_params());

    std::vector<NegativeCycle> batch;
    while (ncb.next(batch)) {
        std::unordered_map<int,int> lnodes_hist;
        int accepted_in_batch = 0;
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
            ++lnodes_hist[L_nodes];
            ++accepted_in_batch;
            const double dens = 1.0 / std::max(1, L_nodes);

            // Bump along positive path edges (FULL-ids) and mirror to POS-domain used mask
            for (const auto& pe : cyc.pos_edges()) {
                auto feid = edge_idx[pe];
                round_.selected_count_full[(size_t)feid] += dens;
            }

            // Also give density credit to the negative anchor (u,v) reconstructed from path endpoints
            auto neg_eid = edge_idx[cyc.neg_edge()];
            round_.selected_count_full[(size_t)neg_eid] += dens;
        }
		// after finishing processing 'batch' (and before batch.clear())
		ncb.flush_emitted_to([&](int full_eid, double d){
		    const int pid = Gt_.full2pos[full_eid];
		    if (pid >= 0 && pid < (int)round_.used_in_batch_pos.size()) {
		        round_.used_in_batch_pos[(size_t)pid] += d;
		    }
		});
        {
            std::vector<std::pair<int,int>> hist_pairs;
            hist_pairs.reserve(lnodes_hist.size());
            for (auto &kv : lnodes_hist) hist_pairs.emplace_back(kv.first, kv.second);
            std::sort(hist_pairs.begin(), hist_pairs.end(), [](const auto& a, const auto& b){ return a.first < b.first; });
            std::ostringstream oss;
            oss << "[SP-BATCH] accepted=" << accepted_in_batch << " emitted=" << emitted_batch << " L_nodes_hist=";
            for (const auto& kv : hist_pairs) oss << " " << kv.first << ":" << kv.second;
            std::cout << oss.str() << "\n";
        }
        batch.clear();
    }

    round_.emitted_cycles = g_sp_cycles_accepted;
}

void TriangleCyclePipeline::commit_stage_(const SignedGraphForMIP& G_) {
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

    const double delta = S_.ema_delta;
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

    // Compute (emitted − accepted) on POS edges via pos mapping
    size_t H_rej_nz = 0;
    double H_rej_sum = 0.0;
	const double rho = C_.ema.rho; // new knob, try 0.05–0.1
    if (Gt_.m_pos > 0 && !round_.used_in_batch_pos.empty()) {
        // Map accepted FULL-edge densities to POS indices
        std::vector<double> acc_pos(Gt_.m_pos, 0.0);
        for (int pid = 0; pid < Gt_.m_pos; ++pid) {
            int eid = Gt_.pos2full[pid];
            if (eid >= 0 && eid < m) {
				acc_pos[pid] = round_.selected_count_full[eid];
				Hnew[eid] += rho * round_.selected_count_full[eid];
			}
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

    // ── PROBE: summarize update magnitudes
    std::cout << "[COMMIT] omega_inc_nz=" << omega_inc_nz
              << " omega_inc_sum=" << omega_inc_sum
              << " clamp_min=" << omega_clamp_min
              << " clamp_max=" << omega_clamp_max
              << " H_rej_nz=" << H_rej_nz
              << " H_rej_sum=" << H_rej_sum
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
}
