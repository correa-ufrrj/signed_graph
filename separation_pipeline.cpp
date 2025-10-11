/*separation_pipeline.cpp */
#include "separation_pipeline.h"
#include <unordered_map>
#include <limits>
#include <cmath>
#include <algorithm>
#include "triangle_bucket_batch.h"
#include "negative_cycle_batch.h"   // for SP stage

// Forward-declare the C hooks with correct linkage so RAII scope can call them.
extern "C" {
void TBB_set_active(void* ctx,
                    void (*emit)(void*, int, double),
                    void (*accept)(void*, int, double),
                    int  (*budget)(void*, int));
void TBB_clear_active();
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

// SP stage: stash only the count here (no persistent updates in Round 4/5)
static thread_local int g_sp_cycles_accepted = 0;

} // namespace

extern "C" {
// Called by TriangleBucketBatch during enumeration/selection
void TBB_on_emit(int edge_id, double used_density) {
    if (g_tbb_vt.emit) g_tbb_vt.emit(g_tbb_vt.ctx, edge_id, used_density);
}
void TBB_on_accept(int edge_id, double density) {
    if (g_tbb_vt.accept) g_tbb_vt.accept(g_tbb_vt.ctx, edge_id, density);
}
int TBB_budget_override(int base) {
    return g_tbb_vt.budget ? g_tbb_vt.budget(g_tbb_vt.ctx, base) : base;
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

TriangleCyclePipeline::TriangleCyclePipeline(SignedGraphForMIP& G,
                                             SeparationPersistent& persistent,
                                             SeparationConfig cfg)
    : G_(G), S_(persistent), C_(cfg)
{
    // Ensure persistent edge-aligned arrays are sized
    S_.init_sizes_if_needed(G_);
    rebuild_pos_maps();
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
    const int m = G_.edge_count();
    round_.selected_count_full.assign(m, 0.0);
    g_sp_cycles_accepted = 0;

    // === 1) Choose switching s and apply it ===
    G_.restore_switched_sign();
    if (phase == Phase::Fractional && x_hat && y_hat) {
        // Load salience if available (λ_LP blending is handled in weights)
        (void)G_.weighting_from_fractional(*x_hat, *y_hat);
    }
    auto s = G_.greedy_switching();
    G_.switching_from_partition(s);

    // === 2) Rebuild positive-subgraph maps under the new switching ===
    rebuild_pos_maps();

    // === 3) Build salience (FULL-edge array) ===
    build_salience_(x_hat, y_hat);

    // === 4) Build working weights ω′ on E⁺ ===
    build_omega_prime_pos_();

    // === Run stages ===
    pre_enumeration_reheat_();
    triangle_stage_();
    sp_stage_();

    // === Commit persistent updates (Round 6) ===
    commit_stage_();

    Result R{};
    R.triangles_accepted = (int)tri_selected_.size();
    R.cycles_accepted    = g_sp_cycles_accepted;
    // Accepted_keys to be assembled in later steps when materialization is wired here.
    return R;
}

void TriangleCyclePipeline::build_salience_(const std::vector<double>* /*x_hat*/,
                                            const std::vector<double>* y_hat)
{
    const int m = G_.edge_count();
    round_.sal_full.assign(m, 0.0);

    // Preferred: graph-prepared salience (populated by weighting_from_fractional()).
    const std::vector<double>& gsal = G_.edge_salience_view();
    if ((int)gsal.size() == m) {
        round_.sal_full = gsal; // copy
        return;
    }

    // Fallback: derive from ŷ if provided.
    if (y_hat) {
        const auto& y = *y_hat;
        const int M = std::min(m, (int)y.size());
        for (int i = 0; i < M; ++i) {
            double d = std::fabs(y[i] - 0.5);
            double s = 1.0 - std::min(1.0, 2.0 * d);
            round_.sal_full[i] = (s < 0.0 ? 0.0 : (s > 1.0 ? 1.0 : s));
        }
    }
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

    for (int pid = 0; pid < m_pos; ++pid) {
        const int eid = Gt_.pos2full[pid];

        double base = (eid >= 0 && eid < (int)omega.size()) ? omega[eid] : 1.0;
        double rep  = (eid >= 0 && eid < (int)H.size())     ? C_.ranking.lambda_hist * std::log1p(std::max(0.0, H[eid])) : 0.0;
        double lp   = (eid >= 0 && eid < (int)sal.size())   ? C_.ranking.lambda_LP   * sal[eid] : 0.0;

        double w = base + rep - lp;
        if (w < eps) w = eps;
        if (w > w_cap) w = w_cap;
        round_.omega_prime_pos[pid] = w;
    }
}

void TriangleCyclePipeline::pre_enumeration_reheat_() {
    // Round 6: still no-op (reheat wiring comes later).
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
    if (P.B_tri <= 0) P.B_tri = std::max(1, (int)neg_edges.size());

    TriangleBucketBatch tbb(neg_edges, pos_adj, key2eid, P);

    // Scorer: primary = α(1/ω′(uw)+1/ω′(wv)) + θ·φ, secondary = 0 (viol reserved for later rounds)
    auto scorer = [&](TriangleBucketBatch::Candidate& c){
        const int pid_uw = (c.pos_eid_uw >= 0 && c.pos_eid_uw < (int)Gt_.full2pos.size()) ? Gt_.full2pos[c.pos_eid_uw] : -1;
        const int pid_wv = (c.pos_eid_wv >= 0 && c.pos_eid_wv < (int)Gt_.full2pos.size()) ? Gt_.full2pos[c.pos_eid_wv] : -1;
        double invw = 0.0;
        if (pid_uw >= 0 && pid_uw < (int)round_.omega_prime_pos.size()) invw += 1.0 / std::max(C_.weights.omega_eps, round_.omega_prime_pos[pid_uw]);
        if (pid_wv >= 0 && pid_wv < (int)round_.omega_prime_pos.size()) invw += 1.0 / std::max(C_.weights.omega_eps, round_.omega_prime_pos[pid_wv]);
        auto sal = [&](int eid){ return (eid >= 0 && eid < (int)round_.sal_full.size()) ? round_.sal_full[eid] : 0.0; };
        c.phi = (sal(c.neg_eid) + sal(c.pos_eid_uw) + sal(c.pos_eid_wv)) / 3.0;
        c.score_primary   = C_.ranking.alpha * invw + C_.ranking.theta * c.phi;
        c.score_secondary = 0.0;
        c.viol = 0.0;
    };

    tbb.build_buckets(scorer);

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

    {
        TBBHookScope scope(static_cast<void*>(this), emit_cb, accept_cb, budget_cb);
        std::vector<int> covered; // full eids of negative anchors with nonempty buckets
        const auto& sel = tbb.select(covered);

        // Persist outputs in round storage
        covered_neg_eids_.assign(covered.begin(), covered.end());
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
    }

    round_.emitted_triangles = (int)tri_selected_.size();
}

void TriangleCyclePipeline::sp_stage_() {
    // === Shortest-path cycle enumeration on uncovered anchors (NCB) ===
    g_sp_cycles_accepted = 0;

    NegativeCycleBatch ncb(G_, /*cover=*/false, /*use_triangle_order=*/false, C_.to_ncb_params());

    std::vector<NegativeCycle> batch;
    while (ncb.next(batch)) {
        g_sp_cycles_accepted += (int)batch.size();
        // (Round 6: persistent drift is handled in commit_stage_ via round_ arrays)
        batch.clear();
    }

    round_.emitted_cycles = g_sp_cycles_accepted;
}

void TriangleCyclePipeline::commit_stage_() {
    // ── Step 6: Persistent updates (ω, H via EMA, pool_count) ───────────
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
    for (int eid = 0; eid < m; ++eid) {
        double inc = round_.selected_count_full[eid];
        if (inc != 0.0) {
            double w = S_.omega[eid] + beta_sel * inc;
            if (w < eps) w = eps;
            if (w > w_cap) w = w_cap;
            S_.omega[eid] = w;
        } else {
            // keep ω within clamps
            if (S_.omega[eid] < eps) S_.omega[eid] = eps;
            if (S_.omega[eid] > w_cap) S_.omega[eid] = w_cap;
        }
    }

    // (B) H EMA: H ← δ·H + ρ·accepted + κ·(emitted−accepted)_pos
    // Start with decay
    std::vector<double> Hnew = S_.H;
    for (double& v : Hnew) v *= delta;

    // + ρ·accepted over ALL edges (neg & pos)
    for (int eid = 0; eid < m; ++eid) {
        double acc = round_.selected_count_full[eid];
        if (acc != 0.0) Hnew[eid] += rho * acc;
    }

    // Compute (emitted − accepted) on POS edges via pos mapping
    if (Gt_.m_pos > 0 && !round_.used_in_batch_pos.empty()) {
        // Map accepted FULL-edge densities to POS indices
        std::vector<double> acc_pos(Gt_.m_pos, 0.0);
        for (int pid = 0; pid < Gt_.m_pos; ++pid) {
            int eid = Gt_.pos2full[pid];
            if (eid >= 0 && eid < m) acc_pos[pid] = round_.selected_count_full[eid];
        }
        for (int pid = 0; pid < Gt_.m_pos; ++pid) {
            double emitted  = (pid < (int)round_.used_in_batch_pos.size()) ? round_.used_in_batch_pos[pid] : 0.0;
            double accepted = acc_pos[pid];
            double rejected = emitted - accepted;
            if (rejected < 0.0) rejected = 0.0;
            int eid = Gt_.pos2full[pid];
            if (eid >= 0 && eid < m) Hnew[eid] += kappa * rejected;
        }
    }
    S_.H.swap(Hnew);

    // (C) pool_count: accumulate accepted density on edges of accepted cuts
    for (int eid = 0; eid < m; ++eid) {
        double inc = round_.selected_count_full[eid];
        if (inc != 0.0) S_.pool_count[eid] += inc;
    }

    // NOTE: Stateless de-dup and reheat remain in SeparationPersistent and are
    // handled elsewhere. No additional reheat behavior is introduced here.
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
    // Pass-through for now; annealing may adjust this later.
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
