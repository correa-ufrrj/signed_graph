/*separation_pipeline.cpp */
#include "separation_pipeline.h"
#include <unordered_map>
#include <limits>
#include <cmath>
#include "triangle_bucket_batch.h"

// TriangleBucketBatch → driver hook bridge (strong defs override TBB weak ones)
namespace { TriangleCyclePipeline* g_tbb_active = nullptr; }
extern "C" {
void TBB_on_emit(int edge_id, double used_density) {
    if (g_tbb_active) g_tbb_active->on_emit_(edge_id, used_density);
}
void TBB_on_accept(int edge_id, double density) {
    if (g_tbb_active) g_tbb_active->on_accept_(edge_id, density);
}
int TBB_budget_override(int base) {
    return g_tbb_active ? g_tbb_active->override_budget_(base) : base;
}
}

// ─────────────────────────────────────────────────────────────
// TriangleCyclePipeline — step 1: scaffolding only (no behavior change)
// ─────────────────────────────────────────────────────────────

TriangleCyclePipeline::TriangleCyclePipeline(SignedGraphForMIP& G,
                                             SeparationPersistent& persistent,
                                             SeparationConfig cfg)
    : G_(G), S_(persistent), C_(cfg)
{
    // Ensure persistent edge-aligned arrays are properly sized.
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
    // We read G_.switched_weights (already used in NegativeCycleBatch) to decide E⁺/E⁻.
    Gt_.pos2full.reserve(m);
    const auto& swv = G_.get_switched_weight();
    for (int eid = 0; eid < m; ++eid) {
        double sw = (eid < (int)swv.size() ? swv[eid] : 0.0);
        if (sw > 0.0) {
            Gt_.full2pos[eid] = (int)Gt_.pos2full.size();
            Gt_.pos2full.push_back(eid);
            Gt_.pos_mask[eid] = 1;
        } else if (sw < 0.0) {
            Gt_.neg_mask[eid] = 1;
        } else {
            // zero-weight (shouldn't normally happen); leave unmapped
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
    // Keep persistent arrays in sync with |E| in case the owning model changed graph.
    S_.init_sizes_if_needed(G_);

    // Reset per-round scratch.
    round_.clear();
    const int m = G_.edge_count();
    round_.selected_count_full.assign(m, 0.0);

    // === 1) Choose switching s and apply it ===
    // Always start from the original signature for determinism.
    G_.restore_switched_sign();

    std::shared_ptr<const std::vector<int>> s_ptr;
    if (phase == Phase::Fractional && x_hat && y_hat) {
        // Inject fractional guidance into the graph (fills internal salience too).
        (void)G_.weighting_from_fractional(*x_hat, *y_hat);
        // Try fractional greedy switching; fall back to integer greedy if unavailable.
        if (auto s_opt = G_.fractional_greedy_switching(); s_opt.has_value()) {
            s_ptr = s_opt.value();
        }
    }
    if (!s_ptr) {
        // Pure greedy switching on the (restored) original signature
        auto s = G_.greedy_switching();
        s_ptr = std::make_shared<const std::vector<int>>(std::move(s));
    }
    // Apply switching to rebuild the working signature σ_s inside G_
    G_.switching_from_partition(*s_ptr);

    // === 2) Rebuild positive-subgraph maps under the new switching ===
    rebuild_pos_maps();

    // === 3) Build salience (FULL-edge array). For step 2 this is passive: λ_LP=0 ===
    build_salience_(x_hat, y_hat);

    // === 4) Build working weights ω′ on E⁺ (λ_LP currently 0.0) ===
    build_omega_prime_pos_();

    // The remaining stages are added in steps 3–7
    pre_enumeration_reheat_();
    triangle_stage_();
    sp_stage_();
    commit_stage_();

    Result R; // empty result for step 2
    R.triangles_accepted = 0;
    R.cycles_accepted    = 0;
    return R;
}

// ─────────────────────────────────────────────────────────────
// Step-2+ stubs (minimal implementations so the shell compiles)
// ─────────────────────────────────────────────────────────────

void TriangleCyclePipeline::build_salience_(const std::vector<double>* /*x_hat*/,
                                            const std::vector<double>* y_hat)
{
    // Size to |E| first
    const int m = G_.edge_count();
    round_.sal_full.assign(m, 0.0);

    // Preferred source when available: the graph's prepared salience (normalized [0,1])
    // which is populated by weighting_from_fractional().
    const std::vector<double>& gsal = G_.edge_salience_view();
    if ((int)gsal.size() == m) {
        round_.sal_full = gsal; // copy
        return;
    }

    // Fallback: compute from ŷ if provided (1 - 2*|y-0.5| clamped to [0,1])
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

    const double eps   = std::max(1e-12, C_.omega_eps);
    const double w_cap = (C_.omega_max > 0.0 ? C_.omega_max : std::numeric_limits<double>::infinity());

    const auto& omega = S_.omega;        // |E|
    const auto& H     = S_.H;            // |E|
    const auto& sal   = round_.sal_full; // |E|

    for (int pid = 0; pid < m_pos; ++pid) {
        const int eid = Gt_.pos2full[pid];

        double base = (eid >= 0 && eid < (int)omega.size()) ? omega[eid] : 1.0;
        double rep  = (eid >= 0 && eid < (int)H.size())     ? C_.lambda_hist * std::log1p(std::max(0.0, H[eid])) : 0.0;
        double lp   = (eid >= 0 && eid < (int)sal.size())    ? C_.lambda_LP * sal[eid] : 0.0; // λLP==0 in step 2

        double w = base + rep - lp;
        if (w < eps) w = eps;
        if (w > w_cap) w = w_cap;
        round_.omega_prime_pos[pid] = w;
    }
}

void TriangleCyclePipeline::pre_enumeration_reheat_() {
    // No-op for step 1. Real logic arrives in step 5.
}

void TriangleCyclePipeline::triangle_stage_() {
    // === Step 3: route triangle-first through TriangleBucketBatch ===
    tri_selected_.clear();
    covered_neg_eids_.clear();

    // Build neg edge list (under current switching) and positive adjacency
    const int n = G_.vertex_count();
    const int m = G_.edge_count();

    // Build edge-index map (min(u,v),max(u,v)) -> full eid
    std::unordered_map<long long, int> key2eid;
    key2eid.reserve((size_t)m * 1.2);
    auto key64 = [](int a, int b) -> long long {
        if (a > b) std::swap(a,b);
        return ( (static_cast<long long>(a) << 32) ^ static_cast<unsigned long long>(b) );
    };
    for (const auto& kv : G_.edge_index()) {
        const Edge& e = kv.first; int eid = kv.second;
        key2eid[key64(e.first, e.second)] = eid;
    }

    // Positive adjacency of the switched graph
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

    // Params from config (ensure a positive B_tri)
    TriangleBucketBatch::Params P;
    P.K_tri_per_neg  = C_.K_tri_per_neg;
    P.cap_per_vertex = C_.tri_cap_per_vertex;
    // If unspecified, aim for at most one per bucket as a first pass
    P.B_tri = (C_.B_tri > 0 ? C_.B_tri : std::max(1, (int)neg_edges.size()));

    TriangleBucketBatch tbb(neg_edges, pos_adj, key2eid, P);

    // Scorer: primary = α(1/ω′(uw)+1/ω′(wv)) + θ·φ, secondary = viol/|C| (set 0 for step 3)
    auto scorer = [&](TriangleBucketBatch::Candidate& c){
        const int pid_uw = (c.pos_eid_uw >= 0 && c.pos_eid_uw < (int)Gt_.full2pos.size()) ? Gt_.full2pos[c.pos_eid_uw] : -1;
        const int pid_wv = (c.pos_eid_wv >= 0 && c.pos_eid_wv < (int)Gt_.full2pos.size()) ? Gt_.full2pos[c.pos_eid_wv] : -1;
        double invw = 0.0;
        if (pid_uw >= 0 && pid_uw < (int)round_.omega_prime_pos.size()) invw += 1.0 / std::max(C_.omega_eps, round_.omega_prime_pos[pid_uw]);
        if (pid_wv >= 0 && pid_wv < (int)round_.omega_prime_pos.size()) invw += 1.0 / std::max(C_.omega_eps, round_.omega_prime_pos[pid_wv]);
        auto sal = [&](int eid){ return (eid >= 0 && eid < (int)round_.sal_full.size()) ? round_.sal_full[eid] : 0.0; };
        c.phi = (sal(c.neg_eid) + sal(c.pos_eid_uw) + sal(c.pos_eid_wv)) / 3.0;
        c.score_primary   = C_.alpha * invw + C_.theta * c.phi;
        c.score_secondary = 0.0; // viol density reserved for step 5+
        c.viol = 0.0;
    };

    // Build buckets (pure 1-neg triangles on G^+)
    tbb.build_buckets(scorer);

    // Selection with hook wiring
    g_tbb_active = this;
    std::vector<int> covered; // full eids of negative anchors with nonempty buckets
    const auto& sel = tbb.select(covered);
    g_tbb_active = nullptr;

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
    round_.emitted_triangles = (int)tri_selected_.size();
}


void TriangleCyclePipeline::sp_stage_() {
    // No-op for step 1. Real logic arrives in step 4.
}

void TriangleCyclePipeline::commit_stage_() {
    // No-op for step 3. Persistent ω/H updates and pool bookkeeping arrive in steps 6–7.
}

// ─────────────────────────────────────────────────────────────
// TBB hook handlers → update round state
// ─────────────────────────────────────────────────────────────
void TriangleCyclePipeline::on_emit_(int full_eid, double used_density) {
    if (full_eid < 0 || full_eid >= (int)Gt_.full2pos.size()) return;
    const int pid = Gt_.full2pos[full_eid];
    if (pid < 0 || pid >= (int)round_.omega_prime_pos.size()) return;
    round_.used_in_batch_pos[pid] = used_density;
    // ω′ ← ω′ + β_emit · used_density (clamped)
    double w = round_.omega_prime_pos[pid] + C_.beta_emit * used_density;
    if (w < C_.omega_eps) w = C_.omega_eps;
    if (C_.omega_max > 0.0 && w > C_.omega_max) w = C_.omega_max;
    round_.omega_prime_pos[pid] = w;
}

void TriangleCyclePipeline::on_accept_(int full_eid, double density) {
    if (full_eid < 0 || full_eid >= (int)round_.selected_count_full.size()) return;
    // Per-round density counter for persistent drift (applied in step 6)
    round_.selected_count_full[full_eid] += density;
}

int TriangleCyclePipeline::override_budget_(int base) const {
    // For now, respect TBB's provided budget; future steps may anneal here.
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
