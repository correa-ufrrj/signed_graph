#include "triangle_bucket_batch.h"
#include <unordered_map>

// Weak, overridable integration hooks:
extern "C" {
#if defined(__GNUC__) || defined(__clang__)
__attribute__((weak))
#endif
void TBB_on_emit(int /*edge_id*/, double /*used_density*/) {}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((weak))
#endif
void TBB_on_accept(int /*edge_id*/, double /*density*/) {}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((weak))
#endif
int TBB_budget_override(int base) { return base; }
}


// Example usage (to be wired later):
//
//   TriangleBucketBatch::PosAdj pos_adj = ...;        // G^+_σ
//   TriangleBucketBatch::EdgeIndex edge_idx = ...;    // (u,v)->eid
//   std::vector<std::pair<int,int>> neg_edges = ...;        // anchors
//   TriangleBucketBatch::Params P{K,B,cap};
//   TriangleBucketBatch stream(neg_edges, pos_adj, edge_idx, P);
//
//   auto scorer = [&](TriangleBucketBatch::Candidate& c){
//       // fill c.score_primary, c.score_secondary, optionally c.viol/c.phi
//   };
//   stream.build_buckets(scorer);
//   std::vector<int> covered;
//   const auto& chosen = stream.select(covered);
//
//
// Hook to compute the primary & secondary score and optional cached metrics.
// The caller supplies a functor:
//
//    void scorer(Candidate& c)
//
// that fills c.score_primary, c.score_secondary (and optionally c.viol, c.phi).
//
// Two-pass selection:
//  - Pass 1: secure at most one strong triangle per bucket (respect per-vertex caps)
//  - Pass 2: fill remaining candidates, cycling buckets, up to B_tri (respect caps)
//
// Returns indices into internal 'selected_' vector; also outputs the set of
// covered anchors (neg edge ids) in 'covered_neg_out'.
const std::vector<TriangleBucketBatch::Candidate>& TriangleBucketBatch::select(std::vector<TriangleBucketBatch::EdgeId>& covered_neg_out) {
    selected_.clear();
    taken_.clear();
    covered_neg_out.clear();

    if (buckets_.empty() || P_.B_tri <= 0) return selected_;

    std::vector<int> used_per_vertex(pos_adj_.size(), 0);

    // Round-robin over buckets in deterministic order
    std::vector<EdgeId> keys;
    keys.reserve(buckets_.size());
    for (auto& kv : buckets_) keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());
    // Spec: covered_neg_out lists all nonempty bucket keys (incl. reheat-seeded)
    for (EdgeId key : keys) {
        if (!buckets_[key].empty()) covered_neg_out.push_back(key);
    }

    // Annealed budget and within-batch density (positive edges only)
    int budget = TBB_budget_override(P_.B_tri);
    std::unordered_map<EdgeId, double> used_in_this_batch;


    // Pass 1: one per bucket
    for (EdgeId key : keys) {
        auto& buck = buckets_[key];
        for (const auto& c : buck) {
            if (!respect_cap_(c, used_per_vertex)) continue;
            commit_(c, used_per_vertex);
            // Cross-batch density: add 1/3 to all three edges of the triangle
            TBB_on_accept(c.neg_eid,    1.0/3.0);
            TBB_on_accept(c.pos_eid_uw, 1.0/3.0);
            TBB_on_accept(c.pos_eid_wv, 1.0/3.0);
            // Within-batch mask bump on positive edges
            used_in_this_batch[c.pos_eid_uw] += 1.0/3.0;
            used_in_this_batch[c.pos_eid_wv] += 1.0/3.0;
            TBB_on_emit(c.pos_eid_uw, used_in_this_batch[c.pos_eid_uw]);
            TBB_on_emit(c.pos_eid_wv, used_in_this_batch[c.pos_eid_wv]);
            break; // at most one in Pass 1
        }
        if ((int)selected_.size() >= budget) return selected_;
    }

    // Pass 2: fill remaining up to B_tri
    bool progressed = true;
    while (progressed && (int)selected_.size() < budget) {
        progressed = false;
        for (EdgeId key : keys) {
            auto& buck = buckets_[key];
            for (const auto& c : buck) {
                // Skip if already chosen the same triangle (u,w,v) earlier
                if (already_taken_(c)) continue;
                if (!respect_cap_(c, used_per_vertex)) continue;
                commit_(c, used_per_vertex);
                // Cross-batch density: add 1/3 to all three edges of the triangle
                TBB_on_accept(c.neg_eid,    1.0/3.0);
                TBB_on_accept(c.pos_eid_uw, 1.0/3.0);
                TBB_on_accept(c.pos_eid_wv, 1.0/3.0);
                // Within-batch mask bump on positive edges
                used_in_this_batch[c.pos_eid_uw] += 1.0/3.0;
                used_in_this_batch[c.pos_eid_wv] += 1.0/3.0;
                TBB_on_emit(c.pos_eid_uw, used_in_this_batch[c.pos_eid_uw]);
                TBB_on_emit(c.pos_eid_wv, used_in_this_batch[c.pos_eid_wv]);
                progressed = true;
                break; // move to next bucket
            }
            if ((int)selected_.size() >= budget) break;
        }
    }
    return selected_;
}

// Access selected triangles
const std::vector<TriangleBucketBatch::Candidate>& TriangleBucketBatch::selected() const { return selected_; }

// Access buckets (read-only)
const std::unordered_map<TriangleBucketBatch::EdgeId, std::vector<TriangleBucketBatch::Candidate>>& TriangleBucketBatch::buckets() const { return buckets_; }

const TriangleBucketBatch::Params& TriangleBucketBatch::params() const { return P_; }
TriangleBucketBatch::Params& TriangleBucketBatch::params() { return P_; }

// Map (min(u,v),max(u,v)) to a 64-bit key for unordered_map
long long TriangleBucketBatch::key_(VertexId a, VertexId b) {
    if (a > b) std::swap(a, b);
    return ( (static_cast<long long>(a) << 32) ^ static_cast<unsigned long long>(b) );
}

int TriangleBucketBatch::eid_(int a, int b) const {
    auto it = edge_index_.find(key_(a,b));
    return (it == edge_index_.end()) ? -1 : it->second;
}

void TriangleBucketBatch::ensure_sorted_adjacency_() {
    if (adj_sorted_) return;
    for (auto& nbrs : const_cast<PosAdj&>(pos_adj_)) {
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }
    adj_sorted_ = true;
}

bool TriangleBucketBatch::respect_cap_(const Candidate& c, const std::vector<int>& used) const {
    auto ok = [&](VertexId x){ return used[x] < P_.cap_per_vertex; };
    return ok(c.u) && ok(c.v) && ok(c.w);
}

bool TriangleBucketBatch::already_taken_(const TriangleBucketBatch::Candidate& c) const {
    // cheap set check by (neg_eid, w) would suffice, but keep full triplet
    auto key = std::tuple<VertexId,VertexId,VertexId>(c.u,c.w,c.v);
    return taken_.find(key) != taken_.end();
}

void TriangleBucketBatch::commit_(const TriangleBucketBatch::Candidate& c, std::vector<int>& used) {
    used[c.u]++; used[c.v]++; used[c.w]++;
    selected_.push_back(c);
    taken_.insert(std::tuple<VertexId,VertexId,VertexId>(c.u,c.w,c.v));
}
