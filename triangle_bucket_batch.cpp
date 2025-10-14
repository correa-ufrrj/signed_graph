// triangle_bucket_batch.cpp
#include "triangle_bucket_batch.h"
#include <unordered_map>
#include <algorithm>

extern "C" {
void TBB_on_emit(int edge_id, double used_density);
void TBB_on_accept(int edge_id, double density);
int  TBB_budget_override(int base);
}

// Returns selected candidates; also outputs the set of covered anchors (neg edge ids).
const std::vector<TriangleBucketBatch::Candidate>&
TriangleBucketBatch::select(std::vector<TriangleBucketBatch::EdgeId>& covered_neg_out) {
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
    for (EdgeId key : keys) if (!buckets_[key].empty()) covered_neg_out.push_back(key);

    // Annealed budget and within-batch density (positive edges only)
    int budget = TBB_budget_override(P_.B_tri);
    std::unordered_map<EdgeId, double> used_in_this_batch;
    int budget_used = 0;  // <-- add this

    // Pass 1: one per bucket
    for (EdgeId key : keys) {
        auto& buck = buckets_[key];
        for (const auto& c : buck) {
            if (!respect_cap_(c, used_per_vertex)) continue;
            commit_(c, used_per_vertex);
            ++budget_used;  // <-- count a committed triangle
            TBB_on_accept(c.neg_eid,    1.0/3.0);
            TBB_on_accept(c.pos_eid_uw, 1.0/3.0);
            TBB_on_accept(c.pos_eid_wv, 1.0/3.0);
            used_in_this_batch[c.pos_eid_uw] += 1.0/3.0;
            used_in_this_batch[c.pos_eid_wv] += 1.0/3.0;
            TBB_on_emit(c.pos_eid_uw, used_in_this_batch[c.pos_eid_uw]);
            TBB_on_emit(c.pos_eid_wv, used_in_this_batch[c.pos_eid_wv]);
            break; // at most one in Pass 1
        }
        if ((int)selected_.size() >= budget) {
            stats_.B_eff    = budget_used;           // <-- use budget_used (or selected_.size())
            stats_.selected = (int)selected_.size();
            return selected_;
        }
    }

    // Pass 2: fill remaining up to B_tri
    bool progressed = true;
    while (progressed && (int)selected_.size() < budget) {
        progressed = false;
        for (EdgeId key : keys) {
            auto& buck = buckets_[key];
            for (const auto& c : buck) {
                if (already_taken_(c)) continue;
                if (!respect_cap_(c, used_per_vertex)) continue;
                commit_(c, used_per_vertex);
                ++budget_used;  // <-- count a committed triangle
                TBB_on_accept(c.neg_eid,    1.0/3.0);
                TBB_on_accept(c.pos_eid_uw, 1.0/3.0);
                TBB_on_accept(c.pos_eid_wv, 1.0/3.0);
                used_in_this_batch[c.pos_eid_uw] += 1.0/3.0;
                used_in_this_batch[c.pos_eid_wv] += 1.0/3.0;
                TBB_on_emit(c.pos_eid_uw, used_in_this_batch[c.pos_eid_uw]);
                TBB_on_emit(c.pos_eid_wv, used_in_this_batch[c.pos_eid_wv]);
                progressed = true;
                break;
            }
            if ((int)selected_.size() >= budget) break;
        }
    }

    stats_.B_eff    = budget_used;                 // <-- replace ‘effective_budget_used’
    stats_.selected = (int)selected_.size();
    return selected_;
}

// Access selected triangles
const std::vector<TriangleBucketBatch::Candidate>& TriangleBucketBatch::selected() const { return selected_; }

// Access buckets (read-only)
const std::unordered_map<TriangleBucketBatch::EdgeId, std::vector<TriangleBucketBatch::Candidate>>&
TriangleBucketBatch::buckets() const { return buckets_; }

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
    if (P_.cap_per_vertex <= 0) return true;
    if (c.u < 0 || c.v < 0 || c.w < 0) return false;
    if (c.u >= (int)used.size() || c.v >= (int)used.size() || c.w >= (int)used.size()) return false;
    if (used[c.u] >= P_.cap_per_vertex) return false;
    if (used[c.v] >= P_.cap_per_vertex) return false;
    if (used[c.w] >= P_.cap_per_vertex) return false;
    return true;
}

bool TriangleBucketBatch::already_taken_(const Candidate& c) const {
    // We store exactly the oriented triple (u,w,v); buckets are keyed by (u,v) anchor anyway
    std::tuple<int,int,int> key{c.u, c.w, c.v};
    return taken_.find(key) != taken_.end();
}

void TriangleBucketBatch::commit_(const Candidate& c, std::vector<int>& used) {
    selected_.push_back(c);
    used[c.u] += 1;
    used[c.v] += 1;
    used[c.w] += 1;
    taken_.insert(std::tuple<int,int,int>{c.u, c.w, c.v});
}


