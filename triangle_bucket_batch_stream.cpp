#include "triangle_bucket_batch_stream.h"

// This translation unit is intentionally minimal.
// All logic is in the header to keep the initial integration tiny and
// avoid touching existing build glue. If you prefer, you can move the
// implementation here later with explicit template instantiations.
//
// Example usage (to be wired later):
//
//   TriangleBucketBatchStream::PosAdj pos_adj = ...;        // G^+_σ
//   TriangleBucketBatchStream::EdgeIndex edge_idx = ...;    // (u,v)->eid
//   std::vector<std::pair<int,int>> neg_edges = ...;        // anchors
//   TriangleBucketBatchStream::Params P{K,B,cap};
//   TriangleBucketBatchStream stream(neg_edges, pos_adj, edge_idx, P);
//
//   auto scorer = [&](TriangleBucketBatchStream::Candidate& c){
//       // fill c.score_primary, c.score_secondary, optionally c.viol/c.phi
//   };
//   stream.build_buckets(scorer);
//   std::vector<int> covered;
//   const auto& chosen = stream.select(covered);


// Hook to compute the primary & secondary score and optional cached metrics.
// The caller supplies a functor:
//
//    void scorer(Candidate& c)
//
// that fills c.score_primary, c.score_secondary (and optionally c.viol, c.phi).
template <class Scorer>
void build_TriangleBucketBatchStream::buckets(Scorer&& scorer) {
    buckets_.clear();
    // Build fast membership for positive adjacency (per u)
    // We'll reuse pos_adj_ directly; to check common neighbors, we
    // two-pointer merge after sorting neighbor lists (we sort once lazily).
    ensure_sorted_adjacency_();

    for (const auto& uv : neg_edges_) {
        VertexId u = uv.first;
        VertexId v = uv.second;
        // anchor must exist in edge index (negative under σ)
        EdgeId neg_eid = eid_(u, v);
        if (neg_eid < 0) continue;

        // common neighbors in G^+_σ: N^+(u) ∩ N^+(v)
        const auto& Nu = pos_adj_[u];
        const auto& Nv = pos_adj_[v];
        std::vector<VertexId> W;
        W.reserve(std::min(Nu.size(), Nv.size()));
        std::set_intersection(Nu.begin(), Nu.end(), Nv.begin(), Nv.end(), std::back_inserter(W));

        auto& buck = buckets_[neg_eid];
        for (VertexId w : W) {
            // Check 1-neg condition strictly: (u,v) negative; (u,w) and (w,v) positive by construction
            Candidate c;
            c.u = u; c.v = v; c.w = w;
            c.neg_eid   = neg_eid;
            c.pos_eid_uw= eid_(u, w);
            c.pos_eid_wv= eid_(w, v);
            if (c.pos_eid_uw < 0 || c.pos_eid_wv < 0) continue; // guard (should not happen)
            scorer(c);     // let the caller compute scores/viol/phi
            buck.push_back(std::move(c));
        }

        // Sort by lexicographic (primary desc, secondary desc)
        auto cmp = [](const Candidate& a, const Candidate& b){
            if (a.score_primary != b.score_primary) return a.score_primary > b.score_primary;
            return a.score_secondary > b.score_secondary;
        };
        std::sort(buck.begin(), buck.end(), cmp);
        if (P_.K_tri_per_neg > 0 && (int)buck.size() > P_.K_tri_per_neg) {
            buck.resize(P_.K_tri_per_neg);
        }
    }
}

// Two-pass selection:
//  - Pass 1: secure at most one strong triangle per bucket (respect per-vertex caps)
//  - Pass 2: fill remaining candidates, cycling buckets, up to B_tri (respect caps)
//
// Returns indices into internal 'selected_' vector; also outputs the set of
// covered anchors (neg edge ids) in 'covered_neg_out'.
const std::vector<Candidate>& TriangleBucketBatchStream::select(std::vector<EdgeId>& covered_neg_out) {
    selected_.clear();
    covered_neg_out.clear();

    if (buckets_.empty() || P_.B_tri <= 0) return selected_;

    std::vector<int> used_per_vertex(pos_adj_.size(), 0);

    // Round-robin over buckets in deterministic order
    std::vector<EdgeId> keys;
    keys.reserve(buckets_.size());
    for (auto& kv : buckets_) keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());

    // Pass 1: one per bucket
    for (EdgeId key : keys) {
        auto& buck = buckets_[key];
        for (const auto& c : buck) {
            if (!respect_cap_(c, used_per_vertex)) continue;
            commit_(c, used_per_vertex);
            covered_neg_out.push_back(key);
            break; // at most one in Pass 1
        }
        if ((int)selected_.size() >= P_.B_tri) return selected_;
    }

    // Pass 2: fill remaining up to B_tri
    bool progressed = true;
    while (progressed && (int)selected_.size() < P_.B_tri) {
        progressed = false;
        for (EdgeId key : keys) {
            auto& buck = buckets_[key];
            for (const auto& c : buck) {
                // Skip if already chosen the same triangle (u,w,v) earlier
                if (already_taken_(c)) continue;
                if (!respect_cap_(c, used_per_vertex)) continue;
                commit_(c, used_per_vertex);
                progressed = true;
                break; // move to next bucket
            }
            if ((int)selected_.size() >= P_.B_tri) break;
        }
    }
    return selected_;
}

// Access selected triangles
const std::vector<Candidate>& TriangleBucketBatchStream::selected() const { return selected_; }

// Access buckets (read-only)
const std::unordered_map<EdgeId, std::vector<Candidate>>& TriangleBucketBatchStream::buckets() const { return buckets_; }

const Params& TriangleBucketBatchStream::params() const { return P_; }
Params& params() { return P_; }

// Map (min(u,v),max(u,v)) to a 64-bit key for unordered_map
static long long TriangleBucketBatchStream::key_(VertexId a, VertexId b) {
    if (a > b) std::swap(a, b);
    return ( (static_cast<long long>(a) << 32) ^ static_cast<unsigned long long>(b) );
}

EdgeId TriangleBucketBatchStream::eid_(VertexId a, VertexId b) const {
    auto it = edge_index_.find(key_(a,b));
    return (it == edge_index_.end()) ? -1 : it->second;
}

void TriangleBucketBatchStream::ensure_sorted_adjacency_() {
    if (adj_sorted_) return;
    for (auto& nbrs : const_cast<PosAdj&>(pos_adj_)) {
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }
    adj_sorted_ = true;
}

bool TriangleBucketBatchStream::respect_cap_(const Candidate& c, const std::vector<int>& used) const {
    auto ok = [&](VertexId x){ return used[x] < P_.cap_per_vertex; };
    return ok(c.u) && ok(c.v) && ok(c.w);
}

bool TriangleBucketBatchStream::already_taken_(const Candidate& c) const {
    // cheap set check by (neg_eid, w) would suffice, but keep full triplet
    auto key = std::tuple<VertexId,VertexId,VertexId>(c.u,c.w,c.v);
    return taken_.find(key) != taken_.end();
}

void TriangleBucketBatchStream::commit_(const Candidate& c, std::vector<int>& used) {
    used[c.u]++; used[c.v]++; used[c.w]++;
    selected_.push_back(c);
    taken_.insert(std::tuple<VertexId,VertexId,VertexId>(c.u,c.w,c.v));
}
