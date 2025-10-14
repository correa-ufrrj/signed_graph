/* triangle_bucket_batch.h */
#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include <limits>
#include <cstddef>
#include <tuple>
#include <functional>
#include <utility>

// Design:
//  - Iterate negative edges (anchors) under the *current switched signs*
//  - For each (u,v) in E^-_σ, scan common positive neighbors w with (u,w),(w,v) in E^+_σ
//  - Form 1-neg triangle candidates and bucket by anchor negative edge
//  - Two-pass choose: at most one strong candidate per bucket, then fill up to B_tri
//
// Notes:
//  - The stream *does not* build cuts; it returns compact triangle
//    descriptors that existing cut builders can consume.
//  - Scoring hooks are exposed so callers can plug their current
//    violation/density score without changing public API later.
//  - No dependency on project headers; integrate using your existing
//    mapping (vertex/edge indices, y-index per edge, etc.) in the caller.
// Custom hasher/equality for std::tuple<int,int,int> so unordered_set works on libstdc++17.
namespace tbb_detail {
struct TupleHash {
    std::size_t operator()(const std::tuple<int,int,int>& t) const noexcept {
        const int a = std::get<0>(t);
        const int b = std::get<1>(t);
        const int c = std::get<2>(t);
        std::size_t h = 0;
        auto mix = [&](std::size_t v) {
            h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        };
        h ^= std::hash<int>{}(a);
        mix(std::hash<int>{}(b));
        mix(std::hash<int>{}(c));
        return h;
    }
};
struct TupleEq {
    bool operator()(const std::tuple<int,int,int>& x,
                    const std::tuple<int,int,int>& y) const noexcept {
        return x == y;
    }
};
} // namespace tbb_detail
//
// Types:
//   VertexId: integer id of a vertex
//   EdgeId  : integer id of an undirected edge (consistent with model y-index)
//
// In triangle_bucket_batch.h (minimal additions)
class TriangleBucketBatch {
public:
    using VertexId = int;
    using EdgeId   = int;

    struct Candidate {
        // Triangle (u,w,v) with anchor negative edge (u,v)
        VertexId u{}, w{}, v{};
        EdgeId   neg_eid{-1};       // y-index / internal id of (u,v)
        EdgeId   pos_eid_uw{-1};    // y-index / internal id of (u,w)
        EdgeId   pos_eid_wv{-1};    // y-index / internal id of (w,v)
        double   score_primary{0.0}; // caller-provided (e.g., 1/ω'+θ*sal)
        double   score_secondary{0.0}; // caller-provided (e.g., viol/|C|)
        // Optional cached values the caller may want to keep:
        double   viol{0.0};
        double   phi{0.0};
    };

    // Simple adjacency for G^+_σ: pos_adj[u] = list of neighbors (v) with (u,v) \in E^+_σ
    // Edge index map: undirected (min(u,v), max(u,v)) -> EdgeId (e.g., y-index)
    using PosAdj = std::vector<std::vector<VertexId>>;
    using EdgeIndex = std::unordered_map<long long, EdgeId>;

    // Settings (caps/truncation)
    struct Params {
        int K_tri_per_neg = 8;      // prefilter top-K per bucket (not an acceptance quota)
        int B_tri         = 64;     // global per-batch budget
        int cap_per_vertex= 3;      // per-vertex cap across the batch (small integer)
    };

	// neg_edges: list of (u,v) that are negative under current switching
	// pos_adj  : adjacency of G^+_σ
	// edge_index: map (min(u,v),max(u,v)) -> y-index / internal edge id
	// params   : see Params
	explicit TriangleBucketBatch(const std::vector<std::pair<VertexId,VertexId>>& neg_edges,
	                                   const PosAdj& pos_adj,
	                                   const EdgeIndex& edge_index,
	                                   TriangleBucketBatch::Params params)
	    : neg_edges_(neg_edges), pos_adj_(pos_adj), edge_index_(edge_index), P_(params) {
			    stats_ = Stats{};
			    stats_.neg_edges = (int)neg_edges.size();
    			stats_.B_tri     = params.B_tri;
		}
	    
	
	// Convenience overload: uses default Params{} (avoids default-arg inside class)
	explicit TriangleBucketBatch(const std::vector<std::pair<VertexId,VertexId>>& neg_edges,
	                                   const PosAdj& pos_adj,
	                                   const EdgeIndex& edge_index)
	    : TriangleBucketBatch(neg_edges, pos_adj, edge_index, Params{}) {}

    // Hook to compute the primary & secondary score and optional cached metrics.
    // The caller supplies a functor:
    //
    //    void scorer(Candidate& c)
    //
    // that fills c.score_primary, c.score_secondary (and optionally c.viol, c.phi).
    template <class Scorer>
    void build_buckets(Scorer&& scorer);
    // Two-pass selection:
    //  - Pass 1: secure at most one strong triangle per bucket (respect per-vertex caps)
    //  - Pass 2: fill remaining candidates, cycling buckets, up to B_tri (respect caps)
    //
    // Returns indices into internal 'selected_' vector; also outputs the set of
    // covered anchors (neg edge ids) in 'covered_neg_out'.
    const std::vector<Candidate>& select(std::vector<EdgeId>& covered_neg_out);

    // Access selected triangles
    const std::vector<Candidate>& selected() const;

    // Access buckets (read-only)
    const std::unordered_map<EdgeId, std::vector<Candidate>>& buckets() const;

    const TriangleBucketBatch::Params& params() const;
    TriangleBucketBatch::Params& params();

    struct Stats {
        int         neg_edges       = 0;      // number of negative anchors provided
        int         buckets_nonempty= 0;      // how many anchors produced >=1 candidate
        long long   candidates      = 0;      // total candidates across all buckets
        int         B_tri           = 0;      // configured budget at build time
        int         B_eff           = 0;      // effective budget actually used in select()
        int         selected        = 0;      // triangles selected
    };

    const Stats& stats() const { return stats_; }

private:
    // Map (min(u,v),max(u,v)) to a 64-bit key for unordered_map
    static long long key_(VertexId a, VertexId b);
    EdgeId eid_(VertexId a, VertexId b) const;
    void ensure_sorted_adjacency_();
    bool respect_cap_(const Candidate& c, const std::vector<int>& used) const;
    bool already_taken_(const Candidate& c) const;
    void commit_(const Candidate& c, std::vector<int>& used);

    const std::vector<std::pair<VertexId,VertexId>>& neg_edges_;
    const PosAdj& pos_adj_;
    const EdgeIndex& edge_index_;
    Params P_;

    bool adj_sorted_{false};
    std::unordered_map<EdgeId, std::vector<Candidate>> buckets_;
    std::vector<Candidate> selected_;
    std::unordered_set<std::tuple<VertexId,VertexId,VertexId>, tbb_detail::TupleHash, tbb_detail::TupleEq> taken_;
    Stats stats_; // <- new lightweight stats holder
};
// --- Template implementation ---
template <class Scorer>
inline void TriangleBucketBatch::build_buckets(Scorer&& scorer) {
    // Do not clear: allow reheat-preseeded buckets (buck^(0)) to participate
    // buckets_.clear();
    ensure_sorted_adjacency_();

    for (const auto& uv : neg_edges_) {
        VertexId u = uv.first;
        VertexId v = uv.second;
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
            Candidate c;
            c.u = u; c.v = v; c.w = w;
            c.neg_eid    = neg_eid;
            c.pos_eid_uw = eid_(u, w);
            c.pos_eid_wv = eid_(w, v);
            if (c.pos_eid_uw < 0 || c.pos_eid_wv < 0) continue;
            scorer(c);
            buck.push_back(std::move(c));
        }

        // Sort by (primary desc, then secondary desc)
        std::sort(buck.begin(), buck.end(),
                  [](const Candidate& a, const Candidate& b){
                      if (a.score_primary != b.score_primary) return a.score_primary > b.score_primary;
                      return a.score_secondary > b.score_secondary;
                  });
        if (P_.K_tri_per_neg > 0 && (int)buck.size() > P_.K_tri_per_neg) {
            buck.resize(P_.K_tri_per_neg);
        }
    }
    // Now compute the light stats without exposing internals
    int nonempty = 0;
    long long total = 0;
    for (const auto& bucket : buckets_) { // assume internal: std::vector<std::vector<Candidate>> buckets_;
	    const auto& vec = bucket.second;
	    if (!vec.empty()) ++nonempty;
	    total += static_cast<long long>(vec.size());
	}
    stats_.buckets_nonempty = nonempty;
    stats_.candidates       = total;
}