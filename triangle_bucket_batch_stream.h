#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include <limits>
#include <cstddef>

// Minimal, header-first triangle-first, bucketed batch stream.
// This class is intentionally self-contained so it can be introduced
// without refactoring existing sources. Wiring into the current
// separation flow can be done in a later step.
//
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
//
// Types:
//   VertexId: integer id of a vertex
//   EdgeId  : integer id of an undirected edge (consistent with model y-index)
//
class TriangleBucketBatchStream {
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

    // Constructor
    //
    // neg_edges: list of (u,v) that are negative under current switching
    // pos_adj  : adjacency of G^+_σ
    // edge_index: map (min(u,v),max(u,v)) -> y-index / internal edge id
    // params   : see Params
    explicit TriangleBucketBatchStream(const std::vector<std::pair<VertexId,VertexId>>& neg_edges,
                                       const PosAdj& pos_adj,
                                       const EdgeIndex& edge_index,
                                       Params params = {})
        : neg_edges_(neg_edges), pos_adj_(pos_adj), edge_index_(edge_index), P_(params) {}

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

    const Params& params() const;
    Params& params();

private:
    // Map (min(u,v),max(u,v)) to a 64-bit key for unordered_map
    static long long key_(VertexId a, VertexId b);
    EdgeId eid_(VertexId a, VertexId b) const;
    void ensure_sorted_adjacency_();
    bool respect_cap_(const Candidate& c, const std::vector<int>& used) const;
    bool already_taken_(const Candidate& c) const;
    void commit_(const Candidate& c, std::vector<int>& used);

private:
    const std::vector<std::pair<VertexId,VertexId>>& neg_edges_;
    const PosAdj& pos_adj_;
    const EdgeIndex& edge_index_;
    Params P_;

    bool adj_sorted_{false};
    std::unordered_map<EdgeId, std::vector<Candidate>> buckets_;
    std::vector<Candidate> selected_;
    std::unordered_set<std::tuple<VertexId,VertexId,VertexId>> taken_;
};
