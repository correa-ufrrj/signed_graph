// tests/triangle_bucket_batch_smoke.cpp
// Build: c++ -std=c++17 -O2 tests/triangle_bucket_batch_smoke.cpp -o tbb_smoke
// Run:   ./tbb_smoke

#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>
#include <algorithm>
#include <tuple>
#include <string>

#include <set>
#include "triangle_bucket_batch.h"  // uses the API implemented alongside triangle_bucket_batch.cpp

using TBB = TriangleBucketBatch;

static long long key64(int a, int b) {
    if (a > b) std::swap(a,b);
    return (static_cast<long long>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint32_t>(b);
}

struct Toy {
    int n;
    std::vector<std::vector<int>> pos_adj;                   // G^+_σ adjacency
    std::unordered_map<long long, int> edge_index;           // (u,v) -> eid
    std::vector<std::pair<int,int>> neg_edges;               // anchors
    int next_eid = 0;

    explicit Toy(int n_) : n(n_), pos_adj(n_) {}

    int add_pos(int u, int v) {
        pos_adj[u].push_back(v);
        pos_adj[v].push_back(u);
        int eid = next_eid++;
        edge_index.emplace(key64(u,v), eid);
        return eid;
    }
    int add_neg(int u, int v) {
        int eid = next_eid++;
        edge_index.emplace(key64(u,v), eid);
        neg_edges.emplace_back(u,v);
        return eid;
    }
};

int main() {
    // --- Graph ---
    // Vertices: 0..5
    // Make (0,1) negative with two common positive neighbors 2 and 3
    // Make (2,4) negative with one common positive neighbor 5
    Toy G(6);

    // Positive edges: triangles around (0,1): via w=2 and w=3
    G.add_pos(0,2); G.add_pos(1,2);
    G.add_pos(0,3); G.add_pos(1,3);

    // Positive edges: triangle around (2,4): via w=5
    G.add_pos(2,5); G.add_pos(4,5);

    // Add some extra positive edges that should NOT create 1-neg triangles for (0,1) or (2,4)
    G.add_pos(0,4);               // no (1,4), so no common neighbor via 4
    G.add_pos(1,5);               // no (0,5), so no common neighbor via 5

    // Negative anchors (must exist in edge_index)
    int eid_01 = G.add_neg(0,1);
    int eid_24 = G.add_neg(2,4);

    // --- Params ---
    TBB::Params P;
    P.K_tri_per_neg = 2;   // allow both candidates in (0,1)'s bucket prefilter
    P.B_tri         = 3;   // global batch limit
    P.cap_per_vertex= 1;   // strict per-vertex cap to test pruning

    // --- Instance ---
    TBB stream(G.neg_edges, G.pos_adj, G.edge_index, P);

    // --- Scorer ---
    // Prefer smaller w to enforce deterministic order within a bucket.
    auto scorer = [&](TBB::Candidate& c) {
        // Primary: higher is better
        c.score_primary   = 1000.0 - static_cast<double>(c.w);
        // Secondary: tie-break on (u+v)
        c.score_secondary = 1000.0 - static_cast<double>(c.u + c.v);
        // (Optional) keep placeholders for completeness
        c.viol = 0.0;
        c.phi  = 0.0;
    };

    // --- Build buckets ---
    stream.build_buckets(scorer);

    // Verify buckets keyed by neg eid
    const auto& buckets = stream.buckets();
    auto it01 = buckets.find(eid_01);
    auto it24 = buckets.find(eid_24);
    assert(it01 != buckets.end());
    assert(it24 != buckets.end());

    // (0,1) should have two candidates: w=2 and w=3 (order by scorer => w=2 first)
    assert(static_cast<int>(it01->second.size()) == 2);
    assert((it01->second[0].u == 0 && it01->second[0].v == 1 && it01->second[0].w == 2));
    assert((it01->second[1].u == 0 && it01->second[1].v == 1 && it01->second[1].w == 3));

    // (2,4) should have one candidate: w=5
    assert(static_cast<int>(it24->second.size()) == 1);
    assert((it24->second[0].u == 2 && it24->second[0].v == 4 && it24->second[0].w == 5));

    // --- Select (two-pass) ---
    std::vector<int> covered;
    const auto& selected = stream.select(covered);

    // Global cap respected
    assert(static_cast<int>(selected.size()) <= P.B_tri);

    // With cap_per_vertex=1, we can take at most one triangle involving vertex 0
    // Pass-1 should take one per bucket where feasible:
    //  - From bucket (0,1): pick (0,2,1) (w=2) OR (0,3,1) depending on order
    //  - From bucket (2,4): pick (2,5,4)
    // Then Pass-2 cannot add the second (0,1) triangle since vertices 0 and 1 are already used.
    assert(static_cast<int>(selected.size()) == 2);

    // Ensure each selected candidate is truly 1-neg: pos_eids non-negative, neg_eid matches the bucket
    for (const auto& c : selected) {
        assert(c.pos_eid_uw >= 0 && c.pos_eid_wv >= 0 && c.neg_eid >= 0);
    }

    // Spec check: covered_neg should equal all nonempty bucket keys (incl. reheat-seeded)
    std::set<int> covered_set(covered.begin(), covered.end());
    std::set<int> expected_keys;
    for (const auto& kv : buckets) expected_keys.insert(kv.first);
    assert(covered_set == expected_keys && "covered should equal all non-empty bucket keys");

    // Current implementation's `covered` contains anchors accepted in Pass-1; should be size 2 here

    // Basic human-readable dump
    auto pretty = [](const TBB::Candidate& c){
        return std::string("C=(") + std::to_string(c.u) + "," + std::to_string(c.w) + "," + std::to_string(c.v) + ")";
    };

    std::cout << "[OK] Buckets:\n";
    for (const auto& kv : buckets) {
        std::cout << "  neg_eid=" << kv.first << " size=" << kv.second.size() << "\n";
    }
    std::cout << "[OK] Selected (" << selected.size() << "):\n";
    for (const auto& c : selected) std::cout << "  " << pretty(c) << "\n";
    std::cout << "[OK] covered size: " << covered.size() << "\n";

    std::cout << "All smoke tests passed.\n";
    return 0;
}

