// tests/triangle_bucket_batch_smoke.cpp
// Build: g++ -std=c++17 -O2 -Wall -Wextra -I.. -Iinclude ../triangle_bucket_batch.cpp tests/triangle_bucket_batch_smoke.cpp -o build/triangle_bucket_batch_smoke
// Run:   ./build/triangle_bucket_batch_smoke

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

#include "../triangle_bucket_batch.h"

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

static std::string pretty(const TBB::Candidate& c) {
    return "C=(" + std::to_string(c.u) + "," + std::to_string(c.w) + "," + std::to_string(c.v) + ")";
}

int main() {
    // Common scorer: prefer smaller w; tie-break on (u+v)
    auto scorer = [&](TBB::Candidate& c) {
        c.score_primary   = 1000.0 - static_cast<double>(c.w);
        c.score_secondary = 1000.0 - static_cast<double>(c.u + c.v);
        c.viol = 0.0; c.phi = 0.0;
    };

    // ----------------------------------------------------------------------
    // Scenario 1: baseline (no accidental extra common neighbor for (2,4))
    // ----------------------------------------------------------------------
    {
        std::cout << "=== Scenario 1: baseline ===\n";
        Toy G(6);

        // Positives around (0,1): common neighbors w=2,3
        G.add_pos(0,2); G.add_pos(1,2);
        G.add_pos(0,3); G.add_pos(1,3);

        // Positives around (2,4): common neighbor w=5
        G.add_pos(2,5); G.add_pos(4,5);

        // Extra edges that DO NOT create extra common neighbors for (0,1) or (2,4)
        G.add_pos(3,5);
        G.add_pos(1,5);

        // Negative anchors
        int eid_01 = G.add_neg(0,1);
        int eid_24 = G.add_neg(2,4);

        TBB::Params P;
        P.K_tri_per_neg = 2;
        P.B_tri         = 3;
        P.cap_per_vertex= 1;

        TBB stream(G.neg_edges, G.pos_adj, G.edge_index, P);
        stream.build_buckets(scorer);

        const auto& buckets = stream.buckets();
        auto it01 = buckets.find(eid_01);
        auto it24 = buckets.find(eid_24);
        assert(it01 != buckets.end());
        assert(it24 != buckets.end());

        // (0,1): two candidates (w=2, w=3) ordered by scorer (2 first)
        assert(static_cast<int>(it01->second.size()) == 2);
        assert((it01->second[0].u == 0 && it01->second[0].v == 1 && it01->second[0].w == 2));
        assert((it01->second[1].u == 0 && it01->second[1].v == 1 && it01->second[1].w == 3));

        // (2,4): exactly one candidate (w=5)
        assert(static_cast<int>(it24->second.size()) == 1);
        assert((it24->second[0].u == 2 && it24->second[0].v == 4 && it24->second[0].w == 5));

        std::cout << "Buckets (Scenario 1):\n";
        for (const auto& kv : buckets) {
            std::cout << "  neg_eid=" << kv.first << " size=" << kv.second.size() << " ->";
            for (const auto& c : kv.second) std::cout << " " << pretty(c);
            std::cout << "\n";
        }

        std::vector<int> covered;
        const auto& selected = stream.select(covered);

        std::cout << "Selected (Scenario 1): " << selected.size() << "\n";
        for (const auto& c : selected) std::cout << "  " << pretty(c) << "\n";

        for (const auto& c : selected) {
            assert(c.pos_eid_uw >= 0 && c.pos_eid_wv >= 0 && c.neg_eid >= 0);
        }

        // covered must equal all nonempty bucket keys
        std::set<int> covered_set(covered.begin(), covered.end());
        std::set<int> expected_keys;
        for (const auto& kv : buckets) expected_keys.insert(kv.first);
        assert(covered_set == expected_keys && "covered should equal all non-empty bucket keys (scenario 1)");

        // With cap_per_vertex=1 and a shared vertex (2) across buckets,
        // only one triangle overall is feasible.
        assert(static_cast<int>(selected.size()) == 1);
    }

    // ----------------------------------------------------------------------
    // Scenario 2: intentionally add an extra common neighbor for (2,4)
    // ----------------------------------------------------------------------
    {
        std::cout << "=== Scenario 2: extra common neighbor for (2,4) ===\n";
        Toy G2(6);

        // Positives around (0,1)
        G2.add_pos(0,2); G2.add_pos(1,2);
        G2.add_pos(0,3); G2.add_pos(1,3);

        // Positives around (2,4)
        G2.add_pos(2,5); G2.add_pos(4,5);

        // This creates a new common neighbor for (2,4): w = 0
        G2.add_pos(0,4);
        // Harmless extra
        G2.add_pos(1,5);

        int eid2_01 = G2.add_neg(0,1);
        int eid2_24 = G2.add_neg(2,4);

        TBB::Params P2;
        P2.K_tri_per_neg = 2;
        P2.B_tri         = 3;
        P2.cap_per_vertex= 1;

        TBB stream2(G2.neg_edges, G2.pos_adj, G2.edge_index, P2);
        stream2.build_buckets(scorer);

        const auto& buckets2 = stream2.buckets();
        auto it2_01 = buckets2.find(eid2_01);
        auto it2_24 = buckets2.find(eid2_24);
        assert(it2_01 != buckets2.end());
        assert(it2_24 != buckets2.end());

        // (2,4) now has two candidates: w=0 and w=5 (w=0 first by scorer)
        assert(static_cast<int>(it2_24->second.size()) == 2);
        assert((it2_24->second[0].u == 2 && it2_24->second[0].v == 4 && it2_24->second[0].w == 0));
        assert((it2_24->second[1].u == 2 && it2_24->second[1].v == 4 && it2_24->second[1].w == 5));

        std::cout << "Buckets (Scenario 2):\n";
        for (const auto& kv : buckets2) {
            std::cout << "  neg_eid=" << kv.first << " size=" << kv.second.size() << " ->";
            for (const auto& c : kv.second) std::cout << " " << pretty(c);
            std::cout << "\n";
        }

        std::vector<int> covered2;
        const auto& selected2 = stream2.select(covered2);

        std::cout << "Selected (Scenario 2): " << selected2.size() << "\n";
        for (const auto& c : selected2) std::cout << "  " << pretty(c) << "\n";

        // covered must equal all nonempty bucket keys
        std::set<int> covered_set2(covered2.begin(), covered2.end());
        std::set<int> expected_keys2;
        for (const auto& kv : buckets2) expected_keys2.insert(kv.first);
        assert(covered_set2 == expected_keys2 && "covered should equal all non-empty bucket keys (scenario 2)");

        // With cap_per_vertex=1 and vertex 2 involved in both buckets,
        // selection can only take one triangle overall (the other conflicts on vertex 2).
        assert(static_cast<int>(selected2.size()) == 1);
    }

    // ----------------------------------------------------------------------
    // Scenario 3: prefilter truncation (K_tri_per_neg = 1)
    // ----------------------------------------------------------------------
    {
        std::cout << "=== Scenario 3: prefilter truncation K_tri_per_neg=1 ===\n";
        Toy G3(6);

        // Same as Scenario 2 to ensure buckets start with >=2 candidates
        G3.add_pos(0,2); G3.add_pos(1,2);
        G3.add_pos(0,3); G3.add_pos(1,3);
        G3.add_pos(2,5); G3.add_pos(4,5);
        G3.add_pos(0,4); // makes an extra candidate for (2,4)
        G3.add_pos(1,5);

        int eid3_01 = G3.add_neg(0,1);
        int eid3_24 = G3.add_neg(2,4);

        TBB::Params P3;
        P3.K_tri_per_neg = 1;  // <-- prefilter should truncate per bucket to 1
        P3.B_tri         = 10; // large budget so truncation is visible in buckets
        P3.cap_per_vertex= 3;  // relaxed cap to avoid capping effects here

        TBB stream3(G3.neg_edges, G3.pos_adj, G3.edge_index, P3);
        stream3.build_buckets(scorer);
        const auto& buckets3 = stream3.buckets();

        auto it3_01 = buckets3.find(eid3_01);
        auto it3_24 = buckets3.find(eid3_24);
        assert(it3_01 != buckets3.end());
        assert(it3_24 != buckets3.end());

        // Truncated to top-1 by scorer:
        //  - (0,1): w=2 wins over w=3
        //  - (2,4): w=0 wins over w=5
        assert(static_cast<int>(it3_01->second.size()) == 1);
        assert(static_cast<int>(it3_24->second.size()) == 1);
        assert((it3_01->second[0].u == 0 && it3_01->second[0].v == 1 && it3_01->second[0].w == 2));
        assert((it3_24->second[0].u == 2 && it3_24->second[0].v == 4 && it3_24->second[0].w == 0));

        std::cout << "Buckets (Scenario 3):\n";
        for (const auto& kv : buckets3) {
            std::cout << "  neg_eid=" << kv.first << " size=" << kv.second.size() << " ->";
            for (const auto& c : kv.second) std::cout << " " << pretty(c);
            std::cout << "\n";
        }

        // Optional: selection should pick both (budget large, cap relaxed)
        std::vector<int> covered3;
        const auto& selected3 = stream3.select(covered3);
        std::cout << "Selected (Scenario 3): " << selected3.size() << "\n";
        for (const auto& c : selected3) std::cout << "  " << pretty(c) << "\n";

        // covered must equal all nonempty bucket keys
        std::set<int> covered_set3(covered3.begin(), covered3.end());
        std::set<int> expected_keys3;
        for (const auto& kv : buckets3) expected_keys3.insert(kv.first);
        assert(covered_set3 == expected_keys3 && "covered should equal all non-empty bucket keys (scenario 3)");

        // With truncation to 1 per bucket and relaxed caps, expect exactly 2 selected
        assert(static_cast<int>(selected3.size()) == 2);
    }

    std::cout << "All smoke tests passed.\n";
    return 0;
}
