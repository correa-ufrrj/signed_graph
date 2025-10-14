// ============================================================================
// File: include/cycle_key.h
// Purpose: switching-invariant, orientation-agnostic keys for 1-neg cycles
// Simplified to leverage existing Edge/NegativeCycle types from signed_graph.h
// ============================================================================
#pragma once
#include <vector>
#include <algorithm>
#include <cstdint>
#include "signed_graph.h"   // Edge, NegativeCycle

namespace fmkey {

// Canonicalize an unordered pair (u,v) to (min,max)
inline std::pair<int,int> mm(int u, int v) {
    if (u > v) std::swap(u, v);
    return {u, v};
}

// A 1-negative-cycle key = (neg edge endpoints) + sorted list of distinct
// positive edges (each canonicalized as (min,max)).
struct CycleKey {
    Edge neg = Edge(0,0);                                // neg edge endpoints (min,max)
    std::vector<Edge> pos;             // sorted, unique (min,max) positive edges
};

// Construct from a triangle (neg uv, pos uw & wv)
inline CycleKey make_from_triangle(int u, int v, int w) {
    CycleKey k; auto ab = Edge(u, v); k.neg = ab;
    k.pos = { Edge(u, w), Edge(v, w) };
    if (k.pos[1] < k.pos[0]) std::swap(k.pos[0], k.pos[1]);
    return k;
}

// Construct from a general 1-neg cycle (neg uv + positive path)
inline CycleKey make_from_cycle(const NegativeCycle& C) {
    CycleKey k; auto ab = C.neg_edge();
    k.neg = ab;
    k.pos.reserve(C.pos_edges().size());
    for (auto const& e : C.pos_edges()) k.pos.push_back(e);
    std::sort(k.pos.begin(), k.pos.end());
    k.pos.erase(std::unique(k.pos.begin(), k.pos.end()), k.pos.end());
    return k;
}

// 64-bit anchor-only key for (u,v) (useful for reheat seed maps)
inline uint64_t anchor64(int u, int v) {
    if (u > v) std::swap(u, v);
    return (static_cast<uint64_t>(static_cast<uint32_t>(u)) << 32) |
           static_cast<uint32_t>(v);
}

// Hash & equality so CycleKey can be used in unordered_{set,map}
struct CycleKeyHash {
    std::size_t operator()(CycleKey const& k) const noexcept {
        uint64_t h = (static_cast<uint64_t>(static_cast<uint32_t>(k.neg.first)) << 32)
                   ^ static_cast<uint32_t>(k.neg.second);
        for (auto const& p : k.pos) {
            uint64_t t = (static_cast<uint64_t>(static_cast<uint32_t>(p.first)) << 32)
                       ^ static_cast<uint32_t>(p.second) ^ 0x9e3779b97f4a7c15ULL;
            // mix a bit
            h ^= t; h = (h << 7) ^ (h >> 3);
        }
        return static_cast<std::size_t>(h);
    }
};

struct CycleKeyEq {
    bool operator()(CycleKey const& x, CycleKey const& y) const noexcept {
        return x.neg == y.neg && x.pos == y.pos;
    }
};

} // namespace fmkey



