// ============================================================================
// File: include/cycle_key.h
// Why: single source of truth for the switching-invariant, orientation-aware key
// ============================================================================
#pragma once
#include <vector>
#include <cstddef>
#include <functional>

namespace fmkey {

struct CycleKey {
    // Canonical y-index sequence around the cycle (undirected edges, canonical rotation + direction)
    std::vector<int> y_idx;
    // Parity/right-hand-side; 0 for 1-neg cycles in your pipeline
    int rhs = 0;
    // Whether canonical direction is reversed vs. vertex order used to build the key
    bool reversed = false;
};

// Stable hashing across processes; keep simple + fast.
struct CycleKeyHash {
    std::size_t operator()(const CycleKey& k) const noexcept {
        std::size_t h = static_cast<std::size_t>(k.rhs) * 0x9e3779b97f4a7c15ULL ^ (k.reversed ? 0x85ebca6b : 0xc2b2ae35);
        for (int v : k.y_idx) {
            // FNV-1a-ish mix
            h ^= static_cast<std::size_t>(v) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
        }
        return h;
    }
};

struct CycleKeyEq {
    bool operator()(const CycleKey& a, const CycleKey& b) const noexcept {
        return a.rhs == b.rhs && a.reversed == b.reversed && a.y_idx == b.y_idx;
    }
};

} // namespace fmkey



