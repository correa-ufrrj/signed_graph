// ============================================================================
// File: include/reheat_pool.h
// Why: decouple storage/policies for reheated cycles from FrustrationModelXY
// ============================================================================
#pragma once
#include <unordered_map>
#include <vector>
#include "cycle_key.h"

struct ReheatItem {
    std::vector<int> cyc_vertices; // cycle as vertex ring; owner canonicalizes to CycleKey
    int ttl = 3;                   // rounds to live; decremented after selection
};

class ReheatPool {
public:
    using Key    = fmkey::CycleKey;
    using Value  = ReheatItem;
    using MapT   = std::unordered_map<Key, Value, fmkey::CycleKeyHash, fmkey::CycleKeyEq>;
    using iterator = MapT::iterator;
    using const_iterator = MapT::const_iterator;

    ReheatPool() = default;

    // --- STL-like surface to minimize call-site changes ---
    bool        empty() const noexcept;
    std::size_t size()  const noexcept;
    void        clear()       noexcept;

    iterator       begin()       noexcept;
    iterator       end()         noexcept;
    const_iterator begin() const noexcept;
    const_iterator end()   const noexcept;
    const_iterator cbegin()const noexcept;
    const_iterator cend()  const noexcept;

    iterator erase(iterator it) noexcept;

    // Insert or assign whole item by key
    void upsert(const Key& key, const Value& val);
    // Convenience: mutating access if you still prefer operator[]
    Value& operator[](const Key& key);

    // Lookups
    bool          contains(const Key& key) const;
    const Value*  find(const Key& key) const;
    Value*        find_mutable(const Key& key);

private:
    MapT pool_;
};



