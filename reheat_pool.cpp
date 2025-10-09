// ============================================================================
// File: src/reheat_pool.cpp
// ============================================================================
#include "reheat_pool.h"

bool ReheatPool::empty() const noexcept { return pool_.empty(); }
std::size_t ReheatPool::size() const noexcept { return pool_.size(); }
void ReheatPool::clear() noexcept { pool_.clear(); }

ReheatPool::iterator ReheatPool::begin() noexcept { return pool_.begin(); }
ReheatPool::iterator ReheatPool::end()   noexcept { return pool_.end(); }
ReheatPool::const_iterator ReheatPool::begin() const noexcept { return pool_.begin(); }
ReheatPool::const_iterator ReheatPool::end()   const noexcept { return pool_.end(); }
ReheatPool::const_iterator ReheatPool::cbegin() const noexcept { return pool_.cbegin(); }
ReheatPool::const_iterator ReheatPool::cend()   const noexcept { return pool_.cend(); }

ReheatPool::iterator ReheatPool::erase(iterator it) noexcept { return pool_.erase(it); }

void ReheatPool::upsert(const Key& key, const Value& val) {
    auto it = pool_.find(key);
    if (it == pool_.end()) pool_.emplace(key, val);
    else it->second = val;
}

ReheatPool::Value& ReheatPool::operator[](const Key& key) { return pool_[key]; }

bool ReheatPool::contains(const Key& key) const { return pool_.find(key) != pool_.end(); }

const ReheatPool::Value* ReheatPool::find(const Key& key) const {
    auto it = pool_.find(key);
    return it == pool_.end() ? nullptr : &it->second;
}

ReheatPool::Value* ReheatPool::find_mutable(const Key& key) {
    auto it = pool_.find(key);
    return it == pool_.end() ? nullptr : &it->second;
}

