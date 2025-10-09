// ============================================================================
// File: src/reheat_pool.cpp
// ============================================================================
#include "reheat_pool.h"
#include <algorithm>

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

std::pair<std::size_t, std::size_t> ReheatPool::prune_by_ttl() {
    std::size_t decremented = 0;
    std::size_t erased = 0;
    for (auto it = pool_.begin(); it != pool_.end(); ) {
        auto& item = it->second;
        if (item.ttl > 0) {
            --item.ttl;
            ++decremented;
            if (item.ttl == 0) {
                it = pool_.erase(it);
                ++erased;
            } else {
                ++it;
            }
        } else {
            it = pool_.erase(it);
            ++erased;
        }
    }
    return {decremented, erased};
}

void ReheatPool::update_ema(const Key& key, double viol, double alpha) {
    // Clamp alpha for safety.
    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;
    auto it = pool_.find(key);
    if (it == pool_.end()) return;
    auto& item = it->second;
    item.ema_viol = alpha * viol + (1.0 - alpha) * item.ema_viol;
    item.last_viol = viol;
}

std::pair<std::size_t, std::size_t> ReheatPool::prune_by_ttl_and_ema(double ema_min) {
    std::size_t decremented = 0;
    std::size_t erased = 0;
    for (auto it = pool_.begin(); it != pool_.end(); ) {
        auto& item = it->second;
        bool erase = false;
        // TTL policy: decrement >0 and erase if hits 0; erase if already <=0.
        if (item.ttl > 0) {
            --item.ttl;
            ++decremented;
            if (item.ttl == 0) {
                erase = true;
            }
        } else {
            erase = true;
        }
        // EMA policy: prune if below threshold (chronically non-violated).
        if (!erase && item.ema_viol < ema_min) {
            erase = true;
        }
        if (erase) {
            it = pool_.erase(it);
            ++erased;
        } else {
            ++it;
        }
    }
    return {decremented, erased};
}
