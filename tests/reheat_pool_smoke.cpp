// tests/reheat_pool_smoke.cpp
// Build: c++ -std=c++17 -O2 tests/reheat_pool_smoke.cpp -o reheat_smoke
// Run:   ./reheat_smoke

#include <cassert>
#include <iostream>
#include <vector>

#include "cycle_key.h"
#include "reheat_pool.h"

static void test_basic_upsert_find() {
    ReheatPool pool;
    fmkey::CycleKey k1; k1.y_idx = {1,2,3}; k1.rhs = 0; k1.reversed = false;
    fmkey::CycleKey k2; k2.y_idx = {3,2,1}; k2.rhs = 0; k2.reversed = true;

    // upsert + size
    pool.upsert(k1, ReheatItem{{0,1,2}, 2});
    assert(!pool.empty() && pool.size() == 1);
    assert(pool.contains(k1));
    assert(!pool.contains(k2));

    // find immutable/mutable
    const ReheatItem* p1 = pool.find(k1);
    assert(p1 && p1->ttl == 2 && p1->cyc_vertices.size() == 3);
    ReheatItem* p1m = pool.find_mutable(k1);
    assert(p1m); p1m->ttl = 3;
    assert(pool.find(k1)->ttl == 3);

    // operator[] creates default
    pool[k2].ttl = 1;
    assert(pool.size() == 2);
}

static void test_ttl_decay_and_erase() {
    ReheatPool pool;
    fmkey::CycleKey k; k.y_idx = {7,3,5}; k.rhs = 0; k.reversed = false;
    pool.upsert(k, ReheatItem{{0,1,2}, 2});

    // two-step decay with erase()
    for (auto it = pool.begin(); it != pool.end(); ) {
        if (it->second.ttl > 0) { --(it->second.ttl); ++it; } else { it = pool.erase(it); }
    }
    assert(!pool.empty());
    for (auto it = pool.begin(); it != pool.end(); ) {
        if (it->second.ttl > 0) { --(it->second.ttl); ++it; } else { it = pool.erase(it); }
    }
    assert(pool.empty());
}

int main() {
    test_basic_upsert_find();
    test_ttl_decay_and_erase();
    std::cout << "[OK] ReheatPool smoke tests passed.\n";
    return 0;
}
