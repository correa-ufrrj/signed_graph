#include "../cycle_key.h"
#include "../reheat_pool.h"

#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

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
    assert(p1m);
    p1m->ttl = 3;
    assert(pool.find(k1)->ttl == 3);

    // operator[] creates default
    pool[k2].ttl = 1;
    pool[k2].cyc_vertices = {5,6,7};
    assert(pool.size() == 2 && pool.contains(k2));
}

static void test_ttl_decay_and_erase() {
    ReheatPool pool;

    fmkey::CycleKey kA; kA.y_idx = {7,3,5}; kA.rhs = 0; kA.reversed = false;
    fmkey::CycleKey kB; kB.y_idx = {8,4,2}; kB.rhs = 0; kB.reversed = false;

    pool.upsert(kA, ReheatItem{{0,1,2}, 2}); // A: ttl=2
    pool.upsert(kB, ReheatItem{{3,4,5}, 1}); // B: ttl=1

    // Pass 1: A -> 1 (kept), B -> 0 (erased)
    auto [dec1, er1] = pool.prune_by_ttl();
    assert(dec1 == 2);
    assert(er1  == 1);
    assert(pool.size() == 1);

    // Pass 2: remaining -> 0 (erased)
    auto [dec2, er2] = pool.prune_by_ttl();
    assert(dec2 == 1);
    assert(er2  == 1);
    assert(pool.empty());
}

int main() {
    test_basic_upsert_find();
    test_ttl_decay_and_erase();
    std::cout << "[OK] ReheatPool::prune_by_ttl() tests passed.\n";
    return 0;
}
