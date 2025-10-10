#ifndef NEG_CYCLE_BATCH_TEST_STUB
#define NEG_CYCLE_BATCH_TEST_STUB
#endif
#include "../negative_cycle_batch.h"
#include <cassert>
#include <iostream>

static inline void add_edge_stub(SignedGraphForMIP& G, int u, int v, double w) {
    igraph_add_edge(&G.g, u, v);
    igraph_integer_t eid;
    igraph_get_eid(&G.g, &eid, u, v, /*directed=*/0, /*error=*/1);
    if ((size_t)eid >= G.switched_weights.size()) G.switched_weights.resize((size_t)eid+1, 0.0);
    G.switched_weights[(size_t)eid] = w;
}

int main() {
    // Build a 7-node undirected signed graph
    SignedGraphForMIP G;
    igraph_empty(&G.g, 7, /*directed=*/IGRAPH_UNDIRECTED);
    G.switched_weights.clear();

    // Anchor A: (0,1) negative with TWO triangle candidates via w=2 and w=3
    add_edge_stub(G, 0,1, -1.0);
    add_edge_stub(G, 0,2,  1.0); add_edge_stub(G, 1,2, 1.0);
    add_edge_stub(G, 0,3,  1.0); add_edge_stub(G, 1,3, 1.0);

    // Anchor B: (4,6) negative, NO triangle; but positive path 4-5-2-6
    add_edge_stub(G, 4,6, -1.0);
    add_edge_stub(G, 4,5,  1.0);
    add_edge_stub(G, 5,2,  1.0);
    add_edge_stub(G, 2,6,  1.0);

    // Some extra positives that don’t create a triangle for (4,6)
    add_edge_stub(G, 0,4, 1.0);

    NegativeCycleBatch stream(G, /*cover=*/true, /*use_triangle_order=*/true);
    std::vector<NegativeCycle> out;
    bool ok = stream.next(out);
    assert(ok && "stream.next() should emit at least one batch");

    int tri_cnt = 0, sp_cnt = 0;
    bool saw_neg_46 = false;
    for (const auto& C : out) {
        // triangle if exactly two positive edges on the path
        if ((C.neg.first == 0 && C.neg.second == 1) || (C.neg.first == 1 && C.neg.second == 0)) {
            if (C.path.size() == 2) ++tri_cnt;
        }
        if ((C.neg.first == 4 && C.neg.second == 6) || (C.neg.first == 6 && C.neg.second == 4)) {
            saw_neg_46 = true;
            if (C.path.size() >= 3) ++sp_cnt; // longer-than-triangle path
        }
    }

    // Expect: at least one triangle for (0,1)
    assert(tri_cnt >= 1);
    // Expect: SP ran for uncovered (4,6)
    assert(saw_neg_46);
    assert(sp_cnt >= 1);

    std::cout << "[OK] triangles=" << tri_cnt
              << " sp=" << sp_cnt
              << " total=" << out.size() << "\n";
    std::cout << "[OK] batches_emitted=" << stream.batches_emitted()
              << " total_cycles=" << stream.total_cycles_emitted() << "\n";

    igraph_destroy(&G.g);
    return 0;
}
