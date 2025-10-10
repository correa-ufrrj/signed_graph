#include <cassert>
#include <iostream>
#include <vector>
#include <utility>
#include <cstdint>

#include "../negative_cycle_batch.h"

static inline void add_edge_stub(SignedGraphForMIP& G, int u, int v, double w) {
    igraph_add_edge(&G.g, u, v);
    igraph_integer_t eid; igraph_get_eid(&G.g, &eid, u, v, /*directed=*/0, /*error=*/1);
    if ((size_t)eid >= G.switched_weights.size()) G.switched_weights.resize((size_t)eid+1, 0.0);
    G.switched_weights[(size_t)eid] = w;
}

int main() {
    // Tiny graph exercising both stages:
    //  - (0,1) negative with two triangle connectors via w=2 and w=3
    //  - (4,6) negative without triangle, but a positive path 4-5-2-6 (SP stage)
    SignedGraphForMIP G; igraph_empty(&G.g, 7, /*directed=*/IGRAPH_UNDIRECTED); G.switched_weights.clear();

    // Triangle bucket anchors (0,1)
    add_edge_stub(G, 0,1, -1.0);
    add_edge_stub(G, 0,2,  1.0); add_edge_stub(G, 1,2, 1.0);
    add_edge_stub(G, 0,3,  1.0); add_edge_stub(G, 1,3, 1.0);

    // SP-only anchor (4,6)
    add_edge_stub(G, 4,6, -1.0);
    add_edge_stub(G, 4,5,  1.0);
    add_edge_stub(G, 5,2,  1.0);
    add_edge_stub(G, 2,6,  1.0);

    // Some extra positives (not forming triangle for 4-6)
    add_edge_stub(G, 0,4, 1.0);

    // Run one batch with triangle-first enabled
    NegativeCycleBatch stream(G, /*cover=*/true, /*use_triangle_order=*/true);
    std::vector<NegativeCycle> out;
    bool ok = stream.next(out);
    assert(ok && "pipeline should emit a non-empty batch");

    int tri_cnt = 0, sp_cnt = 0; bool saw_01 = false, saw_46 = false;
    for (const auto& C : out) {
        const bool is01 = ( (C.neg.first==0 && C.neg.second==1) || (C.neg.first==1 && C.neg.second==0) );
        const bool is46 = ( (C.neg.first==4 && C.neg.second==6) || (C.neg.first==6 && C.neg.second==4) );
        if (is01) { saw_01 = true; if (C.path.size()==2) ++tri_cnt; }
        if (is46) { saw_46 = true; if (C.path.size()>=3) ++sp_cnt; }
    }

    // Expectations: at least one triangle on (0,1) and at least one SP cycle on (4,6)
    assert(saw_01 && tri_cnt >= 1);
    assert(saw_46 && sp_cnt >= 1);

    std::cout << "[PIPE] triangles=" << tri_cnt
              << " sp=" << sp_cnt
              << " total=" << out.size() << "\n";
    std::cout << "[PIPE] batches_emitted=" << stream.batches_emitted()
              << " total_cycles=" << stream.total_cycles_emitted() << "\n";

    igraph_destroy(&G.g);
    return 0;
}
