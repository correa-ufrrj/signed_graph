// File: tests/sg_rand_stats.cpp
#include <iostream>
#include <iomanip>
#include <vector>
#include <queue>
#include <cstdint>
#include <string>
#include <stdexcept>

#include "signed_graph.h"

namespace {
int count_components_any(const SignedGraph& G) {
    const int n = G.vertex_count();
    std::vector<char> vis(n, 0);
    int comps = 0;
    std::vector<int> q;
    q.reserve(n);

    for (int s = 0; s < n; ++s) {
        if (vis[s]) continue;
        ++comps;
        q.clear();
        q.push_back(s);
        vis[s] = 1;
        for (size_t qi = 0; qi < q.size(); ++qi) {
            const int u = q[qi];
            for (int v : G.neighbors(u)) {       // sign-agnostic adjacency
                if (!vis[v]) { vis[v] = 1; q.push_back(v); }
            }
        }
    }
    return comps;
}
}

int main(int argc, char** argv) {
    try {
        if (argc < 3 || argc > 4) {
            std::cerr << "usage: " << argv[0] << " <n:int>=vertices <p:double in [0,1]> [seed:uint64]\n";
            return 1;
        }
        const int n = std::stoi(argv[1]);
        const double p = std::stod(argv[2]);
        const uint64_t seed = (argc == 4) ? static_cast<uint64_t>(std::stoull(argv[3])) : std::random_device{}();

        if (n < 3) { std::cerr << "error: n must be >= 3\n"; return 2; }
        if (!(p >= 0.0 && p <= 1.0)) { std::cerr << "error: p must be in [0,1]\n"; return 3; }

        auto G = SignedGraph::make_random(n, p, seed);

        const int nV = G->vertex_count();
        const int m  = G->edge_count();
        const int mpos = G->edge_count(+1);
        const int mneg = G->edge_count(-1);

        const double denom = (nV > 1) ? (static_cast<double>(nV) * (nV - 1) / 2.0) : 1.0;
        const double dens      = (denom > 0.0) ? (static_cast<double>(m)    / denom) : 0.0;
        const double dens_pos  = (denom > 0.0) ? (static_cast<double>(mpos) / denom) : 0.0;
        const double dens_neg  = (denom > 0.0) ? (static_cast<double>(mneg) / denom) : 0.0;

        const int ncc = count_components_any(*G);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "vertices: " << nV << "\n";
        std::cout << "edges:    " << m  << "\n";
        std::cout << "edges(+): " << mpos << "\n";
        std::cout << "edges(-): " << mneg << "\n";
        std::cout << "density:      " << dens     << "\n";
        std::cout << "density(+):   " << dens_pos << "\n";
        std::cout << "density(-):   " << dens_neg << "\n";
        std::cout << "components:   " << ncc << "\n";

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "fatal: " << ex.what() << "\n";
        return 99;
    }
}
