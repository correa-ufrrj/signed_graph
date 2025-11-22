// edge.h
#pragma once

#define IGRAPH_ENABLE_LGPL
#include <igraph/igraph.h>
#include <igraph/igraph_vector.h>
#include <igraph/igraph_attributes.h>

#include <vector>
#include <unordered_map>
#include <iostream>
#include <utility>

struct Edge {
    int first, second;
    Edge() : first(0), second(0) {}
    Edge(int a, int b) : first(std::min(a, b)), second(std::max(a, b)) {}

    bool operator==(const Edge& other) const {
        return first == other.first && second == other.second;
    }

    bool operator<(const Edge& other) const {
        return (first < other.first) || (first == other.first && second < other.second);
    }

    bool is_adjacent_to(const Edge& other) const {
        return first == other.first || first == other.second ||
               second == other.first || second == other.second;
    }
};

std::ostream& operator<<(std::ostream& os, const Edge& e);

struct SignedEdge {
    Edge points;
    int  sign;

    bool operator==(const SignedEdge& other) const {
        // Equality by undirected endpoints only (ignore sign)
        return (points.first == other.points.first && points.second == other.points.second) ||
               (points.first == other.points.second && points.second == other.points.first);
    }
};

struct EdgePolarity {
    Edge points;
    double  polarity;

    bool operator==(const EdgePolarity& other) const {
        // Equality by undirected endpoints only (ignore polarity)
        return (points.first == other.points.first && points.second == other.points.second) ||
               (points.first == other.points.second && points.second == other.points.first);
    }
};

struct EdgeHash {
    size_t operator()(const Edge& edge) const {
        return std::hash<int>()(edge.first) ^ std::hash<int>()(edge.second << 1);
    }
};

using EdgeMapType = std::unordered_map<Edge, int, EdgeHash>;

std::ostream& operator<<(std::ostream& os, const SignedEdge& e);

class EdgeIndexesView {
private:
    const igraph_t* g; const EdgeMapType& map;
public:
    explicit EdgeIndexesView(const EdgeMapType& m, const igraph_t* graph) : g(graph), map(m) {}
    class const_iterator {
    private: const igraph_t* g; igraph_integer_t eid; 
    public:
        using value_type = std::pair<Edge, int>;
        const_iterator(const igraph_t* g, igraph_integer_t eid) : g(g), eid(eid) {}
        value_type operator*() const { igraph_integer_t u,v; igraph_edge(g, eid, &u, &v); return {{(int)u,(int)v}, (int)eid}; }
        const_iterator& operator++() { ++eid; return *this; }
        bool operator!=(const const_iterator& other) const { return eid != other.eid; }
    };
    const_iterator begin() const { return const_iterator{g, 0}; }
    const_iterator end()   const { return const_iterator{g, igraph_ecount(g)}; }
    int operator[](const Edge& e) const { auto it = map.find(e); if (it == map.end()) throw std::out_of_range("Edge not found"); return it->second; }
};

class NegativeCycle {
public:
    NegativeCycle(const Edge& neg_edge, std::vector<Edge>&& pos_edges);

    const Edge& neg_edge() const;
    const std::vector<Edge>& pos_edges() const;

private:
    Edge neg_edge_;
    std::vector<Edge> pos_edges_;
    
};

// Inline definitions to satisfy linker
inline NegativeCycle::NegativeCycle(const Edge& neg_edge, std::vector<Edge>&& pos_edges)
    : neg_edge_(neg_edge), pos_edges_(std::move(pos_edges)) {}

inline const Edge& NegativeCycle::neg_edge() const { return neg_edge_; }
inline const std::vector<Edge>& NegativeCycle::pos_edges() const { return pos_edges_; }

std::ostream& operator<<(std::ostream& os, const NegativeCycle& nc);
