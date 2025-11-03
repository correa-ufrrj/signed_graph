// edge.h
#pragma once

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
    double salience; // in [0,1]

    bool operator==(const SignedEdge& other) const {
        // Equality by undirected endpoints only (ignore sign/salience)
        return (points.first == other.points.first && points.second == other.points.second) ||
               (points.first == other.points.second && points.second == other.points.first);
    }
};

struct EdgeHash {
    size_t operator()(const Edge& edge) const {
        int a = std::min(edge.first, edge.second);
        int b = std::max(edge.first, edge.second);
        return std::hash<int>()(a) ^ std::hash<int>()(b << 1);
    }
};

std::ostream& operator<<(std::ostream& os, const SignedEdge& e);

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
