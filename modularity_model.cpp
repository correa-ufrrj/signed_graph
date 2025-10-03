// File: modularity_model.cpp

#include "modularity_model.h"
#include "signed_graph_mip.h"
#include <stdexcept>
#include <cmath>
#include <iostream>

ModularityModel::ModularityModel(SignedGraphForMIP& g, int cut_flags)
    : FrustrationModelXY(*(new SignedGraphForMIP(&g, [&]() {
        std::cerr << "[ModularityModel] Computing modularity weights..." << std::endl;

        const auto& d_plus = g.plus_degrees_view();
        const auto& d_minus = g.minus_degrees_view();

        double m_plus = 0.0, m_minus = 0.0;
        for (int d : d_plus) m_plus += d;
        for (int d : d_minus) m_minus += d;
        m_plus /= 2.0;
        m_minus /= 2.0;

        std::vector<double> weights(g.edge_count());
        for (const auto& [edge, sign] : g.signs_view()) {
            int u = edge.first;
            int v = edge.second;

            double w = 0.0;
            if (sign > 0) {
                double term1 = (d_plus[u] * d_plus[v]) / (2.0 * m_plus);
                double term2 = (d_minus[u] * d_minus[v]) / (2.0 * m_minus);
                w = 1.0 - term1 + term2;
            } else {
                double term1 = (d_minus[u] * d_minus[v]) / (2.0 * m_minus);
                double term2 = (d_plus[u] * d_plus[v]) / (2.0 * m_plus);
                w = - 1.0 - term1 + term2;
            }

            int eid = g.edge_index()[edge];

            weights[eid] = w;
        }
        std::cerr << "[ModularityModel] Weight computation completed." << std::endl;
        return weights;
    }())), cut_flags), cloned_graph_ptr(&dynamic_cast<SignedGraphForMIP&>(this->graph)) {

    std::cerr << "[ModularityModel] Constructor entered." << std::endl;

    if (cut_flags & NET_DEGREE_CUTS) {
        throw std::invalid_argument("ModularityModel does not support NET_DEGREE_CUTS as a cut flag.");
    }

    if (!cloned_graph_ptr) {
        throw std::runtime_error("[ModularityModel] dynamic_cast to SignedGraphForMIP failed.");
    }
}

ModularityModel::~ModularityModel() {
    std::cerr << "[ModularityModel] Destructor called. Deleting cloned graph pointer." << std::endl;
    delete cloned_graph_ptr;
}

void ModularityModel::build() {
    std::cerr << "[ModularityModel] Build called." << std::endl;
    FrustrationModelXY::build();
}

void ModularityModel::solve() {
    std::cerr << "[ModularityModel] Solve called." << std::endl;
    FrustrationModelXY::solve();
    modularity_value = get_modularity_value();
    std::cerr << "[ModularityModel] Modularity value computed: " << modularity_value << std::endl;
}

double ModularityModel::get_modularity_value() const {
    const auto& weights = cloned_graph_ptr->weights_view();

    double sum_abs_weights = 0.0;
    double sum_signed_weights = 0.0;

    for (const auto& [edge, weight] : weights) {
        sum_abs_weights += std::abs(weight);
        sum_signed_weights += weight;
    }

    if (sum_abs_weights == 0.0) {
        throw std::runtime_error("Total absolute weight is zero, cannot compute modularity value.");
    }

    double fr_value = FrustrationModel::get_frustration_index();
    double q_prime = (sum_abs_weights / (2.0*cloned_graph_ptr->edge_count())) * (0.5 * (1.0 + (sum_signed_weights / sum_abs_weights)) - fr_value);

    return q_prime;
}

double ModularityModel::get_frustration_index() const {
	return get_modularity_value();
}
