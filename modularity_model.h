// File: modularity_model.h
#pragma once

#include "frustration_model_xy.h"

class ModularityModel : public FrustrationModelXY {
public:
    // Constructor: enforces no NET_DEGREE_CUTS in cut_flags
    explicit ModularityModel(SignedGraphForMIP& g, int cut_flags = 0);
	~ModularityModel();

    // Builds the MIP model with custom modularity weights
    void build() override;

    // Solves the optimization and stores the modularity-related result
    void solve() override;

	// To solution export
	double get_frustration_index() const override;

    // Getter for the modularity score derived from the FrustrationModel result
    double get_modularity_value() const;

private:
    SignedGraphForMIP* cloned_graph_ptr = nullptr;
    double modularity_value = 0.0;
};
