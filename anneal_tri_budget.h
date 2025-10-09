#pragma once

#include <cmath>
#include <algorithm>

namespace tri_anneal {

struct Params {
  double gamma_min = 0.75;  // lower bound (0<gamma_min<1)
  double gamma_max = 0.90;  // upper bound (gamma_min<=gamma_max<1)
  double v0        = 0.0;   // target violation
  double tau       = 0.25;  // softness > 0
  int    B_floor   = 1;     // minimum cap
};

inline double logistic(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

inline int anneal_B_tri(int B_prev, double mean_violation_prev, const Params& p) {
  const double T = logistic((mean_violation_prev - p.v0) / p.tau);
  const double gamma = p.gamma_min + (p.gamma_max - p.gamma_min) * T;
  const double B_real = std::ceil(gamma * static_cast<double>(B_prev));
  const int B_next = std::max(p.B_floor, static_cast<int>(B_real));
  return B_next;
}

} // namespace tri_anneal

extern "C" {
// Set the mean triangle violation from the previous round.
void TBB_set_mean_tri_violation(double mean_violation_prev);
// Set annealing parameters (see tri_anneal::Params fields).
void TBB_set_anneal_params(double gamma_min, double gamma_max, double v0, double tau, int B_floor);
// Budget override hook used by TriangleBucketBatch (weak in the batch TU, strong here).
int TBB_budget_override(int base_B_tri_prev);
}
