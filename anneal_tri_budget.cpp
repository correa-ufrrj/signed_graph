#include "anneal_tri_budget.h"

namespace {
  // State updated by the owner before each batch.
  tri_anneal::Params g_params{};
  double g_mean_violation_prev = 0.0;
}

extern "C" void TBB_set_mean_tri_violation(double mean_violation_prev) {
  g_mean_violation_prev = mean_violation_prev;
}

extern "C" void TBB_set_anneal_params(double gamma_min, double gamma_max, double v0, double tau, int B_floor) {
  g_params.gamma_min = gamma_min;
  g_params.gamma_max = gamma_max;
  g_params.v0        = v0;
  g_params.tau       = tau;
  g_params.B_floor   = B_floor;
}

extern "C" int TBB_budget_override(int base_B_tri_prev) {
  if (base_B_tri_prev <= 0) {
    return g_params.B_floor;
  }
  return tri_anneal::anneal_B_tri(base_B_tri_prev, g_mean_violation_prev, g_params);
}
