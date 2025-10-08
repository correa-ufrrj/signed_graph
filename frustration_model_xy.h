// File: frustration_model_xy.h
#pragma once

#include "frustration_model.h"
// Lightweight headers used by the reheat pool / stateless dedup key.
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>

// ---------------------------------------------------------------------
// Switching-invariant, orientation-aware key for cycle identification.
// Header-only so callbacks and owner share the same type.
// ---------------------------------------------------------------------
namespace fmkey {
  struct CycleKey {
    std::vector<int> y_idx;  // canonical y-indices of the cycle
    int rhs{0};
    bool reversed{false};
  };
  struct CycleKeyHash {
    size_t operator()(const CycleKey& k) const noexcept {
      std::uint64_t h = 1469598103934665603ull;                 // FNV-1a
      auto mix = [&](std::uint64_t v){ h ^= v; h *= 1099511628211ull; };
      for (int v : k.y_idx) mix(static_cast<std::uint64_t>(static_cast<std::uint32_t>(v)));
      mix(static_cast<std::uint64_t>(static_cast<std::uint32_t>(k.rhs)));
      mix(static_cast<std::uint64_t>(k.reversed ? 1u : 0u));
      return static_cast<size_t>(h);
    }
  };
  struct CycleKeyEq { bool operator()(const CycleKey& a, const CycleKey& b) const noexcept { return a.rhs==b.rhs && a.reversed==b.reversed && a.y_idx==b.y_idx; } };
} // namespace fmkey

class FrustrationModelXY : public FrustrationModel {
public:
    explicit FrustrationModelXY(SignedGraphForMIP& g, int cut_flags = 0);

    void build() override;
    void solve() override;
    void export_solution(const std::string& file_prefix, bool with_svg) const;

private:
    std::vector<IloBoolVar> x;
    std::vector<IloNumVar> y;
    bool triangle_cut_added_last_round = false;

    // ===== Reheat pool: switching-agnostic storage of previously seen cycles =====
	struct ReheatItem {
	  // canonical, switching-agnostic representation
	  std::vector<int> cyc_vertices;   // ordered cycle
	  int ttl = 3;
	  double last_viol = 0.0;
	};
    // Lives on the owner so it persists across callback duplications/rounds
    std::unordered_map<fmkey::CycleKey, ReheatItem, fmkey::CycleKeyHash, fmkey::CycleKeyEq> reheat_pool_;
    std::unordered_set<fmkey::CycleKey, fmkey::CycleKeyHash, fmkey::CycleKeyEq> in_model_keys_;   // already added
    std::unordered_set<fmkey::CycleKey, fmkey::CycleKeyHash, fmkey::CycleKeyEq> recent_keys_;     // de-dupe within recent rounds
    struct ReheatParams {
        int max_ttl = 3;
        int keep_cap = 5000;          // safety cap for pool size
        int seed_cap_per_edge = 8;    // at most this many reheated items per neg edge
    } reheat_;

    // Seed buckets from reheated items that are currently violated; defined in .cpp
    void reheat_seed_buckets_(
        /*in/out*/ std::unordered_map<int, std::vector<int>>& by_neg_eid,
        const std::vector<char>& neg_edge_mask,
        const IloNumArray& xhat,
        const IloNumArray& yhat);
    // After selection, decay TTLs and prune the pool; defined in .cpp
    void reheat_after_selection_(
        const std::vector<int>& chosen_cids,
        const std::vector<std::vector<int>>& cand_pos_edges,
        const IloNumArray& xhat,
        const IloNumArray& yhat);

    // Custom cut logic
	IloRange generate_cycle_cut_standard(IloEnv& env, const std::vector<Edge>& all_edges);
	std::vector<std::pair<IloRange, std::string>> generate_cycle_cuts(IloEnv& env, const std::vector<Edge>& all_edges) override;
    std::vector<std::pair<IloRange, std::string>> generate_pending_triangle_cuts(IloEnv& env, const TriangleInequalities& t) override;
	std::vector<std::pair<IloRange, std::string>> generate_positive_triangle_cuts(IloEnv& env, const std::vector<Edge>& all_edges);

    class NegativeTriangleCutGenerator : public IloCplex::UserCutCallbackI {
    public:
        FrustrationModelXY& owner;

        NegativeTriangleCutGenerator(IloEnv env,
                                     FrustrationModelXY& owner);

        ~NegativeTriangleCutGenerator();

        void main() override;

        IloCplex::CallbackI* duplicateCallback() const override {
            return new (getEnv()) NegativeTriangleCutGenerator(getEnv(), owner);
        }
    };

    class NegativeCycleCutGenerator : public IloCplex::UserCutCallbackI {
    public:
        FrustrationModelXY& owner;

        NegativeCycleCutGenerator(IloEnv env,
                                  FrustrationModelXY& owner);
		~NegativeCycleCutGenerator();

        void main() override;

        IloCplex::CallbackI* duplicateCallback() const override {
            return (new (getEnv()) NegativeCycleCutGenerator(getEnv(), owner));
        }
    };

    class ConditionalCycleCutGenerator : public NegativeCycleCutGenerator {
    public:
        ConditionalCycleCutGenerator(IloEnv env, FrustrationModelXY& owner)
            : NegativeCycleCutGenerator(env, owner) {}

        void main() override {
            if (!owner.triangle_cut_added_last_round) {
                NegativeCycleCutGenerator::main();
            }
        }

        IloCplex::CallbackI* duplicateCallback() const override {
            return new (getEnv()) ConditionalCycleCutGenerator(getEnv(), owner);
        }
    };
    
    class SwitchingHeuristicCallback : public IloCplex::HeuristicCallbackI {
    private:
        FrustrationModelXY& owner;

    public:
        SwitchingHeuristicCallback(IloEnv env,
                                   FrustrationModelXY& owner)
            : IloCplex::HeuristicCallbackI(env), owner(owner) {}

        void main() override;

        IloCplex::CallbackI* duplicateCallback() const override {
            return (new (getEnv()) SwitchingHeuristicCallback(getEnv(), owner));
        }
    };
};
