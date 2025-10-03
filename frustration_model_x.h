// File: frustration_model_x.h
#pragma once

#include "frustration_model.h"

class FrustrationModelX : public FrustrationModel {
public:
    explicit FrustrationModelX(SignedGraphForMIP& g, int cut_flags = 0);

    void build() override;
    void solve() override;
    void export_solution(const std::string& file_prefix, bool with_svg) const;

private:
    std::vector<IloBoolVar> x;
    IloNumVar z;

	inline IloRange build_edge_partition_cut(const std::vector<std::pair<int, int>>& ordered_pairs);
    bool vertex_partition_compare(int u, int v, const IloNumArray& x_vals) const;
	inline std::vector<std::pair<int, int>> lift_edge_partition_cuts(const std::vector<std::pair<int, int>>& undecided, const IloNumArray& x_vals);
	inline std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>> generate_edge_partition_ineq(const IloNumArray& x_vals);
    inline std::vector<std::pair<IloRange, std::vector<std::pair<int, int>>>>  generate_cycle_ineq(const std::vector<Edge>& all_edges, const IloNumArray& x_vals);
    std::vector<std::pair<IloRange, std::string>> generate_cycle_cuts(IloEnv& env, const std::vector<Edge>& all_edges) override;
    std::vector<std::pair<IloRange, std::string>> generate_pending_triangle_cuts(IloEnv& env, const TriangleInequalities& t);

    class EdgePartitionCutGenerator : public IloCplex::UserCutCallbackI {
    public:
        FrustrationModelX& owner;

        EdgePartitionCutGenerator(IloEnv env,
                                  FrustrationModelX& owner);
        ~EdgePartitionCutGenerator();

        void main() override;

        IloCplex::CallbackI* duplicateCallback() const override {
            return (new (getEnv()) EdgePartitionCutGenerator(getEnv(), owner));
        }
    };

    class NegativeCycleCutGenerator : public IloCplex::UserCutCallbackI {
    public:
        FrustrationModelX& owner;

        NegativeCycleCutGenerator(IloEnv env,
                                  FrustrationModelX& owner);
        ~NegativeCycleCutGenerator();

        void main() override;

        IloCplex::CallbackI* duplicateCallback() const override {
            return (new (getEnv()) NegativeCycleCutGenerator(getEnv(), owner));
        }
    };

    class LazyEdgePartitionSeparator : public IloCplex::LazyConstraintCallbackI {
    public:
        FrustrationModelX& owner;

        LazyEdgePartitionSeparator(IloEnv env,
                           FrustrationModelX& owner);
        ~LazyEdgePartitionSeparator();

        void main() override;

        IloCplex::CallbackI* duplicateCallback() const override {
            return (new (getEnv()) LazyEdgePartitionSeparator(getEnv(), owner));
        }
    };

    class SwitchingHeuristicCallback : public IloCplex::HeuristicCallbackI {
    public:
        FrustrationModelX& owner;

        SwitchingHeuristicCallback(IloEnv env,
                                   FrustrationModelX& owner)
            : IloCplex::HeuristicCallbackI(env), owner(owner) {}

        void main() override;

        IloCplex::CallbackI* duplicateCallback() const override {
            return (new (getEnv()) SwitchingHeuristicCallback(getEnv(), owner));
        }
    };
};
