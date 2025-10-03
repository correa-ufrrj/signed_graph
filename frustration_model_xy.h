// File: frustration_model_xy.h
#pragma once

#include "frustration_model.h"

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
