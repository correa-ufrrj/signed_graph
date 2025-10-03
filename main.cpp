// File: main.cpp
#include "frustration_model_xy.h"
#include "frustration_model_x.h"
#include "modularity_model.h"
#include <filesystem>
#include <iostream>
#include <vector>
#include <string>
#include <memory>

void print_help() {
    std::cout << "Usage: ./frustration_index <graph.csv> [graph2.csv ...] [--model=xy|x|mod] [--netdeg] [--negcycles] [--negcyclescover] [--triangles] [--svg] [--help]\n";
    std::cout << "\nOptions:\n";
    std::cout << "  --model=xy|x|mod   Choose model variant (default: xy)\n";
    std::cout << "  --netdeg       Enable net degree inequalities\n";
    std::cout << "  --negcycles    Enable negative cycle inequalities\n";
    std::cout << "  --negcyclescover Enable negative cycle inequalities as a cover\n";
    std::cout << "  --triangles    Enable all negative triangle inequalities\n";
    std::cout << "  --svg          Export graph solution as SVG visualization\n";
    std::cout << "  --help         Show this help message and exit\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    int cut_flags = FrustrationModel::NO_CUTS;
    bool with_svg = false;
    std::vector<std::string> files;
    std::string model_type = "xy";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help") {
            print_help();
            return 0;
        } else if (arg == "--netdeg") {
            cut_flags |= FrustrationModel::NET_DEGREE_CUTS;
        } else if (arg == "--negcycles") {
            cut_flags |= FrustrationModel::NEGATIVE_CYCLE_CUTS;
        } else if (arg == "--negcyclescover") {
            cut_flags |= FrustrationModel::NEGATIVE_CYCLE_COVER_CUTS;
        } else if (arg == "--triangles") {
            cut_flags |= FrustrationModel::TRIANGLE_CUTS;
        } else if (arg == "--svg") {
            with_svg = true;
        } else if (arg.rfind("--model=", 0) == 0) {
		    model_type = arg.substr(8);
		    if (model_type != "xy" && model_type != "x" && model_type != "mod") {
		        std::cerr << "Invalid model type: " << model_type << "\n";
		        return 1;
		    }
		    std::cout << "Model type set to: " << model_type << "\n";
		} else {
            files.push_back(arg);
        }
    }
    
    if (model_type == "mod" && (cut_flags & FrustrationModel::NET_DEGREE_CUTS)) {
    	std::cerr << "Error: --netdeg is not supported with --model=mod.\n";
    	return 1;
	}

    for (const auto& filename : files) {
        try {
            std::cout << "\n[Processing] " << filename;
            if (cut_flags == FrustrationModel::NO_CUTS) std::cout << " (no cuts)";
            else std::cout << " (cuts=" << cut_flags << ")";
            std::cout << std::endl;

            SignedGraphForMIP g(filename);
            g.print_info();

            std::unique_ptr<FrustrationModel> model;
			if (model_type == "xy") {
			    model = std::make_unique<FrustrationModelXY>(g, cut_flags);
			} else if (model_type == "x") {
			    model = std::make_unique<FrustrationModelX>(g, cut_flags);
			} else if (model_type == "mod") {
			    model = std::make_unique<ModularityModel>(g, cut_flags);
			}

            std::cout << "[Model] Created model." << std::endl;
            model->build();
            std::cout << "[Model] Built model." << std::endl;
            model->solve();
            std::cout << "[Model] Solved model." << std::endl;
	
			if (model_type == "mod") {
			    auto* mod_ptr = dynamic_cast<ModularityModel*>(model.get());
			    std::cout << "[Modularity] Q'(G) = " << mod_ptr->get_modularity_value() << std::endl;
			} else {
			    std::cout << "[Frustration Index] = " << model->get_frustration_index() << std::endl;
			}
            model->print_solution();

            std::string prefix = std::filesystem::path(filename).stem().string();
            std::string cut_suffix;
            if (cut_flags == FrustrationModel::NO_CUTS) {
                cut_suffix = "_nocuts";
            } else {
                if (cut_flags & FrustrationModel::NET_DEGREE_CUTS) cut_suffix += "_netdeg";
                if (cut_flags & FrustrationModel::NEGATIVE_CYCLE_CUTS) cut_suffix += "_negcycles";
                if (cut_flags & FrustrationModel::NEGATIVE_CYCLE_COVER_CUTS) cut_suffix += "_negcyclescover";
                if (cut_flags & FrustrationModel::TRIANGLE_CUTS) cut_suffix += "_triangles";
            }

			prefix += cut_suffix;
			if (model_type != "xy") {
			    prefix += "_" + model_type;
			}
            model->export_solution(prefix, with_svg);

        } catch (const std::exception& e) {
            std::cerr << "Error processing " << filename << ": " << e.what() << std::endl;
        }
    }

    return 0;
}
