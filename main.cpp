// File: src/main.cpp
#include "frustration_model_xy.h"
#include "frustration_model_x.h"
#include "modularity_model.h"
#include "signed_graph_mip.h"
#include <filesystem>
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <optional>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include "fm_debug.h"

static void print_build_flags() {
  std::cout << "[BuildFlags] EXP=" << FM_EXPERIMENT_MODE
            << " PIPE=" << FM_PROBE_PIPELINE
            << " SWITCH=" << FM_PROBE_SWITCH
            << " PHASE=" << FM_PROBE_PHASE
            << " TBB=" << FM_PROBE_TBB
            << " SP=" << FM_PROBE_SP
            << " CB=" << FM_PROBE_CALLBACK
            << " GREEDY=" << FM_PROBE_GREEDY
            << " ROUND=" << FM_PROBE_ROUND << "\n";
}

static void print_help() {
    std::cout << "Usage:\n"
              << "  ./frustration_index <graph.csv> [graph2.csv ...] [opts]\n"
              << "  ./frustration_index --rand-mip n p q [seed] [opts]\n\n"
              << "Options:\n"
              << "  --model=xy|x|mod    Choose model variant (default: xy)\n"
              << "  --netdeg            Enable net degree inequalities\n"
              << "  --negcycles         Enable negative cycle inequalities\n"
              << "  --triangles         Enable all negative triangle inequalities\n"
              << "  --svg               Export graph solution as SVG visualization\n"
              << "  --help              Show this help message and exit\n\n"
              << "Notes:\n"
              << "  Passing both --negcycles and --triangles enables the TriangleCyclePipeline.\n"
              << "  Stats are always printed for every run (files and random graphs).\n";
}

namespace {
struct Job {
    enum class Kind { File, RandMIP } kind;
    std::string filename;  // if File
    int n = 0;
    double p = 0.0;
    double q = 1.0;
    uint64_t seed = 0;
    bool has_seed = false;
};

static bool parse_int(const std::string& s, int& out) {
    char* end=nullptr; long v = std::strtol(s.c_str(), &end, 10);
    if (!end || *end) return false; out = static_cast<int>(v); return true;
}
static bool parse_u64(const std::string& s, uint64_t& out) {
    char* end=nullptr; unsigned long long v = std::strtoull(s.c_str(), &end, 10);
    if (!end || *end) return false; out = static_cast<uint64_t>(v); return true;
}
static bool parse_double(const std::string& s, double& out) {
    char* end=nullptr; double v = std::strtod(s.c_str(), &end);
    if (!end || *end) return false; out = v; return true;
}

static std::string make_prefix_for_random(const Job& j) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);
    oss << "rand_mip_n" << j.n << "_p" << j.p << "_q" << j.q;
    if (j.has_seed) oss << "_seed" << j.seed;
    return oss.str();
}

static int count_components_any(const SignedGraphForMIP& G) {
    const int n = G.vertex_count();
    std::vector<char> vis(n, 0);
    std::vector<int> q; q.reserve(n);
    int comps = 0;
    for (int s = 0; s < n; ++s) {
        if (vis[s]) continue;
        ++comps; q.clear(); q.push_back(s); vis[s] = 1;
        for (size_t qi = 0; qi < q.size(); ++qi) {
            const int u = q[qi];
            for (int v : G.neighbors(u)) {
                if (!vis[v]) { vis[v] = 1; q.push_back(v); }
            }
        }
    }
    return comps;
}

// File: src/main.cpp  (only the print_stats() body changes)
static void print_stats(const SignedGraphForMIP& G) {
    const int nV   = G.vertex_count();
    const int m    = G.edge_count();
    const int mpos = G.edge_count(+1);
    const int mneg = G.edge_count(-1);

    const double denom_all = (nV > 1) ? (double)nV * (nV - 1) / 2.0 : 1.0;
    const double dens      = (denom_all > 0.0) ? (m    / denom_all) : 0.0;
    const double dens_pos  = (denom_all > 0.0) ? (mpos / denom_all) : 0.0;
    const double dens_neg  = (denom_all > 0.0) ? (mneg / denom_all) : 0.0;

    const int components   = count_components_any(G);

    // Frustration density: ratio (# frustrated edges) / (# edges)
    const int frustrated   = G.frustrated_edges();                // inherited from SignedGraph
    const double frust_den = (m > 0) ? (double)frustrated / m : 0.0;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "[STATS] vertices="    << nV
              << " edges="              << m
              << " pos="                << mpos
              << " neg="                << mneg
              << " density="            << dens
              << " density(+)="         << dens_pos
              << " density(-)="         << dens_neg
              << " components="         << components
              << " frustrated="         << frustrated
              << " frust_density="      << frust_den
              << "\n";
}

static std::unique_ptr<SignedGraphForMIP> build_from_job(const Job& job) {
    if (job.kind == Job::Kind::File) {
        return std::make_unique<SignedGraphForMIP>(job.filename);
    }
    const uint64_t seed = job.has_seed ? job.seed : std::random_device{}();
    return SignedGraphForMIP::make_random(job.n, job.p, job.q, seed);
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    int cut_flags = FrustrationModel::NO_CUTS;
    bool with_svg = false;
    std::string model_type = "xy";
    std::vector<Job> jobs;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help") {
            print_help();
            return 0;
        } else if (arg == "--netdeg") {
            cut_flags |= FrustrationModel::NET_DEGREE_CUTS;
        } else if (arg == "--negcycles") {
            cut_flags |= FrustrationModel::NEGATIVE_CYCLE_CUTS;
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
        } else if (arg == "--rand-mip") {
            if (i + 3 > argc) { std::cerr << "error: --rand-mip needs n p q [seed]\n"; return 2; }
            Job j; j.kind = Job::Kind::RandMIP;
            if (!parse_int(argv[++i], j.n)) { std::cerr << "bad n\n"; return 2; }
            if (!parse_double(argv[++i], j.p)) { std::cerr << "bad p\n"; return 2; }
            if (!parse_double(argv[++i], j.q)) { std::cerr << "bad q\n"; return 2; }
            if (i + 1 < argc && argv[i+1][0] != '-') { j.has_seed = parse_u64(argv[++i], j.seed); }
            jobs.push_back(j);
        } else {
            jobs.push_back(Job{ Job::Kind::File, arg });
        }
    }

    if (model_type == "mod" && (cut_flags & FrustrationModel::NET_DEGREE_CUTS)) {
        std::cerr << "Error: --netdeg is not supported with --model=mod.\n";
        return 1;
    }

    print_build_flags();

    if ((cut_flags & FrustrationModel::NEGATIVE_CYCLE_CUTS) &&
        (cut_flags & FrustrationModel::TRIANGLE_CUTS)) {
        std::cout << "[Flags] Using TriangleCyclePipeline (negcycles + triangles)\n";
    }

    if (jobs.empty()) {
        std::cerr << "Nothing to process. Provide CSVs or a random spec.\n";
        return 1;
    }

    for (const auto& job : jobs) {
        try {
            std::string label;
            std::string prefix;
            if (job.kind == Job::Kind::File) {
                label = job.filename;
                prefix = std::filesystem::path(job.filename).stem().string();
            } else {
                label = "[Random MIP]";
                prefix = make_prefix_for_random(job);
            }

            std::cout << "\n[Processing] " << label;
            if (cut_flags == FrustrationModel::NO_CUTS) std::cout << " (no cuts)";
            else std::cout << " (cuts=" << cut_flags << ")";
            std::cout << std::endl;

            std::unique_ptr<SignedGraphForMIP> g = build_from_job(job);
            g->print_info();
            print_stats(*g); // ALWAYS print stats

            std::unique_ptr<FrustrationModel> model;
            if (model_type == "xy") {
                model = std::make_unique<FrustrationModelXY>(*g, cut_flags);
            } else if (model_type == "x") {
                model = std::make_unique<FrustrationModelX>(*g, cut_flags);
            } else { // mod
                model = std::make_unique<ModularityModel>(*g, cut_flags);
            }

            if (auto* xy = dynamic_cast<FrustrationModelXY*>(model.get())) {
                SeparationConfig C;
                C.modes.use_triangles = (cut_flags & FrustrationModel::TRIANGLE_CUTS) != 0;
                C.modes.use_negcycles = (cut_flags & FrustrationModel::NEGATIVE_CYCLE_CUTS) != 0;
                xy->configure_separation(C);
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

            std::string cut_suffix;
            if (cut_flags == FrustrationModel::NO_CUTS) {
                cut_suffix = "_nocuts";
            } else {
                if (cut_flags & FrustrationModel::NET_DEGREE_CUTS) cut_suffix += "_netdeg";
                if (cut_flags & FrustrationModel::NEGATIVE_CYCLE_CUTS) cut_suffix += "_negcycles";
                if (cut_flags & FrustrationModel::TRIANGLE_CUTS)      cut_suffix += "_triangles";
            }

            std::string out_prefix = prefix + cut_suffix;
            if (model_type != "xy") out_prefix += "_" + model_type;
            model->export_solution(out_prefix, with_svg);
            print_stats(*g); // ALWAYS print stats

        } catch (const std::exception& e) {
            if (job.kind == Job::Kind::File)
                std::cerr << "Error processing " << job.filename << ": " << e.what() << std::endl;
            else
                std::cerr << "Error in random job: " << e.what() << std::endl;
        }
    }

    return 0;
}
