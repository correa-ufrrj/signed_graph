# README.md

## Frustration Index Solver

This project provides a complete C++ implementation to compute the **frustration index** of signed graphs using:
- [igraph](https://igraph.org/c/) for graph handling
- [IBM CPLEX 22.1.2](https://www.ibm.com/products/ilog-cplex-optimization-studio) for MIP solving with branch-and-cut

### 📁 Project Structure
```
├── CMakeLists.txt
├── config.h.in
├── main.cpp
├── signed_graph.cpp / .h
├── signed_graph_mip.cpp / .h
├── frustration_model.cpp / .h
├── frustration_model_x.cpp / .h
├── frustration_model_xy.cpp / .h
├── modularity_model.cpp / .h
├── scripts/
│   └── run_batch.sh
│   └── run_all_experiments.sh
```

---

## 🧱 Build Instructions

```bash
mkdir build && cd build
cmake ..
make -j
```

This compiles:
- `frustration_index`: command line executable

> CPLEX 22.1.2 must be installed and the `CPLEX_ROOT` updated in `CMakeLists.txt`.

---

## 🚀 Running the Solver
```bash
./frustration_index <graph.csv> [graph2.csv ...] [--model=xy|x|mod] [--netdeg] [--negcycles] [--negcyclescover] [--triangles] [--svg] [--help]
```
- `--model=xy|x|mod`: Choose model variant (default: xy)
- `--netdeg`: Enable net degree inequalities
- `--negcycles`: Enable negative cycle inequalities
- `--negcyclescover`: Enable negative cycle inequalities as a cover
- `--triangles`: Enable all negative triangle inequalities
- `--svg`: Export graph solution as SVG visualization
- `--help`: Show this help message and exit

Outputs:
- Graph statistics
- Objective value (frustration index)

---

## 🔎 Batch Evaluation
Use the script to run the solver on many `.csv` files recursively:

```bash
./scripts/run_batch.sh path/to/dataset_root --no-build --svg
```
- `--no-build`: skip building the project

Adds SVG diagrams to results folders if enabled.

- Supports nested subfolders
- Logs:
  - `results_summary.txt`: all runs
  - `logs/<relative>.log`: per-instance output
  - `logs/<folder>/summary.csv`: per-folder stats

Skips already-processed files with prompt. Use `Y` to skip all.
---

## 🔁 Batch Cut Strategy Evaluation

This script automates batch evaluation of different cut strategies for a dataset folder structure.

### ✅ Usage

```bash
./scripts/run_all_experiments.sh /path/to/dataset_root
```

Runs all selected solver variants on datasets in subfolders using `run_batch.sh`.

### 📁 Input

* `DATASET_DIR`: Root folder containing input `.csv` files (may contain nested folders).

### 🧠 Enabled Strategy (Active):

```bash
[Batch] Negative-cycle cuts only...
```

### 💤 Skipped Variants (Commented Out)

Uncomment to enable other strategy combinations:

```bash
[Batch] No cuts
./run_batch.sh "$DATASET_DIR" --no-build

[Batch] Net-degree cuts only
./run_batch.sh "$DATASET_DIR" --netdeg --no-build

[Batch] All triangles cuts only
./run_batch.sh "$DATASET_DIR" --triangles --no-build

[Batch] Negative-cycle cover cuts only
./run_batch.sh "$DATASET_DIR" --negcyclescover --no-build

[Batch] Net-degree + triangles cuts
./run_batch.sh "$DATASET_DIR" --netdeg --triangles --no-build

[Batch] Net-degree + negative-cycle cuts
./run_batch.sh "$DATASET_DIR" --netdeg --negcycles --no-build

[Batch] Net-degree + negative-cycle cover cuts
./run_batch.sh "$DATASET_DIR" --netdeg --negcyclescover --no-build

[Batch] Triangles + negative-cycle cuts
./run_batch.sh "$DATASET_DIR" --triangles --negcycles --no-build

[Batch] Triangles + negative-cycle cover cuts
./run_batch.sh "$DATASET_DIR" --triangles --negcyclescover --no-build

[Batch] Net-degree + triangles + negative-cycle cuts
./run_batch.sh "$DATASET_DIR" --netdeg --triangles --negcycles --no-build
```

### 🛠 Underlying Script

This script invokes:

```bash
run_batch.sh
```

### 📊 Output Overview

Each batch execution generates:

* `results_summary.txt`: All run results across datasets.
* `logs/<relative>.log`: Individual instance logs.
* `logs/<folder>/summary.csv`: Folder-wise aggregated stats.

### ✨ Features

* Automatically detects nested folders
* Optional build skip (`--no-build`)
* Easily extendable by uncommenting cut strategies
* Prompt before re-processing existing outputs

---

## 📊 Output Format
Example row in `results_summary.txt`:
```
graph1.csv,realistic/graph1.csv,9,0,0.2451
```

---

## 📂 CSV Input Format
Each line of the input file should be:
```
<node1>,<node2>,<sign>
```
Where `sign` is `+1` or `-1`.

Example:
```
0,1,1
1,2,-1
2,0,1
```
---

## 🧪 Testing Tips
- Include both small and dense graphs
- Test graphs with known frustration indices

---

## 🛠️ Troubleshooting
- Make sure CPLEX shared libraries are visible in `LD_LIBRARY_PATH`
- Ensure `.csv` inputs are well-formatted (no blank lines)

