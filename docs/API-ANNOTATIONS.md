# docs/API-ANNOTATIONS.md
# API & Complexity Annotations (No Code Changes)

Scope: public interfaces in:
- `signed_graph.h` (previous section above)
- `modularity_model.h` (previous section above)
- `frustration_model_xy.h` (new here)
- `signed_graph_mip.h` + `NegativeCycleBatchStream` (new here)

---

## Core value types

### `struct Edge { int first, second; }`
- **Semantics:** Undirected edge between two vertices. Constructor normalizes `(a,b)` to `(min,max)`.
- **API:**
  - `Edge(int a,int b)`: O(1). **Pre:** `a` and `b` in `[0, n)`. **Post:** `first<=second`.
  - `operator==(const Edge&)`: O(1).
  - `is_adjacent_to(const Edge&)`: O(1). True if the two edges share an endpoint.
- **Hotspots:** Ensure all callers pass normalized vertex ids; constructor already normalizes.

### `struct SignedEdge { Edge points; int sign; }`
- **Semantics:** An edge plus its sign.
- **API:** `operator==(const SignedEdge&)` compares by endpoints only (ignores `sign`).
- **Hotspot:** Potential confusion—equality does not include sign. If used as a key or in sets, two edges with different signs compare equal (intentional?).

### `struct EdgeHash`
- **Semantics:** Hash for `Edge` keyed by endpoints independent of order.
- **Complexity:** O(1).
- **Hotspot:** `std::hash<int>(b << 1)` hashes the *shifted* `b`; okay, but shifting then hashing is unusual. Collisions are still fine for `unordered_map`.

---

## `class NegativeCycle`
- **Semantics:** A negative cycle represented by exactly one negative edge and a list of positive edges closing the cycle.
- **API:**
  - `NegativeCycle(const Edge& neg_edge, std::vector<Edge>&& pos_edges)`: O(k).
  - Accessors: `neg_edge()`, `pos_edges()`: O(1).
- **Use:** Validation/diagnostics (e.g., cycle correctness, edge-disjointness checks).

---

## `class SignedGraph`
- **Ownership/Lifetime:** Wraps an `igraph_t g` plus `weights` and derived degree stats. Destructor must `igraph_destroy(&g)`.
- **Copying:** Copy constructors are **protected** and pointer-based; external copying is via `clone()` → `unique_ptr<SignedGraph>`. **Hotspot:** Enforce Rule-of-Five; avoid accidental shallow copies of `igraph_t`.

### Construction & lifetime
- `SignedGraph(const std::string& file_path)`
  - **Purpose:** Build graph + weights from file.
  - **Complexity:** Depends on parser; typically O(n + m).
  - **Pre:** File exists, format matches expectations. **Fail:** Throws on IO/parse errors.

- `std::unique_ptr<SignedGraph> clone() const`
  - **Purpose:** Deep clone, optionally using protected ctors.
  - **Complexity:** O(n + m).
  - **Post:** Independent `igraph_t`, `weights`, degrees, maps.

- `virtual ~SignedGraph()`
  - **Purpose:** Release `igraph_t` and buffers.
  - **Hotspot:** Double-destroy if copies exist; ensure deep copy only via `clone()`.

### Topology/degree caches
- Private: `compute_degrees()`
  - **Purpose:** Fill `d_plus`, `d_minus`, `switched_d_plus`, `switched_d_minus`.
  - **Complexity:** O(m).
  - **Hotspot:** Definition of sign for `w=0.0` (see SignedEdgesView); negative zero treated as negative.

- `const std::vector<double>& plus_degrees_view() const` / `minus_degrees_view() const`
  - **Purpose:** Expose cached signed degrees.
  - **Complexity:** O(1) access; consumers iterate O(n).
  - **Pre:** Must be in sync with `weights` (ensure recomputed after switching).

- `const long double max_degree_vertex() const`
  - **Purpose:** Return **what exactly?** Name suggests a *vertex id*, type is `long double`.  
    **Hotspot:** Return type mismatches likely intent; if it returns a degree value, name should reflect it; if it returns a vertex id, type should be integer. Verify at call sites.

### Neighborhood/iteration
- `const std::vector<int> neighbors(int u) const`
  - **Purpose:** Vertex adjacency.
  - **Complexity:** O(deg(u)) to collect neighbors.
  - **Pre:** `0 ≤ u < n`.

- **Views** (zero-alloc iteration over edges + weights):
  - `SignedEdgesView signs_view() const`
    - **Iter::operator***: Builds `(Edge, sign)` using `igraph_edge(g,eid)` then sign from `weights[eid]`.  
      **Sign rule:** `-1` if `w<0` or `w==0 && signbit(w)`, else `+1`.
    - **Complexity:** `begin/end`: O(1). Iteration: O(1) per edge → O(m) total.
    - `operator[](eid)`: O(1). `operator[](Edge)`: expected O(1) average if uses index map (verify impl).
    - **Hotspot:** Negative zero is rare; ensure weight producers don’t emit `-0.0` accidentally.

  - `WeightsView weights_view() const`
    - **Iter::operator***: Returns `(Edge, weight)`; O(1)/edge.
    - `operator[](Edge)`: Expected O(1) avg if uses map.
    - **Hotspot:** Keep `weights.size() == igraph_ecount(g)` invariant.

  - `EdgeIndexesView edge_index() const`
    - **Purpose:** Map `Edge → eid`.
    - **Complexity:** Lookup O(1) average; building map O(m).
    - **Hotspot:** Ensure construction happens once and stays in sync after switching or graph edits.

### Frustration/triangles/switching
- `std::vector<int> negative_triangle_count_per_vertex() const`
  - **Purpose:** Per-vertex count of negative triangles.
  - **Complexity:** Depends on algorithm:  
    - With adjacency lists & node ordering: ~O(∑_{(u,v)} min(deg(u),deg(v))) ≈ O(m·√m) typical sparse.  
    - Worst-case dense: O(n³).  
  - **Hotspot:** Avoid recomputing across passes; cache if used in greedy tiebreaks.

- `int negative_triangle_count_of_vertex(int u) const`
  - **Purpose:** Count negative triangles incident to `u`.
  - **Complexity:** O(∑ over neighbors v of min(deg(u),deg(v))) ≈ O(m_u · √m) typical.
  - **Pre:** `0 ≤ u < n`.

- `int frustrated_edges(const std::vector<int>& partition) const`
  - **Purpose:** Count edges whose sign disagrees with `partition` (e.g., positive across cut or negative inside).
  - **Complexity:** O(m).
  - **Pre:** `partition.size()==n` with labels in `{0,1}` (or {-1,+1}—confirm).  
  - **Hotspot:** Label convention must match `switching_from_partition`.

- `std::vector<Edge> frustrated_edges_keys(const std::vector<int>& partition) const`
  - **Purpose:** Return the set of frustrated edges.
  - **Complexity:** O(m) + O(k) to output.
  - **Hotspot:** Consider memory churn if called in loops; prefer counting unless keys truly needed.

- `void switching_from_partition(const std::vector<int>& s)`
  - **Purpose:** Flip signs/weights according to vertex labels into `switched_weights` and update `switched_*` degrees.
  - **Complexity:** O(m).
  - **Pre:** `s.size()==n`.
  - **Post:** `switched_*` views reflect switched graph; `restore_switched_sign()` reverts.
  - **Hotspot:** Keep a clear invariant: base `weights` never mutated; all ops read from the correct active buffer.

- `void single_switching(int u)` **(public)** and `void single_switching(int u, igraph_vector_int_t* incident)` **(private)**
  - **Purpose:** Flip incident edges w.r.t. a single vertex.
  - **Complexity:** O(deg(u)).
  - **Hotspot:** Overload shadowing; ensure callers don’t accidentally call the private variant. Incident list reuse avoids repeated allocations—good for hot loops.

- `const std::vector<int> greedy_switching()`
  - **Purpose:** Run default greedy switching heuristic (no kicks).
  - **Complexity:** Typically O(m) per sweep; #sweeps small in practice.
  - **Post:** Returns final partition/labels.
  - **Hotspot:** Determinism (tie breaks). Ensure consistent seed or stable ordering.

- `std::vector<int> greedy_switching_base(const std::function<int(int,int)>& cmp_fn, const GreedyKickOptions& opts)`
  - **Purpose:** Configurable greedy with kick options.
  - **Complexity:** Baseline ~O(m) per pass; kicks add bounded extra work (`neighbor_cap`, `triangle_cap_per_u` limit worst-case).
  - **Key options (effects):**
    - `neg_edge_threshold_abs/frac`: Skip kicks unless sufficiently many negative edges exist; avoids noise on near-balanced graphs.
    - `use_weighted_degree`: Primary score from signed weighted degree; O(1) per candidate if cached.
    - `use_triangle_tiebreak` + `triangle_beta`: Adds O(triangles around u) work per tiebreak; capped.
    - `delta_m_minus_cap` / `delta_m_minus_penalty`: Enforce/penalize increases in negative edges; stabilizes convergence.
    - `edge_salience` + `kick_salience_bias`: Bias move scoring by exogenous edge importance in [0,1].
    - Caps (`neighbor_cap`, `triangle_cap_per_u`): Prevents quadratic blowups on hubs.
  - **Hotspots:**  
    - Avoid recomputing degrees/triangle counts from scratch per move—maintain deltas.  
    - `cmp_fn` must be pure/stable for reproducibility.

- `void restore_switched_sign()`
  - **Purpose:** Revert to original `weights` / degree caches.
  - **Complexity:** O(m).
  - **Hotspot:** Ensure all views after restore read base buffers.

### Graph info, validation, IO
- `int vertex_count() const`, `int edge_count() const`
  - **Complexity:** O(1).

- `void print_info() const`
  - **Purpose:** Human-readable summary (n, m, degree stats).
  - **Hotspot:** Avoid heavy computations/logging in large runs.

- `bool are_cycles_edge_disjoint(const std::vector<NegativeCycle>& cycles) const`
  - **Purpose:** Check pairwise edge disjointness.
  - **Complexity:** O(total_edges_in_cycles) average with a set; O(k²) naive.
  - **Hotspot:** Use `unordered_set<Edge,EdgeHash>` to keep O(1) inserts/lookups.

- `bool are_cycles_sign_correct(const std::vector<NegativeCycle>& cycles, bool expect_negative = true) const`
  - **Purpose:** Verify each cycle sign product matches `expect_negative`.
  - **Complexity:** O(total_edges_in_cycles).
  - **Hotspot:** Sign of zero weights; ensure the same sign convention as elsewhere.

- `void save_partition_svg(const std::vector<int>& partition, const std::string& filename, bool custom_layout) const`
  - **Purpose:** Persist SVG visualization.
  - **Complexity:** Layout dominates; IO O(n+m).
  - **Hotspot:** Non-deterministic layouts reduce reproducibility unless `custom_layout` fixes coordinates.

---

## Iteration helpers (operator<<)
- `std::ostream& operator<<(std::ostream&, const Edge&)` / `(const SignedEdge&)` / `(const NegativeCycle&)`
  - **Purpose:** Debug printing.
  - **Hotspot:** Avoid in tight loops; guard by verbosity.

---

## `class ModularityModel : public FrustrationModelXY`
- **Role:** MIP-based solver for modularity tailored to signed graphs, built atop the XY frustration model base.

### Lifecycle & interface
- `explicit ModularityModel(SignedGraphForMIP& g, int cut_flags = 0)`
  - **Purpose:** Initialize model for graph `g`; **enforces** that `NET_DEGREE_CUTS` is not in `cut_flags`.
  - **Complexity:** O(1) aside from shallow setup.
  - **Hotspot:** Validate `cut_flags`; document overrides to caller expectations.

- `~ModularityModel()`
  - **Purpose:** Cleanup resources; likely owns `cloned_graph_ptr` if a clone was created in `build()`.
  - **Hotspot:** Clear ownership: who allocates/frees `cloned_graph_ptr`? Avoid leaks/double-free; prefer smart pointers internally.

- `void build() override`
  - **Purpose:** Construct MIP variables/constraints/objective for modularity on (possibly cloned/switched) graph.
  - **Complexity:** Typically O(n + m) to build; memory O(n + m) variables/constraints (formulation-dependent).
  - **Pre:** Not already built or handle idempotently. Graph immutable during/after build.
  - **Hotspots:**  
    - Avoid rebuilding on every `solve()` call.  
    - Use integer vars when possible for stability.  
    - If cloning the graph (see `cloned_graph_ptr`), document when/why (e.g., sign switching).

- `void solve() override`
  - **Purpose:** Run solver, store labels/partition and objective → sets `modularity_value` and base class outputs.
  - **Complexity:** Solver-dependent; worst-case exponential; practical time driven by heuristics and cuts.
  - **Pre:** `build()` done. Valid solver environment.
  - **Hotspots:** Timeouts, mip gap targets, warm starts from greedy/XY solution to stabilize and speed up.

- `double get_frustration_index() const override`
  - **Purpose:** Expose base objective value (for compatibility with XY base).
  - **Complexity:** O(1).

- `double get_modularity_value() const`
  - **Purpose:** Return modularity score derived from solved labels/partition.
  - **Complexity:** O(1).
  - **Hotspot:** Define scale/range (e.g., [-1,1]) and whether it’s normalized; keep consistent with publications.

- **Private:** `SignedGraphForMIP* cloned_graph_ptr = nullptr;`
  - **Purpose:** If building on a transformed graph (e.g., switched), hold ownership.
  - **Hotspot:** Ownership semantics; prefer `std::unique_ptr` to avoid manual delete in dtor.

---

## Cross-cutting hotspots & contracts

1. **Sign convention**
   - Throughout, sign is derived from `weights[eid]` using:  
     `sign = (w == 0.0 && signbit(w)) ? -1 : (w >= 0.0 ? 1 : -1)`.  
   - **Contract:** Producers must not emit `-0.0` unless truly negative; test with inputs having zeros.

2. **Indexing invariants**
   - `weights.size() == igraph_ecount(g)`; edge id `eid` stable across views.  
   - `edge_index()` must be rebuilt if edges are reordered/added.

3. **Thread safety**
   - `igraph` is not generally thread-safe for concurrent writes.  
   - Guard shared state when running greedy or MIP preprocessing in parallel.

4. **Performance**
   - Prefer delta updates for local search (`O(deg(u))` per move) over recomputing `O(m)`.  
   - Reserve capacity for adjacency/containers; avoid heap churn inside inner loops.  
   - Use `int64_t` for counters/accumulators to avoid overflow on large graphs.

5. **Determinism**
   - Tie-breaks and traversal order affect outputs. Document default `cmp_fn` and RNG seeds.

6. **I/O and large runs**
   - Avoid `operator<<` in tight loops.  
   - When writing SVG, gate by size thresholds to prevent huge files.

---

## FrustrationModelXY (frustration_model_xy.h)

**Role:** MIP formulation using XY variables: `x` (binary labels) and `y` (edge auxiliaries). Adds triangle/cycle user cuts; may seed heuristics via switching.

### Lifecycle
- `explicit FrustrationModelXY(SignedGraphForMIP& g, int cut_flags = 0)`
  - **Purpose:** Initialize with graph and cut configuration.
  - **Pre:** `g` alive for object lifetime.
  - **Complexity:** O(1).

- `void build() override`
  - **Purpose:** Create variables/constraints/objective in CPLEX.
  - **Side effects:** Initializes `x` (size `n`) and `y` (size `m`), registers callbacks if enabled.
  - **Complexity:** O(n + m) model build; memory proportional to created vars/constraints.
  - **Failure modes:** CPLEX env/alloc errors; invalid cut flags.

- `void solve() override`
  - **Purpose:** Optimize model; may attach callbacks/heuristics.
  - **Complexity:** Solver-dependent (worst-case exponential). Practical dominated by cut density and heuristics.
  - **Pre:** `build()` done.

- `void export_solution(const std::string& file_prefix, bool with_svg) const`
  - **Purpose:** Persist solution artifacts (labels, y-values, metrics). Optional SVG.
  - **Complexity:** O(n + m) IO.

### Variables
- `std::vector<IloBoolVar> x` — node labels (binary).
- `std::vector<IloNumVar>  y` — edge auxiliaries (continuous/binary per formulation).
  - **Invariant:** `x.size()==n`, `y.size()==m`.

### Custom cuts
- `IloRange generate_cycle_cut_standard(IloEnv&, const std::vector<Edge>& all_edges)`
  - **Purpose:** Single cycle inequality (standard form).
  - **Complexity:** O(|cycle|) to build one cut.
- `std::vector<std::pair<IloRange, std::string>> generate_cycle_cuts(...) override`
  - **Purpose:** Batch generation for detected negative cycles.
  - **Complexity:** Sum O(|C|) over cycles in batch.
- `std::vector<std::pair<IloRange, std::string>> generate_pending_triangle_cuts(...) override`
  - **Purpose:** Add triangle inequalities not yet posted.
  - **Complexity:** O(#triangles posted this round).
- `std::vector<std::pair<IloRange, std::string>> generate_positive_triangle_cuts(...)`
  - **Purpose:** Extra tightening for positive triangles.
  - **Note:** Guard against duplication; track posted IDs/tags.

### Callbacks (CPLEX)
- `NegativeTriangleCutGenerator : IloCplex::UserCutCallbackI`
- `NegativeCycleCutGenerator   : IloCplex::UserCutCallbackI`
- `ConditionalCycleCutGenerator: NegativeCycleCutGenerator`
  - **Behavior:** Skips cycle cuts if `triangle_cut_added_last_round==false`.
- `SwitchingHeuristicCallback  : IloCplex::HeuristicCallbackI`
  - **Purpose:** Inject a feasible solution via greedy switching to warm start.
  - **Hotspot:** Avoid frequent/expensive rebuilds during callback; cache graph stats.

**Callback contracts & pitfalls**
- `duplicateCallback()` must allocate via `new (getEnv()) ...` (already done).
- Threading: CPLEX may run callbacks in parallel contexts—use only thread-local or callback-safe state. No shared mutation without sync.
- Performance: Bound number of cuts per round; track whether a triangle cut was actually added to avoid fruitless rounds (`triangle_cut_added_last_round` already present).

**Testing ideas**
- Tiny graphs with known negative cycles/triangles: check that cuts are generated once and are valid (don’t cut off known feasible integral solutions).
- Heuristic injection: ensure objective improves and solution is accepted (check incumbent logs).

---

## SignedGraphForMIP (signed_graph_mip.h)

**Role:** Graph adapter for MIP: maintains fractional/derived weights to guide cuts and heuristics; exposes negative-cycle stream.

### Data members
- `edge_to_eid : unordered_map<Edge, igraph_integer_t>` — O(1) avg edge lookup.
- `frac_weights : vector<double>` — last fractional `y` (or related) per edge.
- `mask_weights : vector<double>` — `|frac - 0.5|` (edge-aligned). **Range:** [0, 0.5].
- `salience_full_ : vector<double>` — `1 - min(1, 2*|frac-0.5|)` in [0,1]; peaks at 1 near 0.5.
  - **Invariant:** All vectors length `m`; recomputed together.

**Accessors**
- `const std::vector<double>& get_mask_weights() const`
- `const std::vector<double>& edge_salience_view() const`
  - **Purpose:** Feature for heuristics/kicks; returns `salience_full_` if non-empty; else `mask_weights`.
  - **Complexity:** O(1).

**Lifecycle**
- Ctors:
  - `(const std::string& file_path)` — build from file.
  - `(const SignedGraph* other)` / `(const SignedGraph* other, std::vector<double> new_weights)`
  - **Complexity:** O(n + m) each.
  - **Ownership:** Derives from `SignedGraph`; dtor tears down added state.

**Key methods**
- `bool weighting_from_fractional(const std::vector<double>& x, const std::vector<double>& y)`
  - **Purpose:** Ingest fractional solution and update `frac_weights`, `mask_weights`, `salience_full_`.
  - **Complexity:** O(m).
  - **Pre:** `y.size()==m`.  
  - **Return:** true if accepted/updated, false if inconsistent.
  - **Hotspot:** Numeric noise near 0.5; consider epsilon when deriving masks/salience.

- `const std::vector<int> greedy_switching()`
  - **Purpose:** Heuristic partition on current (possibly fractional-informed) weights.
  - **Complexity:** ~O(m) per pass.

- `std::optional<std::shared_ptr<const std::vector<int>>> fractional_greedy_switching(const SignedGraph::GreedyKickOptions& opts)`
- `std::optional<std::shared_ptr<const std::vector<int>>> fractional_greedy_switching()`
  - **Purpose:** Optional heuristic using fractional-guided salience and kick options.
  - **Complexity:** Baseline O(m) per pass; capped by options.
  - **Return:** `nullopt` if unsuitable/no improvement.

- **Negative cycle discovery (public wrappers)**
  - `std::vector<NegativeCycle> find_switched_lower_bound(bool cover = false)`
  - `std::vector<std::vector<NegativeCycle>> find_switched_lower_bound_grouped(bool cover) const`
  - **Purpose:** Generate sets of negative cycles for cuts or lower bounds.
  - **Complexity:** Depends on stream settings; typically near-linear per batch.

- **Stream factory**
  - `NegativeCycleBatchStream open_negative_cycle_stream(bool cover, bool use_triangle_order = false) const`
  - **Purpose:** On-demand, batched cycle enumeration with optional triangle-aware ordering.

**Pitfalls & contracts**
- Keep `edge_to_eid` synchronized with `g`; rebuild if edges are reordered.
- When using fractional data, define a stable epsilon range around 0.5 to avoid flip-flopping salience across iterations.

**Tests**
- Inject synthetic `y` around 0.5; verify `mask_weights` and `salience_full_` match formulas and bounds.
- Verify `edge_salience_view()` fallback works (empty `salience_full_` → returns `mask_weights`).

---

## NegativeCycleBatchStream (signed_graph_mip.h)

**Role:** Incremental emitter of negative cycles in **batches**. Supports coverage mode and triangle-aware ordering; blends LP scores to prioritize edges.

### Public API
- `explicit NegativeCycleBatchStream(const SignedGraphForMIP& G, bool cover, bool use_triangle_order = false)`
  - **Purpose:** Initialize stream over `G`.
  - **Complexity:** O(n + m) to prepare caches/mappings (see internal state).

- `bool next(std::vector<NegativeCycle>& out)`
  - **Purpose:** Produce next batch of cycles into `out`.
  - **Return:** `false` when no more cycles; `true` otherwise.
  - **Complexity:** Amortized near-linear in size of emitted batch; depends on strategy/caps.

- Stats:
  - `int total_cycles_emitted() const`
  - `int batches_emitted() const`
  - **Complexity:** O(1).

- `void set_lp_scores_full_edges(const std::vector<double>& s, double alpha=1.0, double beta=0.3)`
  - **Purpose:** Blend LP scores into edge priorities for selection.
  - **Pre:** `s.size()==m`; `alpha,beta` within documented ranges.
  - **Complexity:** O(m).

### Internal state (performance-relevant)
- Graph sizes: `vcount_`, `ecount_`.
- Per-edge weights: `base_pos_` (selection bias), `neg_hard_` (neg edges), `reuse_accum_` (cross-batch penalty).
- Negative edge list: `neg_edges_`; disconnected nodes: `disconnected_`.
- Saved weights: `saved_weights_` (+ `_pos_` for positive-only subgraph); flags to avoid repeated copies.
- Positive-only graph: `g_pos_` with `full2pos_eid_`, `pos2full_eid_`.
- Triangle counts & caps: `neg_tri_vert_`, `tri_used_per_vertex_`, `tri_cap_per_vertex_` / `tri_cap_per_v_`.
- Degrees: `neg_deg_`, `pos_deg_`.
- Coverage knobs: `K_tri_per_neg_`, `overlap_penalty_gamma_`, `cross_batch_penalty_scale_`.
- LP blending: `lp_score_full_`, `use_lp_weights_`, `lp_alpha_`, `lp_beta_`.

### Complexity notes
- Initial build (`build_initial_state_()`): O(n + m).
- Each batch:
  - Build/refresh masks (`build_mask_for_batch_()`): O(m) or O(m_pos).
  - Triangle-aware ordering: computing/refreshing counts ~O(∑ min(deg(u),deg(v))) for affected edges; capped by `tri_cap_per_vertex_`.
- Overall: Typically **near-linear per batch** with caps preventing quadratic behavior on hubs.

### Hotspots & pitfalls
- **Memory churn:** Reuse vectors; avoid reallocations across batches.
- **Edge ID mappings:** Keep `full2pos_eid_`/`pos2full_eid_` correct; negative edges must map to `-1` in `full2pos_eid_`.
- **Saved weights flags:** Ensure RAII-like restore to prevent leaking modified weights back to `G`.
- **Reentrancy/Thread safety:** Stream is **not** thread-safe; one consumer at a time.
- **Scoring stability:** Penalization (`overlap_penalty_gamma_`, `cross_batch_penalty_scale_`) should be monotone to avoid cycling.
- **Triangle caps:** `tri_cap_per_vertex_` small by default; document impact on coverage/optimality of cuts.

### Tests
- Small graphs with known negative cycles: assert batches are disjoint/covering depending on `cover`.
- With LP scores: monotone relationship—edges with higher LP weight should be prioritized (within randomness caps).
- Positive-only graph build: check counts and mapping invariants.

---

## Cross-cutting recommendations (docs only)
- **Ownership:** Prefer smart pointers (`std::unique_ptr`) for `cloned_graph_ptr` to clarify lifetime; avoid manual `delete` in dtors.
- **Determinism:** Document randomness/tie-breaks and provide a seed. Make callbacks respect determinism (or clearly state they don’t).
- **Numerics:** Use consistent tolerance for “fractional 0/1” decisions (e.g., `1e-6` / `1e-8`) across salience, masks, and cut separation.
- **Big-O guardrails:** Keep per-iteration work to O(m) with explicit caps; never scan O(n²).

