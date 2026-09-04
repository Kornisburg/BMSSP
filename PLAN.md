# BMSSP — Implementation Plan

Implements the Duan–Mao–Mao–Shu–Yin algorithm *"Breaking the Sorting Barrier for Directed
Single-Source Shortest Paths"* (arXiv:2504.17033, STOC 2025) in Rust, driven by Test-Driven
Development and a custom (zero-extra-dep) benchmark harness.

Source spec: `spec.md`. Paper: <https://arxiv.org/abs/2504.17033>.

## Locked-in decisions

- **Write from scratch.** Existing Rust ports (`quantumrag/BMSSP`, `danalec/DMMSY-SSSP-rs`,
  `lucas-montes/bmssp`, `johannhipp/bmssp`) and the honest-benchmark study
  `primaryaesthetics/logtwothirds` are consulted only to disambiguate paper pseudocode.
- **Goal: balanced.** Correctness first; then performance as far as feasible; report honestly
  either way. Expectation set by `logtwothirds`: tuned Dijkstra likely beats BMSSP at every
  runnable size (extrapolated crossover ~ n = 2^400000). We will measure, not assert.
- **Custom minimal bench harness.** No criterion/proptest (not in cargo cache); hand-rolled
  timing + operation counters + `rand`/`rand_chacha` (cached) for deterministic graphs.
- **First session: Phases 0–1.** Scaffold + Dijkstra baseline + correct-but-slow BMSSP
  verified vs Dijkstra on n ≤ 5000 across all graph families.

## Project layout

```
Cargo.toml            edition 2021; deps: rand 0.8, rand_chacha 0.3
src/lib.rs            module re-exports
src/graph.rs          CSR Graph {offsets,to,weight} + generators + WeightDist
src/dijkstra.rs       binary-heap Dijkstra baseline (MinState, INF) + counters
src/params.rs         k = floor((log2 n)^(1/3)), t = floor((log2 n)^(2/3)), l = ceil(log2 n / t)
src/counters.rs       operation counters (relaxations, heap ops, recursion, queue ops, pivots)
src/queue.rs          partial-order queue (BTreeMap v0: Insert/Pull/BatchPrepend)
src/bmssp.rs          FindPivots / BaseCase / BMSSP recursion / barrier_breaker_sssp driver
src/bin/bench_sssp.rs custom bench CLI: median wall-clock + counters + markdown table
benches/sssp.rs       `cargo bench` target (harness = false) reusing bench core
tests/                oracle + property + handcrafted tests
PLAN.md, PROGRESS.md, BENCHMARKS.md, ALGORITHM.md, AUDIT.md   (living docs)
```


## Algorithm parameters (paper)

- `k = floor((log2 n)^(1/3))`, `t = floor((log2 n)^(2/3))`, `l = ceil(log2 n / t)`, each `>= 1`.
- Base-2 logs guarantee the top-level loop bound `k * 2^(l*t) >= n` (all reachable vertices
  settle at the top level => successful execution, B' = infinity).
- Top level: `BMSSP(l, B = inf, S = {src})` with `dhat[src] = 0`.

## Recursion semantics (faithful to paper, as disambiguated by reference ports)

- `FindPivots(B, S)`: up to `k` rounds of bounded Bellman-Ford relaxation from seeds S (only
  vertices with `dhat < B`), building predecessor forest F (parent of first best improver).
  `W` = visited set; early stop when `|W| > k|S|` => `P = S`; else `P` = roots in S whose
  F-subtree has `>= k` vertices. Returns `(P, W)`.
- `BaseCase(B, S)`: bounded **multi-source** mini-Dijkstra from all seeds with `dhat < B`,
  settles up to `k+1` vertices. If `|U0| <= k`: return `(B, U0)`. Else return
  `(B' = max dhat in U0, U = {v in U0 : dhat[v] < B'})`.
- `BMSSP(l, B, S)`: l=0 => BaseCase. Else: FindPivots; seed partial queue D with pivots;
  loop `while !D.is_empty()`: Pull -> `(S_i, B_i)` (whole smallest bucket + separation bound
  = next key or B); recurse `BMSSP(l-1, B_i, S_i)` -> `(B'_i, U_i)`; extend U; re-prepend
  S_i leftovers with `dhat in [B'_i, B_i)`; relax from U_i; for cand: if `cand in [B_i, B)`
  Insert into D, elif `cand in [B'_i, B_i)` collect into K; BatchPrepend K. After loop,
  sweep W: `U += {x in W : dhat[x] < B'}`. Return `(B, U)` (always runs to D-empty in the
  Phase-1 reference variant; no `k*2^(l*t)` cap yet — ablation later).
- Relaxation uses `<=` (equal-key relax is legal; the interval logic prevents re-inserting
  already-settled vertices), heap/base-case pushes use strict `<` (avoids zero-weight churn).
- `EPS = 1e-12` for interval membership; final distances are exact shortest-path sums
  (dist updates use exact comparison), so integer-weight tests can match Dijkstra exactly.

## TDD loop (each step red -> green -> refactor)

1. Scaffold crate + stubs. `cargo build` clean, dijkstra tests red.
2. Implement dijkstra; tests vs Floyd–Warshall oracle on small graphs (green).
3. Implement queue; model-based randomized test (green).
4. Implement bmssp stubs with failing property tests (red).
5. Implement BaseCase -> FindPivots -> BMSSP -> driver, iterate until property tests green
   on n <= 5000 across all families (including ties / zero weights / unreachable).
6. Bench harness; record tables in BENCHMARKS.md; update PROGRESS.md.

## Roadmap (later sessions)

- **Phase 2** — Constant-degree transform (vertex splitting, directed), re-verify vs Dijkstra.
- **Phase 3** — Lemma 3.3 block-based queue (D0/D1 block lists + BST of block upper bounds),
  differentially tested vs BTreeMap model before drop-in.
- **Phase 4** — Engineering polish: CSR locality, arena/pre-sized buffers, iterative stack
  recursion, optional rayon (note: recursion is sequential over shared dhat).
- **Phase 5** — Scaling to n = 10^4..10^7; ablation matrix (no-pivots, k/t schedules, BTreeMap
  vs block queue, +-transform, recursion vs iterative); perf counters; asymptotics
  t/(m log^(2/3) n) vs t/(m + n log n).
- **Phase 6** — Reproducibility docs: ALGORITHM.md, AUDIT.md, BENCHMARKS.md (auto-updated),
  pinned seeds, exact repro commands.

## Guardrails

- Correctness is non-negotiable: distances must match Dijkstra within 1e-9 (or exactly for
  integer weights); every change re-runs the full test suite.
- Report honestly where BMSSP wins/loses; never tune results to hide a loss.
- Deterministic: pinned seeds (SplitMix64/ChaCha8), recorded in every report.
