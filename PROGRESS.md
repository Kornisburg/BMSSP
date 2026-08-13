# PROGRESS — living log

## Session 1 (2026-08-05) — Phases 0–1

Goal: full scaffold + Dijkstra baseline + correct (slow) BMSSP, verified vs Dijkstra on
n <= 5000 across all graph families, with custom bench harness.

### Status

- [x] PLAN.md written
- [x] Crate scaffold (Cargo.toml, modules, stubs)
- [x] graph.rs CSR + generators
- [x] params.rs / counters.rs
- [x] dijkstra.rs baseline + oracle tests
- [x] queue.rs BTreeMap PartialQueue + model test
- [x] bmssp.rs (BaseCase, FindPivots, BMSSP, driver) — CORRECT
- [x] Property tests bmssp vs dijkstra
- [x] Bench harness
- [x] Results recorded below

### Bench snapshot (release, `BMSSP_BENCH_ITERS=3`, n=10000 default / 100000 for the big rows)

| family | n | m | dijk(ms) | bmssp(ms) | speedup | verified |
|---|---|---|---|---|---|---|
| er_c2 | 10000 | 20000 | 1.066 | 3.178 | 0.34 | true |
| er_c4 | 10000 | 40000 | 1.901 | 3.834 | 0.50 | true |
| er_c8 | 10000 | 80000 | 2.783 | 7.312 | 0.38 | true |
| er_c4_1e5 | 100000 | 400000 | 42.821 | 105.101 | 0.41 | true |
| grid_100 | 10000 | 19800 | 1.145 | 7.041 | 0.16 | true |
| grid_316 | 99856 | 199080 | 15.620 | 217.226 | 0.07 | true |
| pl_c2 | 10000 | 20000 | 0.002 | 0.010 | 0.17 | true |
| layered | 100000 | 398000 | 22.146 | 1760.904 | 0.01 | true |
| er_c4_real | 10000 | 40000 | 1.948 | 20.735 | 0.09 | true |

As expected for a Phase-1 reference (no partial-execution cap, BTreeMap queue, full recursion),
BMSSP is 2–100x slower than binary-heap Dijkstra; worst case is `layered` (long chains defeat the
interval structure). This matches the honest-benchmark literature (BMSSP loses to tuned Dijkstra at
every runnable size). Phase-2 work targets the asymptotic upgrades and a block-based queue.

### Key findings (updated as we go)

1. **Algorithm design, phase-1 reference variant.** Implemented BMSSP(l, B, S) with the
   partial-order queue (Insert / Pull / BatchPrepend), FindPivots (bounded Bellman-Ford
   sweep + shortest-path-forest pivot rule), BaseCase (bounded multi-source mini-Dijkstra
   settling k+1 vertices), and the two-phase interval loop. Params k, t, l from `params(n)`.
   Runs to D-empty (B' = B always); no k·2^(lt) cap yet (Phase 2).

2. **Correctness bugs found & fixed (via exhaustive counterexample search + trace):**
   - `base_case` must return EVERY vertex it touches (settled AND still-in-heap), not just
     the k+1 settled ones — otherwise un-settled discovered vertices are lost.
   - The main loop must route discoveries on strict improvement into pq / batch. The first
     version checked `cand < dist[v]` AFTER relaxing, which was always false → nothing routed.
   - A level's relax loop must gate the distance update itself by `cand < b - EPS`. A child
     with bound bi=1 was "stealing" an improvement (cand=2, above its bound) without routing
     it, so the parent could never re-route it. Fix: above-bound discoveries are left for the
     parent to re-relax and route. (This replaces the reference's route-on-every-edge, which
     requires an iteration cap to terminate.)
   - Boundary/tie buckets: base_case returning all settled vertices (incl. boundary) avoids
     the infinite re-add loop and the tied-bucket loss; routing on strict improvement only
     avoids re-inserting completed vertices.
   - `find_pivots` marks/propagates vertices reached within the bound even on equal keys, so
     the W-sweep completes; pivots fall back to S on early-stop / no-qualifier.

3. **Verification.** `cargo test`: 9 lib + 7 property + 8 handcrafted = 24 green (vs
   binary-heap Dijkstra; Floyd-Warshall oracle on small). Exhaustive brute force over all
   simple digraphs n=2..4 (weights {0,1,2}, all sources) + 200k random each for n=5,6:
   no counterexample. Counterexamples that were found and fixed: n=3 triangle with a zero
   edge; n=6 graph where a child level stole an above-bound improvement.

4. **Debug aids.** `BMSSP_TRACE=1` env var prints relax events (R/B lines). Brute force in
   `examples/bruteforce.rs`; repro harness in `examples/dbg.rs`.

## Session 2 (2026-08-05) — Phase 2: constant-degree transform

### Status

- [x] `src/transform.rs`: `to_constant_degree` — splits every vertex of out-degree > 2 into a
      balanced binary tree of zero-weight aux nodes (root = original vertex). Result has max
      out-degree <= 2; original vertices keep ids `0..n`; distance-preserving. Size
      n' = n + sum_{deg>2}(2·deg − 2), m' ≈ 3m (O(n+m)).
- [x] `tests/transform.rs` (6 tests): max out-degree <= 2; distances preserved across all
      families (er/dense/grid/chain/power_law/layered, incl. zero weights); weighted-edge
      (head, weight) multiset preserved; BMSSP on transformed matches Dijkstra on original;
      exhaustive n=2..4 all simple digraphs. All green; full suite now 30 tests.
- [x] Bench harness extended with `bench_one_transformed` rows (suffix `_tr`).

### Findings

1. Transform works and is exact (verified=true everywhere, projected back to original ids).
2. Cost is real: high-degree hubs blow up the graph (~7x for er_c8: n 10k -> 150k; power-law
   ~70x: n 100k -> 700k). BMSSP-on-transformed is always slower than BMSSP-on-original and
   slower than Dijkstra on any variant (speedup 0.00–0.58). Worst case `layered_tr` (n=697k)
   is ~70s vs Dijkstra 66ms — the aux trees create long zero-weight chains that the interval
   structure cannot skip.
3. Conclusion (matches `logtwothirds`): the constant-degree transform is a theoretical
   device, not a practical speedup here; BMSSP remains behind tuned Dijkstra at every
   runnable size. Grid (`right`+`down`, out-degree <= 2) transforms to itself.

### Open / next

- [ ] **Review follow-ups (latent, not live)**
  - [ ] `base_case` marks seeds at init, so `u0` counts only *non-seed* settlements and the
        returned B' (`max dist over u0`) is not the paper's `max dhat in U0`. Harmless now
        (B' is unconsumed; `touched` still returns seeds), but must be settled before the
        Phase-3 D-empty/early-stop B' logic lands.
  - [ ] `find_pivots` forest parents are only reset for seeds; a W vertex with no incoming
        equality edge can keep a stale parent, and zero-weight ties can break the layer order
        for the reverse-order `sub` accumulation. Both only skew pivot *selection* (fallback
        `P = S` preserves correctness) — perf concern, revisit in the pivot-quality ablation.
- [ ] Phase 3: complexity-relevant upgrades (bound `k·2^(l·t)`, partial execution, D-empty
      vs early stop) + ablation matrix.
- [ ] Phase 3b: Lemma 3.3 block-based queue, differentially tested vs BTreeMap model.
- [ ] Phase 4: engineering polish (CSR locality, arena buffers, iterative recursion).
- [ ] Phase 5: scaling n = 10^4..10^7; ablation (no-pivots, k/t schedules, BTreeMap vs block,
      ±transform, recursion vs iterative); asymptotics t/(m log^(2/3) n) vs t/(m + n log n).
- [ ] Phase 6: reproducibility docs (ALGORITHM.md, AUDIT.md, auto BENCHMARKS.md).
