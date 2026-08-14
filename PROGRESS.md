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
- [x] **Phase 3: complexity-relevant upgrades (bound `k·2^(l·t)`, partial execution, D-empty
      vs early stop) + ablation matrix** — see Session 3.
- [x] **Phase 3b: Lemma 3.3 block-based queue, differentially tested vs BTreeMap model.**
- [ ] Phase 4: engineering polish (CSR locality, arena buffers, iterative recursion).
- [ ] Phase 5: scaling n = 10^4..10^7; ablation (no-pivots, k/t schedules, BTreeMap vs block,
      ±transform, recursion vs iterative); asymptotics t/(m log^(2/3) n) vs t/(m + n log n).
- [ ] Phase 6: reproducibility docs (ALGORITHM.md, AUDIT.md, auto BENCHMARKS.md).

## Session 3 (2026-08-14) — Phase 3 + 3b: partial execution, block queue, ablation

Goal: bound the per-call workload by `k·2^(l·t)` (halt + hand up unprocessed work), implement
the Lemma 3.3 block-based queue with real Pull, and expose the full ablation matrix in the bench.

### Status

- [x] `BmsspConfig { k, t, l, use_pivots, partial_execution, queue_impl }`; bench knobs
      `BMSSP_NO_PIVOTS`, `BMSSP_QUEUE=map|block`, `BMSSP_PARTIAL`, `BMSSP_BENCH_ITERS`.
- [x] Partial execution: halting at `|U| > k·2^(l·t)` (distinct vertices, epoch-deduped via
      dedicated `marked_u` array + `u_epoch` — the shared `marked` is clobbered by children's
      interior epochs); on halt the call drains its queue into U so the caller relaxes the
      dropped items' edges.
- [x] `drain()` on `PartialQueue`/`BlockQueue`/`QueueOps`; halt branch drains queue into u_total
      with a comment citing Remark 3.4.
- [x] Routing deviation kept (documented): every strict improvement is routed
      (`cand >= bi - EPS -> pq`, else `to_batch`), not the paper's interval skip; W-vertices
      discovered by FindPivots are inserted into the call's own queue (dist < b − EPS) instead
      of relying on the paper's ≤-equality re-insertion (Remark 3.4) — this is what keeps
      halting calls from losing them.
- [x] BlockQueue (Lemma 3.3): fixed-block `F = 2^((l-1)·t)`, per-block min-buckets, max-boundary
      split, gap placement; differential test `block_queue_matches_model` vs the BTreeMap model.
- [x] `tests/config_variants.rs` (8 variants: pivots ±, queue ×, partial ×) verified vs Dijkstra
      on all sized families.
- [x] Bench ablation wired: per-graph `from_n(g.n)` for k/t/l (the top level must stay
      successful: its bound `k·2^(l·t)` exceeds `|U| ≤ n`), env knobs for pivots/queue/partial.

### Bugs found & fixed (partial execution)

1. **Queue dropped on halt.** The halting call returned only `u_total` (its touched set), so
   discovered-but-unprocessed queue items vanished. Fixed by draining the queue into U.
2. **FindPivots-set dists were never queued.** `find_pivots` improves `dhat` for W vertices
   without inserting them; on halt the W-sweep excludes them (`dhat >= B'_i`), so their edges
   are never relaxed at the parent level. Case: `layered(20,15,4,5,Int{0,9})`, pivots=true,
   vertex 202 (path `0→…→172→186→202`, w=0 edge `186→202`): dist[186]=21 set by the pivot
   sweep, level-2 call halts with `B'_2=18`, vertex 202 stays unreachable. Fixed by queuing the
   W-vertices (dist < b − EPS) directly — verified via `config_variants` (3/3 green).

### Bench snapshot (release, `BMSSP_BENCH_ITERS=1`, partial execution ON, Block queue, pivots ON)

| family | n | m | dijk(ms) | bmssp(ms) | speedup | verified |
|---|---|---|---|---|---|---|
| er_c2 | 10000 | 20000 | 0.986 | 9.128 | 0.11 | true |
| er_c4 | 10000 | 40000 | 1.824 | 14.935 | 0.12 | true |
| er_c8 | 10000 | 80000 | 2.900 | 18.343 | 0.16 | true |
| er_c4_1e5 | 100000 | 400000 | 41.211 | 369.900 | 0.11 | true |
| grid_100 | 10000 | 19800 | 1.310 | 11.289 | 0.12 | true |
| grid_316 | 99856 | 199080 | 13.557 | 149.582 | 0.09 | true |
| pl_c2 | 10000 | 20000 | 0.002 | 0.019 | 0.09 | true |
| pl_c4_1e5 | 100000 | 400000 | 0.029 | 0.395 | 0.07 | true |
| layered | 100000 | 398000 | 19.989 | 190.546 | 0.10 | true |
| er_c4_real | 10000 | 40000 | 1.753 | 23.015 | 0.08 | true |

### Key findings

1. Every ablation variant is exact (all `verified=true`, incl. `_tr` rows) on the full bench
   set. Partial execution + Block queue cut `layered` 1.76s -> 0.19s (Phase-1 snapshot) while
   staying exact — the interval structure + halt bound now cap the pathological chains.
2. BMSSP remains 6–15x slower than tuned Dijkstra at every runnable size; the block queue and
   partial execution shrink the constant but don't change the asymptotic picture — consistent
   with the honest-benchmark literature.

### Open / next

- [ ] Phase 4: engineering polish (CSR locality, arena buffers, iterative recursion).
- [ ] Phase 5: scaling n = 10^4..10^7; ablation (no-pivots, k/t schedules, BTreeMap vs block,
      ±transform, recursion vs iterative); asymptotics t/(m log^(2/3) n) vs t/(m + n log n).
- [ ] Phase 6: reproducibility docs (ALGORITHM.md, AUDIT.md, auto BENCHMARKS.md).
