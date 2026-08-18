# BENCHMARKS

Auto-oriented living table. Re-run:

```bash
BMSSP_BENCH_ITERS=1 BMSSP_PARTIAL=1 BMSSP_QUEUE=block \
  cargo run --release --bin bench_sssp
# optional: BMSSP_NO_PIVOTS=1, BMSSP_QUEUE=map, BMSSP_SCALE=1 (larger n)
```

Pinned generator seeds live in `src/bin/bench_sssp.rs` (`0xB0555EED`, …).

## Session 5 snapshot (Phase 4 polish)

Machine: local macOS release (`lto=thin`). Ablation: **pivots=on, Block queue, partial=on**, `BMSSP_BENCH_ITERS=1`.

| family | n | m | dijk(ms) | bmssp(ms) | speedup | verified |
|---|---|---|---|---|---|---|
| er_c2 | 10000 | 20000 | 0.648 | 6.136 | 0.11 | true |
| er_c4 | 10000 | 40000 | 1.163 | 8.703 | 0.13 | true |
| er_c8 | 10000 | 80000 | 1.598 | 10.090 | 0.16 | true |
| er_c4_1e5 | 100000 | 400000 | 14.281 | 143.581 | 0.10 | true |
| grid_100 | 10000 | 19800 | 0.673 | 6.011 | 0.11 | true |
| grid_316 | 99856 | 199080 | 8.068 | 61.010 | 0.13 | true |
| pl_c2 | 10000 | 20000 | 0.001 | 0.013 | 0.09 | true |
| pl_c4_1e5 | 100000 | 400000 | 0.013 | 0.063 | 0.20 | true |
| layered | 100000 | 398000 | 11.764 | 85.582 | 0.14 | true |
| er_c4_real | 10000 | 40000 | 1.100 | 12.885 | 0.09 | true |

Constant-degree transform rows (suffix `_tr`): all `verified=true`; BMSSP still behind Dijkstra (speedup 0.06–0.57). Worst practical blow-up remains high-degree → aux-tree expansion (`layered_tr` n≈697k).

### vs Session 3

Same ablation shape. `layered` 190ms → **86ms**; `er_c4_1e5` 370ms → **144ms**; `grid_316` 150ms → **61ms**. Depth-stamped U + arena buffers + CSR adjacency sort cut constants without changing the “Dijkstra still wins” conclusion.

## Honest takeaway

At every runnable size here, tuned binary-heap Dijkstra is ~6–12× faster. The asymptotic crossover claimed by theory is not visible yet; matching `logtwothirds` / public ports. Correctness is solid (`verified=true` on the full matrix).
