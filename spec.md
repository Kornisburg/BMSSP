**Implementing the Duan–Mao et al. (BMSSP / “Breaking the Sorting Barrier”) algorithm in Rust via agentic coding, with a solid Dijkstra baseline and a rigorous comparison methodology.**

The paper gives a deterministic \(O(m \log^{2/3} n)\) SSSP for directed graphs with non-negative real weights in the comparison-addition model. A production-quality version that *actually* beats highly-tuned Dijkstra on realistic sparse graphs is a non-trivial engineering project (most public implementations still use ordinary heaps and are slower or only win on very large instances). Agentic coding is an excellent fit because the algorithm has clear modular pieces (FindPivots, BaseCase, BMSSP recursion, partial-order data structure) that can be built, tested, and optimized iteratively.

### 1. Baseline: Standard Dijkstra

Start with a clean, fast reference implementation. Use adjacency lists (or CSR for better locality) and `BinaryHeap`.

```rust
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(Copy, Clone, Eq, PartialEq)]
struct State {
    cost: f64,
    node: usize,
}

impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        // min-heap via reverse
        other.cost.partial_cmp(&self.cost).unwrap_or(Ordering::Equal)
            .then_with(|| self.node.cmp(&other.node))
    }
}
impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub fn dijkstra(adj: &[Vec<(usize, f64)>], src: usize) -> Vec<f64> {
    let n = adj.len();
    let mut dist = vec![f64::INFINITY; n];
    dist[src] = 0.0;
    let mut heap = BinaryHeap::new();
    heap.push(State { cost: 0.0, node: src });

    while let Some(State { cost, node }) = heap.pop() {
        if cost > dist[node] { continue; }
        for &(v, w) in &adj[node] {
            let next = cost + w;
            if next < dist[v] {
                dist[v] = next;
                heap.push(State { cost: next, node: v });
            }
        }
    }
    dist
}
```

Optional stronger baseline: `petgraph::algo::dijkstra` or a dial/radix-heap variant for integer weights. Measure both wall-clock time and number of decrease-key / extract-min operations.

### 2. Agentic Coding Setup (Harness + Guidance + Loop)

**Harness (the test/bench infrastructure that the agent lives inside)**

```toml
# Cargo.toml
[package]
name = "bmssp-rs"
edition = "2021"

[dependencies]
# optional: petgraph, rand, rayon

[dev-dependencies]
criterion = { version = "0.5", features = ["html_reports"] }
proptest = "1.4"
rand = "0.8"
# for real graphs later: snap, or download DIMACS/SNAP datasets

[[bench]]
name = "sssp"
harness = false
```

Essential harness components:
- Unit + integration tests that assert distances match Dijkstra (within `1e-9` or exact for rationals).
- Property-based tests (`proptest`) generating random directed graphs (Erdős–Rényi sparse, preferential attachment, grids).
- Criterion benchmarks with multiple graph families and sizes.
- Counters: number of relaxations, heap operations, recursive calls, pivot selections.
- Optional: `flamegraph`, `cargo-instruments`, `heaptrack`, or Linux `perf` for cache misses / branch prediction.
- Graph generators + loaders for SNAP / DIMACS / synthetic sparse graphs (`m ≈ 2–10 n`).

**Guidance (system / developer prompts you feed the agent)**

Provide the paper’s pseudocode (FindPivots, BaseCase, BMSSP) almost verbatim, plus invariants:
- “Every incomplete vertex with \(d(v) < B\) has a shortest path that passes through a complete vertex in the current \(S\).”
- Parameters: \(k = \lfloor\log^{1/3} n\rfloor\), \(t = \lfloor\log^{2/3} n\rfloor\), \(l = \lceil\log n / t\rceil\).
- Data-structure contract (Insert / BatchPrepend / Pull of the \(M\) smallest).
- “First make it correct and match Dijkstra exactly; only then optimise. Never sacrifice correctness for speed.”
- Explicit phases and success criteria (see loop below).

Keep a living `SPEC.md` and `PROGRESS.md` that the agent updates.

**The agentic loop (typical multi-day / multi-session workflow)**

1. **Scaffold** – Agent creates the crate, graph type, Dijkstra baseline, and empty BMSSP stubs. Harness already runs and shows “all tests fail / benches exist”.
2. **Correctness core** – Implement BaseCase (limited Dijkstra) + FindPivots. Test on tiny graphs. Agent iterates until distances match.
3. **Recursion** – Wire BMSSP recursive structure (still using a normal `BinaryHeap` or `BTreeMap` for the frontier). Keep matching Dijkstra on graphs up to a few thousand vertices.
4. **Data structure** – Replace the frontier with a closer approximation of Lemma 3.3 (blocks + BST of upper bounds, or at least a structure that supports cheap BatchPrepend of smaller keys). Re-verify correctness.
5. **Engineering polish** – Degree reduction (vertex splitting to constant degree), arena allocation / pre-sized vectors to avoid recursion allocation, cache-friendly CSR layout, careful float handling / tie-breaking, optional parallelisation of independent subproblems.
6. **Measurement & ablation** – Agent runs the full criterion suite after every meaningful change, records tables, and decides the next experiment (e.g. “remove pivots → how much slower?”, “different \(k,t\) schedules”, “true block DS vs BTreeMap”).
7. **Scaling & profiling** – Move to larger graphs (\(n = 10^5 \dots 10^7\)), real networks, and low-level counters. Agent proposes and applies micro-optimisations guided by profiles.
8. **Documentation & reproducibility** – Final report with tables, plots, and exact reproduction commands.

Use any coding agent (Claude Code, Cursor, Aider, Continue, custom harness with tool calling, etc.). The key is the closed loop: *generate → compile/test/bench → observe metrics → refine*, with the harness enforcing correctness at every step.

### 3. High-level Structure of the New Algorithm in Rust

Mirror the paper:

```rust
struct BmsspState {
    dhat: Vec<f64>,
    // optional: parent / forest for FindPivots
    // the partial-order frontier structure
}

fn find_pivots(/* B, S, &mut dhat, k */) -> (Vec<usize>, Vec<usize>); // P, W
fn base_case(/* B, S, &mut dhat, k */) -> (f64, Vec<usize>);
fn bmssp(/* l, B, S, &mut state, k, t */) -> (f64, Vec<usize>);

pub fn barrier_breaker_sssp(adj: &[Vec<(usize, f64)>], src: usize) -> Vec<f64> {
    // compute k, t, l from n
    // initialise dhat[src] = 0
    // call bmssp(l, f64::INFINITY, vec![src], ...)
    // return dhat
}
```

Start with a simple frontier (`BinaryHeap` or `BTreeMap<f64, Vec<usize>>`). Later replace it with a block-based structure that supports the three operations with the claimed amortised costs. Several public Rust crates already exist (search “bmssp”, “DunMaoSSSP”, “duan_sssp”, “DMMSY”)—use them as reference or starting points, then improve the data structure and constant factors.

### 4. How to Compare Them – Full Engineering Exploration

**Correctness (non-negotiable)**
- Exact distance vectors must match Dijkstra (or differ only on unreachable nodes / within floating-point tolerance).
- Property tests + a suite of hand-crafted graphs (chains, dense clusters, grids, graphs with many equal weights).

**Performance metrics**
- Wall-clock time (Criterion, multiple runs, warm-up).
- Operation counts: extract-min / insert / decrease-key / edge relaxations / recursive calls / pivots selected.
- Memory high-water mark and allocation volume.
- Hardware counters (cache misses, branch misses) via `perf` or Instruments.
- Asymptotic behaviour: plot time / \((m \log^{2/3} n)\) vs time / \((m + n \log n)\) across increasing \(n\).

**Graph families (essential for honest comparison)**
- Synthetic sparse: \(m = c\cdot n\) for \(c = 2,4,8,16\).
- Scale-free / power-law.
- Road networks / social graphs (SNAP, DIMACS).
- Grids and layered graphs (where the theory is most favourable).
- Both small (\(n\le 10^4\), where Dijkstra wins) and large (\(n\ge 10^6\)).

**Ablation & sensitivity**
- With vs without FindPivots.
- Different parameter schedules for \(k,t\).
- Ordinary heap vs better partial-order structure.
- With vs without constant-degree reduction.
- Recursion vs iterative version (stack simulation).
- Single-threaded vs limited parallelism of independent BMSSP subproblems.

**Expected practical picture (from existing public implementations)**
- Correctness is achievable.
- On modest sizes Dijkstra (especially a well-tuned one) is usually faster because of lower constants and better locality.
- On sufficiently large sparse graphs, and once the frontier data structure is closer to the paper’s design, the new algorithm can pull ahead (some Rust ports already report substantial speed-ups on particular instances such as LiveJournal-scale graphs). The engineering goal is to close the constant-factor gap and make the asymptotic advantage visible earlier.

**Reporting template**
Keep a markdown table (auto-updated by the agent) with columns: graph, \(n\), \(m\), Dijkstra ms, BMSSP ms, speed-up, #heap ops Dijkstra, #heap ops BMSSP, #relaxations, notes.

### Practical Starting Advice

1. Clone or study an existing Rust port for the skeleton, then rebuild the critical data structure yourself under the agentic loop.
2. Get correctness first on \(n\le 5\,000\); only then scale.
3. Instrument *everything*—operation counts are more informative than wall time while you are still iterating.
4. Treat the block-based partial-order queue (Lemma 3.3) as a separate mini-project; a correct but slower `BTreeMap` version already lets you explore the rest of the algorithm.
5. Reproducibility: pin seeds, record exact graph generation parameters, and keep Criterion HTML reports + raw counters.

This workflow turns a complex theoretical algorithm into a systematic engineering exploration: the harness keeps you honest, the guidance keeps the agent aligned with the paper, and the measurement loop tells you exactly where the remaining constant-factor work lies. You will end up with both a working implementation and a clear understanding of *when* (and *why*) the new algorithm beats classic Dijkstra in practice.
