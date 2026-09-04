use std::time::Instant;

use bmssp_rs::bmssp::{BmsspConfig, BmsspEngine, QueueKind};
use bmssp_rs::counters::Counters;
use bmssp_rs::dijkstra::dijkstra;
use bmssp_rs::graph::{Graph, WeightDist};
use bmssp_rs::transform::to_constant_degree;

fn close(a: &[f64], b: &[f64]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    for (x, y) in a.iter().zip(b) {
        if x.is_infinite() && y.is_infinite() {
            continue;
        }
        if x.is_infinite() || y.is_infinite() {
            return false;
        }
        let scale = 1.0 + x.abs().max(y.abs());
        if (x - y).abs() > 1e-9 * scale {
            return false;
        }
    }
    true
}

pub struct Row {
    pub family: String,
    pub n: usize,
    pub m: usize,
    pub dijk_ms: f64,
    pub bmssp_ms: f64,
    pub verified: bool,
    pub d_extract: u64,
    pub d_relax: u64,
    pub b_relax: u64,
    pub b_pulls: u64,
    pub b_recursive: u64,
    pub b_pivots: u64,
}

pub fn bench_one(g: &Graph, src: u32, family: &str, iters: usize, ab: &Ablation) -> Row {
    let cfg = BmsspConfig {
        use_pivots: ab.use_pivots,
        queue_impl: ab.queue_impl,
        partial_execution: ab.partial_execution,
        ..BmsspConfig::from_n(g.n)
    };
    let mut dc = Counters::new();
    let d = dijkstra(g, src, &mut dc);
    let mut bc = Counters::new();
    let b = BmsspEngine::new(g, cfg.clone(), &mut bc).run(src);
    let verified = close(&d, &b);

    let mut dijk_ms = f64::MAX;
    for _ in 0..iters {
        let mut c = Counters::new();
        let t = Instant::now();
        dijkstra(g, src, &mut c);
        dijk_ms = dijk_ms.min(t.elapsed().as_secs_f64() * 1e3);
    }
    let mut bmssp_ms = f64::MAX;
    for _ in 0..iters {
        let mut c = Counters::new();
        let t = Instant::now();
        BmsspEngine::new(g, cfg.clone(), &mut c).run(src);
        bmssp_ms = bmssp_ms.min(t.elapsed().as_secs_f64() * 1e3);
    }

    Row {
        family: family.to_string(),
        n: g.n,
        m: g.m(),
        dijk_ms,
        bmssp_ms,
        verified,
        d_extract: dc.heap_extract_min,
        d_relax: dc.relaxations,
        b_relax: bc.relaxations,
        b_pulls: bc.queue_pull,
        b_recursive: bc.recursive_calls,
        b_pivots: bc.pivots_found,
    }
}

/// Bench the constant-degree transform: run both baselines on the transformed
/// graph; `verified` means BMSSP(transformed) projected to original ids matches
/// Dijkstra on the *original* graph.
pub fn bench_one_transformed(
    orig: &Graph,
    src: u32,
    family: &str,
    iters: usize,
    ab: &Ablation,
) -> Row {
    let t = to_constant_degree(orig);
    let cfg = BmsspConfig {
        use_pivots: ab.use_pivots,
        queue_impl: ab.queue_impl,
        partial_execution: ab.partial_execution,
        ..BmsspConfig::from_n(t.n)
    };
    let d_orig = dijkstra(orig, src, &mut Counters::new());

    let mut dc = Counters::new();
    let _d = dijkstra(&t, src, &mut dc);
    let mut bc = Counters::new();
    let b = BmsspEngine::new(&t, cfg.clone(), &mut bc).run(src);
    let verified = close(&d_orig, &b[..orig.n]);

    let mut dijk_ms = f64::MAX;
    for _ in 0..iters {
        let mut c = Counters::new();
        let tt = Instant::now();
        dijkstra(&t, src, &mut c);
        dijk_ms = dijk_ms.min(tt.elapsed().as_secs_f64() * 1e3);
    }
    let mut bmssp_ms = f64::MAX;
    for _ in 0..iters {
        let mut c = Counters::new();
        let tt = Instant::now();
        BmsspEngine::new(&t, cfg.clone(), &mut c).run(src);
        bmssp_ms = bmssp_ms.min(tt.elapsed().as_secs_f64() * 1e3);
    }

    Row {
        family: format!("{family}_tr"),
        n: t.n,
        m: t.m(),
        dijk_ms,
        bmssp_ms,
        verified,
        d_extract: dc.heap_extract_min,
        d_relax: dc.relaxations,
        b_relax: bc.relaxations,
        b_pulls: bc.queue_pull,
        b_recursive: bc.recursive_calls,
        b_pivots: bc.pivots_found,
    }
}

pub fn header() -> String {
    format!(
        "{:<10} {:>10} {:>12} {:>10} {:>10} {:>7} {:>8} {:>10} {:>10} {:>10} {:>10} {:>9} {:>6}",
        "family",
        "n",
        "m",
        "dijk(ms)",
        "bmssp(ms)",
        "speedup",
        "verified",
        "d#ext",
        "d#rel",
        "b#rel",
        "b#pull",
        "b#rec",
        "b#piv"
    )
}

pub fn fmt_row(r: &Row) -> String {
    let speedup = if r.bmssp_ms > 0.0 {
        r.dijk_ms / r.bmssp_ms
    } else {
        0.0
    };
    format!(
        "{:<10} {:>10} {:>12} {:>10.3} {:>10.3} {:>7.2} {:>8} {:>10} {:>10} {:>10} {:>10} {:>9} {:>6}",
        r.family, r.n, r.m, r.dijk_ms, r.bmssp_ms, speedup, r.verified,
        r.d_extract, r.d_relax, r.b_relax, r.b_pulls, r.b_recursive, r.b_pivots
    )
}

fn build_set() -> Vec<(String, Graph)> {
    let int = WeightDist::Int { min: 1, max: 100 };
    let real = WeightDist::Real { min: 0.0, max: 1.0 };
    let scale = std::env::var("BMSSP_SCALE").is_ok();
    let mut v = Vec::new();
    let base: Vec<(&str, Graph)> = if scale {
        // Phase 5 scaling set (n up to ~1e6). Skip transform pass for these.
        vec![
            ("er_c4_1e4", er(10_000, 4, &int)),
            ("er_c4_1e5", er(100_000, 4, &int)),
            ("er_c4_3e5", er(300_000, 4, &int)),
            ("er_c4_1e6", er(1_000_000, 4, &int)),
            ("grid_316", grid(316, 316, &int)),
            ("grid_500", grid(500, 500, &int)),
            ("layered_1e5", layered(200, 500, 4, &real)),
            ("layered_3e5", layered(300, 1000, 4, &real)),
        ]
    } else {
        vec![
            ("er_c2", er(10_000, 2, &int)),
            ("er_c4", er(10_000, 4, &int)),
            ("er_c8", er(10_000, 8, &int)),
            ("er_c4_1e5", er(100_000, 4, &int)),
            ("grid_100", grid(100, 100, &int)),
            ("grid_316", grid(316, 316, &int)),
            ("pl_c2", power_law(10_000, 2, &int)),
            ("pl_c4_1e5", power_law(100_000, 4, &int)),
            ("layered", layered(200, 500, 4, &real)),
            ("er_c4_real", er(10_000, 4, &real)),
        ]
    };
    for (label, g) in base {
        v.push((label.to_string(), g));
    }
    v
}

fn er(n: usize, c: usize, wd: &WeightDist) -> Graph {
    bmssp_rs::graph::er_random(n, c, 0xB0555EED, wd)
}
fn grid(nx: usize, ny: usize, wd: &WeightDist) -> Graph {
    bmssp_rs::graph::grid(nx, ny, 0x9A11D, wd)
}
fn power_law(n: usize, c: usize, wd: &WeightDist) -> Graph {
    bmssp_rs::graph::power_law(n, c, 0x50AE1A, wd)
}
fn layered(layers: usize, width: usize, out: usize, wd: &WeightDist) -> Graph {
    bmssp_rs::graph::layered(layers, width, out, 0x1A1A, wd)
}

pub fn run_all(iters: usize, ab: &Ablation) -> String {
    let mut out = String::new();
    out.push_str(&header());
    out.push('\n');
    let set = build_set();
    let scale = std::env::var("BMSSP_SCALE").is_ok();
    for (family, g) in &set {
        let src = 0u32;
        let row = bench_one(g, src, family, iters, ab);
        out.push_str(&fmt_row(&row));
        out.push('\n');
        // Asymptotic normalizers (Phase 5): time / (m log^{2/3} n) vs time / (m + n log n).
        if scale && row.n > 1 && row.m > 0 {
            let logn = (row.n as f64).log2();
            let denom_bm = (row.m as f64) * logn.powf(2.0 / 3.0);
            let denom_dj = row.m as f64 + (row.n as f64) * logn;
            out.push_str(&format!(
                "  asymptotics: bmssp/(m log^(2/3)n)={:.3e}  dijk/(m+n log n)={:.3e}\n",
                row.bmssp_ms / denom_bm,
                row.dijk_ms / denom_dj
            ));
        }
    }
    if !scale {
        out.push_str("\nconstant-degree transform (BMSSP run on the transformed graph):\n");
        for (family, g) in &set {
            let row = bench_one_transformed(g, 0u32, family, iters, ab);
            out.push_str(&fmt_row(&row));
            out.push('\n');
        }
    }
    out
}

/// Ablation knobs from the environment: BMSSP_NO_PIVOTS=1 disables pivots,
/// BMSSP_QUEUE=map|block selects the queue, BMSSP_PARTIAL=1 enables the
/// k*2^(l*t) workload bound. k/t/l are always derived from the *graph's* n so
/// the top-level call stays successful (its bound exceeds |U| <= n).
pub struct Ablation {
    pub use_pivots: bool,
    pub queue_impl: QueueKind,
    pub partial_execution: bool,
}

impl std::fmt::Debug for Ablation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "pivots={} queue={:?} partial={}",
            self.use_pivots, self.queue_impl, self.partial_execution
        )
    }
}

pub fn ablation_from_env() -> Ablation {
    let queue_impl = match std::env::var("BMSSP_QUEUE").as_deref() {
        Ok("block") => QueueKind::Block,
        _ => QueueKind::BTreeMap,
    };
    Ablation {
        use_pivots: std::env::var("BMSSP_NO_PIVOTS").is_err(),
        queue_impl,
        partial_execution: std::env::var("BMSSP_PARTIAL").is_ok(),
    }
}

#[allow(dead_code)]
fn main() {
    let iters = std::env::var("BMSSP_BENCH_ITERS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);
    let ab = ablation_from_env();
    println!("iters={iters} ablation={ab:?}\n");
    print!("{}", run_all(iters, &ab));
}
