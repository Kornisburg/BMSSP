use std::time::Instant;

use bmssp_rs::bmssp::{barrier_breaker_sssp, BmsspConfig, BmsspEngine};
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

pub fn bench_one(g: &Graph, src: u32, family: &str, iters: usize, use_pivots: bool) -> Row {
    let mut dc = Counters::new();
    let d = dijkstra(g, src, &mut dc);
    let mut bc = Counters::new();
    let b = if use_pivots {
        barrier_breaker_sssp(g, src, &mut bc)
    } else {
        let cfg = BmsspConfig {
            use_pivots: false,
            ..BmsspConfig::from_n(g.n)
        };
        BmsspEngine::new(g, cfg, &mut bc).run(src)
    };
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
        if use_pivots {
            barrier_breaker_sssp(g, src, &mut c);
        } else {
            let cfg = BmsspConfig {
                use_pivots: false,
                ..BmsspConfig::from_n(g.n)
            };
            BmsspEngine::new(g, cfg, &mut c).run(src);
        }
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
    use_pivots: bool,
) -> Row {
    let t = to_constant_degree(orig);
    let d_orig = dijkstra(orig, src, &mut Counters::new());

    let mut dc = Counters::new();
    let _d = dijkstra(&t, src, &mut dc);
    let mut bc = Counters::new();
    let b = if use_pivots {
        barrier_breaker_sssp(&t, src, &mut bc)
    } else {
        let cfg = BmsspConfig {
            use_pivots: false,
            ..BmsspConfig::from_n(t.n)
        };
        BmsspEngine::new(&t, cfg, &mut bc).run(src)
    };
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
        if use_pivots {
            barrier_breaker_sssp(&t, src, &mut c);
        } else {
            let cfg = BmsspConfig {
                use_pivots: false,
                ..BmsspConfig::from_n(t.n)
            };
            BmsspEngine::new(&t, cfg, &mut c).run(src);
        }
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
        "family", "n", "m", "dijk(ms)", "bmssp(ms)", "speedup", "verified",
        "d#ext", "d#rel", "b#rel", "b#pull", "b#rec", "b#piv"
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
    let mut v = Vec::new();
    for (label, g) in [
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
    ] {
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

pub fn run_all(iters: usize, use_pivots: bool) -> String {
    let mut out = String::new();
    out.push_str(&header());
    out.push('\n');
    let set = build_set();
    for (family, g) in &set {
        let src = 0u32;
        let row = bench_one(g, src, family, iters, use_pivots);
        out.push_str(&fmt_row(&row));
        out.push('\n');
    }
    out.push_str("\nconstant-degree transform (BMSSP run on the transformed graph):\n");
    for (family, g) in &set {
        let row = bench_one_transformed(g, 0u32, family, iters, use_pivots);
        out.push_str(&fmt_row(&row));
        out.push('\n');
    }
    out
}

#[allow(dead_code)]
fn main() {
    let iters = std::env::var("BMSSP_BENCH_ITERS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);
    let use_pivots = std::env::var("BMSSP_NO_PIVOTS").is_err();
    println!("iters={iters} use_pivots={use_pivots}\n");
    print!("{}", run_all(iters, use_pivots));
}
