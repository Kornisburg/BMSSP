//! Every BmsspConfig variant (queue backend x partial execution x pivots) must
//! match the Dijkstra baseline on the same graph families as the default
//! configuration.

mod common;

use bmssp_rs::bmssp::{barrier_breaker_sssp, BmsspConfig, BmsspEngine, QueueKind};
use bmssp_rs::counters::Counters;
use bmssp_rs::graph::{er_random, WeightDist};
use common::assert_close;
use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

fn variants_for(n: usize) -> Vec<BmsspConfig> {
    // Use graph-native k/t/l so the top-level partial bound k·2^(l·t) dominates
    // n (successful top-level execution) instead of a fixed from_n(200).
    let base = BmsspConfig::from_n(n);
    let mut v = Vec::new();
    for partial in [false, true] {
        for queue in [QueueKind::BTreeMap, QueueKind::Block] {
            for pivots in [true, false] {
                v.push(BmsspConfig {
                    partial_execution: partial,
                    queue_impl: queue,
                    use_pivots: pivots,
                    ..base.clone()
                });
            }
        }
    }
    v
}

fn check_variants(g: &bmssp_rs::graph::Graph, src: u32, eps: f64) {
    let d = bmssp_rs::dijkstra::dijkstra(g, src, &mut Counters::new());
    for cfg in variants_for(g.n) {
        let mut c = Counters::new();
        let b = BmsspEngine::new(g, cfg.clone(), &mut c).run(src);
        if let Err(e) = common::try_close(&d, &b, eps) {
            let mut edges: Vec<String> = Vec::new();
            for u in 0..g.n {
                for i in g.edge_range(u) {
                    edges.push(format!("{} {} {}", u, g.to[i], g.weight[i]));
                }
            }
            panic!("n={} src={src} cfg={cfg:?} {e}\nedges: {}", g.n, edges.join(" "));
        }
        let _ = c; // counters alive: exercises drop paths too
    }
}

#[test]
fn variants_match_on_small_random_graphs() {
    let mut rng = ChaCha8Rng::seed_from_u64(0xC0FFEE);
    let dists = [
        WeightDist::Unit,
        WeightDist::Int { min: 0, max: 3 },
        WeightDist::Real { min: 0.0, max: 1.0 },
    ];
    for trial in 0..200u64 {
        let n = rng.gen_range(1..=20usize);
        let wd = &dists[(trial % dists.len() as u64) as usize];
        let p = rng.gen_range(0.0..0.4);
        let mut edges = Vec::new();
        for u in 0..n as u32 {
            for v in 0..n as u32 {
                if rng.gen_bool(p) {
                    edges.push((u, v, wd.sample(&mut rng)));
                }
            }
        }
        let g = bmssp_rs::graph::Graph::from_edges(n, &edges);
        let src = rng.gen_range(0..n as u32);
        let eps = if matches!(wd, WeightDist::Real { .. }) { 1e-9 } else { 0.0 };
        check_variants(&g, src, eps);
    }
}

#[test]
fn variants_match_on_sized_families() {
    let wd = WeightDist::Int { min: 1, max: 10 };
    for n in [2usize, 5, 20, 100, 400] {
        for seed in 0..3u64 {
            let g = er_random(n, 4, seed, &wd);
            check_variants(&g, 0, 0.0);
        }
    }    let wd2 = WeightDist::Int { min: 0, max: 9 };
    let g = bmssp_rs::graph::grid(10, 12, 5, &wd2);
    check_variants(&g, 0, 0.0);
    check_variants(&g, 60, 0.0);
    let g = bmssp_rs::graph::chain(300, 3, &wd2);
    check_variants(&g, 0, 0.0);
    let g = bmssp_rs::graph::layered(20, 15, 4, 5, &wd2);
    check_variants(&g, 0, 0.0);
    let g = bmssp_rs::graph::power_law(300, 3, 9, &wd2);
    check_variants(&g, 0, 0.0);
}

#[test]
fn variants_agree_with_default_on_dense_graphs() {
    let wd = WeightDist::Int { min: 1, max: 5 };
    for n in [4usize, 10, 30] {
        for p in [0.2, 0.8] {
            let g = bmssp_rs::graph::dense(n, p, 3, &wd);
            for src in 0..n as u32 {
                let d = barrier_breaker_sssp(&g, src, &mut Counters::new());
                for cfg in variants_for(g.n) {
                    let mut c = Counters::new();
                    let b = BmsspEngine::new(&g, cfg, &mut c).run(src);
                    assert_close(&d, &b, 0.0);
                }
            }
        }
    }
}