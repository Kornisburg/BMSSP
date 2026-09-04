//! Randomized property tests: BMSSP must match the Dijkstra baseline on all
//! graph families, sizes, weight distributions, and several sources.

mod common;

use bmssp_rs::graph::{er_random, WeightDist};
use common::{assert_bmssp_matches_dijkstra, floyd_warshall};
use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

fn exact(wd: &WeightDist) -> bool {
    matches!(wd, WeightDist::Int { .. } | WeightDist::Unit)
}

#[test]
fn many_small_random_graphs_match_dijkstra() {
    let mut rng = ChaCha8Rng::seed_from_u64(0x5EED);
    let dists = [
        WeightDist::Unit,
        WeightDist::Int { min: 0, max: 3 },
        WeightDist::Int { min: 1, max: 100 },
        WeightDist::Real { min: 0.0, max: 1.0 },
    ];
    for trial in 0..3000u64 {
        let n = rng.gen_range(1..=25usize);
        let wd = &dists[(trial % dists.len() as u64) as usize];
        // random edge probability (sparse..dense), plus occasional weird shapes
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
        let eps = if exact(wd) { 0.0 } else { 1e-9 };
        assert_bmssp_matches_dijkstra(&g, src, eps);
    }
}

#[test]
fn er_sparse_across_sizes_and_densities() {
    let wd = WeightDist::Int { min: 1, max: 10 };
    for n in [2usize, 3, 5, 10, 50, 200, 1000, 5000] {
        for c in [1usize, 2, 4, 8, 16] {
            for seed in 0..4u64 {
                let g = er_random(
                    n,
                    c,
                    seed.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(1),
                    &wd,
                );
                for src in [0u32, (n / 2) as u32] {
                    assert_bmssp_matches_dijkstra(&g, src, 0.0);
                }
            }
        }
    }
}

#[test]
fn grids_match() {
    let wd = WeightDist::Int { min: 0, max: 9 };
    for (nx, ny) in [(2, 2), (3, 5), (8, 8), (20, 20), (45, 45)] {
        let g = bmssp_rs::graph::grid(nx, ny, 0xA11, &wd);
        for src in [0u32, (nx * ny / 2) as u32] {
            assert_bmssp_matches_dijkstra(&g, src, 0.0);
        }
    }
}

#[test]
fn dense_graphs_match() {
    let wd = WeightDist::Int { min: 1, max: 5 };
    for n in [4usize, 8, 20] {
        for p in [0.2, 0.5, 0.9] {
            let g = bmssp_rs::graph::dense(n, p, 7, &wd);
            for src in 0..n as u32 {
                assert_bmssp_matches_dijkstra(&g, src, 0.0);
            }
        }
    }
}

#[test]
fn chains_and_layered_and_powerlaw_match() {
    let wd = WeightDist::Int { min: 0, max: 9 };
    let g = bmssp_rs::graph::chain(2000, 3, &wd);
    assert_bmssp_matches_dijkstra(&g, 0, 0.0);
    assert_bmssp_matches_dijkstra(&g, 1000, 0.0);
    let g = bmssp_rs::graph::layered(50, 40, 4, 5, &wd);
    assert_bmssp_matches_dijkstra(&g, 0, 0.0);
    let g = bmssp_rs::graph::power_law(2000, 3, 9, &wd);
    assert_bmssp_matches_dijkstra(&g, 0, 0.0);
}

#[test]
fn tie_heavy_graphs_match() {
    // unit weights create huge numbers of equal-distance ties
    let g = er_random(500, 8, 11, &WeightDist::Unit);
    assert_bmssp_matches_dijkstra(&g, 0, 0.0);
    let g = bmssp_rs::graph::grid(25, 25, 13, &WeightDist::Unit);
    assert_bmssp_matches_dijkstra(&g, 0, 0.0);
    let g = er_random(500, 4, 17, &WeightDist::Int { min: 1, max: 2 });
    assert_bmssp_matches_dijkstra(&g, 0, 0.0);
}

#[test]
fn medium_graphs_match_floyd_oracle() {
    // cross-check against the Floyd-Warshall oracle (not just Dijkstra)
    for seed in 0..8u64 {
        let g = er_random(12, 3, seed, &WeightDist::Int { min: 0, max: 5 });
        let oracle = floyd_warshall(&g);
        for src in 0..12u32 {
            let b = bmssp_rs::bmssp::barrier_breaker_sssp(
                &g,
                src,
                &mut bmssp_rs::counters::Counters::new(),
            );
            common::assert_close(&b, &oracle[src as usize * 12..(src as usize + 1) * 12], 0.0);
        }
    }
}
