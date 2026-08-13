//! Hand-crafted graphs with known shortest-path structure.

mod common;

use bmssp_rs::graph::Graph;
use common::assert_bmssp_matches_dijkstra;

fn edges(n: usize, e: &[(u32, u32, f64)]) -> Graph {
    Graph::from_edges(n, e)
}

#[test]
fn chain_exact_distances() {
    let g = edges(5, &[(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 4, 1.0)]);
    let d = bmssp_rs::bmssp::barrier_breaker_sssp(&g, 0, &mut bmssp_rs::counters::Counters::new());
    assert_eq!(d, vec![0.0, 1.0, 2.0, 3.0, 4.0]);
    // source in the middle: right side unreachable, left side reachable
    let d = bmssp_rs::bmssp::barrier_breaker_sssp(&g, 2, &mut bmssp_rs::counters::Counters::new());
    assert_eq!(d[0], f64::INFINITY);
    assert_eq!(d[1], f64::INFINITY);
    assert_eq!(d[2], 0.0);
    assert_eq!(d[4], 2.0);
}

#[test]
fn two_paths_with_equal_length() {
    // 0->1 (1), 0->2 (1), 1->3 (1), 2->3 (1): dist[3] = 2 via either path
    let g = edges(4, &[(0, 1, 1.0), (0, 2, 1.0), (1, 3, 1.0), (2, 3, 1.0)]);
    let d = bmssp_rs::bmssp::barrier_breaker_sssp(&g, 0, &mut bmssp_rs::counters::Counters::new());
    assert_eq!(d[3], 2.0);
    assert_eq!(d, vec![0.0, 1.0, 1.0, 2.0]);
}

#[test]
fn zero_weight_cycle() {
    // 0 -> 1 (0), 1 -> 0 (0), 1 -> 2 (5): dist[2] = 5, dist[0] = 0
    let g = edges(3, &[(0, 1, 0.0), (1, 0, 0.0), (1, 2, 5.0)]);
    assert_bmssp_matches_dijkstra(&g, 0, 0.0);
    let d = bmssp_rs::bmssp::barrier_breaker_sssp(&g, 0, &mut bmssp_rs::counters::Counters::new());
    assert_eq!(d[2], 5.0);
}

#[test]
fn star_with_hub() {
    // hub 0 connected to 100 spokes with distinct weights; spokes link pairwise.
    let n = 101;
    let mut e = Vec::new();
    for i in 1..n {
        e.push((0, i as u32, i as f64));
        if i < n - 1 {
            e.push((i as u32, (i + 1) as u32, 1000.0));
        }
    }
    let g = edges(n, &e);
    assert_bmssp_matches_dijkstra(&g, 0, 0.0);
}

#[test]
fn dense_complete_graph() {
    let n = 6;
    let mut e = Vec::new();
    for u in 0..n {
        for v in 0..n {
            if u != v {
                e.push((u as u32, v as u32, ((u + v) % 5 + 1) as f64));
            }
        }
    }
    let g = edges(n, &e);
    for src in 0..n as u32 {
        assert_bmssp_matches_dijkstra(&g, src, 0.0);
    }
}

#[test]
fn unreachable_vertices() {
    // 0 -> 1 -> 2, and isolated 3,4,5
    let g = edges(6, &[(0, 1, 1.0), (1, 2, 1.0)]);
    let d = bmssp_rs::bmssp::barrier_breaker_sssp(&g, 0, &mut bmssp_rs::counters::Counters::new());
    assert_eq!(d[3], f64::INFINITY);
    assert_eq!(d[2], 2.0);
}

#[test]
fn single_vertex_graph() {
    let g = edges(1, &[]);
    let d = bmssp_rs::bmssp::barrier_breaker_sssp(&g, 0, &mut bmssp_rs::counters::Counters::new());
    assert_eq!(d, vec![0.0]);
}

#[test]
fn parallel_edges() {
    // 0 -> 1 via weight 5 and weight 1: min is 1
    let g = edges(2, &[(0, 1, 5.0), (0, 1, 1.0)]);
    let d = bmssp_rs::bmssp::barrier_breaker_sssp(&g, 0, &mut bmssp_rs::counters::Counters::new());
    assert_eq!(d[1], 1.0);
}
