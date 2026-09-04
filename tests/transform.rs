mod common;

use bmssp_rs::counters::Counters;
use bmssp_rs::dijkstra::dijkstra;
use bmssp_rs::graph::{chain, dense, er_random, grid, layered, power_law, Graph, WeightDist};
use bmssp_rs::transform::{max_out_degree, to_constant_degree};
use common::{assert_bmssp_matches_dijkstra, assert_close};

#[test]
fn transform_max_out_degree_is_at_most_two() {
    let wd = WeightDist::Int { min: 1, max: 10 };
    let graphs = [
        er_random(50, 3, 1, &wd),
        er_random(200, 8, 2, &wd),
        dense(30, 0.3, 3, &wd),
        dense(60, 0.1, 4, &wd),
        power_law(100, 3, 5, &wd),
        grid(15, 15, 6, &wd),
        chain(40, 7, &wd),
    ];
    for g in graphs {
        let t = to_constant_degree(&g);
        assert!(
            max_out_degree(&t) <= 2,
            "max out-degree {} > 2",
            max_out_degree(&t)
        );
    }
}

#[test]
fn transform_preserves_distances_across_families() {
    let wd = WeightDist::Int { min: 1, max: 20 };
    let mut cases = Vec::new();
    cases.push((er_random(80, 4, 11, &wd), 0u32));
    cases.push((er_random(80, 4, 11, &wd), 37u32));
    cases.push((dense(40, 0.25, 12, &wd), 13u32));
    cases.push((grid(20, 20, 13, &wd), 0u32));
    cases.push((chain(50, 14, &wd), 25u32));
    cases.push((power_law(90, 3, 15, &wd), 7u32));
    cases.push((layered(6, 10, 3, 16, &wd), 0u32));
    for (g, src) in cases {
        let d_orig = dijkstra(&g, src, &mut Counters::new());
        let t = to_constant_degree(&g);
        let d_trans = dijkstra(&t, src, &mut Counters::new());
        assert_close(&d_orig, &d_trans[..g.n], 1e-9);
    }
}

#[test]
fn transform_handles_zero_and_parallel_edges() {
    let wd = WeightDist::Int { min: 0, max: 5 };
    let g = er_random(60, 4, 21, &wd);
    let d_orig = dijkstra(&g, 0, &mut Counters::new());
    let t = to_constant_degree(&g);
    let d_trans = dijkstra(&t, 0, &mut Counters::new());
    assert_close(&d_orig, &d_trans[..g.n], 1e-9);
}

#[test]
fn bmssp_on_transformed_graph_matches_dijkstra_original() {
    let wd = WeightDist::Int { min: 1, max: 15 };
    for g in [
        er_random(70, 4, 31, &wd),
        dense(35, 0.3, 32, &wd),
        grid(18, 18, 33, &wd),
        power_law(80, 3, 34, &wd),
        layered(5, 12, 4, 35, &wd),
    ] {
        let t = to_constant_degree(&g);
        let d = dijkstra(&g, 0, &mut Counters::new());
        let b = bmssp_rs::bmssp::barrier_breaker_sssp(&t, 0, &mut Counters::new());
        assert_close(&d, &b[..g.n], 1e-9);
        assert_bmssp_matches_dijkstra(&t, 0, 1e-9);
    }
}

#[test]
fn transform_preserves_weighted_edge_multiset() {
    let g = dense(25, 0.4, 41, &WeightDist::Int { min: 1, max: 7 });
    let t = to_constant_degree(&g);
    // Every weighted edge in the transformed graph is exactly one original edge
    // (leaf carries it unchanged); every other edge is a zero-weight aux edge.
    let mut orig: Vec<(u32, u64)> = Vec::new();
    for u in 0..g.n {
        for i in g.edge_range(u) {
            orig.push((g.to[i], g.weight[i] as u64));
        }
    }
    let mut trans: Vec<(u32, u64)> = Vec::new();
    for u in 0..t.n {
        for i in t.edge_range(u) {
            let w = t.weight[i];
            if w > 0.0 {
                trans.push((t.to[i], w as u64));
            }
        }
    }
    orig.sort_unstable();
    trans.sort_unstable();
    assert_eq!(orig, trans, "weighted-edge multiset changed by transform");
}

#[test]
fn transform_exhaustive_small_graphs() {
    // Exhaustively check n = 2..=4 over all simple digraphs with weights {0,1,2}.
    let weights = [0.0f64, 1.0, 2.0];
    for n in 2..=4usize {
        for mask in 0u32..(1u32 << (n * n)) {
            let mut edges: Vec<(u32, u32, f64)> = Vec::new();
            for u in 0..n {
                for v in 0..n {
                    if u == v {
                        continue;
                    }
                    if mask & (1 << (u * n + v)) != 0 {
                        let w = weights[((mask >> 16) as usize + u + v) % weights.len()];
                        edges.push((u as u32, v as u32, w));
                    }
                }
            }
            let g = Graph::from_edges(n, &edges);
            let t = to_constant_degree(&g);
            assert!(max_out_degree(&t) <= 2);
            for src in 0..n as u32 {
                let d = dijkstra(&g, src, &mut Counters::new());
                let dt = dijkstra(&t, src, &mut Counters::new());
                assert_close(&d, &dt[..n], 1e-9);
            }
        }
    }
}
