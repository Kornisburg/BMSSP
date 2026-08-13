//! Debug trace for a BMSSP-vs-Dijkstra counterexample.
use bmssp_rs::bmssp::barrier_breaker_sssp;
use bmssp_rs::counters::Counters;
use bmssp_rs::dijkstra::dijkstra;
use bmssp_rs::graph::Graph;

fn main() {
    let edges = [
        (0u32, 1u32, 0.0),
        (0, 5, 1.0),
        (1, 2, 2.0),
        (1, 5, 1.0),
        (2, 3, 0.0),
        (3, 5, 0.0),
        (4, 0, 0.0),
        (5, 2, 2.0),
        (5, 3, 2.0),
    ];
    let n = 6usize;
    let src = 4u32;
    let g = Graph::from_edges(n, &edges);
    let mut dc = Counters::new();
    let d = dijkstra(&g, src, &mut dc);
    let mut bc = Counters::new();
    let b = barrier_breaker_sssp(&g, src, &mut bc);
    println!("dijkstra: {d:?}");
    println!("bmssp:    {b:?}");
    println!("dijkstra counters: {dc:?}");
    println!("bmssp counters:    {bc:?}");
    for (i, (x, y)) in d.iter().zip(&b).enumerate() {
        if (x - y).abs() > 1e-9 {
            println!("MISMATCH at vertex {i}: dijkstra={x:?} bmssp={y:?}");
        }
    }
}
