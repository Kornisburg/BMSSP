//! Brute-force search for a minimal graph where BMSSP disagrees with Dijkstra.
use bmssp_rs::bmssp::barrier_breaker_sssp;
use bmssp_rs::counters::Counters;
use bmssp_rs::dijkstra::dijkstra;
use bmssp_rs::graph::Graph;

fn check(n: usize, edges: &[(u32, u32, f64)], src: u32, eps: f64) -> bool {
    let g = Graph::from_edges(n, edges);
    let d = dijkstra(&g, src, &mut Counters::new());
    let b = barrier_breaker_sssp(&g, src, &mut Counters::new());
    d.iter()
        .zip(&b)
        .all(|(x, y)| {
            if x.is_infinite() && y.is_infinite() {
                true
            } else if x.is_infinite() || y.is_infinite() {
                false
            } else {
                (x - y).abs() <= eps * (1.0 + x.abs().max(y.abs()))
            }
        })
}

fn main() {
    let weights = [0.0f64, 1.0, 2.0];
    for n in 2..=4usize {
        // enumerate all simple digraphs on n vertices with weights from {0,1,2}
        let mut found = false;
        let mut edges: Vec<(u32, u32, f64)> = Vec::new();
        for mask in 0u32..(1u32 << (n * n)) {
            edges.clear();
            for u in 0..n {
                for v in 0..n {
                    if u == v {
                        continue;
                    }
                    let bit = u * n + v;
                    if mask & (1 << bit) != 0 {
                        let w = weights[((mask >> 16) as usize + u + v) % weights.len()];
                        edges.push((u as u32, v as u32, w));
                    }
                }
            }
            for src in 0..n as u32 {
                if !check(n, &edges, src, 1e-12) {
                    println!("counterexample n={n} src={src} edges={edges:?}");
                    found = true;
                    break;
                }
            }
            if found {
                break;
            }
        }
        if found {
            println!("(first counterexample in n={n})");
            return;
        }
        println!("n={n}: all graphs ok");
    }
    println!("no counterexample found up to n=4 (exhaustive)");

    // random sampling for larger n
    use rand::Rng;
    let mut rng = rand::thread_rng();
    let mut edges: Vec<(u32, u32, f64)> = Vec::new();
    for n in 5..=6usize {
        for _ in 0..200_000 {
            edges.clear();
            let p = match n {
                5 => rng.gen_range(0.15..0.45),
                _ => rng.gen_range(0.08..0.3),
            };
            for u in 0..n {
                for v in 0..n {
                    if u != v && rng.gen_bool(p) {
                        let w = weights[rng.gen_range(0..weights.len())];
                        edges.push((u as u32, v as u32, w));
                    }
                }
            }
            for src in 0..n as u32 {
                if !check(n, &edges, src, 1e-12) {
                    println!("counterexample n={n} src={src} edges={edges:?}");
                    return;
                }
            }
        }
        println!("n={n}: 200k random graphs ok");
    }
    println!("no counterexample found");
}
