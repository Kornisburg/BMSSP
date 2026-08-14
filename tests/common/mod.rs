use bmssp_rs::counters::Counters;
use bmssp_rs::dijkstra::dijkstra;
use bmssp_rs::graph::Graph;

/// Floyd-Warshall all-pairs oracle for tiny graphs (weights non-negative).
#[allow(dead_code)]
pub fn floyd_warshall(g: &Graph) -> Vec<f64> {
    let n = g.n;
    let inf = f64::INFINITY;
    let mut d = vec![inf; n * n];
    for u in 0..n {
        d[u * n + u] = 0.0;
    }
    for u in 0..n {
        for i in g.edge_range(u) {
            let v = g.to[i] as usize;
            d[u * n + v] = d[u * n + v].min(g.weight[i]);
        }
    }
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                let nd = d[i * n + k] + d[k * n + j];
                if nd < d[i * n + j] {
                    d[i * n + j] = nd;
                }
            }
        }
    }
    d
}

pub fn try_close(a: &[f64], b: &[f64], eps: f64) -> Result<(), String> {
    if a.len() != b.len() {
        return Err("length mismatch".into());
    }
    for (x, y) in a.iter().zip(b) {
        if x.is_infinite() && y.is_infinite() {
            continue;
        }
        if x.is_infinite() || y.is_infinite() {
            return Err(format!("reachability mismatch: {x} vs {y}"));
        }
        let scale = 1.0 + x.abs().max(y.abs());
        if (x - y).abs() > eps * scale {
            return Err(format!("dist mismatch: {x} vs {y}"));
        }
    }
    Ok(())
}

pub fn assert_close(a: &[f64], b: &[f64], eps: f64) {
    if let Err(e) = try_close(a, b, eps) {
        panic!("{e}");
    }
}

/// Assert bmssp distances match Dijkstra (tolerance 0 for exact integer sums).
#[allow(dead_code)]
pub fn assert_bmssp_matches_dijkstra(g: &Graph, src: u32, eps: f64) {
    let d = dijkstra(g, src, &mut Counters::new());
    let b = bmssp_rs::bmssp::barrier_breaker_sssp(g, src, &mut Counters::new());
    assert_close(&d, &b, eps);
}
