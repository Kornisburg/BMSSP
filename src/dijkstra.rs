use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::counters::Counters;
use crate::graph::Graph;

pub const INF: f64 = f64::INFINITY;

#[derive(Copy, Clone)]
pub struct MinState {
    pub cost: f64,
    pub node: u32,
}

impl PartialEq for MinState {
    fn eq(&self, other: &Self) -> bool {
        self.cost == other.cost && self.node == other.node
    }
}

impl Eq for MinState {}

impl Ord for MinState {
    fn cmp(&self, other: &Self) -> Ordering {
        // min-heap via reverse
        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(Ordering::Equal)
            .then_with(|| self.node.cmp(&other.node))
    }
}

impl PartialOrd for MinState {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Standard binary-heap Dijkstra on a CSR graph (lazy-deletion).
pub fn dijkstra(g: &Graph, src: u32, counters: &mut Counters) -> Vec<f64> {
    let n = g.n;
    let mut dist = vec![INF; n];
    dist[src as usize] = 0.0;
    let mut heap: BinaryHeap<MinState> = BinaryHeap::new();
    heap.push(MinState {
        cost: 0.0,
        node: src,
    });
    counters.heap_insert += 1;

    while let Some(MinState { cost, node }) = heap.pop() {
        counters.heap_extract_min += 1;
        if cost > dist[node as usize] {
            continue;
        }
        for i in g.edge_range(node as usize) {
            let v = g.to[i];
            let w = g.weight[i];
            counters.relaxations += 1;
            let next = cost + w;
            if next < dist[v as usize] {
                dist[v as usize] = next;
                heap.push(MinState {
                    cost: next,
                    node: v,
                });
                counters.heap_insert += 1;
            }
        }
    }
    dist
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::WeightDist;

    pub fn floyd_warshall(g: &Graph) -> Vec<f64> {
        let n = g.n;
        let mut d = vec![INF; n * n];
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

    pub fn assert_close(a: &[f64], b: &[f64], eps: f64) {
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b) {
            if x.is_infinite() && y.is_infinite() {
                continue;
            }
            assert!(
                (x - y).abs() <= eps,
                "dist mismatch: {x} vs {y} (src-dist vectors differ)"
            );
        }
    }

    #[test]
    fn dijkstra_matches_floyd_on_tiny_graphs() {
        let wd = WeightDist::Int { min: 1, max: 10 };
        for n in [1usize, 2, 3, 5, 8] {
            for c in [1usize, 2, 4] {
                for seed in 0..20u64 {
                    let g = crate::graph::er_random(n, c, seed, &wd);
                    let oracle = floyd_warshall(&g);
                    for src in 0..n as u32 {
                        let d = dijkstra(&g, src, &mut Counters::new());
                        assert_close(&d, &oracle[src as usize * n..(src as usize + 1) * n], 1e-9);
                    }
                }
            }
        }
    }

    #[test]
    fn dijkstra_handles_unreachable() {
        let g = crate::graph::chain(5, 1, &WeightDist::Unit);
        let d = dijkstra(&g, 3, &mut Counters::new());
        assert_eq!(d[0], INF);
        assert_eq!(d[4], 1.0);
    }
}
