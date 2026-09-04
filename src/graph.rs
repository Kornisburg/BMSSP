use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::collections::HashSet;
use std::ops::Range;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WeightDist {
    /// Every edge has weight exactly 1.
    Unit,
    /// Uniform integer weights in `min..=max` (exact in f64 for small sums).
    Int { min: i64, max: i64 },
    /// Uniform real weights in `min..max`.
    Real { min: f64, max: f64 },
}

impl WeightDist {
    pub fn sample(&self, rng: &mut impl Rng) -> f64 {
        match *self {
            WeightDist::Unit => 1.0,
            WeightDist::Int { min, max } => rng.gen_range(min..=max) as f64,
            WeightDist::Real { min, max } => rng.gen_range(min..max),
        }
    }
}

/// Compressed-Sparse-Row directed graph.
#[derive(Debug, Clone)]
pub struct Graph {
    pub n: usize,
    pub offsets: Vec<usize>,
    pub to: Vec<u32>,
    pub weight: Vec<f64>,
}

impl Graph {
    pub fn empty(n: usize) -> Self {
        Graph {
            n,
            offsets: vec![0; n + 1],
            to: Vec::new(),
            weight: Vec::new(),
        }
    }

    pub fn m(&self) -> usize {
        self.to.len()
    }

    pub fn out_degree(&self, u: usize) -> usize {
        self.offsets[u + 1] - self.offsets[u]
    }

    pub fn edge_range(&self, u: usize) -> Range<usize> {
        self.offsets[u]..self.offsets[u + 1]
    }

    pub fn from_edges(n: usize, edges: &[(u32, u32, f64)]) -> Self {
        let mut out: Vec<Vec<(u32, f64)>> = vec![Vec::new(); n];
        for &(u, v, w) in edges {
            assert!(u < n as u32 && v < n as u32, "vertex out of range");
            assert!(
                w.is_finite() && w >= 0.0,
                "weight must be finite, non-negative"
            );
            out[u as usize].push((v, w));
        }
        // Sort each adjacency list by (weight, target) for scan locality in
        // Dijkstra / BaseCase / FindPivots relax loops.
        for list in &mut out {
            list.sort_by(|a, b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
        }
        let mut offsets = vec![0usize; n + 1];
        for u in 0..n {
            offsets[u + 1] = offsets[u] + out[u].len();
        }
        let m = offsets[n];
        let mut to = Vec::with_capacity(m);
        let mut weight = Vec::with_capacity(m);
        for list in &out {
            for &(v, w) in list {
                to.push(v);
                weight.push(w);
            }
        }
        Graph {
            n,
            offsets,
            to,
            weight,
        }
    }

    pub fn to_adjacency(&self) -> Vec<Vec<(usize, f64)>> {
        let mut adj = vec![Vec::new(); self.n];
        for (u, list) in adj.iter_mut().enumerate() {
            for i in self.edge_range(u) {
                list.push((self.to[i] as usize, self.weight[i]));
            }
        }
        adj
    }
}

fn rng_for(seed: u64) -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(seed)
}

fn random_target(rng: &mut impl Rng, n: usize, u: u32) -> u32 {
    if n <= 1 {
        return 0;
    }
    let mut v = rng.gen_range(0..n as u32);
    while v == u {
        v = rng.gen_range(0..n as u32);
    }
    v
}

/// Erdos-Renyi sparse graph with m = c*n directed edges (self-loops avoided).
pub fn er_random(n: usize, c: usize, seed: u64, wd: &WeightDist) -> Graph {
    let mut rng = rng_for(seed);
    let m = n.saturating_mul(c);
    let mut edges = Vec::with_capacity(m);
    for _ in 0..m {
        let u = rng.gen_range(0..n as u32);
        let v = random_target(&mut rng, n, u);
        edges.push((u, v, wd.sample(&mut rng)));
    }
    Graph::from_edges(n, &edges)
}

/// Erdos-Renyi dense graph with edge probability p.
pub fn dense(n: usize, p: f64, seed: u64, wd: &WeightDist) -> Graph {
    let mut rng = rng_for(seed);
    let mut edges = Vec::new();
    for u in 0..n as u32 {
        for v in 0..n as u32 {
            if u != v && rng.gen_bool(p) {
                edges.push((u, v, wd.sample(&mut rng)));
            }
        }
    }
    Graph::from_edges(n, &edges)
}

/// Directed grid: east + south edges on an nx x ny lattice.
pub fn grid(nx: usize, ny: usize, seed: u64, wd: &WeightDist) -> Graph {
    let mut rng = rng_for(seed);
    let id = |x: usize, y: usize| (y * nx + x) as u32;
    let mut edges = Vec::new();
    for y in 0..ny {
        for x in 0..nx {
            if x + 1 < nx {
                edges.push((id(x, y), id(x + 1, y), wd.sample(&mut rng)));
            }
            if y + 1 < ny {
                edges.push((id(x, y), id(x, y + 1), wd.sample(&mut rng)));
            }
        }
    }
    Graph::from_edges(nx * ny, &edges)
}

/// Forward chain 0 -> 1 -> ... -> n-1.
pub fn chain(n: usize, seed: u64, wd: &WeightDist) -> Graph {
    let mut rng = rng_for(seed);
    let mut edges = Vec::new();
    for u in 0..n.saturating_sub(1) {
        edges.push((u as u32, (u + 1) as u32, wd.sample(&mut rng)));
    }
    Graph::from_edges(n, &edges)
}

/// Layered graph: `layers` layers of `width` vertices; each vertex sends
/// `out` edges to random vertices in the next layer (theory-friendly).
pub fn layered(layers: usize, width: usize, out: usize, seed: u64, wd: &WeightDist) -> Graph {
    let mut rng = rng_for(seed);
    let mut edges = Vec::new();
    let n = layers * width;
    for layer in 0..layers.saturating_sub(1) {
        for pos in 0..width {
            let u = (layer * width + pos) as u32;
            for _ in 0..out {
                let v = ((layer + 1) * width + rng.gen_range(0..width)) as u32;
                edges.push((u, v, wd.sample(&mut rng)));
            }
        }
    }
    Graph::from_edges(n, &edges)
}

/// Barabasi-Albert preferential attachment: each new vertex i sends `c` edges
/// to earlier vertices chosen proportional to total degree.
pub fn power_law(n: usize, c: usize, seed: u64, wd: &WeightDist) -> Graph {
    let mut rng = rng_for(seed);
    let mut edges: Vec<(u32, u32, f64)> = Vec::new();
    let mut deg = vec![0usize; n];
    let m0 = (c + 1).min(n);
    for u in 0..m0 {
        for v in 0..m0 {
            if u != v {
                edges.push((u as u32, v as u32, wd.sample(&mut rng)));
                deg[u] += 1;
                deg[v] += 1;
            }
        }
    }
    for i in m0..n {
        let c_targets = c.min(i);
        let mut chosen: HashSet<u32> = HashSet::with_capacity(c_targets);
        while chosen.len() < c_targets {
            let total: usize = deg[..i].iter().map(|d| d + 1).sum();
            let mut r = rng.gen_range(0..total);
            let mut t = 0;
            for (j, d) in deg[..i].iter().enumerate() {
                let w = d + 1;
                if r < w {
                    t = j;
                    break;
                }
                r -= w;
            }
            if !chosen.contains(&(t as u32)) {
                chosen.insert(t as u32);
                edges.push((i as u32, t as u32, wd.sample(&mut rng)));
                deg[i] += 1;
                deg[t] += 1;
            }
        }
    }
    Graph::from_edges(n, &edges)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csr_invariants() {
        let g = er_random(100, 4, 42, &WeightDist::Unit);
        assert_eq!(g.offsets.len(), 101);
        assert_eq!(g.offsets[0], 0);
        assert_eq!(g.offsets[100], g.to.len());
        assert_eq!(g.m(), 400);
        for u in 0..100 {
            let r = g.edge_range(u);
            assert_eq!(r.len(), g.out_degree(u));
            for i in r {
                assert!(g.to[i] < 100);
                assert!(g.weight[i] >= 0.0);
            }
        }
    }

    #[test]
    fn generators_produce_valid_graphs() {
        let wd = WeightDist::Int { min: 1, max: 10 };
        let graphs = [
            er_random(50, 4, 1, &wd),
            dense(20, 0.3, 2, &wd),
            grid(8, 8, 3, &wd),
            chain(50, 4, &wd),
            layered(5, 10, 3, 5, &wd),
            power_law(50, 2, 6, &wd),
        ];
        for g in &graphs {
            for u in 0..g.n {
                for i in g.edge_range(u) {
                    assert!(g.to[i] < g.n as u32);
                    assert!(g.weight[i].is_finite() && g.weight[i] >= 0.0);
                }
            }
        }
    }
}
