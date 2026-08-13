use std::collections::BinaryHeap;

use crate::counters::Counters;
use crate::dijkstra::{MinState, INF};
use crate::graph::Graph;
use crate::params::params;
use crate::queue::PartialQueue;

const EPS: f64 = 1e-12;

#[derive(Debug, Clone)]
pub struct BmsspConfig {
    pub k: usize,
    pub t: usize,
    pub l: usize,
    /// ablation: set false to skip FindPivots (P = S, W = empty).
    pub use_pivots: bool,
}

impl BmsspConfig {
    pub fn from_n(n: usize) -> Self {
        let (k, t, l) = params(n);
        BmsspConfig {
            k,
            t,
            l,
            use_pivots: true,
        }
    }
}

pub struct BmsspEngine<'a> {
    g: &'a Graph,
    dist: Vec<f64>,
    cfg: BmsspConfig,
    counters: &'a mut Counters,
    marked: Vec<u64>,
    parent: Vec<u32>,
    sub: Vec<usize>,
    epoch: u64,
    trace: bool,
}

/// Top-level driver: BMSSP(l, B = infinity, S = {src}).
pub fn barrier_breaker_sssp(g: &Graph, src: u32, counters: &mut Counters) -> Vec<f64> {
    let cfg = BmsspConfig::from_n(g.n);
    BmsspEngine::new(g, cfg, counters).run(src)
}

impl<'a> BmsspEngine<'a> {
    pub fn new(g: &'a Graph, cfg: BmsspConfig, counters: &'a mut Counters) -> Self {
        BmsspEngine {
            g,
            dist: vec![INF; g.n],
            cfg,
            counters,
            marked: vec![0; g.n],
            parent: vec![u32::MAX; g.n],
            sub: vec![0; g.n],
            epoch: 0,
            trace: std::env::var("BMSSP_TRACE").is_ok(),
        }
    }

    pub fn run(&mut self, src: u32) -> Vec<f64> {
        assert!(src < self.g.n as u32, "source out of range");
        self.dist.fill(INF);
        self.dist[src as usize] = 0.0;
        let l = self.cfg.l;
        let (_, _u) = self.bmssp(l, INF, &[src]);
        self.dist.clone()
    }

    /// BMSSP(l, B, S) -> (B', U). Phase-1 reference variant: always runs to
    /// D-empty, so B' = B (successful execution). Returns the set of vertices
    /// completed below B that were touched by this execution.
    fn bmssp(&mut self, l: usize, b: f64, s: &[u32]) -> (f64, Vec<u32>) {
        debug_assert!(!s.is_empty());
        self.counters.recursive_calls += 1;
        if l == 0 {
            return self.base_case(b, s);
        }

        let (pivots, wset) = if self.cfg.use_pivots {
            self.find_pivots(b, s)
        } else {
            (s.to_vec(), Vec::new())
        };

        let mut pq = PartialQueue::new(b);
        for &x in &pivots {
            if self.dist[x as usize] < b - EPS {
                pq.insert(x, self.dist[x as usize]);
                self.counters.queue_insert += 1;
            }
        }

        let mut u_total: Vec<u32> = Vec::new();
        while !pq.is_empty() {
            let (si, bi) = pq.pull();
            self.counters.queue_pull += 1;
            if si.is_empty() {
                continue;
            }
            // Child execution bounded by this bucket's separation bound.
            let (_, ui) = self.bmssp(l - 1, bi, &si);
            u_total.extend_from_slice(&ui);

            // Relax from every vertex completed by the child. Discoveries are
            // routed on strict improvement only: within the current bound but
            // below this bucket's separation -> front-loaded batch (processed
            // next); at/above it -> main structure (processed in order).
            let mut to_batch: Vec<(u32, f64)> = Vec::new();
            for &u in &ui {
                let du = self.dist[u as usize];
                for i in self.g.edge_range(u as usize) {
                    let v = self.g.to[i];
                    let w = self.g.weight[i];
                    self.counters.relaxations += 1;
                    let cand = du + w;
                    if cand < b - EPS && cand < self.dist[v as usize] - EPS {
                        self.dist[v as usize] = cand;
                        self.parent[v as usize] = u;
                        if self.trace {
                            eprintln!("R: {u} -> {v} cand={cand} b={b} bi={bi}");
                        }
                        if cand < bi - EPS {
                            to_batch.push((v, cand));
                        } else {
                            pq.insert(v, cand);
                            self.counters.queue_insert += 1;
                        }
                    }
                }
            }
            if !to_batch.is_empty() {
                pq.batch_prepend(&to_batch);
                self.counters.queue_batch_prepend += 1;
            }
        }

        for &x in &wset {
            if self.dist[x as usize] < b - EPS {
                u_total.push(x);
            }
        }
        (b, u_total)
    }

    /// BaseCase, generalised to a multi-source seed set. Bounded mini-Dijkstra
    /// settling at most k+1 vertices total (seeds included). Returns every
    /// vertex whose distance this call touched (settled plus still-in-heap),
    /// so no reachable vertex's out-edges are ever skipped by a caller.
    fn base_case(&mut self, b: f64, s: &[u32]) -> (f64, Vec<u32>) {
        self.counters.base_case_calls += 1;
        self.epoch += 1;
        let e = self.epoch;

        let mut u0: Vec<u32> = Vec::new();
        let mut touched: Vec<u32> = Vec::new();
        let mut heap: BinaryHeap<MinState> = BinaryHeap::new();
        for &x in s {
            if self.dist[x as usize] < b - EPS {
                if self.marked[x as usize] != e {
                    self.marked[x as usize] = e;
                    touched.push(x);
                }
                heap.push(MinState {
                    cost: self.dist[x as usize],
                    node: x,
                });
                self.counters.heap_insert += 1;
            }
        }

        while let Some(MinState { cost, node }) = heap.pop() {
            self.counters.heap_extract_min += 1;
            if cost > self.dist[node as usize] + EPS {
                continue;
            }
            if cost >= b - EPS {
                break;
            }
            if self.marked[node as usize] != e {
                self.marked[node as usize] = e;
                u0.push(node);
                touched.push(node);
            }
            if u0.len() > self.cfg.k {
                break;
            }
            for i in self.g.edge_range(node as usize) {
                let v = self.g.to[i];
                let w = self.g.weight[i];
                self.counters.relaxations += 1;
                let cand = cost + w;
                if cand < self.dist[v as usize] - EPS && cand < b - EPS {
                    self.dist[v as usize] = cand;
                    self.parent[v as usize] = node;
                    if self.trace {
                        eprintln!("B: {node} -> {v} cand={cand} b={b}");
                    }
                    if self.marked[v as usize] != e {
                        self.marked[v as usize] = e;
                        touched.push(v);
                    }
                    heap.push(MinState { cost: cand, node: v });
                    self.counters.heap_insert += 1;
                }
            }
        }

        let boundary = if u0.len() <= self.cfg.k {
            b
        } else {
            u0.iter()
                .map(|&x| self.dist[x as usize])
                .fold(0.0f64, f64::max)
        };
        (boundary, touched)
    }

    /// FindPivots (Algorithm 1 / Lemma 3.2): k rounds of bounded Bellman-Ford
    /// from seeds S building a predecessor forest F. W = visited set; P = roots
    /// in S whose F-subtree has >= k vertices (or P = S on early stop / fallback).
    fn find_pivots(&mut self, b: f64, s: &[u32]) -> (Vec<u32>, Vec<u32>) {
        self.counters.find_pivots_calls += 1;
        self.epoch += 1;
        let e = self.epoch;

        let mut w: Vec<u32> = Vec::new();
        let mut frontier: Vec<u32> = Vec::new();
        for &x in s {
            if self.dist[x as usize] < b - EPS {
                if self.marked[x as usize] != e {
                    self.marked[x as usize] = e;
                    w.push(x);
                }
                frontier.push(x);
            }
        }

        let mut early_stop = false;
        for _ in 0..self.cfg.k {
            if frontier.is_empty() {
                break;
            }
            let mut next: Vec<u32> = Vec::new();
            for &u in &frontier {
                let du = self.dist[u as usize];
                for i in self.g.edge_range(u as usize) {
                    let v = self.g.to[i];
                    let wt = self.g.weight[i];
                    self.counters.relaxations += 1;
                    let cand = du + wt;
                    if cand < b - EPS && cand <= self.dist[v as usize] + EPS {
                        if cand < self.dist[v as usize] - EPS {
                            self.dist[v as usize] = cand;
                        }
                        if self.marked[v as usize] != e {
                            self.marked[v as usize] = e;
                            w.push(v);
                            next.push(v);
                        }
                    }
                }
            }
            frontier = next;
            if w.len() > self.cfg.k * s.len().max(1) {
                early_stop = true;
                break;
            }
        }

        let pivots = if early_stop || w.is_empty() {
            s.to_vec()
        } else {
            // Forest F: for each v in W, pick the smallest u with an equality
            // edge u -> v (cand == dist[v]) as its parent.
            for &x in s {
                self.parent[x as usize] = u32::MAX;
            }
            for &u in &w {
                let du = self.dist[u as usize];
                for i in self.g.edge_range(u as usize) {
                    let v = self.g.to[i];
                    let wt = self.g.weight[i];
                    let cand = du + wt;
                    if self.marked[v as usize] == e
                        && (cand - self.dist[v as usize]).abs() <= EPS
                        && (self.parent[v as usize] == u32::MAX || u < self.parent[v as usize])
                    {
                        self.parent[v as usize] = u;
                    }
                }
            }
            self.sub.fill(0);
            for &v in &w {
                self.sub[v as usize] = 1;
            }
            for &v in w.iter().rev() {
                let p = self.parent[v as usize];
                if p != u32::MAX {
                    self.sub[p as usize] += self.sub[v as usize];
                }
            }
            let mut pivots: Vec<u32> = Vec::new();
            for &x in s {
                if self.marked[x as usize] == e
                    && self.parent[x as usize] == u32::MAX
                    && self.sub[x as usize] >= self.cfg.k
                {
                    pivots.push(x);
                }
            }
            if pivots.is_empty() {
                s.to_vec()
            } else {
                pivots
            }
        };
        self.counters.pivots_found += pivots.len() as u64;
        (pivots, w)
    }
}
