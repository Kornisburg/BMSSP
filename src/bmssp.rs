use std::collections::BinaryHeap;

use crate::counters::Counters;
use crate::dijkstra::{MinState, INF};
use crate::graph::Graph;
use crate::params::params;
use crate::queue::{BlockQueue, PartialQueue, QueueOps};

const EPS: f64 = 1e-12;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QueueKind {
    BTreeMap,
    Block,
}

#[derive(Debug, Clone)]
pub struct BmsspConfig {
    pub k: usize,
    pub t: usize,
    pub l: usize,
    /// ablation: set false to skip FindPivots (P = S, W = empty).
    pub use_pivots: bool,
    /// ablation: halt a level once it has completed > k*2^(l*t) vertices
    /// (Lemma 3.1 partial execution). false = run to D-empty (B' = B always).
    pub partial_execution: bool,
    /// queue implementation behind the partial-order structure.
    pub queue_impl: QueueKind,
}

impl BmsspConfig {
    pub fn from_n(n: usize) -> Self {
        let (k, t, l) = params(n);
        BmsspConfig {
            k,
            t,
            l,
            use_pivots: true,
            partial_execution: false,
            queue_impl: QueueKind::BTreeMap,
        }
    }

    /// k * 2^(l*t): the Lemma 3.1 partial-execution workload bound for a call
    /// at level `l`. Saturates to usize::MAX instead of overflowing.
    pub fn partial_limit(&self, l: usize) -> usize {
        self.k.saturating_mul(
            1usize
                .checked_shl((l * self.t).min(63) as u32)
                .unwrap_or(usize::MAX),
        )
    }

    /// Lemma 3.3 block size M = 2^((l-1)*t) for a call at level `l >= 1`.
    pub fn block_size(&self, l: usize) -> usize {
        if l == 0 {
            return 1;
        }
        1usize
            .checked_shl(((l - 1) * self.t).min(63) as u32)
            .unwrap_or(usize::MAX)
    }
}

pub struct BmsspEngine<'a> {
    g: &'a Graph,
    dist: Vec<f64>,
    cfg: BmsspConfig,
    counters: &'a mut Counters,
    marked: Vec<u64>,
    marked_u: Vec<u64>,
    parent: Vec<u32>,
    sub: Vec<usize>,
    epoch: u64,
    u_epoch: u64,
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
            marked_u: vec![0; g.n],
            u_epoch: 0,
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

    /// BMSSP(l, B, S) -> (B', U). With `partial_execution` off (default) it runs
    /// to D-empty and B' = B. With it on, a level halts once it has completed
    /// more than k*2^(l*t) vertices and returns B' = the last child's boundary.
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

        let mut pq = self.make_queue(l, b);
        for &x in &pivots {
            if self.dist[x as usize] < b - EPS {
                pq.insert(x, self.dist[x as usize]);
                self.counters.queue_insert += 1;
            }
        }
        // FindPivots' bounded Bellman-Ford improves dhat for W vertices without
        // routing them, and (k rounds) never relaxes their own edges. The paper
        // rescues these on upper levels via the <= equality re-insertion
        // (Remark 3.4); we instead queue W directly so this call processes them
        // (edges relaxed, out-neighbours routed) instead of losing them.
        for &x in &wset {
            if self.dist[x as usize] < b - EPS {
                pq.insert(x, self.dist[x as usize]);
                self.counters.queue_insert += 1;
            }
        }

        let partial_limit = self.cfg.partial_limit(l);
        let mut u_total: Vec<u32> = Vec::new();
        let mut halted = false;
        let mut last_bp = b;
        // The paper's |U| is a set: count each completed vertex once (children
        // may touch the same vertex many times across calls), so the partial
        // workload bound k*2^(l*t) is compared against a distinct count.
        self.u_epoch += 1;
        let ue = self.u_epoch;
        while !pq.is_empty() {
            if self.cfg.partial_execution && u_total.len() > partial_limit {
                halted = true;
                if self.trace {
                    eprintln!("HALT l={l} b={b} |U|={} limit={partial_limit} last_bp={last_bp}", u_total.len());
                }
                // Unprocessed queue items were improved by this call but never
                // completed, so their edges were never relaxed. Hand them up so
                // the caller relaxes from them; the paper instead re-inserts on
                // equality (Remark 3.4) to re-use lower-level relaxations.
                for (v, _) in pq.drain() {
                    if self.marked_u[v as usize] != ue {
                        self.marked_u[v as usize] = ue;
                        u_total.push(v);
                    }
                }
                break;
            }
            let (si, bi) = pq.pull();
            self.counters.queue_pull += 1;
            if si.is_empty() {
                continue;
            }
            // Child execution bounded by this bucket's separation bound.
            let (bp_i, ui) = self.bmssp(l - 1, bi, &si);
            last_bp = bp_i;
            for &x in &ui {
                if self.marked_u[x as usize] != ue {
                    self.marked_u[x as usize] = ue;
                    u_total.push(x);
                }
            }

            // Relax from every vertex completed by the child, routing on strict
            // improvement only. Every strict improvement is re-queued: the
            // paper's lower bound B'_i (cand < B'_i -> nothing) is only sound
            // when completed vertices keep exact distances, which our
            // touched-based relaxation (it relaxes from still-in-heap vertices)
            // does not guarantee. So we route on the safe interval: [0, B_i) ->
            // front-loaded batch, [B_i, B) -> main structure.
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
                            eprintln!("R: {u} -> {v} cand={cand} b={b} bi={bi} bp={bp_i}");
                        }
                        if cand >= bi - EPS {
                            pq.insert(v, cand);
                            self.counters.queue_insert += 1;
                        } else {
                            to_batch.push((v, cand));
                        }
                    }
                }
            }
            if !to_batch.is_empty() {
                pq.batch_prepend(&to_batch);
                self.counters.queue_batch_prepend += 1;
            }
        }

        let b_prime = if halted {
            self.counters.partial_executions += 1;
            last_bp
        } else {
            b
        };
        for &x in &wset {
            if self.dist[x as usize] < b_prime - EPS && self.marked_u[x as usize] != ue {
                self.marked_u[x as usize] = ue;
                u_total.push(x);
            }
        }
        (b_prime, u_total)
    }

    /// BaseCase, generalised to a multi-source seed set. Bounded mini-Dijkstra
    /// settling at most k+1 vertices total (seeds included); returns B' = the
    /// paper's `max dhat in U0` when a boundary was hit (else B), and every
    /// vertex whose distance this call touched (settled plus still-in-heap), so
    /// no reachable vertex's out-edges are ever skipped by a caller.
    fn base_case(&mut self, b: f64, s: &[u32]) -> (f64, Vec<u32>) {
        self.counters.base_case_calls += 1;
        self.epoch += 1;
        let e = self.epoch;

        let mut u0: Vec<u32> = Vec::new();
        let mut touched: Vec<u32> = Vec::new();
        let mut heap: BinaryHeap<MinState> = BinaryHeap::new();
        let mut settled = 0usize;
        let mut boundary = 0.0f64;
        for &x in s {
            if self.dist[x as usize] < b - EPS {
                if self.marked[x as usize] != e {
                    self.marked[x as usize] = e;
                    touched.push(x);
                }
                settled += 1;
                boundary = boundary.max(self.dist[x as usize]);
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
                settled += 1;
                boundary = boundary.max(cost);
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
            if settled > self.cfg.k {
                break;
            }
        }

        let bprime = if settled <= self.cfg.k { b } else { boundary };
        (bprime, touched)
    }

    /// Construct the partial-order queue for a call at level `l` (Lemma 3.3
    /// block size M = 2^((l-1)t); ignored by the BTreeMap implementation).
    fn make_queue(&self, l: usize, b: f64) -> QueueOps {
        match self.cfg.queue_impl {
            QueueKind::BTreeMap => QueueOps::Map(PartialQueue::new(b)),
            QueueKind::Block => QueueOps::Block(BlockQueue::new(b, self.cfg.block_size(l))),
        }
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
            // edge u -> v (cand == dist[v]) as its parent. Reset parents for
            // every W vertex first so no stale parent survives across calls.
            for &x in &w {
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
