use crate::graph::Graph;

/// Convert `g` into an equivalent graph with maximum out-degree 2 by splitting
/// every vertex of out-degree > 2 into a binary tree of zero-weight auxiliary
/// nodes (vertex splitting, directed). Original vertices keep ids `0..g.n`;
/// distances to them are unchanged, so SSSP on the transformed graph projected
/// back to `0..g.n` is exact.
///
/// Size: total transformed vertices = n + sum over deg>2 of (2*deg - 2) = O(m).
pub fn to_constant_degree(g: &Graph) -> Graph {
    let n = g.n;
    let mut aux_count = 0usize;
    for u in 0..n {
        let d = g.out_degree(u);
        if d > 2 {
            aux_count += 2 * d - 2;
        }
    }

    let total = n + aux_count;
    let mut adj: Vec<Vec<(u32, f64)>> = vec![Vec::new(); total];
    let mut next = n as u32;
    for u in 0..n {
        let d = g.out_degree(u);
        if d <= 2 {
            for i in g.edge_range(u) {
                adj[u].push((g.to[i], g.weight[i]));
            }
        } else {
            let edges: Vec<(u32, f64)> = g
                .edge_range(u)
                .map(|i| (g.to[i], g.weight[i]))
                .collect();
            build_tree(&mut adj, &mut next, u as u32, &edges, 0, d);
        }
    }
    debug_assert_eq!(next as usize, total);

    let mut offsets = vec![0usize; total + 1];
    for (v, list) in adj.iter().enumerate() {
        offsets[v + 1] = offsets[v] + list.len();
    }
    let mut to = vec![0u32; offsets[total]];
    let mut weight = vec![0.0f64; offsets[total]];
    for v in 0..total {
        for (k, &(t, w)) in adj[v].iter().enumerate() {
            to[offsets[v] + k] = t;
            weight[offsets[v] + k] = w;
        }
    }
    Graph {
        n: total,
        offsets,
        to,
        weight,
    }
}

/// Build a binary tree over `edges[lo..hi]` rooted at node `id`, allocating
/// fresh aux node ids from `next` upward. Every internal node has exactly two
/// children joined by zero-weight edges; every leaf carries one original edge.
fn build_tree(
    adj: &mut Vec<Vec<(u32, f64)>>,
    next: &mut u32,
    id: u32,
    edges: &[(u32, f64)],
    lo: usize,
    hi: usize,
) {
    if hi - lo == 1 {
        adj[id as usize].push(edges[lo]);
        return;
    }
    let mid = lo + (hi - lo).div_ceil(2);
    let left = *next;
    *next += 1;
    let right = *next;
    *next += 1;
    adj[id as usize].push((left, 0.0));
    adj[id as usize].push((right, 0.0));
    build_tree(adj, next, left, edges, lo, mid);
    build_tree(adj, next, right, edges, mid, hi);
}

/// Maximum out-degree over all vertices (for tests).
pub fn max_out_degree(g: &Graph) -> usize {
    (0..g.n).map(|u| g.out_degree(u)).max().unwrap_or(0)
}
