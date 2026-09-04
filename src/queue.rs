use std::cmp::Ordering;
use std::collections::BTreeMap;

/// Total-order wrapper for f64 keys (no NaN; `total_cmp` gives a total order).
#[derive(Debug, Clone, Copy)]
struct K(f64);

impl PartialEq for K {
    fn eq(&self, other: &Self) -> bool {
        self.0.to_bits() == other.0.to_bits()
    }
}

impl Eq for K {}

impl Ord for K {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

impl PartialOrd for K {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Phase-1 partial-order queue (Lemma 3.3 semantics, BTreeMap-backed).
///
/// Values are key/value pairs (value = current `dhat`). Operations:
/// - `insert(v, key)`: add a pair.
/// - `pull()`: remove the smallest bucket (all pairs with the smallest key) and
///   return it together with the separation bound `B_i` = smallest remaining key
///   (or the queue's upper bound `B` if empty).
/// - `batch_prepend(items)`: front-load a batch of smaller keys (via insert into
///   the map, which sorts them ahead of existing larger keys).
///
/// A vertex may appear at several keys; buckets are drained whole to preserve
/// the separation invariant. Asymptotically this is O(log Q) per op; the block
/// structure of the paper (`BlockQueue`) is the Phase-3 upgrade.
pub struct PartialQueue {
    map: BTreeMap<K, Vec<u32>>,
    bound: f64,
}

impl PartialQueue {
    pub fn new(bound: f64) -> Self {
        PartialQueue {
            map: BTreeMap::new(),
            bound,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    pub fn insert(&mut self, v: u32, key: f64) {
        debug_assert!(key.is_finite());
        self.map.entry(K(key)).or_default().push(v);
    }

    /// Returns (bucket, B_i): the smallest bucket and the bound of the rest.
    pub fn pull(&mut self) -> (Vec<u32>, f64) {
        if let Some((_, mut bucket)) = self.map.pop_first() {
            let bi = self.map.keys().next().map(|k| k.0).unwrap_or(self.bound);
            (std::mem::take(&mut bucket), bi)
        } else {
            (Vec::new(), self.bound)
        }
    }

    pub fn batch_prepend(&mut self, items: &[(u32, f64)]) {
        for &(v, k) in items {
            debug_assert!(k.is_finite());
            self.map.entry(K(k)).or_default().push(v);
        }
    }

    /// Remove and return everything still queued.
    pub fn drain(&mut self) -> Vec<(u32, f64)> {
        let mut out = Vec::new();
        for (k, vs) in std::mem::take(&mut self.map) {
            for v in vs {
                out.push((v, k.0));
            }
        }
        out
    }
}

/// Identity key for a block: `(upper bound, sequence)` so two blocks sharing an
/// upper bound remain distinct.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct UbKey {
    ub: K,
    seq: u64,
}

impl Ord for UbKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.ub
            .cmp(&other.ub)
            .then_with(|| self.seq.cmp(&other.seq))
    }
}

impl PartialOrd for UbKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// A block: at most `m` key/value pairs, unordered inside, with `ub` an upper
/// bound on all values (we keep it exact at insertion time; it may go stale
/// high after Pull deletions, which is safe for the ordering invariant).
#[derive(Debug)]
struct Block {
    items: Vec<(u32, f64)>,
    ub: f64,
}

/// Partition an overflowing block into two blocks that keep the inter-block
/// value ordering: everything strictly below the block's max goes to the lower
/// half and the max-valued pairs to the upper half, so the upper half never
/// carries a value smaller than the max of any block that precedes it. If the
/// block's values are all equal the median split is used instead (either half
/// is valid then).
fn split_block(mut blk: Block) -> (Block, Block) {
    blk.items.sort_by(|a, b| a.1.total_cmp(&b.1));
    let below = blk.items.partition_point(|it| it.1 < blk.ub);
    let mid = if below == 0 || below == blk.items.len() {
        blk.items.len() / 2
    } else {
        below
    };
    let right = blk.items.split_off(mid);
    let ub_l = blk
        .items
        .iter()
        .map(|it| it.1)
        .fold(f64::NEG_INFINITY, f64::max);
    let ub_r = right
        .iter()
        .map(|it| it.1)
        .fold(f64::NEG_INFINITY, f64::max);
    (
        Block {
            items: blk.items,
            ub: ub_l,
        },
        Block {
            items: right,
            ub: ub_r,
        },
    )
}

/// Block-based partial-order queue (Lemma 3.3, `D1` structure). Items live in
/// value-ordered blocks of at most `m` pairs, keyed in a BST by block upper
/// bound, so an Insert locates its block in O(log(#blocks)) = O(log(N/m)).
///
/// `Pull` returns the `m` smallest values with ties at the boundary taken whole
/// (so the separation bound strictly exceeds every returned value, preserving
/// the algorithm's interval invariant even on equal keys). `BatchPrepend` is
/// routed through the ordered insert path rather than the paper's separate O(1)
/// front-list `D0`, because our routing does not guarantee prepends are smaller
/// than everything already queued; this is the documented simplification vs
/// Lemma 3.3's amortized O(1) prepend.
#[derive(Debug)]
pub struct BlockQueue {
    blocks: BTreeMap<UbKey, Block>,
    bound: f64,
    m: usize,
    seq: u64,
}

impl BlockQueue {
    pub fn new(bound: f64, m: usize) -> Self {
        BlockQueue {
            blocks: BTreeMap::new(),
            bound,
            m: m.max(1),
            seq: 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }

    pub fn insert(&mut self, v: u32, key: f64) {
        debug_assert!(key.is_finite());
        let m = self.m;

        if self.blocks.is_empty() {
            self.seq += 1;
            self.blocks.insert(
                UbKey {
                    ub: K(key),
                    seq: self.seq,
                },
                Block {
                    items: vec![(v, key)],
                    ub: key,
                },
            );
            return;
        }

        // First block whose upper bound is >= key.
        let cand = self
            .blocks
            .range(UbKey { ub: K(key), seq: 0 }..)
            .next()
            .map(|(k, _)| (k.ub, k.seq));
        let Some((ub, seq)) = cand else {
            // key exceeds every block's max: extend the last block if it has
            // room, else start a fresh trailing block.
            let (lub, lseq) = {
                let (k, _) = self.blocks.iter().next_back().unwrap();
                (k.ub, k.seq)
            };
            if self.blocks[&UbKey { ub: lub, seq: lseq }].items.len() < m {
                let mut blk = self.blocks.remove(&UbKey { ub: lub, seq: lseq }).unwrap();
                blk.items.push((v, key));
                blk.ub = key;
                self.blocks.insert(
                    UbKey {
                        ub: K(key),
                        seq: lseq,
                    },
                    blk,
                );
            } else {
                self.seq += 1;
                self.blocks.insert(
                    UbKey {
                        ub: K(key),
                        seq: self.seq,
                    },
                    Block {
                        items: vec![(v, key)],
                        ub: key,
                    },
                );
            }
            return;
        };

        let target = UbKey { ub, seq };
        let min = self.blocks[&target]
            .items
            .iter()
            .map(|it| it.1)
            .fold(f64::INFINITY, f64::min);

        if key < min {
            // key belongs in the gap before `target` (every block before it has
            // max strictly below key). Extend the predecessor if it has room,
            // else start a fresh block between them.
            let prev = self.blocks.range(..target).next_back().map(|(k, _)| *k);
            if let Some(pk) = prev {
                if self.blocks[&pk].items.len() < m {
                    let mut blk = self.blocks.remove(&pk).unwrap();
                    blk.items.push((v, key));
                    blk.ub = key;
                    self.blocks.insert(
                        UbKey {
                            ub: K(key),
                            seq: pk.seq,
                        },
                        blk,
                    );
                    return;
                }
            }
            self.seq += 1;
            self.blocks.insert(
                UbKey {
                    ub: K(key),
                    seq: self.seq,
                },
                Block {
                    items: vec![(v, key)],
                    ub: key,
                },
            );
            return;
        }

        // key in [min, ub]: append into the found block; split when it grows
        // past m (the split keeps the ordering, see split_block).
        let blk = self.blocks.get_mut(&target).unwrap();
        blk.items.push((v, key));
        if blk.items.len() > m {
            let blk = self.blocks.remove(&target).unwrap();
            let (b1, b2) = split_block(blk);
            self.seq += 1;
            self.blocks.insert(
                UbKey {
                    ub: K(b1.ub),
                    seq: self.seq,
                },
                b1,
            );
            self.seq += 1;
            self.blocks.insert(
                UbKey {
                    ub: K(b2.ub),
                    seq: self.seq,
                },
                b2,
            );
        }
    }

    pub fn batch_prepend(&mut self, items: &[(u32, f64)]) {
        for &(v, k) in items {
            self.insert(v, k);
        }
    }

    /// Returns (bucket, B_i): the `m` smallest values (ties taken whole) and the
    /// smallest remaining value, or `bound` when the structure is empty.
    pub fn pull(&mut self) -> (Vec<u32>, f64) {
        if self.blocks.is_empty() {
            return (Vec::new(), self.bound);
        }
        let m = self.m;
        let keys: Vec<UbKey> = self.blocks.keys().copied().collect();

        // Tie-aware front collection: gather blocks from the lowest ub until we
        // have >= m pairs and no equal-value pair straddles the next block.
        let mut collected: Vec<(u32, f64)> = Vec::new();
        let mut scanned: Vec<UbKey> = Vec::new();
        let mut cur_max = f64::NEG_INFINITY;
        for i in 0..keys.len() {
            let blk = &self.blocks[&keys[i]];
            scanned.push(keys[i]);
            for &it in &blk.items {
                collected.push(it);
                if it.1 > cur_max {
                    cur_max = it.1;
                }
            }
            if collected.len() >= m {
                let next_ties = i + 1 < keys.len()
                    && self.blocks[&keys[i + 1]]
                        .items
                        .iter()
                        .any(|it| it.1 <= cur_max);
                if !next_ties {
                    break;
                }
            }
        }

        // Everything drained: return all of it (the collection loop only stops
        // early with >= m pairs, so len < m means every block was scanned).
        if collected.len() < m {
            for k in &scanned {
                self.blocks.remove(k);
            }
            let vs = collected.iter().map(|&(v, _)| v).collect();
            return (vs, self.bound);
        }

        // S' = the m smallest values, with ties at the m-th taken whole; the
        // leftovers (which may outlive a fully-scanned structure) stay queued.
        collected.sort_by(|a, b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
        let vm = collected[m - 1].1;
        let take = collected.partition_point(|it| it.1 <= vm);
        let s: Vec<u32> = collected[..take].iter().map(|&(v, _)| v).collect();
        let remaining = collected[take..].to_vec();

        for k in &scanned {
            self.blocks.remove(k);
        }
        self.rebuild_blocks(&remaining);
        let x = self.min_value().unwrap_or(self.bound);
        (s, x)
    }

    fn min_value(&self) -> Option<f64> {
        self.blocks.first_key_value().map(|(_, blk)| {
            blk.items
                .iter()
                .map(|it| it.1)
                .fold(f64::INFINITY, f64::min)
        })
    }

    /// Re-insert sorted `items` as fresh blocks of at most `m` pairs each,
    /// preserving the inter-block value ordering.
    fn rebuild_blocks(&mut self, items: &[(u32, f64)]) {
        for chunk in items.chunks(self.m) {
            let ub = chunk
                .iter()
                .map(|it| it.1)
                .fold(f64::NEG_INFINITY, f64::max);
            self.seq += 1;
            self.blocks.insert(
                UbKey {
                    ub: K(ub),
                    seq: self.seq,
                },
                Block {
                    items: chunk.to_vec(),
                    ub,
                },
            );
        }
    }

    /// Remove and return everything still queued (in block order).
    pub fn drain(&mut self) -> Vec<(u32, f64)> {
        let mut out = Vec::new();
        for (_, blk) in std::mem::take(&mut self.blocks) {
            out.extend(blk.items);
        }
        out
    }
}

/// Uniform interface used by the BMSSP engine for both queue backends.
pub enum QueueOps {
    Map(PartialQueue),
    Block(BlockQueue),
}

impl QueueOps {
    pub fn is_empty(&self) -> bool {
        match self {
            QueueOps::Map(q) => q.is_empty(),
            QueueOps::Block(q) => q.is_empty(),
        }
    }

    pub fn insert(&mut self, v: u32, key: f64) {
        match self {
            QueueOps::Map(q) => q.insert(v, key),
            QueueOps::Block(q) => q.insert(v, key),
        }
    }

    pub fn pull(&mut self) -> (Vec<u32>, f64) {
        match self {
            QueueOps::Map(q) => q.pull(),
            QueueOps::Block(q) => q.pull(),
        }
    }

    pub fn batch_prepend(&mut self, items: &[(u32, f64)]) {
        match self {
            QueueOps::Map(q) => q.batch_prepend(items),
            QueueOps::Block(q) => q.batch_prepend(items),
        }
    }

    /// Remove and return everything still queued.
    pub fn drain(&mut self) -> Vec<(u32, f64)> {
        match self {
            QueueOps::Map(q) => q.drain(),
            QueueOps::Block(q) => q.drain(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pull_returns_smallest_bucket_and_bound() {
        let mut q = PartialQueue::new(100.0);
        q.insert(1, 5.0);
        q.insert(2, 5.0);
        q.insert(3, 9.0);
        q.insert(4, 1.0);
        let (b0, bi) = q.pull();
        assert_eq!(b0, vec![4]);
        assert_eq!(bi, 5.0);
        let (b1, bi) = q.pull();
        assert_eq!(b1, vec![1, 2]);
        assert_eq!(bi, 9.0);
        let (b2, bi) = q.pull();
        assert_eq!(b2, vec![3]);
        assert_eq!(bi, 100.0);
        assert!(q.is_empty());
    }

    #[test]
    fn batch_prepend_sorts_ahead() {
        let mut q = PartialQueue::new(100.0);
        q.insert(1, 9.0);
        q.batch_prepend(&[(2, 3.0), (3, 3.0)]);
        let (b, _) = q.pull();
        assert_eq!(b, vec![2, 3]);
    }

    #[test]
    fn empty_pull_returns_bound() {
        let mut q = PartialQueue::new(42.0);
        let (b, bi) = q.pull();
        assert!(b.is_empty());
        assert_eq!(bi, 42.0);
    }

    #[test]
    fn block_queue_basic() {
        let mut q = BlockQueue::new(100.0, 2);
        q.insert(1, 5.0);
        q.insert(2, 5.0);
        q.insert(3, 9.0);
        q.insert(4, 1.0);
        let (b, bi) = q.pull();
        // m=2, values {1,5,5,9}: the 2 smallest are {1,5}; ties at 5 taken whole.
        let mut b = b;
        b.sort_unstable();
        assert_eq!(b, vec![1, 2, 4]);
        assert_eq!(bi, 9.0);
        let (b, bi) = q.pull();
        assert_eq!(b, vec![3]);
        assert_eq!(bi, 100.0);
        assert!(q.is_empty());
    }

    /// Model-based differential test: BlockQueue vs a simple sorted-vector
    /// model under randomized ops, including heavy ties.
    #[test]
    fn block_queue_matches_model() {
        use rand::Rng;
        use rand::SeedableRng;
        use rand_chacha::ChaCha8Rng;
        let mut rng = ChaCha8Rng::seed_from_u64(0xB10C);
        for &m in &[1usize, 2, 3, 8] {
            for bound in [50.0f64, 1000.0] {
                let mut q = BlockQueue::new(bound, m);
                let mut model: Vec<(u32, f64)> = Vec::new();
                let mut log: Vec<String> = Vec::new();
                for _ in 0..4000 {
                    match rng.gen_range(0..4u32) {
                        0 => {
                            let v = rng.gen_range(0..12u32);
                            let k = (rng.gen_range(0..40) as f64) * 0.25;
                            q.insert(v, k);
                            model.push((v, k));
                            log.push(format!("I {v} {k}"));
                        }
                        1 => {
                            let n = rng.gen_range(0..=5usize);
                            let items: Vec<(u32, f64)> = (0..n)
                                .map(|_| {
                                    let v = rng.gen_range(0..12u32);
                                    let k = (rng.gen_range(0..40) as f64) * 0.25;
                                    (v, k)
                                })
                                .collect();
                            q.batch_prepend(&items);
                            model.extend_from_slice(&items);
                            log.push(format!("B {items:?}"));
                        }
                        2 => {
                            type Snapshot = Vec<(UbKey, Vec<(u32, f64)>, f64)>;
                            let q_pre: Snapshot = q
                                .blocks
                                .iter()
                                .map(|(k, b)| (*k, b.items.clone(), b.ub))
                                .collect();
                            let m_pre = model.clone();
                            let (sb, xb) = q.pull();
                            let (sm, xm) = model_pull(&mut model, m, bound);
                            let mut sb = sb;
                            sb.sort_unstable();
                            let mut sm = sm;
                            sm.sort_unstable();
                            assert_eq!(
                                sb, sm,
                                "bucket mismatch m={m}\nq_pre={q_pre:?}\nmodel_pre={m_pre:?}\nops={}",
                                log.join(" ")
                            );
                            assert_eq!(xb, xm, "separation mismatch m={m}\nqueue={:?}", q);
                            assert_eq!(q.is_empty(), model.is_empty());
                            log.push(format!("P {:?}", sm));
                        }
                        _ => {
                            assert_eq!(q.is_empty(), model.is_empty());
                        }
                    }
                    // invariant: all items in block i <= all items in block i+1
                    {
                        let q_count: usize = q.blocks.values().map(|b| b.items.len()).sum();
                        assert_eq!(
                            q_count,
                            model.len(),
                            "item count mismatch m={m}\nqueue={q:?}\nmodel={model:?}\nops={}",
                            log.join(" ")
                        );
                        let keys: Vec<UbKey> = q.blocks.keys().copied().collect();
                        for w in keys.windows(2) {
                            let lo_max = q.blocks[&w[0]]
                                .items
                                .iter()
                                .map(|it| it.1)
                                .fold(f64::NEG_INFINITY, f64::max);
                            let hi_min = q.blocks[&w[1]]
                                .items
                                .iter()
                                .map(|it| it.1)
                                .fold(f64::INFINITY, f64::min);
                            assert!(
                                lo_max <= hi_min,
                                "order invariant broken m={m}\nqueue={:?}\nops={}",
                                q,
                                log.join(" ")
                            );
                        }
                    }
                }
            }
        }
    }

    fn model_pull(items: &mut Vec<(u32, f64)>, m: usize, bound: f64) -> (Vec<u32>, f64) {
        items.sort_by(|a, b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
        if items.len() <= m {
            let vs = items.iter().map(|&(v, _)| v).collect();
            items.clear();
            return (vs, bound);
        }
        let vm = items[m - 1].1;
        let take = items.partition_point(|it| it.1 <= vm);
        let s: Vec<u32> = items[..take].iter().map(|&(v, _)| v).collect();
        items.drain(..take);
        let x = if items.is_empty() { bound } else { items[0].1 };
        (s, x)
    }
}
