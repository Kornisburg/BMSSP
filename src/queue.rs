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
/// structure of the paper is the Phase-3 upgrade.
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

    pub fn len(&self) -> usize {
        self.map.values().map(Vec::len).sum()
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
}
