#[derive(Debug, Clone, Copy, Default)]
pub struct Counters {
    pub relaxations: u64,
    pub heap_insert: u64,
    pub heap_extract_min: u64,
    pub recursive_calls: u64,
    pub find_pivots_calls: u64,
    pub base_case_calls: u64,
    pub queue_insert: u64,
    pub queue_pull: u64,
    pub queue_batch_prepend: u64,
    pub pivots_found: u64,
}

impl Counters {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn reset(&mut self) {
        *self = Self::default();
    }
}
