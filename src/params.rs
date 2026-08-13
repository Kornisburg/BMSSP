pub const INF: f64 = f64::INFINITY;

/// k = floor((log2 n)^(1/3)), t = floor((log2 n)^(2/3)), l = ceil(log2 n / t).
///
/// Base-2 logs ensure `k * 2^(l*t) >= n`, so the top-level BMSSP call always
/// runs to completion (B' = B = infinity) and settles every reachable vertex.
pub fn params(n: usize) -> (usize, usize, usize) {
    let n = n.max(2);
    let log2n = (n as f64).log2().max(1.0);
    let mut k = log2n.powf(1.0 / 3.0).floor() as usize;
    let mut t = log2n.powf(2.0 / 3.0).floor() as usize;
    let mut l = (log2n / t as f64).ceil() as usize;
    if k < 1 {
        k = 1;
    }
    if t < 1 {
        t = 1;
    }
    if l < 1 {
        l = 1;
    }
    (k, t, l)
}

#[cfg(test)]
mod tests {
    use super::params;

    #[test]
    fn params_are_sane_for_small_n() {
        for n in [1, 2, 3, 5, 10, 100, 5000, 100_000, 10_000_000usize] {
            let (k, t, l) = params(n);
            assert!(k >= 1 && t >= 1 && l >= 1, "n={n}: k={k} t={t} l={l}");
            // top-level bound must dominate n
            let cap = k.saturating_mul(1usize.checked_shl((l * t).min(63) as u32).unwrap_or(usize::MAX));
            assert!(cap >= n, "n={n}: cap {cap} < n");
        }
    }

    #[test]
    fn params_grow_weakly() {
        let a = params(10);
        let b = params(10_000_000);
        assert!(b.0 >= a.0 && b.1 >= a.1 && b.2 >= a.2);
    }
}
