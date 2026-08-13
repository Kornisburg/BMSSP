// cargo bench harness (custom, no criterion): `cargo bench`.
#[path = "../src/bin/bench_sssp.rs"]
mod bench_sssp;

fn main() {
    let iters = std::env::var("BMSSP_BENCH_ITERS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(2);
    let use_pivots = std::env::var("BMSSP_NO_PIVOTS").is_err();
    println!("cargo bench: iters={iters} use_pivots={use_pivots}\n");
    print!("{}", bench_sssp::run_all(iters, use_pivots));
}
