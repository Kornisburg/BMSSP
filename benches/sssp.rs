// cargo bench harness (custom, no criterion): `cargo bench`.
#[path = "../src/bin/bench_sssp.rs"]
mod bench_sssp;

fn main() {
    let iters = std::env::var("BMSSP_BENCH_ITERS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(2);
    let ab = bench_sssp::ablation_from_env();
    println!("cargo bench: iters={iters} ablation={ab:?}\n");
    print!("{}", bench_sssp::run_all(iters, &ab));
}
