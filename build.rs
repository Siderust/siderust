#[path = "scripts/vsop87/mod.rs"]
mod vsop87_build;

use std::{env, path::PathBuf};

fn main() {
    let out_dir = PathBuf::from(
        env::var("OUT_DIR").expect("OUT_DIR not set by Cargo"),
    );

    let data_dir = out_dir.join("vsop87_dataset");

    // Download the vosp87 data and auto-generate the Rust source code.
    vsop87_build::run(data_dir.as_path())
        .expect("VSOP87 codegen failed");
}
