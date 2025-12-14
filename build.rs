#[path = "scripts/vsop87/mod.rs"]
mod vsop87_build;

#[path = "scripts/elp2000/mod.rs"]
mod elp2000_build;

use std::{env, path::PathBuf};

fn main() {
    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR not set by Cargo"));

    // VSOP87
    let data_dir = out_dir.join("vsop87_dataset");
    vsop87_build::run(data_dir.as_path()).expect("VSOP87 codegen failed");

    // ELP2000
    let elp_dir = out_dir.join("elp2000_dataset");
    elp2000_build::run(elp_dir.as_path()).expect("ELP2000 codegen failed");
}
