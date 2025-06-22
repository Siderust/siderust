#[path = "scripts/vsop87/mod.rs"]
mod vsop87_build;

use std::{env, path::PathBuf};

fn main() {
    // OUT_DIR siempre existe y lo da Cargo:
    let out_dir = PathBuf::from(
        env::var("OUT_DIR").expect("OUT_DIR not set by Cargo"),
    );

    // Subcarpeta dentro de OUT_DIR para el dataset
    let data_dir = out_dir.join("vsop87_dataset");

    // run necesita &Path
    vsop87_build::run(data_dir.as_path())
        .expect("VSOP87 codegen failed");
}
