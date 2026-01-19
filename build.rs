// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

#[path = "scripts/vsop87/mod.rs"]
mod vsop87_build;

#[path = "scripts/elp2000/mod.rs"]
mod elp2000_build;

use std::{env, path::PathBuf};

fn main() {
    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR not set by Cargo"));

    // VSOP87
    eprintln!("Building VSOP87 data...");
    let data_dir = out_dir.join("vsop87_dataset");
    vsop87_build::run(data_dir.as_path()).unwrap_or_else(|e| {
        panic!("VSOP87 codegen failed: {}", e);
    });
    eprintln!("VSOP87 data generation complete");

    // ELP2000
    eprintln!("Building ELP2000 data...");
    let elp_dir = out_dir.join("elp2000_dataset");
    elp2000_build::run(elp_dir.as_path()).unwrap_or_else(|e| {
        panic!("ELP2000 codegen failed: {}", e);
    });
    eprintln!("ELP2000 data generation complete");
}
