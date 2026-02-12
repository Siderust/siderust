// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

#[path = "scripts/vsop87/mod.rs"]
mod vsop87_build;

#[path = "scripts/elp2000/mod.rs"]
mod elp2000_build;

#[cfg(feature = "de440")]
#[path = "scripts/jpl/de440/mod.rs"]
mod de440_build;

#[cfg(feature = "de441")]
#[path = "scripts/jpl/de441/mod.rs"]
mod de441_build;

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

    // DE440 (only with feature)
    #[cfg(feature = "de440")]
    {
        eprintln!("Building DE440 data...");
        let de440_dir = out_dir.join("de440_dataset");
        de440_build::run(de440_dir.as_path()).unwrap_or_else(|e| {
            panic!("DE440 codegen failed: {}", e);
        });
        eprintln!("DE440 data generation complete");
    }

    // DE441 (only with feature)
    #[cfg(feature = "de441")]
    {
        eprintln!("Building DE441 data...");
        let de441_dir = out_dir.join("de441_dataset");
        de441_build::run(de441_dir.as_path()).unwrap_or_else(|e| {
            panic!("DE441 codegen failed: {}", e);
        });
        eprintln!("DE441 data generation complete");
    }
}
