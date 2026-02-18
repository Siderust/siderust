// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

#[path = "scripts/vsop87/mod.rs"]
mod vsop87_build;

#[path = "scripts/elp2000/mod.rs"]
mod elp2000_build;

#[path = "scripts/iers/mod.rs"]
mod iers_build;

#[cfg(feature = "de440")]
#[path = "scripts/jpl/de440/mod.rs"]
mod de440_build;

#[cfg(feature = "de441")]
#[path = "scripts/jpl/de441/mod.rs"]
mod de441_build;

use std::{env, path::PathBuf};

fn stub_enabled_for(prefix: &str) -> bool {
    let Ok(raw) = env::var("SIDERUST_JPL_STUB") else {
        return false;
    };
    let raw = raw.trim();
    if raw.is_empty() {
        return false;
    }

    let lower = raw.to_ascii_lowercase();
    if lower == "all" || lower == "1" || lower == "true" || lower == "yes" || lower == "on" {
        return true;
    }

    lower
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|s| !s.is_empty())
        .any(|tok| tok == prefix)
}

fn main() {
    println!("cargo:rerun-if-env-changed=SIDERUST_JPL_STUB");
    println!("cargo:rerun-if-env-changed=SIDERUST_DATASETS_DIR");
    println!("cargo:rustc-check-cfg=cfg(siderust_mock_de441)");

    if stub_enabled_for("de441") {
        println!("cargo:rustc-cfg=siderust_mock_de441");
        println!(
            "cargo:warning=Using DE441 mock ephemeris backend (VSOP87/ELP2000) because SIDERUST_JPL_STUB includes de441; real-DE441 tests are skipped"
        );
    }

    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR not set by Cargo"));
    let datasets_base_dir = env::var_os("SIDERUST_DATASETS_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| out_dir.clone());

    // VSOP87
    eprintln!("Building VSOP87 data...");
    let data_dir = datasets_base_dir.join("vsop87_dataset");
    vsop87_build::run(data_dir.as_path()).unwrap_or_else(|e| {
        panic!("VSOP87 codegen failed: {}", e);
    });
    eprintln!("VSOP87 data generation complete");

    // ELP2000
    eprintln!("Building ELP2000 data...");
    let elp_dir = datasets_base_dir.join("elp2000_dataset");
    elp2000_build::run(elp_dir.as_path()).unwrap_or_else(|e| {
        panic!("ELP2000 codegen failed: {}", e);
    });
    eprintln!("ELP2000 data generation complete");

    // IERS EOP (finals2000A.all)
    eprintln!("Building IERS EOP data...");
    let iers_dir = datasets_base_dir.join("iers_dataset");
    iers_build::run(iers_dir.as_path()).unwrap_or_else(|e| {
        panic!("IERS EOP codegen failed: {}", e);
    });
    eprintln!("IERS EOP data generation complete");

    // DE440 (only with feature)
    #[cfg(feature = "de440")]
    {
        eprintln!("Building DE440 data...");
        let de440_dir = datasets_base_dir.join("de440_dataset");
        de440_build::run(de440_dir.as_path()).unwrap_or_else(|e| {
            panic!("DE440 codegen failed: {}", e);
        });
        eprintln!("DE440 data generation complete");
    }

    // DE441 (only with feature)
    #[cfg(feature = "de441")]
    {
        eprintln!("Building DE441 data...");
        let de441_dir = datasets_base_dir.join("de441_dataset");
        de441_build::run(de441_dir.as_path()).unwrap_or_else(|e| {
            panic!("DE441 codegen failed: {}", e);
        });
        eprintln!("DE441 data generation complete");
    }
}
