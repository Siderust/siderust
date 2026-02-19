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
    println!("cargo:rerun-if-env-changed=SIDERUST_REGEN");
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

    // VSOP87, ELP2000 and IERS tables are pre-generated and committed under
    // src/generated/.  Normal builds (including docs.rs) use those files
    // directly via `include!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/generated/…"))`.
    //
    // Set SIDERUST_REGEN=1 to re-run the full download + codegen pipeline and
    // overwrite the committed files in src/generated/.  The preferred way to
    // do this is via the helper script:
    //
    //   scripts/update_generated_tables.sh
    //
    let regen = env::var("SIDERUST_REGEN")
        .map(|v| {
            let v = v.trim().to_ascii_lowercase();
            v == "1" || v == "true" || v == "yes"
        })
        .unwrap_or(false);

    if regen {
        // Compute the path to src/generated/ relative to CARGO_MANIFEST_DIR.
        let manifest_dir =
            PathBuf::from(env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR not set"));
        let gen_dir = manifest_dir.join("src/generated");

        // VSOP87
        eprintln!("Regenerating VSOP87 data...");
        let data_dir = datasets_base_dir.join("vsop87_dataset");
        vsop87_build::run_regen(data_dir.as_path(), &gen_dir).unwrap_or_else(|e| {
            panic!("VSOP87 codegen failed: {}", e);
        });
        eprintln!("VSOP87 regeneration complete");

        // ELP2000
        eprintln!("Regenerating ELP2000 data...");
        let elp_dir = datasets_base_dir.join("elp2000_dataset");
        elp2000_build::run_regen(elp_dir.as_path(), &gen_dir).unwrap_or_else(|e| {
            panic!("ELP2000 codegen failed: {}", e);
        });
        eprintln!("ELP2000 regeneration complete");

        // IERS EOP (finals2000A.all)
        eprintln!("Regenerating IERS EOP data...");
        let iers_dir = datasets_base_dir.join("iers_dataset");
        iers_build::run_regen(iers_dir.as_path(), &gen_dir).unwrap_or_else(|e| {
            panic!("IERS EOP codegen failed: {}", e);
        });
        eprintln!("IERS EOP regeneration complete");
    } else {
        eprintln!(
            "siderust build: using committed generated tables in src/generated/. \
             Set SIDERUST_REGEN=1 to refresh them."
        );
    }

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
