// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

#[cfg(any(feature = "regen-data", feature = "de440"))]
use std::env;
#[cfg(any(feature = "regen-data", feature = "de440"))]
use std::path::PathBuf;

#[cfg(feature = "regen-data")]
#[path = "scripts/vsop87/mod.rs"]
mod vsop87_build;

#[cfg(feature = "regen-data")]
#[path = "scripts/elp2000/mod.rs"]
mod elp2000_build;

#[cfg(feature = "regen-data")]
#[path = "scripts/iers/mod.rs"]
mod iers_build;

#[cfg(feature = "de440")]
#[path = "scripts/jpl/de440/mod.rs"]
mod de440_build;

fn main() {
    println!("cargo:rerun-if-env-changed=SIDERUST_DATASETS_DIR");
    println!("cargo:rerun-if-env-changed=SIDERUST_REGEN");
    println!("cargo:rerun-if-env-changed=SIDERUST_JPL_STUB");
    println!("cargo:rustc-check-cfg=cfg(siderust_mock_de440)");

    #[cfg(feature = "regen-data")]
    regen_tables();

    #[cfg(not(feature = "regen-data"))]
    eprintln!(
        "siderust build: using committed generated tables in src/generated/. \
         Enable the `regen-data` feature and set SIDERUST_REGEN=1 to refresh them."
    );

    #[cfg(feature = "de440")]
    build_de440();
}

// ── Shared helpers ────────────────────────────────────────────────────────────

/// Returns the base directory for dataset storage.
/// Prefers `SIDERUST_DATASETS_DIR`; falls back to `OUT_DIR`.
#[cfg(any(feature = "regen-data", feature = "de440"))]
fn datasets_base_dir() -> PathBuf {
    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR not set by Cargo"));
    env::var_os("SIDERUST_DATASETS_DIR")
        .map(PathBuf::from)
        .unwrap_or(out_dir)
}

// ── `regen-data` feature ──────────────────────────────────────────────────────

/// Regenerates the committed VSOP87, ELP2000 and IERS tables from source data.
///
/// Only runs when `SIDERUST_REGEN=1` (or `true`/`yes`) is set. Requires the
/// `regen-data` build feature so that `reqwest` is compiled only when needed.
///
/// The preferred way to invoke this is via the helper script:
///   scripts/update_generated_tables.sh
#[cfg(feature = "regen-data")]
fn regen_tables() {
    let regen = env::var("SIDERUST_REGEN")
        .map(|v| {
            let v = v.trim().to_ascii_lowercase();
            v == "1" || v == "true" || v == "yes"
        })
        .unwrap_or(false);

    if !regen {
        eprintln!(
            "siderust build (regen-data): SIDERUST_REGEN not set — skipping table regeneration."
        );
        return;
    }

    let base = datasets_base_dir();
    let manifest_dir =
        PathBuf::from(env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR not set"));
    let gen_dir = manifest_dir.join("src/generated");

    eprintln!("Regenerating VSOP87 data...");
    vsop87_build::run_regen(base.join("vsop87_dataset").as_path(), &gen_dir)
        .unwrap_or_else(|e| panic!("VSOP87 codegen failed: {e}"));
    eprintln!("VSOP87 regeneration complete");

    eprintln!("Regenerating ELP2000 data...");
    elp2000_build::run_regen(base.join("elp2000_dataset").as_path(), &gen_dir)
        .unwrap_or_else(|e| panic!("ELP2000 codegen failed: {e}"));
    eprintln!("ELP2000 regeneration complete");

    eprintln!("Regenerating IERS EOP data...");
    iers_build::run_regen(base.join("iers_dataset").as_path(), &gen_dir)
        .unwrap_or_else(|e| panic!("IERS EOP codegen failed: {e}"));
    eprintln!("IERS EOP regeneration complete");
}

// ── `de440` feature ───────────────────────────────────────────────────────────

/// Builds DE440 coefficient data from the JPL BSP file.
///
/// When `SIDERUST_JPL_STUB=all` (or `de440`), emits `siderust_mock_de440` so
/// the library falls back to `Vsop87Ephemeris` without downloading the BSP.
#[cfg(feature = "de440")]
fn build_de440() {
    if jpl_stub_active() {
        println!("cargo:rustc-cfg=siderust_mock_de440");
    }
    let de440_dir = datasets_base_dir().join("de440_dataset");
    eprintln!("Building DE440 data...");
    de440_build::run(de440_dir.as_path()).unwrap_or_else(|e| panic!("DE440 codegen failed: {e}"));
    eprintln!("DE440 data generation complete");
}

/// Returns `true` when `SIDERUST_JPL_STUB` requests stubbing for DE440.
#[cfg(feature = "de440")]
fn jpl_stub_active() -> bool {
    let Ok(raw) = env::var("SIDERUST_JPL_STUB") else {
        return false;
    };
    let lower = raw.trim().to_ascii_lowercase();
    if lower.is_empty() {
        return false;
    }
    if matches!(lower.as_str(), "all" | "1" | "true" | "yes" | "on") {
        return true;
    }
    lower
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|s| !s.is_empty())
        .any(|tok| tok == "de440")
}
