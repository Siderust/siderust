//! **VSOP87 build‑pipeline entry point**
//!
//! This is the public façade that `build.rs` sees.  All heavy lifting is
//! delegated to the three sub‑modules:
//!
//! | module   | responsibility |
//! |----------|---------------|
//! | [`collect`] | Walk data files, parse them, build the in‑memory [`VersionMap`]. |
//! | [`codegen`] | Turn the map into Rust source code (`String`). |
//! | [`io`]      | Write each `(version → source)` pair into `OUT_DIR`. |
//!
//! The **only** symbol meant to be used externally is [`run`]; everything else
//! stays private to keep the build script’s surface minimal.
//!
//! ```rust
//! // build.rs
//! #[path = "scripts/vsop87/mod.rs"]
//! mod vsop87_build;
//!
//! fn main() {
//!     vsop87_build::run("dataset").unwrap();
//! }
//! ```
//!
//! ---
//! ## Data model in memory
//!
//! The VSOP87 coefficients form a 4‑level hierarchy:
//!
//! ```text
//! version  (char: 'A', 'E', …)
//! └── planet  (String: "EARTH", "MARS", …)
//!     └── coord   (u8: 1 = X, 2 = Y, 3 = Z)
//!         └── T‑power (u8: 0…5) → Vec<Term>
//! ```
//!
//! A single [`Term`] is the triple `(a, b, c)` representing
//! `a · cos(b + c·T)`.
//! The aliases below (`VersionMap`, `PlanetMap`, …) make the nested
//! `BTreeMap` layers bearable in signatures.
//!
//! ---
//! ## Compile‑time vs run‑time
//!
//! All of this runs **at compile time** inside the Cargo build script.  The
//! generated `.rs` files are then **included** by normal code via `include!` or
//! by the compiler picking them up as modules.

use std::{
    collections::BTreeMap,
    env,
    path::PathBuf,
};

use anyhow::Context;

mod collect;
mod codegen;
mod io;

// ---------------------------------------------------------------------------
// Core domain types
// ---------------------------------------------------------------------------

/// A single VSOP87 term – coefficient triplet for `a · cos(b + c·T)`.
#[derive(Clone, Copy, Debug)]
struct Term {
    a: f64,
    b: f64,
    c: f64,
}

// Readability aliases for the deeply‑nested map structure.
/// coordinate → T‑power → Vec<Term>
type CoordMap   = BTreeMap<u8, TPowerMap>;
/// T‑power → Vec<Term>
type TPowerMap  = BTreeMap<u8, Vec<Term>>;
/// planet → CoordMap
type PlanetMap  = BTreeMap<String, CoordMap>;
/// version → PlanetMap
type VersionMap = BTreeMap<char, PlanetMap>;

// ===========================================================================
// Public API – the single entry for build.rs
// ===========================================================================

/// Assemble the **entire** VSOP87 build pipeline.
///
/// 1. **Tell Cargo** to rerun if anything under `data_dir` changes.
/// 2. **Locate** `OUT_DIR` (where build scripts must emit artefacts).
/// 3. **collect → codegen → io** pipeline.
///
/// # Errors
/// Propagates any I/O or parsing failure wrapped in `anyhow::Error`.
pub fn run(data_dir: &str) -> anyhow::Result<()> {
    // 1) Watch data directory for changes so Cargo triggers a rebuild.
    let data_dir = PathBuf::from(data_dir);
    println!("cargo:rerun-if-changed={}", data_dir.display());

    // 2) Resolve OUT_DIR – provided by Cargo at build‑time.
    let out_dir = PathBuf::from(
        env::var("OUT_DIR").context("OUT_DIR not set (missing Cargo build context)")?,
    );

    // 3) Pipeline: parse → generate code → write files.
    let versions = collect::collect_terms(&data_dir)?;
    let modules  = codegen::generate_modules(&versions)?;
    io::write_modules(&modules, &out_dir)?;

    Ok(())
}
