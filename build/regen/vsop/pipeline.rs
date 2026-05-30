// SPDX-License-Identifier: BSD-3-Clause
// Copyright (C) 2026 Vallés Puig, Ramon

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
//! #[path = "scripts/vsop87/mod.rs"]
//! mod vsop87_build;
//! vsop87_build::run("dataset").unwrap();
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

use std::{collections::BTreeMap, env, fs, path::Path, path::PathBuf};

use anyhow::Context;

#[path = "codegen.rs"]
mod codegen;
#[path = "collect.rs"]
mod collect;
#[path = "io.rs"]
mod io;

/// Compute a stable SHA-256 fingerprint of the VSOP87 input dataset
/// directory by hashing every `VSOP87X.xxx` file in sorted name order.
///
/// The resulting hex digest is logged during regeneration and can be
/// pinned via the `SIDERUST_VSOP87_SHA256` environment variable for
/// reproducible builds.
fn dataset_sha256(data_dir: &Path) -> anyhow::Result<(String, u64)> {
    use sha2::{Digest, Sha256};
    use std::fs;
    use std::io::Read;

    let mut paths: Vec<PathBuf> = fs::read_dir(data_dir)
        .with_context(|| format!("read-dir {data_dir:?}"))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .map(|s| s.starts_with("VSOP87"))
                .unwrap_or(false)
        })
        .collect();
    paths.sort();

    let mut hasher = Sha256::new();
    let mut total: u64 = 0;
    for path in &paths {
        // Include the file name in the hash so that file rearrangements
        // can't masquerade as identical content.
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            hasher.update(name.as_bytes());
            hasher.update(b"\0");
        }
        let mut f = std::fs::File::open(path).with_context(|| format!("open {path:?}"))?;
        let mut buf = [0u8; 64 * 1024];
        loop {
            let n = f.read(&mut buf).with_context(|| format!("read {path:?}"))?;
            if n == 0 {
                break;
            }
            total += n as u64;
            hasher.update(&buf[..n]);
        }
    }
    let digest = hasher.finalize();
    let mut s = String::with_capacity(64);
    use std::fmt::Write;
    for b in digest {
        let _ = write!(s, "{:02x}", b);
    }
    Ok((s, total))
}

fn check_pinned_sha256(actual: &str, env_var: &str) -> anyhow::Result<()> {
    if let Ok(expected) = env::var(env_var) {
        let expected = expected.trim().to_ascii_lowercase();
        if !expected.is_empty() && expected != actual {
            anyhow::bail!(
                "VSOP87: SHA-256 mismatch.\n  expected: {expected}\n  actual:   {actual}\n\
                 Refusing to regenerate tables from unverified data."
            );
        }
    }
    Ok(())
}

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
type CoordMap = BTreeMap<u8, TPowerMap>;
/// T‑power → Vec<Term>
type TPowerMap = BTreeMap<u8, Vec<Term>>;
/// planet → CoordMap
type PlanetMap = BTreeMap<String, CoordMap>;
/// version → PlanetMap
type VersionMap = BTreeMap<char, PlanetMap>;

/// Runs the complete VSOP87 build pipeline.
///
/// Steps:
/// 1. Instructs Cargo to rerun the build script if anything under `data_dir` changes.
/// 2. Locates `OUT_DIR` (the directory where build artifacts must be emitted).
/// 3. Parses all VSOP87 data files and builds the in-memory `VersionMap`.
/// 4. Generates Rust source code from the parsed data.
/// 5. Writes the generated code to files in `OUT_DIR`.
///
/// If the raw data directory is absent or empty, writes empty stub files and
/// emits a `cargo:warning`. Download the raw VSOP87 files into `src/vsop/raw/`
/// to regenerate the coefficient tables.
///
/// # Errors
/// Returns any I/O or parsing error wrapped in `anyhow::Error`.
#[allow(dead_code)]
pub(crate) fn run(data_dir: &Path) -> anyhow::Result<()> {
    let out_dir = PathBuf::from(
        env::var("OUT_DIR").context("OUT_DIR not set (missing Cargo build context)")?,
    );

    if !data_dir.exists() || is_empty_dir(data_dir)? {
        println!(
            "cargo:warning=VSOP87: raw data absent at {}; writing empty stubs. \
             Place VSOP87 data files in src/vsop/raw/ to enable table generation.",
            data_dir.display()
        );
        return write_empty_stubs(&out_dir);
    }

    let versions = collect::collect_terms(data_dir)?;
    let modules = codegen::generate_modules(&versions)?;
    io::write_modules(&modules, &out_dir)?;

    Ok(())
}

/// Like [`run`] but writes generated files to `gen_dir` instead of `OUT_DIR`.
///
/// Used by `build.rs` to emit Rust coefficient tables into `OUT_DIR`.
/// Verifies the input dataset SHA-256 against the `SIDERUST_VSOP87_SHA256`
/// environment variable when set, and embeds the SHA-256 + retrieval
/// timestamp as a comment header in every emitted `vsop87X.rs` file.
///
/// If the raw data directory is absent or empty, writes empty stub files and
/// emits a `cargo:warning`. No network access is performed; download the raw
/// VSOP87 files into `src/vsop/raw/` first.
pub(crate) fn run_regen(data_dir: &Path, gen_dir: &Path) -> anyhow::Result<()> {
    if !data_dir.exists() || is_empty_dir(data_dir)? {
        println!(
            "cargo:warning=VSOP87: raw data absent at {}; writing empty stubs. \
             Place VSOP87 data files in src/vsop/raw/ to enable table generation.",
            data_dir.display()
        );
        return write_empty_stubs(gen_dir);
    }

    let (sha, bytes) = dataset_sha256(data_dir)?;
    eprintln!("VSOP87: dataset SHA-256 = {sha} ({bytes} bytes across all files)");
    check_pinned_sha256(&sha, "SIDERUST_VSOP87_SHA256")?;

    let prov = codegen::DatasetProvenance {
        source_url: "https://ftp.imcce.fr/pub/ephem/planets/vsop87/".to_string(),
        sha256_hex: sha,
        retrieved_at: iso8601_now(),
        byte_count: bytes,
    };

    let versions = collect::collect_terms(data_dir)?;
    let modules = codegen::generate_modules_with_provenance(&versions, &prov)?;
    io::write_modules(&modules, gen_dir)?;
    Ok(())
}

/// Returns `true` when `dir` contains no files.
fn is_empty_dir(dir: &Path) -> anyhow::Result<bool> {
    Ok(fs::read_dir(dir)
        .with_context(|| format!("read-dir {dir:?}"))?
        .next()
        .is_none())
}

/// Write empty stub files for all modules that `include!` expects to find in
/// `OUT_DIR`. An empty file is valid Rust and produces an empty module body.
fn write_empty_stubs(gen_dir: &Path) -> anyhow::Result<()> {
    fs::create_dir_all(gen_dir).with_context(|| format!("create-dir {gen_dir:?}"))?;
    for name in ["vsop87a.rs", "vsop87e.rs"] {
        fs::write(
            gen_dir.join(name),
            "// No VSOP87 data: raw coefficient files absent from src/vsop/raw/.\n",
        )
        .with_context(|| format!("write stub {name}"))?;
    }
    Ok(())
}

fn iso8601_now() -> String {
    // Same lightweight formatter as in scripts/iers — duplicated here to
    // avoid threading an extra module through `#[path]` includes.
    use std::time::SystemTime;
    let secs = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = (secs / 86_400) as i64;
    let rem = (secs % 86_400) as u32;
    let hour = (rem / 3_600) as u8;
    let minute = ((rem % 3_600) / 60) as u8;
    let second = (rem % 60) as u8;
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u8;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u8;
    let year = (y + if m <= 2 { 1 } else { 0 }) as i32;
    format!("{year:04}-{m:02}-{d:02}T{hour:02}:{minute:02}:{second:02}Z")
}
