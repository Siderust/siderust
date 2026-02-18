// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Shared build-time pipeline for JPL DE4xx ephemerides.
//!
//! Both DE440 and DE441 share identical processing logic. This module
//! parameterizes the pipeline over a [`DeConfig`] so each version is
//! a thin configuration wrapper.

use super::{daf, spk};
use std::env;
use std::path::Path;
use std::process::Command;

/// NAIF body IDs for the segments we extract.
const SUN_ID: i32 = 10;
const EMB_ID: i32 = 3;
const MOON_ID: i32 = 301;

/// Bodies to extract: (target_id, center_id, output_name)
const BODIES: &[(i32, i32, &str)] = &[
    (SUN_ID, 0, "sun"),
    (EMB_ID, 0, "emb"),
    (MOON_ID, EMB_ID, "moon"),
];

/// Configuration for a specific DE version (DE440, DE441, etc.).
pub struct DeConfig {
    /// File prefix for binary output, e.g. `"de440"` or `"de441"`.
    pub prefix: &'static str,
    /// Human-readable label, e.g. `"DE440"` or `"DE441"`.
    pub label: &'static str,
    /// URL to download the BSP file from NAIF.
    pub bsp_url: &'static str,
    /// BSP filename on disk, e.g. `"de440.bsp"`.
    pub bsp_filename: &'static str,
    /// Minimum plausible BSP file size in bytes (for validation).
    pub min_bsp_size: u64,
    /// Download timeout in seconds (curl/wget).
    pub download_timeout: u64,
    /// Human-readable size hint for log messages, e.g. `"~120 MB"`.
    pub size_hint: &'static str,
    /// Relative repo directory for Git LFS fallback, e.g. `"scripts/de440/dataset"`.
    pub repo_subdir: &'static str,
}

fn jpl_stub_enabled(cfg: &DeConfig) -> bool {
    // `SIDERUST_JPL_STUB`:
    // - unset/empty/0/false: normal behavior (download + parse + codegen)
    // - "1"/"true"/"all": stub all JPL datasets
    // - "de441" or "de440,de441": stub only those prefixes (case-insensitive)
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
        .any(|tok| tok == cfg.prefix.to_ascii_lowercase())
}

/// Run the full DE4xx build pipeline: fetch → parse DAF → extract → codegen.
pub fn run(cfg: &DeConfig, data_dir: &Path) -> anyhow::Result<()> {
    let out_dir =
        std::path::PathBuf::from(std::env::var("OUT_DIR").expect("OUT_DIR not set by Cargo"));

    if jpl_stub_enabled(cfg) {
        println!(
            "cargo:warning=Stubbed {} build (SIDERUST_JPL_STUB set): skipping BSP download/codegen; runtime use will panic",
            cfg.label
        );
        let rs_path = out_dir.join(format!("{}_data.rs", cfg.prefix));
        generate_stub_rust_module(cfg, &rs_path)?;
        eprintln!("  → {} (stub)", rs_path.display());
        return Ok(());
    }

    // 1. Fetch
    let bsp_path = ensure_bsp(cfg, data_dir)?;
    eprintln!("  {} BSP at: {}", cfg.label, bsp_path.display());

    // 2. Parse DAF
    let file_data = std::fs::read(&bsp_path)?;
    let daf = daf::Daf::parse(&file_data)?;
    eprintln!(
        "  DAF parsed: {} summaries, ND={}, NI={}",
        daf.summaries.len(),
        daf.nd,
        daf.ni
    );

    // 3. Extract bodies → binary files + codegen
    let mut generated_bodies: Vec<(&str, spk::SegmentMeta)> = Vec::new();

    for &(target, center, name) in BODIES {
        let summary = daf
            .summaries
            .iter()
            .find(|s| s.target_id == target && s.center_id == center)
            .unwrap_or_else(|| {
                panic!(
                    "{}: segment target={} center={} not found in BSP",
                    cfg.label, target, center
                )
            });

        eprintln!(
            "  Segment {}: target={}, center={}, type={}, records at words {}..{}",
            name,
            summary.target_id,
            summary.center_id,
            summary.data_type,
            summary.start_word,
            summary.end_word
        );

        assert!(
            summary.data_type == 2 || summary.data_type == 3,
            "{}: only SPK Type 2/3 supported, got Type {}",
            cfg.label,
            summary.data_type
        );

        let meta = spk::read_type2_segment(&file_data, &daf, summary)?;
        eprintln!(
            "    ncoeff={}, n_records={}, intlen={:.1}s ({:.1} days), rsize={}",
            meta.ncoeff,
            meta.n_records,
            meta.intlen,
            meta.intlen / 86400.0,
            meta.rsize
        );

        // Write binary data file
        let bin_path = out_dir.join(format!("{}_{}.bin", cfg.prefix, name));
        spk::write_binary(&meta, &bin_path)?;
        eprintln!(
            "    → {} ({} bytes)",
            bin_path.display(),
            std::fs::metadata(&bin_path)?.len()
        );

        generated_bodies.push((name, meta));
    }

    // 4. Generate Rust accessor module
    let rs_path = out_dir.join(format!("{}_data.rs", cfg.prefix));
    generate_rust_module(cfg, &generated_bodies, &rs_path)?;
    eprintln!("  → {}", rs_path.display());

    Ok(())
}

fn generate_stub_rust_module(cfg: &DeConfig, path: &Path) -> anyhow::Result<()> {
    use std::io::Write;

    let mut f = std::fs::File::create(path)?;
    writeln!(f, "// AUTOGENERATED by build.rs — do not edit")?;
    writeln!(
        f,
        "// STUBBED {} dataset: SIDERUST_JPL_STUB requested no-download build.",
        cfg.label
    )?;
    writeln!(f, "// Any attempt to evaluate ephemerides will panic.")?;
    writeln!(f)?;

    writeln!(
        f,
        "#[inline(never)]\nfn unavailable_record(_: usize) -> &'static [f64] {{\n    panic!(\"{} dataset is stubbed (SIDERUST_JPL_STUB). Disable stubbing or provide the BSP to build real data.\")\n}}",
        cfg.label
    )?;
    writeln!(f)?;

    for name in ["sun", "emb", "moon"] {
        writeln!(f, "pub mod {} {{", name)?;
        writeln!(f, "    pub const INIT: f64 = 0.0;")?;
        writeln!(f, "    pub const INTLEN: f64 = 1.0;")?;
        writeln!(f, "    pub const NCOEFF: usize = 0;")?;
        writeln!(f, "    pub const RSIZE: usize = 2;")?;
        writeln!(f, "    pub const N_RECORDS: usize = 1;")?;
        writeln!(f)?;
        writeln!(
            f,
            "    #[inline(never)]\n    pub fn record(i: usize) -> &'static [f64] {{\n        let _ = i;\n        super::unavailable_record(i)\n    }}",
        )?;
        writeln!(f, "}}")?;
        writeln!(f)?;
    }

    writeln!(
        f,
        "/// Pre-built segment descriptors for the {} bodies.",
        cfg.label
    )?;
    writeln!(
        f,
        "/// These require `SegmentDescriptor` to be in scope via `use`."
    )?;
    writeln!(f)?;

    writeln!(f, "pub const SUN: SegmentDescriptor = SegmentDescriptor {{")?;
    writeln!(f, "    init: qtty::Seconds::new(sun::INIT),")?;
    writeln!(f, "    intlen: qtty::Seconds::new(sun::INTLEN),")?;
    writeln!(f, "    ncoeff: sun::NCOEFF,")?;
    writeln!(f, "    n_records: sun::N_RECORDS,")?;
    writeln!(f, "    record_fn: sun::record,")?;
    writeln!(f, "}};")?;
    writeln!(f)?;

    writeln!(f, "pub const EMB: SegmentDescriptor = SegmentDescriptor {{")?;
    writeln!(f, "    init: qtty::Seconds::new(emb::INIT),")?;
    writeln!(f, "    intlen: qtty::Seconds::new(emb::INTLEN),")?;
    writeln!(f, "    ncoeff: emb::NCOEFF,")?;
    writeln!(f, "    n_records: emb::N_RECORDS,")?;
    writeln!(f, "    record_fn: emb::record,")?;
    writeln!(f, "}};")?;
    writeln!(f)?;

    writeln!(
        f,
        "pub const MOON: SegmentDescriptor = SegmentDescriptor {{"
    )?;
    writeln!(f, "    init: qtty::Seconds::new(moon::INIT),")?;
    writeln!(f, "    intlen: qtty::Seconds::new(moon::INTLEN),")?;
    writeln!(f, "    ncoeff: moon::NCOEFF,")?;
    writeln!(f, "    n_records: moon::N_RECORDS,")?;
    writeln!(f, "    record_fn: moon::record,")?;
    writeln!(f, "}};")?;
    writeln!(f)?;

    Ok(())
}

/// Generate the Rust source that includes the binary data and exposes typed accessors.
fn generate_rust_module(
    cfg: &DeConfig,
    bodies: &[(&str, spk::SegmentMeta)],
    path: &Path,
) -> anyhow::Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "// AUTOGENERATED by build.rs — do not edit")?;
    writeln!(f, "// {} embedded Chebyshev coefficient data", cfg.label)?;
    writeln!(f)?;

    for (name, meta) in bodies {
        let upper = name.to_uppercase();
        let byte_count = meta.n_records * meta.rsize * 8;

        writeln!(f, "/// {} Chebyshev data for body `{}`.", cfg.label, name)?;
        writeln!(f, "pub mod {} {{", name)?;
        writeln!(f, "    /// Initial epoch (TDB seconds past J2000).")?;
        writeln!(f, "    pub const INIT: f64 = {:?};", meta.init)?;
        writeln!(f, "    /// Interval length (seconds).")?;
        writeln!(f, "    pub const INTLEN: f64 = {:?};", meta.intlen)?;
        writeln!(
            f,
            "    /// Number of Chebyshev coefficients per coordinate."
        )?;
        writeln!(f, "    pub const NCOEFF: usize = {};", meta.ncoeff)?;
        writeln!(f, "    /// Doubles per record (2 + 3 * NCOEFF).")?;
        writeln!(f, "    pub const RSIZE: usize = {};", meta.rsize)?;
        writeln!(f, "    /// Number of records.")?;
        writeln!(f, "    pub const N_RECORDS: usize = {};", meta.n_records)?;
        writeln!(f)?;
        writeln!(f, "    /// Total byte count of the binary data.")?;
        writeln!(f, "    const BYTE_COUNT: usize = {};", byte_count)?;
        writeln!(f)?;
        writeln!(
            f,
            "    /// 8-byte aligned wrapper for `include_bytes!` data."
        )?;
        writeln!(f, "    #[repr(C, align(8))]")?;
        writeln!(f, "    struct Aligned([u8; BYTE_COUNT]);")?;
        writeln!(f)?;
        writeln!(
            f,
            "    /// Raw coefficient data, embedded with 8-byte alignment."
        )?;
        writeln!(
            f,
            "    static DATA: &Aligned = &Aligned(*include_bytes!(concat!(env!(\"OUT_DIR\"), \"/{}_{}.bin\")));",
            cfg.prefix, name
        )?;
        writeln!(f)?;
        writeln!(
            f,
            "    /// Return the `i`-th record as a slice of `RSIZE` f64 values."
        )?;
        writeln!(f, "    #[inline]")?;
        writeln!(f, "    pub fn record(i: usize) -> &'static [f64] {{")?;
        writeln!(
            f,
            "        debug_assert!(i < N_RECORDS, \"{} {}: record index {{}} out of range (max {{}})\", i, N_RECORDS - 1);",
            cfg.label, upper
        )?;
        writeln!(f, "        let byte_offset = i * RSIZE * 8;")?;
        writeln!(f, "        let ptr = DATA.0.as_ptr();")?;
        writeln!(
            f,
            "        // SAFETY: `DATA` is #[repr(align(8))] so `ptr` is 8-byte aligned."
        )?;
        writeln!(
            f,
            "        // Each record starts at a multiple of RSIZE*8 bytes (also aligned)."
        )?;
        writeln!(
            f,
            "        // The binary data was written as native-endian f64 by the build script."
        )?;
        writeln!(f, "        unsafe {{")?;
        writeln!(
            f,
            "            std::slice::from_raw_parts(ptr.add(byte_offset) as *const f64, RSIZE)"
        )?;
        writeln!(f, "        }}")?;
        writeln!(f, "    }}")?;
        writeln!(f, "}}")?;
        writeln!(f)?;
    }

    // Emit SegmentDescriptor constants so the runtime data.rs needs no macro
    writeln!(
        f,
        "/// Pre-built segment descriptors for the {} bodies.",
        cfg.label
    )?;
    writeln!(
        f,
        "/// These require `SegmentDescriptor` to be in scope via `use`."
    )?;
    writeln!(f)?;
    writeln!(f, "/// Segment descriptor for the Sun (NAIF 10 → SSB).")?;
    writeln!(f, "pub const SUN: SegmentDescriptor = SegmentDescriptor {{")?;
    writeln!(f, "    init: qtty::Seconds::new(sun::INIT),")?;
    writeln!(f, "    intlen: qtty::Seconds::new(sun::INTLEN),")?;
    writeln!(f, "    ncoeff: sun::NCOEFF,")?;
    writeln!(f, "    n_records: sun::N_RECORDS,")?;
    writeln!(f, "    record_fn: sun::record,")?;
    writeln!(f, "}};")?;
    writeln!(f)?;
    writeln!(
        f,
        "/// Segment descriptor for the Earth-Moon Barycenter (NAIF 3 → SSB)."
    )?;
    writeln!(f, "pub const EMB: SegmentDescriptor = SegmentDescriptor {{")?;
    writeln!(f, "    init: qtty::Seconds::new(emb::INIT),")?;
    writeln!(f, "    intlen: qtty::Seconds::new(emb::INTLEN),")?;
    writeln!(f, "    ncoeff: emb::NCOEFF,")?;
    writeln!(f, "    n_records: emb::N_RECORDS,")?;
    writeln!(f, "    record_fn: emb::record,")?;
    writeln!(f, "}};")?;
    writeln!(f)?;
    writeln!(f, "/// Segment descriptor for the Moon (NAIF 301 → EMB).")?;
    writeln!(
        f,
        "pub const MOON: SegmentDescriptor = SegmentDescriptor {{"
    )?;
    writeln!(f, "    init: qtty::Seconds::new(moon::INIT),")?;
    writeln!(f, "    intlen: qtty::Seconds::new(moon::INTLEN),")?;
    writeln!(f, "    ncoeff: moon::NCOEFF,")?;
    writeln!(f, "    n_records: moon::N_RECORDS,")?;
    writeln!(f, "    record_fn: moon::record,")?;
    writeln!(f, "}};")?;
    writeln!(f)?;

    Ok(())
}

// ── BSP acquisition ─────────────────────────────────────────────────────

/// Ensure the BSP file exists in `data_dir`, downloading if necessary.
fn ensure_bsp(cfg: &DeConfig, data_dir: &Path) -> anyhow::Result<std::path::PathBuf> {
    std::fs::create_dir_all(data_dir)?;
    let bsp_path = data_dir.join(cfg.bsp_filename);

    // 1) Early return if already present with valid size
    if bsp_path.exists() {
        let meta = std::fs::metadata(&bsp_path)?;
        if meta.len() > cfg.min_bsp_size {
            eprintln!(
                "  {} BSP already cached ({:.1} MB)",
                cfg.label,
                meta.len() as f64 / 1_000_000.0
            );
            return Ok(bsp_path);
        }
        eprintln!(
            "  {} BSP exists but too small ({} B), re-acquiring...",
            cfg.label,
            meta.len()
        );
        std::fs::remove_file(&bsp_path)?;
    }

    // 2) Try local checkout (Git LFS)
    if let Ok(true) = copy_from_repo(cfg, &bsp_path) {
        return Ok(bsp_path);
    }

    // 3) Download from NAIF as fallback
    eprintln!(
        "  Downloading {} BSP from NAIF ({})...",
        cfg.label, cfg.size_hint
    );
    eprintln!("  URL: {}", cfg.bsp_url);

    let result = download_with_curl(cfg, &bsp_path).or_else(|_| download_with_wget(cfg, &bsp_path));

    match result {
        Ok(()) => {
            let size = std::fs::metadata(&bsp_path)?.len();
            if size < cfg.min_bsp_size {
                anyhow::bail!(
                    "Downloaded file too small ({} bytes), expected >{}",
                    size,
                    cfg.min_bsp_size
                );
            }
            eprintln!("  Downloaded {:.1} MB", size as f64 / 1_000_000.0);
            Ok(bsp_path)
        }
        Err(e) => {
            anyhow::bail!(
                "Failed to download {} BSP. Ensure `curl` or `wget` is installed.\n\
                 You can also manually download:\n  {}\n\
                 and place it at:\n  {}\n\
                 Error: {}",
                cfg.label,
                cfg.bsp_url,
                bsp_path.display(),
                e
            );
        }
    }
}

/// Copy BSP from the repository if available (Git LFS fallback).
fn copy_from_repo(cfg: &DeConfig, dst: &Path) -> anyhow::Result<bool> {
    let src = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(cfg.repo_subdir)
        .join(cfg.bsp_filename);

    eprintln!(
        "  Attempting to copy {} BSP from Git LFS: {:?}",
        cfg.label, src
    );

    if !src.exists() {
        eprintln!("  No local BSP found in source tree");
        return Ok(false);
    }

    let meta = std::fs::metadata(&src)?;
    if meta.len() < cfg.min_bsp_size {
        eprintln!("  Local BSP file too small ({} B), skipping", meta.len());
        return Ok(false);
    }

    eprintln!(
        "  Found valid BSP in source ({:.1} MB), copying...",
        meta.len() as f64 / 1_000_000.0
    );
    std::fs::copy(&src, dst)?;
    eprintln!("  Successfully copied from Git LFS");
    Ok(true)
}

fn download_with_curl(cfg: &DeConfig, dest: &Path) -> anyhow::Result<()> {
    eprintln!("  Trying curl...");
    let timeout_str = cfg.download_timeout.to_string();
    let status = Command::new("curl")
        .args([
            "--fail",
            "--silent",
            "--show-error",
            "--location",
            "--max-time",
            &timeout_str,
            "--output",
        ])
        .arg(dest.as_os_str())
        .arg(cfg.bsp_url)
        .status()?;

    if !status.success() {
        anyhow::bail!("curl exited with status {}", status);
    }
    Ok(())
}

fn download_with_wget(cfg: &DeConfig, dest: &Path) -> anyhow::Result<()> {
    eprintln!("  Trying wget...");
    let timeout_arg = format!("--timeout={}", cfg.download_timeout);
    let status = Command::new("wget")
        .args(["--quiet", &timeout_arg, "--output-document"])
        .arg(dest.as_os_str())
        .arg(cfg.bsp_url)
        .status()?;

    if !status.success() {
        anyhow::bail!("wget exited with status {}", status);
    }
    Ok(())
}
