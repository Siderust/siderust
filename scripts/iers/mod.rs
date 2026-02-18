// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Build-time pipeline for IERS Earth Orientation Parameters (EOP).
//!
//! Downloads and parses `finals2000A.all` from the IERS data center,
//! then generates a Rust source file with a static table of EOP entries
//! that gets baked into the binary at compile time.
//!
//! ## Data source
//!
//! The file `finals2000A.all` is the combined Bulletin A + B EOP series
//! maintained by the IERS. It contains daily values of:
//!
//! | Column  | Quantity                 | Unit          |
//! |---------|--------------------------|---------------|
//! | MJD     | Modified Julian Date     | days (UTC)    |
//! | xp, yp  | Pole coordinates         | arcseconds    |
//! | UT1−UTC | UT1 minus UTC            | seconds       |
//! | dX, dY  | CIP offsets (wrt 2000A)  | milliarcsec   |
//!
//! ## Pipeline
//!
//! 1. **Fetch** — try local `scripts/iers/dataset/`, then download from IERS
//! 2. **Parse** — fixed-width columns per USNO readme.finals2000A
//! 3. **Codegen** — emit `iers_eop_data.rs` with `pub static EOP_TABLE`
//!
//! ## References
//!
//! * Format: <https://maia.usno.navy.mil/ser7/readme.finals2000A>
//! * Data:   <https://datacenter.iers.org/data/csv/finals2000A.all>

use std::{
    env,
    fs::{self, File},
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result};
use std::fmt::Write;

// ═══════════════════════════════════════════════════════════════════════════
// Data types
// ═══════════════════════════════════════════════════════════════════════════

/// A single parsed row from `finals2000A.all`.
#[derive(Debug)]
#[allow(dead_code)]
struct EopRow {
    /// Modified Julian Date (UTC).
    mjd: f64,
    /// Pole x-coordinate (arcseconds).
    xp: f64,
    /// Pole y-coordinate (arcseconds).
    yp: f64,
    /// UT1 − UTC (seconds).
    dut1: f64,
    /// Celestial pole offset dX wrt IAU 2000A (milliarcseconds).
    dx: f64,
    /// Celestial pole offset dY wrt IAU 2000A (milliarcseconds).
    dy: f64,
    /// Whether this row is from Bulletin B (final) or A (preliminary/prediction).
    is_bulletin_b: bool,
}

// ═══════════════════════════════════════════════════════════════════════════
// Parser
// ═══════════════════════════════════════════════════════════════════════════

/// Parse a single fixed-width line from `finals2000A.all`.
///
/// Format (1-indexed columns):
/// ```text
///   8-15    F8.2    MJD (UTC)
///  17       A1      I/P flag for polar motion
///  19-27    F9.6    Bull. A PM-x (arcsec)
///  38-46    F9.6    Bull. A PM-y (arcsec)
///  58       A1      I/P flag for UT1-UTC
///  59-68    F10.7   Bull. A UT1-UTC (sec)
///  96       A1      I/P flag for nutation
///  98-106   F9.3    Bull. A dX (mas)
/// 117-125   F9.3    Bull. A dY (mas)
/// 135-144   F10.6   Bull. B PM-x (arcsec)
/// 145-154   F10.6   Bull. B PM-y (arcsec)
/// 155-165   F11.7   Bull. B UT1-UTC (sec)
/// 166-175   F10.3   Bull. B dX (mas)
/// 176-185   F10.3   Bull. B dY (mas)
/// ```
fn parse_line(line: &str) -> Option<EopRow> {
    // Lines shorter than 15 chars can't hold MJD
    if line.len() < 15 {
        return None;
    }

    // MJD: columns 8–15 (0-indexed: 7..15)
    let mjd: f64 = line.get(7..15)?.trim().parse().ok()?;

    // Skip rows with MJD < 37665.0 (before 1962-01-01, no UT1-UTC data)
    if mjd < 37665.0 {
        return None;
    }

    // Try Bulletin B first (final values, available for historical epochs):
    let (xp, yp, dut1, dx, dy, is_b) = if line.len() >= 185 {
        let b_xp: Option<f64> = line.get(134..144).and_then(|s| s.trim().parse().ok());
        let b_yp: Option<f64> = line.get(144..154).and_then(|s| s.trim().parse().ok());
        let b_dut1: Option<f64> = line.get(154..165).and_then(|s| s.trim().parse().ok());
        let b_dx: Option<f64> = line.get(165..175).and_then(|s| s.trim().parse().ok());
        let b_dy: Option<f64> = line.get(175..185).and_then(|s| s.trim().parse().ok());

        if let (Some(xp), Some(yp), Some(dut1)) = (b_xp, b_yp, b_dut1) {
            (xp, yp, dut1, b_dx.unwrap_or(0.0), b_dy.unwrap_or(0.0), true)
        } else {
            // Fall through to Bulletin A
            parse_bulletin_a(line)?
        }
    } else {
        parse_bulletin_a(line)?
    };

    Some(EopRow {
        mjd,
        xp,
        yp,
        dut1,
        dx,
        dy,
        is_bulletin_b: is_b,
    })
}

/// Parse Bulletin A columns for polar motion, UT1-UTC, and CIP offsets.
fn parse_bulletin_a(line: &str) -> Option<(f64, f64, f64, f64, f64, bool)> {
    // PM-x: cols 19–27 (0-indexed: 18..27)
    let xp: f64 = line.get(18..27)?.trim().parse().ok()?;
    // PM-y: cols 38–46 (0-indexed: 37..46)
    let yp: f64 = line.get(37..46)?.trim().parse().ok()?;

    // UT1-UTC: cols 59–68 (0-indexed: 58..68)
    let dut1: f64 = if line.len() >= 68 {
        line.get(58..68)?.trim().parse().ok()?
    } else {
        return None;
    };

    // dX: cols 98–106 (0-indexed: 97..106)
    let dx: f64 = if line.len() >= 106 {
        line.get(97..106)
            .and_then(|s| s.trim().parse().ok())
            .unwrap_or(0.0)
    } else {
        0.0
    };

    // dY: cols 117–125 (0-indexed: 116..125)
    let dy: f64 = if line.len() >= 125 {
        line.get(116..125)
            .and_then(|s| s.trim().parse().ok())
            .unwrap_or(0.0)
    } else {
        0.0
    };

    Some((xp, yp, dut1, dx, dy, false))
}

/// Parse all rows from the given `finals2000A.all` file.
fn parse_finals(path: &Path) -> Result<Vec<EopRow>> {
    let file = File::open(path).with_context(|| format!("open {path:?}"))?;
    let reader = BufReader::new(file);

    let rows: Vec<EopRow> = reader
        .lines()
        .map_while(Result::ok)
        .filter_map(|line| parse_line(&line))
        .collect();

    anyhow::ensure!(!rows.is_empty(), "No valid EOP rows parsed from {path:?}");
    Ok(rows)
}

// ═══════════════════════════════════════════════════════════════════════════
// Code generator
// ═══════════════════════════════════════════════════════════════════════════

fn generate_rust(rows: &[EopRow]) -> Result<String> {
    let mut out = String::with_capacity(rows.len() * 120);

    writeln!(
        out,
        "// ---------------------------------------------------"
    )?;
    writeln!(
        out,
        "// **AUTOGENERATED** by build.rs – DO NOT EDIT BY HAND"
    )?;
    writeln!(
        out,
        "// ---------------------------------------------------"
    )?;
    writeln!(out)?;
    writeln!(out, "/// A single daily EOP record from IERS finals2000A.")?;
    writeln!(out, "///")?;
    writeln!(
        out,
        "/// All angular quantities are in their native IERS units:"
    )?;
    writeln!(
        out,
        "/// pole coordinates in arcseconds, UT1−UTC in seconds,"
    )?;
    writeln!(out, "/// CIP offsets dX/dY in milliarcseconds.")?;
    writeln!(out, "#[derive(Debug, Clone, Copy)]")?;
    writeln!(out, "pub struct EopEntry {{")?;
    writeln!(out, "    /// Modified Julian Date (UTC).")?;
    writeln!(out, "    pub mjd: f64,")?;
    writeln!(out, "    /// Pole x-coordinate (arcseconds).")?;
    writeln!(out, "    pub xp: f64,")?;
    writeln!(out, "    /// Pole y-coordinate (arcseconds).")?;
    writeln!(out, "    pub yp: f64,")?;
    writeln!(out, "    /// UT1 − UTC (seconds).")?;
    writeln!(out, "    pub dut1: f64,")?;
    writeln!(
        out,
        "    /// Celestial pole offset dX wrt IAU 2000A (milliarcseconds)."
    )?;
    writeln!(out, "    pub dx: f64,")?;
    writeln!(
        out,
        "    /// Celestial pole offset dY wrt IAU 2000A (milliarcseconds)."
    )?;
    writeln!(out, "    pub dy: f64,")?;
    writeln!(out, "}}")?;
    writeln!(out)?;
    writeln!(
        out,
        "/// Embedded EOP table from IERS finals2000A.all ({} entries).",
        rows.len()
    )?;
    writeln!(
        out,
        "/// MJD range: {:.1} – {:.1}.",
        rows.first().unwrap().mjd,
        rows.last().unwrap().mjd,
    )?;
    writeln!(out, "pub static EOP_TABLE: &[EopEntry] = &[")?;

    for r in rows {
        writeln!(
            out,
            "    EopEntry {{ mjd: {:.2}, xp: {:.6}, yp: {:.6}, dut1: {:.7}, dx: {:.3}, dy: {:.3} }},",
            r.mjd, r.xp, r.yp, r.dut1, r.dx, r.dy,
        )?;
    }

    writeln!(out, "];")?;
    Ok(out)
}

// ═══════════════════════════════════════════════════════════════════════════
// Dataset acquisition
// ═══════════════════════════════════════════════════════════════════════════

const FINALS_FILENAME: &str = "finals2000A.all";

/// Primary IERS data source URL.
const IERS_URL: &str = "https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all";

/// Fallback (USNO mirror).
const USNO_URL: &str = "https://maia.usno.navy.mil/ser7/finals2000A.all";

/// Ensure `finals2000A.all` is present in `dir`.
///
/// Tries, in order:
/// 1. Re-use existing file in `dir`
/// 2. Download from IERS (primary) then USNO (fallback)
fn ensure_dataset(dir: &Path) -> Result<PathBuf> {
    let target = dir.join(FINALS_FILENAME);

    // Already present from a previous build?
    if target.exists() {
        eprintln!("IERS EOP: reusing cached {}", target.display());
        return Ok(target);
    }

    // Download
    fs::create_dir_all(dir)?;
    download_finals(&target)?;
    Ok(target)
}

fn download_finals(dst: &Path) -> Result<()> {
    use reqwest::blocking::Client;

    let client = Client::builder()
        .user_agent("siderust-build (rust)")
        .build()?;

    // Try primary, then fallback
    for url in &[IERS_URL, USNO_URL] {
        eprintln!("IERS EOP: downloading from {url}");
        match client.get(*url).send() {
            Ok(resp) if resp.status().is_success() => {
                let bytes = resp.bytes()?;
                if bytes.len() < 1000 {
                    eprintln!(
                        "IERS EOP: response too small ({} bytes), skipping",
                        bytes.len()
                    );
                    continue;
                }
                fs::write(dst, &bytes).with_context(|| format!("write {dst:?}"))?;
                eprintln!("IERS EOP: saved {} bytes to {:?}", bytes.len(), dst);
                return Ok(());
            }
            Ok(resp) => {
                eprintln!("IERS EOP: HTTP {} from {url}", resp.status());
            }
            Err(e) => {
                eprintln!("IERS EOP: download failed from {url}: {e}");
            }
        }
    }
    anyhow::bail!(
        "Could not download finals2000A.all from IERS or USNO. \
         You can manually download it and place it at: {}",
        dst.display()
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// Entrypoint
// ═══════════════════════════════════════════════════════════════════════════

/// Build-time entrypoint: fetch, parse, and generate EOP Rust source.
pub fn run(data_dir: &Path) -> Result<()> {
    let finals_path = ensure_dataset(data_dir)?;

    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed={}", data_dir.display());

    let rows = parse_finals(&finals_path)?;
    eprintln!(
        "IERS EOP: parsed {} rows (MJD {:.1}–{:.1})",
        rows.len(),
        rows.first().unwrap().mjd,
        rows.last().unwrap().mjd,
    );

    let code = generate_rust(&rows)?;

    let out_dir = PathBuf::from(env::var("OUT_DIR")?);
    fs::write(out_dir.join("iers_eop_data.rs"), code.as_bytes())?;
    eprintln!("IERS EOP: iers_eop_data.rs generated");

    Ok(())
}
