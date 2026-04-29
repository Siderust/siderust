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
//! 1. **Fetch**, try local `scripts/iers/dataset/`, then download from IERS
//! 2. **Parse**, fixed-width columns per USNO readme.finals2000A
//! 3. **Codegen**, emit `iers_eop_data.rs` with `pub static EOP_TABLE`
//!
//! ## References
//!
//! * Format: <https://maia.usno.navy.mil/ser7/readme.finals2000A>
//! * Data:   <https://datacenter.iers.org/data/csv/finals2000A.all>

use std::{
    env,
    fs::{self, File},
    io::{BufRead, BufReader, Read},
    path::{Path, PathBuf},
    time::SystemTime,
};

use anyhow::{Context, Result};
use sha2::{Digest, Sha256};
use std::fmt::Write;

/// Provenance record for a single downloaded source dataset.
///
/// Embedded as a comment header in every generated file so that the
/// committed Rust output is reproducible from the original source.
#[derive(Debug, Clone)]
pub(crate) struct SourceProvenance {
    /// Human-readable name of the dataset (e.g. `finals2000A.all`).
    pub source_name: String,
    /// Canonical fetch URL (or `"local cache"` when reused from disk).
    pub source_url: String,
    /// Number of bytes in the input artifact.
    pub byte_count: u64,
    /// Lowercase hex-encoded SHA-256 of the input bytes.
    pub sha256_hex: String,
    /// ISO-8601 timestamp at which the artifact was retrieved or hashed.
    pub retrieved_at: String,
}

impl SourceProvenance {
    fn from_file(name: &str, url: &str, path: &Path) -> Result<Self> {
        let mut file = File::open(path).with_context(|| format!("open {path:?}"))?;
        let mut hasher = Sha256::new();
        let mut buf = [0u8; 64 * 1024];
        let mut count = 0u64;
        loop {
            let n = file.read(&mut buf).with_context(|| format!("read {path:?}"))?;
            if n == 0 {
                break;
            }
            count += n as u64;
            hasher.update(&buf[..n]);
        }
        let digest = hasher.finalize();
        Ok(Self {
            source_name: name.to_string(),
            source_url: url.to_string(),
            byte_count: count,
            sha256_hex: hex_encode(&digest),
            retrieved_at: iso8601_now(),
        })
    }
}

fn hex_encode(bytes: &[u8]) -> String {
    let mut s = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        let _ = write!(s, "{:02x}", b);
    }
    s
}

fn iso8601_now() -> String {
    // Minimal ISO-8601 / RFC-3339 formatter that avoids pulling chrono into
    // the build script. Falls back to `1970-01-01T00:00:00Z` if the system
    // clock predates the Unix epoch (shouldn't happen, but cheap to handle).
    let secs = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let (year, month, day, hour, minute, second) = unix_to_civil(secs);
    format!(
        "{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}Z"
    )
}

/// Convert a Unix timestamp (seconds since 1970-01-01 UTC) to civil
/// `(year, month, day, hour, minute, second)`. Algorithm from Howard
/// Hinnant's date library; valid for `secs >= 0`.
fn unix_to_civil(secs: u64) -> (i32, u8, u8, u8, u8, u8) {
    let days = (secs / 86_400) as i64;
    let rem = (secs % 86_400) as u32;
    let hour = (rem / 3_600) as u8;
    let minute = ((rem % 3_600) / 60) as u8;
    let second = (rem % 60) as u8;

    // Hinnant's days_from_civil inverse.
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
    (year, m, d, hour, minute, second)
}

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

fn generate_rust(rows: &[EopRow], provenance: &SourceProvenance) -> Result<String> {
    let mut out = String::with_capacity(rows.len() * 120);

    writeln!(
        out,
        "// ───────────────────────────────────────────────────────────────────"
    )?;
    writeln!(
        out,
        "// **AUTOGENERATED** by build.rs – DO NOT EDIT BY HAND"
    )?;
    writeln!(out, "//")?;
    writeln!(out, "// Source name : {}", provenance.source_name)?;
    writeln!(out, "// Source URL  : {}", provenance.source_url)?;
    writeln!(out, "// Source bytes: {}", provenance.byte_count)?;
    writeln!(out, "// Source SHA256: {}", provenance.sha256_hex)?;
    writeln!(out, "// Retrieved at: {}", provenance.retrieved_at)?;
    writeln!(out, "// Generator   : siderust scripts/iers/mod.rs")?;
    writeln!(
        out,
        "// ───────────────────────────────────────────────────────────────────"
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
#[allow(dead_code)]
pub fn run(data_dir: &Path) -> Result<()> {
    let finals_path = ensure_dataset(data_dir)?;

    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed={}", data_dir.display());

    let provenance = SourceProvenance::from_file(FINALS_FILENAME, IERS_URL, &finals_path)?;
    eprintln!(
        "IERS EOP: input SHA-256 = {} ({} bytes)",
        provenance.sha256_hex, provenance.byte_count
    );

    let rows = parse_finals(&finals_path)?;
    eprintln!(
        "IERS EOP: parsed {} rows (MJD {:.1}–{:.1})",
        rows.len(),
        rows.first().unwrap().mjd,
        rows.last().unwrap().mjd,
    );

    let code = generate_rust(&rows, &provenance)?;

    let out_dir = PathBuf::from(env::var("OUT_DIR")?);
    fs::write(out_dir.join("iers_eop_data.rs"), code.as_bytes())?;
    eprintln!("IERS EOP: iers_eop_data.rs generated");

    Ok(())
}

/// Like [`run`] but writes `iers_eop_data.rs` to `gen_dir` instead of `OUT_DIR`.
///
/// Used by `build.rs` when `SIDERUST_REGEN=1` to overwrite the committed
/// table in `src/generated/`. Verifies the input SHA-256 against the
/// `SIDERUST_IERS_SHA256` environment variable when set, so CI/release
/// builds can pin a specific Bulletin A/B vintage.
pub fn run_regen(data_dir: &Path, gen_dir: &Path) -> Result<()> {
    let finals_path = ensure_dataset(data_dir)?;

    let provenance = SourceProvenance::from_file(FINALS_FILENAME, IERS_URL, &finals_path)?;
    eprintln!(
        "IERS EOP: input SHA-256 = {} ({} bytes)",
        provenance.sha256_hex, provenance.byte_count
    );

    if let Ok(expected) = env::var("SIDERUST_IERS_SHA256") {
        let expected = expected.trim().to_ascii_lowercase();
        if !expected.is_empty() && expected != provenance.sha256_hex {
            anyhow::bail!(
                "IERS EOP: SHA-256 mismatch.\n  expected: {expected}\n  actual:   {}\n\
                 The downloaded finals2000A.all does not match the pinned hash. \
                 Refusing to overwrite committed tables with unverified data.",
                provenance.sha256_hex
            );
        }
    }

    let rows = parse_finals(&finals_path)?;
    eprintln!(
        "IERS EOP: parsed {} rows (MJD {:.1}–{:.1})",
        rows.len(),
        rows.first().unwrap().mjd,
        rows.last().unwrap().mjd,
    );

    let code = generate_rust(&rows, &provenance)?;

    fs::create_dir_all(gen_dir)?;
    fs::write(gen_dir.join("iers_eop_data.rs"), code.as_bytes())?;
    eprintln!("IERS EOP: iers_eop_data.rs regenerated in {:?}", gen_dir);

    Ok(())
}
