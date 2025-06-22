//! VSOP87 data collector
//!
//! Walks over every text file under `data_dir`, recognises the VSOP87 data format
//! (grouped by version, planet, Cartesian coordinate and T‑exponent) and packs
//! everything into a nested [`VersionMap`].
//!
//! The file format comes from the original VSOP87 model published by Bretagnon & Simon.
//! Each file name encodes the version (A, B, C, E…) and the extension encodes the
//! planet. Inside every file there are *blocks* introduced by a header line like:
//!
//! ```text
//! VSOP87 VERSION A12   VENUS VARIABLE  2  (XYZ) *T**3     342 TERMS
//! ```
//!
//! The header tells us:
//! 1. **Planet**    (`VENUS`)
//! 2. **Coordinate** (1 = X, 2 = Y, 3 = Z)
//! 3. **T‑power**   (3 → terms multiplied by _T_^3)
//! 4. **Count**     (how many lines / terms follow the header)
//!
//! Each data line that follows is either
//!
//! * `S  K  A  B  C` (five floats) if both sine and cosine components are present, or
//! * `0  0  A  B  C` if the coefficient `A` is already the resultant cosine amplitude.
//!
//! The two first numbers, *S* and *K*, correspond to `S·sin()` and `K·cos()`
//! components respectively.  We fold them into a single cosine term using
//! trigonometric identities (see [`sk_to_term`]).

use super::{Term, VersionMap};

use anyhow::Context;
use regex::Regex;
use walkdir::WalkDir;
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

// ---------------------------------------------------------------------------
// Regexes
// ---------------------------------------------------------------------------

/// Recognises VSOP87 data files; captures the **version letter** (A, E, …)
const REGEX_FILE: &str = r"^VSOP87([A-Z])\.[A-Za-z]{3}$";

/// Recognises block headers inside each file.
///
/// Capture groups:
/// 1. planet  (e.g. "VENUS")
/// 2. coord   (1/2/3 → X/Y/Z)
/// 3. T power (0…5)
/// 4. count   (number of terms that follow)
const HEADER_REGEX: &str =
    r"VSOP87 VERSION [A-Z]?\d+\s+(\S+)\s+VARIABLE\s+(\d+)\s+\(XYZ\)\s+\*T\*\*(\d+)\s+(\d+)\s+TERMS";


/// Fold `S·sin(x) + K·cos(x)` into `R·cos(x − α)` where
/// * `R = √(S² + K²)` and `α = atan2(S, K)`.
///
/// If `A` is already non‑zero, we assume the line already stores the folded
/// amplitude and simply return it.  Lines where all coefficients evaluate to
/// ~0 are ignored (`None`).
fn sk_to_term(s: f64, k: f64, a: f64, b: f64, c: f64) -> Option<Term> {
    if a.abs() > 1e-15 {
        // The amplitude is already given directly.
        Some(Term { a, b, c })
    } else if s.abs() > 1e-15 || k.abs() > 1e-15 {
        // Combine sine+cosine into a single cosine term.
        let r = (s * s + k * k).sqrt();
        let alpha = s.atan2(k); // Angle shift when folding S·sin+K·cos.
        Some(Term {
            a: r,
            b: b - alpha,
            c,
        })
    } else {
        // Zero term, can be discarded.
        None
    }
}

/// Split a raw data line into its five floating‑point components `(S, K, A, B, C)`.
/// Returns `None` if the line has fewer than five columns or if parsing fails.
fn parse_data_line(line: &str) -> Option<(f64, f64, f64, f64, f64)> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 5 {
        return None;
    }
    let len = parts.len();
    let c = parts[len - 1].parse().ok()?;
    let b = parts[len - 2].parse().ok()?;
    let a = parts[len - 3].parse().ok()?;
    let k = parts[len - 4].parse().ok()?;
    let s = parts[len - 5].parse().ok()?;
    Some((s, k, a, b, c))
}


/// Build a [`VersionMap`] from the dataset found under `data_dir`.
///
/// * Walks every file matching [`REGEX_FILE`].
/// * For each file, iterates over every line, detects header blocks via
///   [`HEADER_REGEX`] and collects the following `count` data lines.
/// * Converts each data line into a [`Term`] via [`sk_to_term`].
/// * Inserts the term at `versions[version][planet][coord][T_power]`.
pub fn collect_terms(data_dir: &Path) -> anyhow::Result<VersionMap> {
    // --- Compile regexes ---------------------------------------------------
    let file_re   = Regex::new(REGEX_FILE)?;
    let header_re = Regex::new(HEADER_REGEX)?;

    let mut versions: VersionMap = BTreeMap::new();

    // --- Walk every file under `data_dir` ----------------------------------
    for entry in WalkDir::new(data_dir)
        .min_depth(1)    // skip the root dir itself
        .into_iter()
        .filter_map(Result::ok) // discard IO errors
    {
        let path  = entry.path();
        let fname = path.file_name().unwrap().to_string_lossy();

        // Only keep files whose names match the VSOP87 convention.
        let caps = match file_re.captures(&fname) {
            Some(c) => c,
            None    => continue,            // not a VSOP87 file → skip
        };
        let v_char = caps[1].chars().next().unwrap(); // e.g. 'A', 'E'

        // --- Open current file --------------------------------------------
        let file   = File::open(path)
            .with_context(|| format!("Could not open {path:?}"))?;
        let reader = BufReader::new(file);

        // Per‑block state as we scan through the file.
        let mut cur_planet = String::new(); // planet currently being parsed
        let mut cur_coord  : u8 = 0;        // 1 =s X, 2 = Y, 3 = Z
        let mut cur_t      : u8 = 0;        // exponent of T
        let mut remaining  : u32 = 0;       // terms left in current block

        // --- Process every line -------------------------------------------
        for line in reader.lines() {
            let line = line?;

            // 1) Header → update block context and term counter
            if let Some(cap) = header_re.captures(&line) {
                cur_planet = cap[1].to_uppercase();
                cur_coord  = cap[2].parse::<u8>().unwrap();
                cur_t      = cap[3].parse::<u8>().unwrap();
                remaining  = cap[4].parse::<u32>().unwrap();
                continue; // header processed, go to next line
            }

            // 2) Data line inside current block
            if remaining > 0 {
                if let Some((s, k, a, b, c)) = parse_data_line(&line) {
                    if let Some(term) = sk_to_term(s, k, a, b, c) {
                        // Navigate (or create) the nested maps and push the term.
                        let planet_map = versions.entry(v_char).or_default();
                        let coord_map  = planet_map.entry(cur_planet.clone()).or_default();
                        let t_map      = coord_map.entry(cur_coord).or_default();
                        let vec_terms  = t_map.entry(cur_t).or_default();
                        vec_terms.push(term);
                    }
                    remaining -= 1; // one less term to read in this block
                }
            }
        }
    }

    Ok(versions)
}
