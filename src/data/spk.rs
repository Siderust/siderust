// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! SPK Type 2 segment reader.
//!
//! Reads Chebyshev polynomial coefficient segments from a parsed DAF container.
//! This is a runtime-capable version of the parser used by the build-time
//! pipeline in `scripts/jpl/`.

use super::daf::{Daf, Summary};
use super::DataError;

/// Metadata and coefficient data for one SPK Type 2 segment.
pub struct SegmentData {
    /// Initial epoch (TDB seconds past J2000).
    pub init: f64,
    /// Interval length (seconds).
    pub intlen: f64,
    /// Doubles per record (2 + 3 × ncoeff).
    pub rsize: usize,
    /// Chebyshev polynomial degree + 1 (coefficients per coordinate).
    pub ncoeff: usize,
    /// Number of records.
    pub n_records: usize,
    /// Flattened record data: `n_records × rsize` f64 values.
    pub records: Vec<f64>,
}

/// NAIF body IDs for the three segments we need.
pub const SUN_TARGET: i32 = 10;
pub const SUN_CENTER: i32 = 0;
pub const EMB_TARGET: i32 = 3;
pub const EMB_CENTER: i32 = 0;
pub const MOON_TARGET: i32 = 301;
pub const MOON_CENTER: i32 = 3; // EMB

/// A set of three SPK segments (Sun, EMB, Moon) extracted from a BSP file.
pub struct BspSegments {
    /// Sun (NAIF 10 → SSB)
    pub sun: SegmentData,
    /// Earth-Moon Barycenter (NAIF 3 → SSB)
    pub emb: SegmentData,
    /// Moon (NAIF 301 → EMB)
    pub moon: SegmentData,
}

/// Read an SPK Type 2 segment from raw file data.
pub fn read_type2_segment(
    file_data: &[u8],
    daf: &Daf,
    summary: &Summary,
) -> Result<SegmentData, DataError> {
    let end = summary.end_word;

    let n_records = daf.read_f64_at_word(file_data, end) as usize;
    let rsize = daf.read_f64_at_word(file_data, end - 1) as usize;
    let intlen = daf.read_f64_at_word(file_data, end - 2);
    let init = daf.read_f64_at_word(file_data, end - 3);

    if !(5..=200).contains(&rsize) {
        return Err(DataError::Parse(format!(
            "Implausible rsize={} for SPK Type 2 segment",
            rsize
        )));
    }

    let ncoeff = (rsize - 2) / 3;
    if 2 + 3 * ncoeff != rsize {
        return Err(DataError::Parse(format!(
            "rsize={} is not 2 + 3k for any k (ncoeff would be {})",
            rsize, ncoeff
        )));
    }

    if n_records == 0 || n_records > 10_000_000 {
        return Err(DataError::Parse(format!(
            "Implausible n_records={}",
            n_records
        )));
    }

    let total_doubles = n_records * rsize;
    let mut records = Vec::with_capacity(total_doubles);

    let data_start_word = summary.start_word;
    for i in 0..total_doubles {
        let word = data_start_word + i;
        records.push(daf.read_f64_at_word(file_data, word));
    }

    Ok(SegmentData {
        init,
        intlen,
        rsize,
        ncoeff,
        n_records,
        records,
    })
}

/// Parse a BSP file and extract the three required segments (Sun, EMB, Moon).
pub fn parse_bsp(file_data: &[u8]) -> Result<BspSegments, DataError> {
    let daf = Daf::parse(file_data)?;

    let find_segment = |target: i32, center: i32, name: &str| -> Result<SegmentData, DataError> {
        let summary = daf
            .summaries
            .iter()
            .find(|s| s.target_id == target && s.center_id == center)
            .ok_or_else(|| {
                DataError::Parse(format!(
                    "BSP: segment target={} center={} ({}) not found",
                    target, center, name
                ))
            })?;

        if summary.data_type != 2 && summary.data_type != 3 {
            return Err(DataError::Parse(format!(
                "BSP: {} segment is Type {} (only Type 2/3 supported)",
                name, summary.data_type
            )));
        }

        read_type2_segment(file_data, &daf, summary)
    };

    let sun = find_segment(SUN_TARGET, SUN_CENTER, "Sun")?;
    let emb = find_segment(EMB_TARGET, EMB_CENTER, "EMB")?;
    let moon = find_segment(MOON_TARGET, MOON_CENTER, "Moon")?;

    Ok(BspSegments { sun, emb, moon })
}
