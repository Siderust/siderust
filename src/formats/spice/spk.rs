// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! SPK Type 2 segment reader.
//!
//! Reads Chebyshev polynomial coefficient segments from a parsed DAF container.
//! This is a runtime-capable version of the parser; build-time JPL extraction
//! lives in `siderust-archive/generators/jpl/`.

use super::daf::{Daf, Summary};
use super::SpiceError;

/// Metadata and coefficient data for one SPK Type 2 segment.
#[derive(Debug)]
pub struct SegmentData {
    /// SPK segment type (2 = position Chebyshev, 3 = state Chebyshev).
    pub data_type: i32,
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

/// A parsed SPK segment with the target-center metadata needed for chaining.
#[derive(Debug)]
pub struct IndexedSegmentData {
    /// NAIF body ID of the segment target.
    pub target_id: i32,
    /// NAIF body ID of the segment center.
    pub center_id: i32,
    /// NAIF frame ID.
    pub frame_id: i32,
    /// Segment start epoch in TDB seconds past J2000.
    pub start_et: f64,
    /// Segment end epoch in TDB seconds past J2000.
    pub end_et: f64,
    /// Chebyshev coefficients and record metadata.
    pub data: SegmentData,
}

/// NAIF body IDs for the three segments we need.
pub const SUN_TARGET: i32 = 10;
/// NAIF center ID for the Sun segment (Solar System Barycenter = 0).
pub const SUN_CENTER: i32 = 0;
/// NAIF target ID for the Earth-Moon Barycenter segment.
pub const EMB_TARGET: i32 = 3;
/// NAIF center ID for the EMB segment (Solar System Barycenter = 0).
pub const EMB_CENTER: i32 = 0;
/// NAIF target ID for the Moon segment (301 = Moon).
pub const MOON_TARGET: i32 = 301;
/// NAIF center ID for the Moon segment (Earth-Moon Barycenter = 3).
pub const MOON_CENTER: i32 = 3; // EMB
/// NAIF target ID for Earth center.
pub const EARTH_TARGET: i32 = 399;
/// NAIF center ID for the Earth segment in DE440-class kernels.
pub const EARTH_CENTER: i32 = 3; // EMB

/// Core SPK segments extracted from a planetary BSP file.
#[derive(Debug)]
pub struct BspSegments {
    /// Sun (NAIF 10 → SSB)
    pub sun: SegmentData,
    /// Earth-Moon Barycenter (NAIF 3 → SSB)
    pub emb: SegmentData,
    /// Moon (NAIF 301 → EMB)
    pub moon: SegmentData,
    /// Earth (NAIF 399 → EMB), present in DE440-class kernels.
    pub earth: Option<SegmentData>,
}

fn file_word_count(file_data: &[u8]) -> usize {
    file_data.len() / 8
}

fn ensure_word_in_file(file_data: &[u8], word: usize) -> Result<(), SpiceError> {
    if word == 0 {
        return Err(SpiceError::FormatParse(
            "SPK segment word index must be ≥ 1".to_string(),
        ));
    }
    let max = file_word_count(file_data);
    if word > max {
        return Err(SpiceError::FormatParse(format!(
            "SPK segment word {word} exceeds file data ({max} words)"
        )));
    }
    Ok(())
}

fn read_trailer_f64(
    file_data: &[u8],
    daf: &Daf,
    word: usize,
    name: &str,
) -> Result<f64, SpiceError> {
    ensure_word_in_file(file_data, word)?;
    let value = daf.read_f64_at_word(file_data, word);
    if !value.is_finite() {
        return Err(SpiceError::FormatParse(format!(
            "SPK segment trailer {name} at word {word} is not finite ({value})"
        )));
    }
    Ok(value)
}

/// Read an SPK Type 2 or Type 3 segment from raw file data.
pub fn read_segment(
    file_data: &[u8],
    daf: &Daf,
    summary: &Summary,
) -> Result<SegmentData, SpiceError> {
    if summary.start_word == 0 || summary.end_word == 0 {
        return Err(SpiceError::FormatParse(
            "SPK summary start_word/end_word must be positive".to_string(),
        ));
    }
    if summary.end_word < summary.start_word {
        return Err(SpiceError::FormatParse(format!(
            "SPK summary end_word={} precedes start_word={}",
            summary.end_word, summary.start_word
        )));
    }
    if summary.end_word < summary.start_word + 3 {
        return Err(SpiceError::FormatParse(format!(
            "SPK summary span too small for trailer (start={}, end={})",
            summary.start_word, summary.end_word
        )));
    }

    let end = summary.end_word;
    for word in [end, end - 1, end - 2, end - 3] {
        ensure_word_in_file(file_data, word)?;
    }

    let n_records_f = read_trailer_f64(file_data, daf, end, "n_records")?;
    let rsize_f = read_trailer_f64(file_data, daf, end - 1, "rsize")?;
    let intlen = read_trailer_f64(file_data, daf, end - 2, "intlen")?;
    let init = read_trailer_f64(file_data, daf, end - 3, "init")?;

    if intlen <= 0.0 {
        return Err(SpiceError::FormatParse(format!(
            "SPK segment intlen must be positive (got {intlen})"
        )));
    }

    if n_records_f.fract() != 0.0 || n_records_f < 0.0 {
        return Err(SpiceError::FormatParse(format!(
            "SPK segment n_records must be a non-negative integer (got {n_records_f})"
        )));
    }
    let n_records = n_records_f as usize;
    if n_records == 0 || n_records > 10_000_000 {
        return Err(SpiceError::FormatParse(format!(
            "Implausible n_records={n_records}"
        )));
    }

    if rsize_f.fract() != 0.0 {
        return Err(SpiceError::FormatParse(format!(
            "SPK segment rsize must be an integer (got {rsize_f})"
        )));
    }
    let rsize = rsize_f as usize;
    if !(5..=200).contains(&rsize) {
        return Err(SpiceError::FormatParse(format!(
            "Implausible rsize={rsize} for SPK Type 2/3 segment"
        )));
    }

    let coeff_axes = match summary.data_type {
        2 => 3,
        3 => 6,
        other => {
            return Err(SpiceError::FormatParse(format!(
                "SPK segment Type {other} is unsupported (only Type 2/3)"
            )));
        }
    };
    let ncoeff = (rsize - 2) / coeff_axes;
    if 2 + coeff_axes * ncoeff != rsize {
        return Err(SpiceError::FormatParse(format!(
            "rsize={rsize} is not 2 + {coeff_axes}k for SPK Type {} (ncoeff would be {ncoeff})",
            summary.data_type
        )));
    }

    let total_doubles = n_records.checked_mul(rsize).ok_or_else(|| {
        SpiceError::FormatParse(format!(
            "SPK record buffer size overflows (n_records={n_records}, rsize={rsize})"
        ))
    })?;

    let data_end_word = summary
        .start_word
        .checked_add(total_doubles)
        .and_then(|w| w.checked_sub(1))
        .ok_or_else(|| {
            SpiceError::FormatParse("SPK segment data extent overflows word index".to_string())
        })?;
    if data_end_word > summary.end_word.saturating_sub(4) {
        return Err(SpiceError::FormatParse(format!(
            "SPK coefficient data ends at word {data_end_word} but trailer starts at {}",
            summary.end_word.saturating_sub(3)
        )));
    }
    ensure_word_in_file(file_data, data_end_word)?;

    let mut records = Vec::with_capacity(total_doubles);
    let data_start_word = summary.start_word;
    for i in 0..total_doubles {
        let word = data_start_word + i;
        ensure_word_in_file(file_data, word)?;
        records.push(daf.read_f64_at_word(file_data, word));
    }

    Ok(SegmentData {
        data_type: summary.data_type,
        init,
        intlen,
        rsize,
        ncoeff,
        n_records,
        records,
    })
}

/// Read an SPK Type 2 segment from raw file data.
///
/// This backwards-compatible helper also accepts Type 3 summaries so callers
/// that already gate on "Type 2 or Type 3" receive the correct record shape.
pub fn read_type2_segment(
    file_data: &[u8],
    daf: &Daf,
    summary: &Summary,
) -> Result<SegmentData, SpiceError> {
    read_segment(file_data, daf, summary)
}

/// Parse every supported J2000 Type 2/3 segment in a BSP file.
pub fn parse_indexed_segments(file_data: &[u8]) -> Result<Vec<IndexedSegmentData>, SpiceError> {
    let daf = Daf::parse(file_data)?;
    let mut segments = Vec::new();

    for summary in &daf.summaries {
        if summary.frame_id != 1 || (summary.data_type != 2 && summary.data_type != 3) {
            continue;
        }
        segments.push(IndexedSegmentData {
            target_id: summary.target_id,
            center_id: summary.center_id,
            frame_id: summary.frame_id,
            start_et: summary.start_et,
            end_et: summary.end_et,
            data: read_segment(file_data, &daf, summary)?,
        });
    }

    if segments.is_empty() {
        return Err(SpiceError::FormatParse(
            "BSP contains no supported J2000 SPK Type 2/3 segments".to_string(),
        ));
    }
    Ok(segments)
}

/// Parse a BSP file and extract the required Sun, EMB, Moon, and optional Earth segments.
pub fn parse_bsp(file_data: &[u8]) -> Result<BspSegments, SpiceError> {
    let daf = Daf::parse(file_data)?;

    let find_summary = |target: i32, center: i32| {
        daf.summaries
            .iter()
            .find(|s| s.target_id == target && s.center_id == center)
    };

    let find_segment = |target: i32, center: i32, name: &str| -> Result<SegmentData, SpiceError> {
        let summary = find_summary(target, center).ok_or_else(|| {
            SpiceError::FormatParse(format!(
                "BSP: segment target={} center={} ({}) not found",
                target, center, name
            ))
        })?;

        if summary.data_type != 2 && summary.data_type != 3 {
            return Err(SpiceError::FormatParse(format!(
                "BSP: {} segment is Type {} (only Type 2/3 supported)",
                name, summary.data_type
            )));
        }

        read_segment(file_data, &daf, summary)
    };

    let sun = find_segment(SUN_TARGET, SUN_CENTER, "Sun")?;
    let emb = find_segment(EMB_TARGET, EMB_CENTER, "EMB")?;
    let moon = find_segment(MOON_TARGET, MOON_CENTER, "Moon")?;
    let earth = match find_summary(EARTH_TARGET, EARTH_CENTER) {
        Some(summary) if summary.data_type == 2 || summary.data_type == 3 => {
            Some(read_segment(file_data, &daf, summary)?)
        }
        Some(summary) => {
            return Err(SpiceError::FormatParse(format!(
                "BSP: Earth segment is Type {} (only Type 2/3 supported)",
                summary.data_type
            )));
        }
        None => None,
    };

    Ok(BspSegments {
        sun,
        emb,
        moon,
        earth,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::formats::spice::daf::{Daf, Summary};

    fn make_daf() -> Daf {
        Daf {
            nd: 2,
            ni: 6,
            summaries: vec![],
        }
    }

    fn make_summary(start_word: usize, end_word: usize, data_type: i32) -> Summary {
        Summary {
            start_et: 0.0,
            end_et: 86400.0,
            target_id: 10,
            center_id: 0,
            frame_id: 1,
            data_type,
            start_word,
            end_word,
        }
    }

    /// Write a little-endian f64 at a given byte offset in the buffer.
    fn w_f64(buf: &mut [u8], offset: usize, val: f64) {
        buf[offset..offset + 8].copy_from_slice(&val.to_le_bytes());
    }

    #[test]
    fn read_type2_segment_bad_rsize_too_small() {
        // Place rsize=1 (too small) at end word
        let end_word = 10;
        let mut buf = vec![0u8; (end_word + 1) * 8];
        w_f64(&mut buf, (end_word - 1) * 8, 2.0); // n_records
        w_f64(&mut buf, (end_word - 2) * 8, 1.0); // rsize = 1 (< 5 → invalid)
        w_f64(&mut buf, (end_word - 3) * 8, 1.0); // intlen
        w_f64(&mut buf, (end_word - 4) * 8, 0.0); // init
        let daf = make_daf();
        let summary = make_summary(1, end_word, 2);
        let result = read_type2_segment(&buf, &daf, &summary);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("rsize") || msg.contains("parse error"));
    }

    #[test]
    fn read_type2_segment_bad_rsize_too_large() {
        let end_word = 10;
        let mut buf = vec![0u8; (end_word + 1) * 8];
        w_f64(&mut buf, (end_word - 1) * 8, 1.0); // n_records
        w_f64(&mut buf, (end_word - 2) * 8, 999.0); // rsize = 999 (> 200 → invalid)
        let daf = make_daf();
        let summary = make_summary(1, end_word, 2);
        let result = read_type2_segment(&buf, &daf, &summary);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("rsize") || msg.contains("parse error"));
    }

    #[test]
    fn read_type2_segment_rsize_not_divisible() {
        // rsize=6: ncoeff = (6-2)/3 = 1, 2 + 3*1 = 5 ≠ 6 → invalid
        let end_word = 10;
        let mut buf = vec![0u8; (end_word + 1) * 8];
        w_f64(&mut buf, (end_word - 1) * 8, 1.0); // n_records
        w_f64(&mut buf, (end_word - 2) * 8, 6.0); // rsize = 6 (not 2+3k)
        let daf = make_daf();
        let summary = make_summary(1, end_word, 2);
        let result = read_type2_segment(&buf, &daf, &summary);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("rsize") || msg.contains("parse error"));
    }

    #[test]
    fn read_type2_segment_zero_n_records() {
        let end_word = 10;
        let mut buf = vec![0u8; (end_word + 1) * 8];
        w_f64(&mut buf, (end_word - 1) * 8, 0.0); // n_records = 0 → invalid
        w_f64(&mut buf, (end_word - 2) * 8, 5.0); // rsize = 5
        let daf = make_daf();
        let summary = make_summary(1, end_word, 2);
        let result = read_type2_segment(&buf, &daf, &summary);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("n_records") || msg.contains("parse error"));
    }

    #[test]
    fn read_type2_segment_valid_minimal() {
        // rsize=5 (ncoeff=1), n_records=1
        // data at word 1..=5, metadata at words 6..=9
        // buf needs at least 9*8 = 72 bytes
        let ncoeff = 1usize;
        let rsize = 2 + 3 * ncoeff; // = 5
        let n_records = 1usize;
        let end_word_idx = 9usize; // 1-based
        let mut buf = vec![0u8; (end_word_idx + 1) * 8];
        // Write coefficient record at words 1..=5 (offsets 0..40)
        for i in 0..rsize {
            w_f64(&mut buf, i * 8, (i as f64) * 10.0);
        }
        // Write metadata
        w_f64(&mut buf, (end_word_idx - 4) * 8, 0.0); // init
        w_f64(&mut buf, (end_word_idx - 3) * 8, 86400.0); // intlen
        w_f64(&mut buf, (end_word_idx - 2) * 8, rsize as f64); // rsize
        w_f64(&mut buf, (end_word_idx - 1) * 8, n_records as f64); // n_records
        let daf = make_daf();
        let summary = make_summary(1, end_word_idx, 2);
        let result = read_type2_segment(&buf, &daf, &summary);
        assert!(
            result.is_ok(),
            "Expected Ok, got: {:?}",
            result.err().map(|e| format!("{}", e))
        );
        let seg = result.unwrap();
        assert_eq!(seg.ncoeff, ncoeff);
        assert_eq!(seg.rsize, rsize);
        assert_eq!(seg.n_records, n_records);
        assert_eq!(seg.records.len(), n_records * rsize);
    }

    #[test]
    fn parse_bsp_on_invalid_data_fails() {
        // A valid DAF with no summaries → parse_bsp fails (no Sun segment)
        let mut buf = vec![0u8; 1024];
        buf[0..8].copy_from_slice(b"DAF/SPK ");
        buf[8..12].copy_from_slice(&2i32.to_le_bytes());
        buf[12..16].copy_from_slice(&6i32.to_le_bytes());
        buf[76..80].copy_from_slice(&0i32.to_le_bytes()); // FWARD = 0
        let result = parse_bsp(&buf);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("Sun") || msg.contains("not found") || msg.contains("parse error"));
    }

    #[test]
    fn parse_bsp_with_wrong_data_type_fails() {
        // A DAF with a Sun summary but data_type=5 (unsupported)
        // This requires building a DAF with a summary
        let mut buf = vec![0u8; 2048];
        buf[0..8].copy_from_slice(b"DAF/SPK ");
        buf[8..12].copy_from_slice(&2i32.to_le_bytes());
        buf[12..16].copy_from_slice(&6i32.to_le_bytes());
        buf[76..80].copy_from_slice(&2i32.to_le_bytes()); // FWARD = 2
        let rec = &mut buf[1024..];
        // next=0, prev=0, nsum=1
        rec[0..8].copy_from_slice(&0.0f64.to_le_bytes());
        rec[8..16].copy_from_slice(&0.0f64.to_le_bytes());
        rec[16..24].copy_from_slice(&1.0f64.to_le_bytes());
        // Summary: start_et=0, end_et=1, target=10(Sun), center=0, frame=1, type=5
        rec[24..32].copy_from_slice(&0.0f64.to_le_bytes()); // start_et
        rec[32..40].copy_from_slice(&1.0f64.to_le_bytes()); // end_et
        rec[40..44].copy_from_slice(&10i32.to_le_bytes()); // target_id=10 (Sun)
        rec[44..48].copy_from_slice(&0i32.to_le_bytes()); // center=0
        rec[48..52].copy_from_slice(&1i32.to_le_bytes()); // frame=1
        rec[52..56].copy_from_slice(&5i32.to_le_bytes()); // data_type=5 (unsupported)
        rec[56..60].copy_from_slice(&1i32.to_le_bytes()); // start_word=1
        rec[60..64].copy_from_slice(&100i32.to_le_bytes()); // end_word=100
        let result = parse_bsp(&buf);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("Type") || msg.contains("supported") || msg.contains("parse error"));
    }

    #[test]
    fn segment_data_constants_are_correct() {
        assert_eq!(SUN_TARGET, 10);
        assert_eq!(SUN_CENTER, 0);
        assert_eq!(EMB_TARGET, 3);
        assert_eq!(EMB_CENTER, 0);
        assert_eq!(MOON_TARGET, 301);
        assert_eq!(MOON_CENTER, 3);
    }

    #[test]
    fn read_segment_rejects_overflowing_record_count() {
        let end_word = 10usize;
        let mut buf = vec![0u8; (end_word + 1) * 8];
        w_f64(&mut buf, (end_word - 4) * 8, 0.0);
        w_f64(&mut buf, (end_word - 3) * 8, 86400.0);
        w_f64(&mut buf, (end_word - 2) * 8, 8.0);
        w_f64(&mut buf, (end_word - 1) * 8, 3_000_000.0); // huge n_records
        let daf = make_daf();
        let summary = make_summary(1, end_word, 2);
        let err = read_type2_segment(&buf, &daf, &summary).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("overflow")
                || msg.contains("exceeds")
                || msg.contains("Implausible")
                || msg.contains("coefficient data ends"),
            "{msg}"
        );
    }

    #[test]
    fn read_segment_rejects_truncated_file() {
        let end_word = 10usize;
        let mut buf = vec![0u8; (end_word + 1) * 8];
        w_f64(&mut buf, (end_word - 4) * 8, 0.0);
        w_f64(&mut buf, (end_word - 3) * 8, 86400.0);
        w_f64(&mut buf, (end_word - 2) * 8, 8.0);
        w_f64(&mut buf, (end_word - 1) * 8, 1.0);
        let daf = make_daf();
        let summary = Summary {
            start_word: 1,
            end_word: 5, // trailer words not present
            ..make_summary(1, end_word, 2)
        };
        assert!(read_type2_segment(&buf, &daf, &summary).is_err());
    }

    #[test]
    fn read_segment_rejects_nan_init() {
        let end_word = 10usize;
        let mut buf = vec![0u8; (end_word + 1) * 8];
        w_f64(&mut buf, (end_word - 4) * 8, f64::NAN);
        w_f64(&mut buf, (end_word - 3) * 8, 86400.0);
        w_f64(&mut buf, (end_word - 2) * 8, 8.0);
        w_f64(&mut buf, (end_word - 1) * 8, 1.0);
        let daf = make_daf();
        let summary = make_summary(1, end_word, 2);
        let msg = read_type2_segment(&buf, &daf, &summary)
            .unwrap_err()
            .to_string();
        assert!(msg.contains("init") || msg.contains("finite"), "{msg}");
    }

    #[test]
    fn read_segment_rejects_non_positive_intlen() {
        let end_word = 10usize;
        let mut buf = vec![0u8; (end_word + 1) * 8];
        w_f64(&mut buf, (end_word - 4) * 8, 0.0);
        w_f64(&mut buf, (end_word - 3) * 8, 0.0);
        w_f64(&mut buf, (end_word - 2) * 8, 8.0);
        w_f64(&mut buf, (end_word - 1) * 8, 1.0);
        let daf = make_daf();
        let summary = make_summary(1, end_word, 2);
        let msg = read_type2_segment(&buf, &daf, &summary)
            .unwrap_err()
            .to_string();
        assert!(msg.contains("intlen"), "{msg}");
    }

    #[test]
    fn read_type3_segment_valid_minimal() {
        let ncoeff = 1usize;
        let rsize = 2 + 6 * ncoeff;
        let n_records = 1usize;
        let start_word = 1usize;
        let end_word_idx = start_word + rsize + 3;
        let mut buf = vec![0u8; (end_word_idx + 1) * 8];
        for i in 0..rsize {
            w_f64(&mut buf, (start_word - 1 + i) * 8, i as f64);
        }
        w_f64(&mut buf, (end_word_idx - 4) * 8, 0.0);
        w_f64(&mut buf, (end_word_idx - 3) * 8, 86400.0);
        w_f64(&mut buf, (end_word_idx - 2) * 8, rsize as f64);
        w_f64(&mut buf, (end_word_idx - 1) * 8, n_records as f64);
        let daf = make_daf();
        let summary = make_summary(start_word, end_word_idx, 3);
        let seg = read_segment(&buf, &daf, &summary).expect("type 3 segment");
        assert_eq!(seg.data_type, 3);
        assert_eq!(seg.rsize, rsize);
    }
}
