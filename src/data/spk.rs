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
#[derive(Debug)]
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
#[derive(Debug)]
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::daf::{Daf, Summary};

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
}
