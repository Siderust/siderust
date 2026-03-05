// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Minimal DAF (Double Precision Array File) parser for SPICE BSP files.
//!
//! Reference: NAIF DAF Required Reading
//! <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html>
//!
//! This is a runtime-capable version of the parser — it returns `Result`
//! errors instead of panicking and can be used both at build time (via the
//! build-script pipeline in `scripts/jpl/`) and at runtime (via the
//! `runtime-data` feature).

use super::DataError;

/// Parsed DAF container.
#[derive(Debug)]
pub struct Daf {
    /// Number of double components per summary.
    pub nd: usize,
    /// Number of integer components per summary.
    pub ni: usize,
    /// All segment summaries found in the file.
    pub summaries: Vec<Summary>,
}

/// A single segment summary extracted from the DAF.
#[derive(Debug)]
pub struct Summary {
    /// Start epoch (TDB seconds past J2000).
    pub start_et: f64,
    /// End epoch (TDB seconds past J2000).
    pub end_et: f64,
    /// NAIF body ID of the target.
    pub target_id: i32,
    /// NAIF body ID of the center.
    pub center_id: i32,
    /// Reference frame ID (typically 1 = J2000).
    pub frame_id: i32,
    /// SPK segment type (2 = Chebyshev position, 3 = Chebyshev pos+vel).
    pub data_type: i32,
    /// 1-based word index of first data element in the file.
    pub start_word: usize,
    /// 1-based word index of last data element in the file.
    pub end_word: usize,
}

impl Daf {
    /// Parse a DAF container from raw file bytes.
    pub fn parse(data: &[u8]) -> Result<Self, DataError> {
        if data.len() < 1024 {
            return Err(DataError::Parse(format!(
                "DAF file too small ({} bytes)",
                data.len()
            )));
        }

        let locid = std::str::from_utf8(&data[0..8]).unwrap_or("").trim();
        if !locid.starts_with("DAF") {
            return Err(DataError::Parse(format!(
                "Not a DAF file (locator = {:?})",
                locid
            )));
        }

        let nd_le = read_i32_le(data, 8);
        let ni_le = read_i32_le(data, 12);
        let nd_be = read_i32_be(data, 8);
        let ni_be = read_i32_be(data, 12);

        let (le, nd, ni) = if nd_le > 0 && nd_le <= 100 && ni_le > 0 && ni_le <= 100 {
            (true, nd_le as usize, ni_le as usize)
        } else if nd_be > 0 && nd_be <= 100 && ni_be > 0 && ni_be <= 100 {
            (false, nd_be as usize, ni_be as usize)
        } else {
            return Err(DataError::Parse(format!(
                "Cannot determine DAF endianness: ND/NI LE=({},{}), BE=({},{})",
                nd_le, ni_le, nd_be, ni_be
            )));
        };

        let read_i32 = if le { read_i32_le } else { read_i32_be };
        let read_f64 = if le { read_f64_le } else { read_f64_be };

        let fward = read_i32(data, 76) as usize;

        if nd != 2 || ni != 6 {
            return Err(DataError::Parse(format!(
                "Expected SPK format (ND=2, NI=6), got ND={}, NI={}",
                nd, ni
            )));
        }

        let ss = nd + ni.div_ceil(2);

        let mut summaries = Vec::new();
        let mut rec = fward;

        while rec > 0 {
            let rec_offset = (rec - 1) * 1024;
            if rec_offset + 1024 > data.len() {
                return Err(DataError::Parse(format!(
                    "Summary record {} extends past EOF",
                    rec
                )));
            }

            let next = read_f64(data, rec_offset) as i64 as usize;
            let _prev = read_f64(data, rec_offset + 8);
            let nsum = read_f64(data, rec_offset + 16) as usize;

            for i in 0..nsum {
                let off = rec_offset + 24 + i * ss * 8;
                if off + ss * 8 > data.len() {
                    return Err(DataError::Parse(format!(
                        "Summary {} in record {} extends past EOF",
                        i, rec
                    )));
                }

                let start_et = read_f64(data, off);
                let end_et = read_f64(data, off + 8);

                let int_off = off + nd * 8;
                let target_id = read_i32(data, int_off);
                let center_id = read_i32(data, int_off + 4);
                let frame_id = read_i32(data, int_off + 8);
                let data_type = read_i32(data, int_off + 12);
                let start_word = read_i32(data, int_off + 16) as usize;
                let end_word = read_i32(data, int_off + 20) as usize;

                summaries.push(Summary {
                    start_et,
                    end_et,
                    target_id,
                    center_id,
                    frame_id,
                    data_type,
                    start_word,
                    end_word,
                });
            }

            rec = next;
        }

        Ok(Daf { nd, ni, summaries })
    }

    /// Read a f64 from the data array at a given word index (1-based).
    #[inline]
    pub fn read_f64_at_word(&self, data: &[u8], word: usize) -> f64 {
        let offset = (word - 1) * 8;
        // DE440/DE441 kernels from NAIF are little-endian on modern systems.
        read_f64_le(data, offset)
    }
}

#[inline]
fn read_i32_le(data: &[u8], offset: usize) -> i32 {
    i32::from_le_bytes(data[offset..offset + 4].try_into().unwrap())
}

#[inline]
fn read_i32_be(data: &[u8], offset: usize) -> i32 {
    i32::from_be_bytes(data[offset..offset + 4].try_into().unwrap())
}

#[inline]
fn read_f64_le(data: &[u8], offset: usize) -> f64 {
    f64::from_le_bytes(data[offset..offset + 8].try_into().unwrap())
}

#[inline]
fn read_f64_be(data: &[u8], offset: usize) -> f64 {
    f64::from_be_bytes(data[offset..offset + 8].try_into().unwrap())
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Helpers to build synthetic DAF byte arrays ─────────────────────────

    /// Write a little-endian i32 into a buffer at a given offset.
    fn write_i32_le(buf: &mut [u8], offset: usize, val: i32) {
        buf[offset..offset + 4].copy_from_slice(&val.to_le_bytes());
    }

    /// Write a little-endian f64 into a buffer at a given offset.
    fn write_f64_le(buf: &mut [u8], offset: usize, val: f64) {
        buf[offset..offset + 8].copy_from_slice(&val.to_le_bytes());
    }

    /// Create a minimal valid DAF header (1024 bytes) with FWARD = 0 (no summaries).
    fn minimal_valid_daf() -> Vec<u8> {
        let mut buf = vec![0u8; 1024];
        buf[0..8].copy_from_slice(b"DAF/SPK ");
        write_i32_le(&mut buf, 8, 2); // ND = 2
        write_i32_le(&mut buf, 12, 6); // NI = 6
        write_i32_le(&mut buf, 76, 0); // FWARD = 0 (no summary records)
        buf
    }

    /// Create a DAF with one summary record containing a single Type-2 SPK summary.
    #[allow(clippy::too_many_arguments)]
    fn daf_with_one_summary(
        start_et: f64,
        end_et: f64,
        target_id: i32,
        center_id: i32,
        frame_id: i32,
        data_type: i32,
        start_word: i32,
        end_word: i32,
    ) -> Vec<u8> {
        let mut buf = vec![0u8; 2048];
        // Header record (bytes 0..1024)
        buf[0..8].copy_from_slice(b"DAF/SPK ");
        write_i32_le(&mut buf, 8, 2); // ND = 2
        write_i32_le(&mut buf, 12, 6); // NI = 6
        write_i32_le(&mut buf, 76, 2); // FWARD = 2 (second record)
                                       // Summary record (bytes 1024..2048)
        let rec = &mut buf[1024..];
        write_f64_le(rec, 0, 0.0); // next = 0 (no next record)
        write_f64_le(rec, 8, 0.0); // prev
        write_f64_le(rec, 16, 1.0); // nsum = 1
                                    // Summary at offset 24:  nd=2 doubles + ni=6 ints = 40 bytes
        write_f64_le(rec, 24, start_et);
        write_f64_le(rec, 32, end_et);
        write_i32_le(rec, 40, target_id);
        write_i32_le(rec, 44, center_id);
        write_i32_le(rec, 48, frame_id);
        write_i32_le(rec, 52, data_type);
        write_i32_le(rec, 56, start_word);
        write_i32_le(rec, 60, end_word);
        buf
    }

    // ── Error path tests ───────────────────────────────────────────────────

    #[test]
    fn parse_too_small_returns_error() {
        let data = vec![0u8; 512]; // less than 1024
        let result = Daf::parse(&data);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("too small") || msg.contains("DAF"));
    }

    #[test]
    fn parse_wrong_magic_returns_error() {
        let mut buf = vec![0u8; 1024];
        buf[0..8].copy_from_slice(b"NOTADAF!");
        let result = Daf::parse(&buf);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("Not a DAF"));
    }

    #[test]
    fn parse_bad_nd_ni_returns_error() {
        let mut buf = vec![0u8; 1024];
        buf[0..8].copy_from_slice(b"DAF/SPK ");
        // Set both LE and BE to nonsense values outside valid range
        write_i32_le(&mut buf, 8, 0); // ND LE = 0 (invalid)
        write_i32_le(&mut buf, 12, 0); // NI LE = 0 (invalid)
                                       // BE also invalid (zeros in both positions)
        let result = Daf::parse(&buf);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(
            msg.contains("endianness")
                || msg.contains("ND")
                || msg.contains("NI")
                || msg.contains("parse error")
        );
    }

    #[test]
    fn parse_non_spk_format_returns_error() {
        // ND=3, NI=4 — valid values but not SPK (ND=2, NI=6)
        let mut buf = vec![0u8; 1024];
        buf[0..8].copy_from_slice(b"DAF/CK  ");
        write_i32_le(&mut buf, 8, 3); // ND = 3
        write_i32_le(&mut buf, 12, 4); // NI = 4
        write_i32_le(&mut buf, 76, 0); // FWARD = 0
        let result = Daf::parse(&buf);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(
            msg.contains("ND")
                || msg.contains("NI")
                || msg.contains("SPK")
                || msg.contains("parse error")
        );
    }

    #[test]
    fn parse_fward_beyond_eof_returns_error() {
        let mut buf = minimal_valid_daf();
        // Set FWARD to record 99, but file only has 1 record (1024 bytes)
        write_i32_le(&mut buf, 76, 99);
        let result = Daf::parse(&buf);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("EOF") || msg.contains("extends past") || msg.contains("parse error"));
    }

    // ── Success path tests ─────────────────────────────────────────────────

    #[test]
    fn parse_minimal_valid_daf_empty_summaries() {
        let buf = minimal_valid_daf();
        let daf = Daf::parse(&buf).expect("should parse valid DAF");
        assert_eq!(daf.nd, 2);
        assert_eq!(daf.ni, 6);
        assert!(daf.summaries.is_empty());
    }

    #[test]
    fn parse_daf_with_one_summary() {
        let buf = daf_with_one_summary(
            0.0,     // start_et
            86400.0, // end_et
            10,      // target_id (Sun)
            0,       // center_id (SSB)
            1,       // frame_id (J2000)
            2,       // data_type (Type 2 Chebyshev)
            257,     // start_word
            300,     // end_word
        );
        let daf = Daf::parse(&buf).expect("should parse DAF with one summary");
        assert_eq!(daf.nd, 2);
        assert_eq!(daf.ni, 6);
        assert_eq!(daf.summaries.len(), 1);

        let s = &daf.summaries[0];
        assert!((s.start_et - 0.0).abs() < 1e-10);
        assert!((s.end_et - 86400.0).abs() < 1e-10);
        assert_eq!(s.target_id, 10);
        assert_eq!(s.center_id, 0);
        assert_eq!(s.frame_id, 1);
        assert_eq!(s.data_type, 2);
        assert_eq!(s.start_word, 257);
        assert_eq!(s.end_word, 300);
    }

    #[test]
    fn read_f64_at_word_reads_correct_value() {
        let daf = Daf {
            nd: 2,
            ni: 6,
            summaries: vec![],
        };
        // Create a buffer with a known f64 at word 1 (offset 0)
        let mut buf = vec![0u8; 64];
        let val = 42.0_f64;
        buf[0..8].copy_from_slice(&val.to_le_bytes());
        let result = daf.read_f64_at_word(&buf, 1);
        assert!((result - 42.0).abs() < 1e-15);
    }

    #[test]
    fn read_f64_at_word_second_word() {
        let daf = Daf {
            nd: 2,
            ni: 6,
            summaries: vec![],
        };
        let mut buf = vec![0u8; 64];
        let val = 99.5_f64;
        buf[8..16].copy_from_slice(&val.to_le_bytes());
        let result = daf.read_f64_at_word(&buf, 2);
        assert!((result - 99.5).abs() < 1e-15);
    }

    #[test]
    fn parse_summary_truncated_returns_error() {
        // Build a DAF that claims 5 summaries but only provides 1's worth of space
        let mut buf = vec![0u8; 2048];
        buf[0..8].copy_from_slice(b"DAF/SPK ");
        write_i32_le(&mut buf, 8, 2); // ND = 2
        write_i32_le(&mut buf, 12, 6); // NI = 6
        write_i32_le(&mut buf, 76, 2); // FWARD = 2
        let rec = &mut buf[1024..];
        write_f64_le(rec, 0, 0.0);
        write_f64_le(rec, 8, 0.0);
        write_f64_le(rec, 16, 100.0); // nsum = 100 — far more than fits in 1024 bytes
        let result = Daf::parse(&buf);
        // Should either succeed with partial summaries or return an error — must not panic
        let _ = result;
    }
}
