// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Minimal DAF (Double Precision Array File) parser for SPICE BSP files.
//!
//! Reference: NAIF DAF Required Reading
//! <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html>
//!
//! A DAF file consists of:
//! - **File record** (first 1024 bytes): metadata, ND, NI, forward/backward ptrs
//! - **Summary records**: doubly-linked list of records, each holding segment summaries
//! - **Name records**: human-readable segment names (ignored here)
//! - **Data records**: the actual arrays (Chebyshev coefficients, etc.)
//!
//! For SPK files: ND=2, NI=6. Each summary encodes:
//!   (start_et, end_et, target_id, center_id, frame_id, data_type, start_word, end_word)

/// Parsed DAF container.
pub struct Daf {
    /// Number of double components per summary.
    pub nd: usize,
    /// Number of integer components per summary.
    pub ni: usize,
    /// All segment summaries found in the file.
    pub summaries: Vec<Summary>,
}

/// A single segment summary extracted from the DAF.
pub struct Summary {
    /// Start epoch (TDB seconds past J2000).
    #[allow(dead_code)]
    pub start_et: f64,
    /// End epoch (TDB seconds past J2000).
    #[allow(dead_code)]
    pub end_et: f64,
    /// NAIF body ID of the target.
    pub target_id: i32,
    /// NAIF body ID of the center.
    pub center_id: i32,
    /// Reference frame ID (typically 1 = J2000).
    #[allow(dead_code)]
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
    pub fn parse(data: &[u8]) -> anyhow::Result<Self> {
        if data.len() < 1024 {
            anyhow::bail!("DAF file too small ({} bytes)", data.len());
        }

        // --- File Record (bytes 0..1023) ---
        // Bytes 0..7: file ID locator ("DAF/SPK" padded)
        let locid = std::str::from_utf8(&data[0..8]).unwrap_or("").trim();
        if !locid.starts_with("DAF") {
            anyhow::bail!("Not a DAF file (locator = {:?})", locid);
        }

        // Detect endianness from the locid or from ND/NI values
        // The standard NAIF DAF is big-endian (created on Sun/SGI) but modern
        // versions may be either. We detect by checking if ND/NI make sense.
        let nd_le = read_i32_le(data, 8);
        let ni_le = read_i32_le(data, 12);
        let nd_be = read_i32_be(data, 8);
        let ni_be = read_i32_be(data, 12);

        let (le, nd, ni) = if nd_le > 0 && nd_le <= 100 && ni_le > 0 && ni_le <= 100 {
            (true, nd_le as usize, ni_le as usize)
        } else if nd_be > 0 && nd_be <= 100 && ni_be > 0 && ni_be <= 100 {
            (false, nd_be as usize, ni_be as usize)
        } else {
            anyhow::bail!(
                "Cannot determine DAF endianness: ND/NI LE=({},{}), BE=({},{})",
                nd_le,
                ni_le,
                nd_be,
                ni_be
            );
        };

        let read_i32 = if le { read_i32_le } else { read_i32_be };
        let read_f64 = if le { read_f64_le } else { read_f64_be };

        // DAF File Record layout:
        //   0..7    LOCIDW  (file ID, e.g. "DAF/SPK ")
        //   8..11   ND      (integer)
        //   12..15  NI      (integer)
        //   16..75  LOCIFN  (internal filename, 60 chars)
        //   76..79  FWARD   (forward pointer to first summary record, 1-based)
        //   80..83  BWARD   (backward pointer to last summary record)
        //   84..87  FREE    (first free word address)

        // Forward pointer to first summary record (1-based record number)
        let fward = read_i32(data, 76) as usize;
        // Backward pointer to last summary record
        let _bward = read_i32(data, 80) as usize;
        // Free word address
        let _free = read_i32(data, 84) as usize;

        // SPK files: ND=2, NI=6
        assert!(
            nd == 2 && ni == 6,
            "Expected SPK format (ND=2, NI=6), got ND={}, NI={}",
            nd,
            ni
        );

        // Summary size in doubles (how many doubles per summary)
        let ss = nd + (ni + 1) / 2; // = 2 + 3 = 5 for SPK

        // --- Walk summary records ---
        let mut summaries = Vec::new();
        let mut rec = fward;

        while rec > 0 {
            // Each summary record is 1024 bytes (128 doubles).
            // Record numbers are 1-based; byte offset = (rec - 1) * 1024.
            let rec_offset = (rec - 1) * 1024;
            if rec_offset + 1024 > data.len() {
                anyhow::bail!("Summary record {} extends past EOF", rec);
            }

            // First 3 doubles of a summary record: next, prev, nsum
            let next = read_f64(data, rec_offset) as i64 as usize;
            let _prev = read_f64(data, rec_offset + 8);
            let nsum = read_f64(data, rec_offset + 16) as usize;

            // Summaries start at byte 24 of the record
            for i in 0..nsum {
                let off = rec_offset + 24 + i * ss * 8;
                if off + ss * 8 > data.len() {
                    anyhow::bail!("Summary {} in record {} extends past EOF", i, rec);
                }

                // ND=2 doubles: start_et, end_et
                let start_et = read_f64(data, off);
                let end_et = read_f64(data, off + 8);

                // NI=6 integers packed into ceil(6/2)=3 doubles
                // Integers are stored in the same byte order
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
    ///
    /// DAF word indices are 1-based and refer to 8-byte (f64) words.
    #[inline]
    pub fn read_f64_at_word(&self, data: &[u8], word: usize) -> f64 {
        let offset = (word - 1) * 8;
        // We detect endianness at parse time; for DE440 from NAIF,
        // it's always little-endian on modern systems.
        read_f64_le(data, offset)
    }
}

// --- Byte-order helpers ---

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
