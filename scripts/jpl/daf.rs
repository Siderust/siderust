// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Minimal DAF (Double Precision Array File) parser for SPICE BSP files.
//!
//! Reference: NAIF DAF Required Reading
//! <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html>
//!
//! This parser is shared by DE440/DE441 build pipelines.

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

        let locid = std::str::from_utf8(&data[0..8]).unwrap_or("").trim();
        if !locid.starts_with("DAF") {
            anyhow::bail!("Not a DAF file (locator = {:?})", locid);
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

        let fward = read_i32(data, 76) as usize;
        let _bward = read_i32(data, 80) as usize;
        let _free = read_i32(data, 84) as usize;

        assert!(
            nd == 2 && ni == 6,
            "Expected SPK format (ND=2, NI=6), got ND={}, NI={}",
            nd,
            ni
        );

        let ss = nd + (ni + 1) / 2;

        let mut summaries = Vec::new();
        let mut rec = fward;

        while rec > 0 {
            let rec_offset = (rec - 1) * 1024;
            if rec_offset + 1024 > data.len() {
                anyhow::bail!("Summary record {} extends past EOF", rec);
            }

            let next = read_f64(data, rec_offset) as i64 as usize;
            let _prev = read_f64(data, rec_offset + 8);
            let nsum = read_f64(data, rec_offset + 16) as usize;

            for i in 0..nsum {
                let off = rec_offset + 24 + i * ss * 8;
                if off + ss * 8 > data.len() {
                    anyhow::bail!("Summary {} in record {} extends past EOF", i, rec);
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
