// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Shared SPK Type 2 segment reader for DE440/DE441 build pipelines.

use super::daf::{Daf, Summary};
use std::path::Path;

/// Metadata and coefficient data for one SPK Type 2 segment.
pub struct SegmentMeta {
    /// Initial epoch (TDB seconds past J2000).
    pub init: f64,
    /// Interval length (seconds).
    pub intlen: f64,
    /// Doubles per record.
    pub rsize: usize,
    /// Chebyshev polynomial degree + 1 (coefficients per coordinate).
    pub ncoeff: usize,
    /// Number of records.
    pub n_records: usize,
    /// Flattened record data: n_records × rsize f64 values.
    pub records: Vec<f64>,
}

/// Read an SPK Type 2 segment from raw file data.
pub fn read_type2_segment(
    file_data: &[u8],
    daf: &Daf,
    summary: &Summary,
) -> anyhow::Result<SegmentMeta> {
    let end = summary.end_word;

    let n_records = daf.read_f64_at_word(file_data, end) as usize;
    let rsize = daf.read_f64_at_word(file_data, end - 1) as usize;
    let intlen = daf.read_f64_at_word(file_data, end - 2);
    let init = daf.read_f64_at_word(file_data, end - 3);

    if rsize < 5 || rsize > 200 {
        anyhow::bail!("Implausible rsize={} for SPK Type 2 segment", rsize);
    }

    let ncoeff = (rsize - 2) / 3;
    if 2 + 3 * ncoeff != rsize {
        anyhow::bail!(
            "rsize={} is not 2 + 3k for any k (ncoeff would be {})",
            rsize,
            ncoeff
        );
    }

    if n_records == 0 || n_records > 10_000_000 {
        anyhow::bail!("Implausible n_records={}", n_records);
    }

    let total_doubles = n_records * rsize;
    let mut records = Vec::with_capacity(total_doubles);

    let data_start_word = summary.start_word;
    for i in 0..total_doubles {
        let word = data_start_word + i;
        records.push(daf.read_f64_at_word(file_data, word));
    }

    Ok(SegmentMeta {
        init,
        intlen,
        rsize,
        ncoeff,
        n_records,
        records,
    })
}

/// Write a segment's coefficient data as raw little-endian f64 bytes.
pub fn write_binary(meta: &SegmentMeta, path: &Path) -> anyhow::Result<()> {
    let mut buf = Vec::with_capacity(meta.records.len() * 8);
    for &val in &meta.records {
        buf.extend_from_slice(&val.to_le_bytes());
    }
    std::fs::write(path, &buf)?;
    Ok(())
}
