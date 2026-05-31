// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # CSV chunk ingestion
//!
//! [`parse_csv_chunk`] consumes a CSV-style ASCII chunk into a
//! `Vec<CatalogRecord>`. The parser is intentionally Siderust-owned and
//! minimal — adapter crates can wrap heavier formats on top.
//!
//! The first non-empty, non-comment line is treated as the header. Recognized
//! column names: `source_id`, `ra`, `dec`, `epoch`, `pmra`, `pmdec`,
//! `parallax`, `radial_velocity`, `g_mag`, `bp_mag`, `rp_mag`, `quality`.
//! Lines beginning with `#` are skipped. Missing optional columns become
//! `None`. Angles in the input are *degrees* and are converted to radians on
//! the way in.

use crate::qtty::{KmPerSeconds, MilliArcseconds, Radians};

use super::record::CatalogRecord;

/// Errors produced by [`parse_csv_chunk`].
#[derive(Debug, thiserror::Error)]
pub enum CatalogIngestError {
    /// A required column was missing or empty for a given record.
    #[error("required column `{column}` missing on record at line {line}")]
    MissingColumn {
        /// Column name.
        column: &'static str,
        /// 1-based line number in the input chunk.
        line: usize,
    },
    /// A column value could not be parsed as a number.
    #[error("invalid number for column `{column}` on line {line}: {raw}")]
    InvalidNumber {
        /// Column name.
        column: &'static str,
        /// 1-based line number.
        line: usize,
        /// The raw column text that failed to parse.
        raw: String,
    },
}

/// Parse a CSV-style ASCII chunk into [`CatalogRecord`]s.
pub fn parse_csv_chunk(input: &str) -> Result<Vec<CatalogRecord>, CatalogIngestError> {
    let mut header: Option<Vec<String>> = None;
    let mut out = Vec::new();
    for (idx, raw_line) in input.lines().enumerate() {
        let line_no = idx + 1;
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
        if header.is_none() {
            header = Some(fields.iter().map(|s| s.to_ascii_lowercase()).collect());
            continue;
        }
        let hdr = header.as_ref().unwrap();
        let col = |name: &'static str| -> Result<Option<&str>, CatalogIngestError> {
            let pos = hdr.iter().position(|h| h == name);
            Ok(pos
                .and_then(|i| fields.get(i))
                .map(|s| s.trim())
                .filter(|s| !s.is_empty()))
        };
        let parse_opt_f64 = |col_name: &'static str| -> Result<Option<f64>, CatalogIngestError> {
            if let Some(s) = col(col_name)? {
                s.parse::<f64>()
                    .map(Some)
                    .map_err(|_| CatalogIngestError::InvalidNumber {
                        column: col_name,
                        line: line_no,
                        raw: s.to_owned(),
                    })
            } else {
                Ok(None)
            }
        };
        let require_f64 = |col_name: &'static str| -> Result<f64, CatalogIngestError> {
            parse_opt_f64(col_name)?.ok_or(CatalogIngestError::MissingColumn {
                column: col_name,
                line: line_no,
            })
        };
        let source_id_s = col("source_id")?.ok_or(CatalogIngestError::MissingColumn {
            column: "source_id",
            line: line_no,
        })?;
        let source_id =
            source_id_s
                .parse::<u64>()
                .map_err(|_| CatalogIngestError::InvalidNumber {
                    column: "source_id",
                    line: line_no,
                    raw: source_id_s.to_owned(),
                })?;
        let ra_deg = require_f64("ra")?;
        let dec_deg = require_f64("dec")?;
        let epoch_jyr = require_f64("epoch")?;
        let quality_ok = match col("quality")? {
            None => true,
            Some(s) => matches!(s.to_ascii_lowercase().as_str(), "1" | "true" | "ok"),
        };
        out.push(CatalogRecord {
            source_id,
            ra: Radians::new(ra_deg.to_radians()),
            dec: Radians::new(dec_deg.to_radians()),
            epoch_jyr,
            pm_ra_cosdec: parse_opt_f64("pmra")?,
            pm_dec: parse_opt_f64("pmdec")?,
            parallax: parse_opt_f64("parallax")?.map(MilliArcseconds::new),
            radial_velocity: parse_opt_f64("radial_velocity")?.map(KmPerSeconds::new),
            g_mag: parse_opt_f64("g_mag")?,
            bp_mag: parse_opt_f64("bp_mag")?,
            rp_mag: parse_opt_f64("rp_mag")?,
            quality_ok,
        });
    }
    Ok(out)
}
