// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Raw Gaia DR3 ingestion structs and lightweight CSV parsing.

use super::error::{GaiaDr3Error, Result};
use super::quality::GaiaDr3QualityFlags;

/// Primitive Gaia DR3 row at the catalogue, CSV, or ADQL boundary.
#[derive(Debug, Clone, PartialEq)]
pub struct GaiaDr3RawSourceRow {
    /// Gaia DR3 `source_id`.
    pub source_id: u64,
    /// Right ascension in degrees.
    pub ra_deg: f64,
    /// Declination in degrees.
    pub dec_deg: f64,
    /// Reference epoch in Julian years.
    pub ref_epoch_jyr: f64,
    /// Proper motion in right ascension times cos(dec), mas/yr.
    pub pmra_mas_per_yr: Option<f64>,
    /// Proper motion in declination, mas/yr.
    pub pmdec_mas_per_yr: Option<f64>,
    /// Annual parallax, mas.
    pub parallax_mas: Option<f64>,
    /// Radial velocity, km/s.
    pub radial_velocity_km_s: Option<f64>,
    /// Gaia G mean magnitude.
    pub phot_g_mean_mag: Option<f64>,
    /// Gaia BP mean magnitude.
    pub phot_bp_mean_mag: Option<f64>,
    /// Gaia RP mean magnitude.
    pub phot_rp_mean_mag: Option<f64>,
    /// Ingestion quality metadata.
    pub quality: GaiaDr3QualityFlags,
}

/// Parse Gaia DR3 CSV text into raw source rows.
///
/// Recognized required column names are `source_id`, `ra`, `dec`, and
/// `ref_epoch`. Optional names are `pmra`, `pmdec`, `parallax`,
/// `radial_velocity`, `phot_g_mean_mag`, `phot_bp_mean_mag`,
/// `phot_rp_mean_mag`, `quality_ok`, and `duplicated_source`.
pub fn parse_gaia_dr3_csv_chunk(input: &str) -> Result<Vec<GaiaDr3RawSourceRow>> {
    let mut header: Option<Vec<String>> = None;
    let mut out = Vec::new();

    for (idx, raw_line) in input.lines().enumerate() {
        let line_no = idx + 1;
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split(',').map(str::trim).collect();
        if header.is_none() {
            header = Some(fields.iter().map(|s| s.to_ascii_lowercase()).collect());
            continue;
        }
        let hdr = header.as_ref().expect("header is initialized");
        let col = |name: &'static str| -> Option<&str> {
            hdr.iter()
                .position(|h| h == name)
                .and_then(|i| fields.get(i))
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
        };
        let parse_opt_f64 = |name: &'static str| -> Result<Option<f64>> {
            if let Some(raw) = col(name) {
                raw.parse::<f64>()
                    .map(Some)
                    .map_err(|_| GaiaDr3Error::InvalidNumber {
                        column: name,
                        line: line_no,
                        raw: raw.to_owned(),
                    })
            } else {
                Ok(None)
            }
        };
        let required_f64 = |name: &'static str| -> Result<f64> {
            parse_opt_f64(name)?.ok_or(GaiaDr3Error::MissingColumn {
                column: name,
                line: line_no,
            })
        };
        let source_id_raw = col("source_id").ok_or(GaiaDr3Error::MissingColumn {
            column: "source_id",
            line: line_no,
        })?;
        let source_id = source_id_raw
            .parse::<u64>()
            .map_err(|_| GaiaDr3Error::InvalidNumber {
                column: "source_id",
                line: line_no,
                raw: source_id_raw.to_owned(),
            })?;

        out.push(GaiaDr3RawSourceRow {
            source_id,
            ra_deg: required_f64("ra")?,
            dec_deg: required_f64("dec")?,
            ref_epoch_jyr: required_f64("ref_epoch")?,
            pmra_mas_per_yr: parse_opt_f64("pmra")?,
            pmdec_mas_per_yr: parse_opt_f64("pmdec")?,
            parallax_mas: parse_opt_f64("parallax")?,
            radial_velocity_km_s: parse_opt_f64("radial_velocity")?,
            phot_g_mean_mag: parse_opt_f64("phot_g_mean_mag")?,
            phot_bp_mean_mag: parse_opt_f64("phot_bp_mean_mag")?,
            phot_rp_mean_mag: parse_opt_f64("phot_rp_mean_mag")?,
            quality: GaiaDr3QualityFlags {
                quality_ok: parse_bool(col("quality_ok")).unwrap_or(true),
                duplicated_source: parse_bool(col("duplicated_source")),
            },
        });
    }

    Ok(out)
}

fn parse_bool(raw: Option<&str>) -> Option<bool> {
    raw.map(|value| {
        matches!(
            value.to_ascii_lowercase().as_str(),
            "1" | "true" | "t" | "yes" | "y" | "ok"
        )
    })
}
