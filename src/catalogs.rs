// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Large catalog ingestion and cone search
//!
//! ## Scientific scope
//!
//! Gaia-style large-catalog records: stable source identifier, ICRS
//! position + epoch, proper motion, parallax, radial velocity (when
//! available), and optional photometric bands and quality flags. The
//! module exposes a small ingestion + cone-search API targeting
//! locally-stored catalog chunks; remote catalog querying belongs in
//! a separate operational layer.
//!
//! ## Technical scope
//!
//! - [`CatalogRecord`] is the column-aligned star record; missing
//!   optional columns are represented as `None`.
//! - [`LargeStarCatalog`] is the in-memory container with cone search
//!   (`cone_search`), magnitude/quality filters, and an epoch-aware
//!   query that propagates positions via proper motion.
//! - [`parse_csv_chunk`] consumes a CSV-style ASCII chunk into a
//!   `Vec<CatalogRecord>`. The parser is intentionally Siderust-owned
//!   and minimal — adapter crates can wrap heavier formats on top.
//!
//! ## References
//!
//! - Lindegren, L. (2018). *Gaia DR2: The catalogue's astrometric
//!   content*. A&A 616 A2.
//! - Górski, K. M. et al. (2005). "HEALPix: a framework for
//!   high-resolution discretization." *ApJ*, 622:759.

#![forbid(unsafe_code)]

use core::f64::consts::TAU;

/// One row of a large stellar catalog.
#[derive(Debug, Clone, PartialEq)]
pub struct CatalogRecord {
    /// Stable source identifier (e.g. Gaia source_id).
    pub source_id: u64,
    /// ICRS right ascension in radians at [`Self::epoch_jyr`].
    pub ra: f64,
    /// ICRS declination in radians at [`Self::epoch_jyr`].
    pub dec: f64,
    /// Catalog epoch in Julian years (e.g. 2016.0 for Gaia DR3).
    pub epoch_jyr: f64,
    /// Proper motion in right ascension * cos(δ), mas/yr.
    pub pm_ra_cosdec: Option<f64>,
    /// Proper motion in declination, mas/yr.
    pub pm_dec: Option<f64>,
    /// Parallax, mas.
    pub parallax: Option<f64>,
    /// Radial velocity, km/s.
    pub radial_velocity: Option<f64>,
    /// G-band mean magnitude.
    pub g_mag: Option<f64>,
    /// BP-band mean magnitude.
    pub bp_mag: Option<f64>,
    /// RP-band mean magnitude.
    pub rp_mag: Option<f64>,
    /// Catalog-supplied quality flag (`true` = good).
    pub quality_ok: bool,
}

impl CatalogRecord {
    /// Propagate the record's ICRS direction by linear proper motion
    /// to the requested Julian-year epoch. Returns `(ra, dec)` in
    /// radians. If proper motion is unavailable the original
    /// `(ra, dec)` is returned unchanged.
    pub fn propagate_to(&self, target_epoch_jyr: f64) -> (f64, f64) {
        let dt_yr = target_epoch_jyr - self.epoch_jyr;
        let (Some(mu_a), Some(mu_d)) = (self.pm_ra_cosdec, self.pm_dec) else {
            return (self.ra, self.dec);
        };
        let mas_to_rad = (core::f64::consts::PI / 180.0) / 3_600_000.0;
        let cos_d = self.dec.cos();
        let d_ra = if cos_d.abs() > 1e-12 {
            (mu_a * mas_to_rad * dt_yr) / cos_d
        } else {
            0.0
        };
        let d_dec = mu_d * mas_to_rad * dt_yr;
        let mut ra = (self.ra + d_ra).rem_euclid(TAU);
        let mut dec = self.dec + d_dec;
        if dec > core::f64::consts::FRAC_PI_2 {
            dec = core::f64::consts::PI - dec;
            ra = (ra + core::f64::consts::PI).rem_euclid(TAU);
        } else if dec < -core::f64::consts::FRAC_PI_2 {
            dec = -core::f64::consts::PI - dec;
            ra = (ra + core::f64::consts::PI).rem_euclid(TAU);
        }
        (ra, dec)
    }
}

/// In-memory large-star catalog.
#[derive(Debug, Default, Clone)]
pub struct LargeStarCatalog {
    records: Vec<CatalogRecord>,
}

/// Constraints for a cone-search query.
#[derive(Debug, Clone, Copy, Default)]
pub struct CatalogFilter {
    /// If set, exclude records with `g_mag` greater than this.
    pub max_g_mag: Option<f64>,
    /// If set, only include records whose `quality_ok` matches.
    pub require_quality: Option<bool>,
}

impl LargeStarCatalog {
    /// Empty catalog.
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of records held.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// True iff the catalog has no records.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Append a single record.
    pub fn push(&mut self, record: CatalogRecord) {
        self.records.push(record);
    }

    /// Extend with an iterator of records.
    pub fn extend(&mut self, iter: impl IntoIterator<Item = CatalogRecord>) {
        self.records.extend(iter);
    }

    /// Borrow the underlying record slice.
    pub fn records(&self) -> &[CatalogRecord] {
        &self.records
    }

    /// Cone search around the supplied ICRS direction. `center_ra`,
    /// `center_dec`, and `radius` are in radians. Records are
    /// evaluated *at their catalog epoch* (no proper-motion shift);
    /// use [`Self::cone_search_at_epoch`] for epoch-aware queries.
    pub fn cone_search(
        &self,
        center_ra: f64,
        center_dec: f64,
        radius: f64,
        filter: CatalogFilter,
    ) -> Vec<&CatalogRecord> {
        self.records
            .iter()
            .filter(|r| {
                inside_cone(r.ra, r.dec, center_ra, center_dec, radius) && passes_filter(r, &filter)
            })
            .collect()
    }

    /// Cone search at a specified epoch. Each record is propagated to
    /// `target_epoch_jyr` via proper motion before the angular-distance
    /// test runs. Records without proper motion are tested at their
    /// catalog epoch.
    pub fn cone_search_at_epoch(
        &self,
        center_ra: f64,
        center_dec: f64,
        radius: f64,
        target_epoch_jyr: f64,
        filter: CatalogFilter,
    ) -> Vec<&CatalogRecord> {
        self.records
            .iter()
            .filter(|r| {
                let (ra, dec) = r.propagate_to(target_epoch_jyr);
                inside_cone(ra, dec, center_ra, center_dec, radius) && passes_filter(r, &filter)
            })
            .collect()
    }

    /// Chunked cone search: splits the underlying record slice into
    /// `chunk_size`-row windows and yields matches per chunk. The
    /// concatenated result is identical to [`Self::cone_search`].
    pub fn cone_search_chunked(
        &self,
        center_ra: f64,
        center_dec: f64,
        radius: f64,
        filter: CatalogFilter,
        chunk_size: usize,
    ) -> impl Iterator<Item = Vec<&CatalogRecord>> {
        let chunks = self.records.chunks(chunk_size.max(1));
        chunks.map(move |chunk| {
            chunk
                .iter()
                .filter(|r| {
                    inside_cone(r.ra, r.dec, center_ra, center_dec, radius)
                        && passes_filter(r, &filter)
                })
                .collect()
        })
    }
}

fn passes_filter(r: &CatalogRecord, f: &CatalogFilter) -> bool {
    if let Some(max_g) = f.max_g_mag {
        if let Some(g) = r.g_mag {
            if g > max_g {
                return false;
            }
        } else {
            return false;
        }
    }
    if let Some(q) = f.require_quality {
        if r.quality_ok != q {
            return false;
        }
    }
    true
}

fn inside_cone(ra: f64, dec: f64, c_ra: f64, c_dec: f64, radius: f64) -> bool {
    // Haversine formula on the unit sphere.
    let d_ra = ra - c_ra;
    let s = (dec - c_dec).sin();
    let s2 = (d_ra / 2.0).sin();
    let a = (s / 2.0).powi(2) + dec.cos() * c_dec.cos() * s2 * s2;
    let c = 2.0 * a.sqrt().min(1.0).asin();
    c <= radius
}

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
///
/// The first non-empty, non-comment line is the header. Recognised
/// column names: `source_id`, `ra`, `dec`, `epoch`, `pmra`, `pmdec`,
/// `parallax`, `radial_velocity`, `g_mag`, `bp_mag`, `rp_mag`,
/// `quality`. Lines beginning with `#` are skipped. Missing optional
/// columns become `None`. Angles in the input are *degrees* and are
/// converted to radians on the way in.
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
            ra: ra_deg.to_radians(),
            dec: dec_deg.to_radians(),
            epoch_jyr,
            pm_ra_cosdec: parse_opt_f64("pmra")?,
            pm_dec: parse_opt_f64("pmdec")?,
            parallax: parse_opt_f64("parallax")?,
            radial_velocity: parse_opt_f64("radial_velocity")?,
            g_mag: parse_opt_f64("g_mag")?,
            bp_mag: parse_opt_f64("bp_mag")?,
            rp_mag: parse_opt_f64("rp_mag")?,
            quality_ok,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    fn rec(id: u64, ra_deg: f64, dec_deg: f64) -> CatalogRecord {
        CatalogRecord {
            source_id: id,
            ra: ra_deg.to_radians(),
            dec: dec_deg.to_radians(),
            epoch_jyr: 2016.0,
            pm_ra_cosdec: None,
            pm_dec: None,
            parallax: None,
            radial_velocity: None,
            g_mag: None,
            bp_mag: None,
            rp_mag: None,
            quality_ok: true,
        }
    }

    #[test]
    fn cone_includes_center() {
        let mut cat = LargeStarCatalog::new();
        cat.push(rec(1, 10.0, 0.0));
        let hits = cat.cone_search(
            10.0_f64.to_radians(),
            0.0,
            1.0_f64.to_radians(),
            CatalogFilter::default(),
        );
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn cone_excludes_far_records() {
        let mut cat = LargeStarCatalog::new();
        cat.push(rec(1, 10.0, 0.0));
        cat.push(rec(2, 80.0, 0.0));
        let hits = cat.cone_search(
            10.0_f64.to_radians(),
            0.0,
            5.0_f64.to_radians(),
            CatalogFilter::default(),
        );
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].source_id, 1);
    }

    #[test]
    fn cone_handles_ra_wraparound() {
        let mut cat = LargeStarCatalog::new();
        cat.push(rec(1, 359.0, 0.0));
        cat.push(rec(2, 1.0, 0.0));
        let hits = cat.cone_search(0.0, 0.0, 2.0_f64.to_radians(), CatalogFilter::default());
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn cone_near_pole_includes_records() {
        let mut cat = LargeStarCatalog::new();
        // All records near the north pole — large RA spread.
        cat.push(rec(1, 0.0, 89.5));
        cat.push(rec(2, 90.0, 89.5));
        cat.push(rec(3, 180.0, 89.5));
        cat.push(rec(4, 270.0, 89.5));
        let hits = cat.cone_search(
            0.0,
            (PI / 2.0) - 0.001,
            2.0_f64.to_radians(),
            CatalogFilter::default(),
        );
        assert_eq!(hits.len(), 4);
    }

    #[test]
    fn cone_with_filter_drops_faint() {
        let mut cat = LargeStarCatalog::new();
        let mut r = rec(1, 10.0, 0.0);
        r.g_mag = Some(8.0);
        cat.push(r);
        let mut r = rec(2, 10.0, 0.0);
        r.g_mag = Some(15.0);
        cat.push(r);
        let hits = cat.cone_search(
            10.0_f64.to_radians(),
            0.0,
            1.0_f64.to_radians(),
            CatalogFilter {
                max_g_mag: Some(10.0),
                require_quality: None,
            },
        );
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].source_id, 1);
    }

    #[test]
    fn epoch_aware_query_shifts_position() {
        // Star at RA=0, Dec=0 with pm_ra=10 deg/yr (deliberately
        // huge) → after 1 yr it moves 10 deg in RA.
        let mut r = rec(1, 0.0, 0.0);
        // Convert 10 deg/yr to mas/yr * cos(dec)=cos(0)=1.
        r.pm_ra_cosdec = Some(10.0 * 3_600_000.0);
        r.pm_dec = Some(0.0);
        let mut cat = LargeStarCatalog::new();
        cat.push(r);
        // Search at RA=10 deg, epoch 2017.0 should hit; at epoch 2016.0 must miss.
        let hits = cat.cone_search_at_epoch(
            10.0_f64.to_radians(),
            0.0,
            0.5_f64.to_radians(),
            2017.0,
            CatalogFilter::default(),
        );
        assert_eq!(hits.len(), 1);
        let hits = cat.cone_search_at_epoch(
            10.0_f64.to_radians(),
            0.0,
            0.5_f64.to_radians(),
            2016.0,
            CatalogFilter::default(),
        );
        assert!(hits.is_empty());
    }

    #[test]
    fn chunked_search_matches_inmemory() {
        let mut cat = LargeStarCatalog::new();
        for i in 0..20 {
            cat.push(rec(i, i as f64, 0.0));
        }
        let full = cat.cone_search(
            5.0_f64.to_radians(),
            0.0,
            3.0_f64.to_radians(),
            CatalogFilter::default(),
        );
        let chunked: Vec<_> = cat
            .cone_search_chunked(
                5.0_f64.to_radians(),
                0.0,
                3.0_f64.to_radians(),
                CatalogFilter::default(),
                4,
            )
            .flat_map(|c| c.into_iter().map(|r| r.source_id))
            .collect();
        let full_ids: Vec<_> = full.iter().map(|r| r.source_id).collect();
        assert_eq!(chunked, full_ids);
    }

    #[test]
    fn csv_parses_required_and_optional_columns() {
        let input = "# comment\nsource_id,ra,dec,epoch,pmra,pmdec,parallax,g_mag,quality\n42,10.0,0.0,2016.0,1.0,2.0,5.0,8.5,1\n43,20.0,1.0,2016.0,,,,,0";
        let recs = parse_csv_chunk(input).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].source_id, 42);
        assert!((recs[0].ra - 10.0_f64.to_radians()).abs() < 1e-12);
        assert_eq!(recs[0].pm_ra_cosdec, Some(1.0));
        assert!(recs[0].quality_ok);
        assert_eq!(recs[1].pm_ra_cosdec, None);
        assert!(!recs[1].quality_ok);
    }

    #[test]
    fn csv_rejects_missing_required_column() {
        let input = "source_id,ra\n1,10.0";
        assert!(parse_csv_chunk(input).is_err());
    }
}
