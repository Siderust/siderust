// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Large catalog ingestion, cone search, and observatory constants
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
//! Ground-based astronomical observations are referenced to a *site*: a
//! geodetic position on the WGS84 ellipsoid plus a typical local
//! atmospheric state (pressure, temperature, relative humidity). These
//! quantities feed every site-dependent computation in the crate —
//! airmass, refraction, atmospheric extinction, night-sky brightness,
//! parallactic correction.
//!
//! ## Technical scope
//!
//! - [`CatalogRecord`] is the column-aligned star record with typed positional
//!   fields ([`crate::qtty::Radians`], [`crate::qtty::MilliArcseconds`], [`crate::qtty::KmPerSeconds`]).
//! - [`LargeStarCatalog`] is the in-memory container with cone search
//!   (`cone_search`), magnitude/quality filters, and an epoch-aware
//!   query that propagates positions via proper motion.
//! - [`parse_csv_chunk`] consumes a CSV-style ASCII chunk into a
//!   `Vec<CatalogRecord>`.
//! - [`Observatory`] carries site name, geodetic position, and reference
//!   atmospheric conditions; named constants for major sites live in
//!   the `observatories` submodule.
//!
//! ## References
//!
//! - Lindegren, L. (2018). *Gaia DR2: The catalogue's astrometric
//!   content*. A&A 616 A2.
//! - National Imagery and Mapping Agency (2000). *WGS84: Its Definition and
//!   Relationships with Local Geodetic Systems*. NIMA TR8350.2, 3rd ed.

#![forbid(unsafe_code)]

pub mod catalog;
pub mod ingest;
pub mod observatories;
pub mod record;

pub use catalog::LargeStarCatalog;
pub use ingest::{parse_csv_chunk, CatalogIngestError};
pub use observatories::Observatory;
pub use record::{CatalogFilter, CatalogRecord};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::Radians;
    use core::f64::consts::PI;

    fn rec(id: u64, ra_deg: f64, dec_deg: f64) -> CatalogRecord {
        CatalogRecord {
            source_id: id,
            ra: Radians::new(ra_deg.to_radians()),
            dec: Radians::new(dec_deg.to_radians()),
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

    fn rad(x: f64) -> Radians {
        Radians::new(x)
    }

    #[test]
    fn cone_includes_center() {
        let mut cat = LargeStarCatalog::new();
        cat.push(rec(1, 10.0, 0.0));
        let hits = cat.cone_search(
            rad(10.0_f64.to_radians()),
            rad(0.0),
            rad(1.0_f64.to_radians()),
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
            rad(10.0_f64.to_radians()),
            rad(0.0),
            rad(5.0_f64.to_radians()),
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
        let hits = cat.cone_search(
            rad(0.0),
            rad(0.0),
            rad(2.0_f64.to_radians()),
            CatalogFilter::default(),
        );
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn cone_near_pole_includes_records() {
        let mut cat = LargeStarCatalog::new();
        cat.push(rec(1, 0.0, 89.5));
        cat.push(rec(2, 90.0, 89.5));
        cat.push(rec(3, 180.0, 89.5));
        cat.push(rec(4, 270.0, 89.5));
        let hits = cat.cone_search(
            rad(0.0),
            rad((PI / 2.0) - 0.001),
            rad(2.0_f64.to_radians()),
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
            rad(10.0_f64.to_radians()),
            rad(0.0),
            rad(1.0_f64.to_radians()),
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
        let mut r = rec(1, 0.0, 0.0);
        r.pm_ra_cosdec = Some(10.0 * 3_600_000.0);
        r.pm_dec = Some(0.0);
        let mut cat = LargeStarCatalog::new();
        cat.push(r);
        let hits = cat.cone_search_at_epoch(
            rad(10.0_f64.to_radians()),
            rad(0.0),
            rad(0.5_f64.to_radians()),
            2017.0,
            CatalogFilter::default(),
        );
        assert_eq!(hits.len(), 1);
        let hits = cat.cone_search_at_epoch(
            rad(10.0_f64.to_radians()),
            rad(0.0),
            rad(0.5_f64.to_radians()),
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
            rad(5.0_f64.to_radians()),
            rad(0.0),
            rad(3.0_f64.to_radians()),
            CatalogFilter::default(),
        );
        let chunked: Vec<_> = cat
            .cone_search_chunked(
                rad(5.0_f64.to_radians()),
                rad(0.0),
                rad(3.0_f64.to_radians()),
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
        assert!((recs[0].ra.value() - 10.0_f64.to_radians()).abs() < 1e-12);
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
