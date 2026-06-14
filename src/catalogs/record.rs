// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Catalog record type and filter
//!
//! [`CatalogRecord`] is the column-aligned row type for large stellar catalogs
//! (Gaia-style). Position is stored as typed [`crate::qtty::Radians`] scalars; parallax and
//! radial velocity use typed quantities. Proper-motion components (mas/yr)
//! remain plain `f64` until a `MilliArcSecondPerYear` unit is added to `qtty`.

use crate::qtty::{KmPerSeconds, MilliArcseconds, Radians};
use core::f64::consts::TAU;

/// One row of a large stellar catalog.
#[derive(Debug, Clone, PartialEq)]
pub struct CatalogRecord {
    /// Stable source identifier (e.g. Gaia source_id).
    pub source_id: u64,
    /// ICRS right ascension at [`Self::epoch_jyr`].
    pub ra: Radians,
    /// ICRS declination at [`Self::epoch_jyr`].
    pub dec: Radians,
    /// Catalog epoch in Julian years (e.g. 2016.0 for Gaia DR3).
    pub epoch_jyr: f64,
    /// Proper motion in right ascension × cos(δ), mas/yr.
    pub pm_ra_cosdec: Option<f64>,
    /// Proper motion in declination, mas/yr.
    pub pm_dec: Option<f64>,
    /// Annual trigonometric parallax.
    pub parallax: Option<MilliArcseconds>,
    /// Radial (line-of-sight) velocity.
    pub radial_velocity: Option<KmPerSeconds>,
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
    /// Propagate the record's ICRS direction by linear proper motion to the
    /// requested Julian-year epoch. Returns `(ra, dec)` in radians. If proper
    /// motion is unavailable the original `(ra, dec)` is returned unchanged.
    pub fn propagate_to(&self, target_epoch_jyr: f64) -> (Radians, Radians) {
        let dt_yr = target_epoch_jyr - self.epoch_jyr;
        let (Some(mu_a), Some(mu_d)) = (self.pm_ra_cosdec, self.pm_dec) else {
            return (self.ra, self.dec);
        };
        let mas_to_rad = (core::f64::consts::PI / 180.0) / 3_600_000.0;
        let cos_d = self.dec.value().cos();
        let d_ra = if cos_d.abs() > 1e-12 {
            (mu_a * mas_to_rad * dt_yr) / cos_d
        } else {
            0.0
        };
        let d_dec = mu_d * mas_to_rad * dt_yr;
        let mut ra = (self.ra.value() + d_ra).rem_euclid(TAU);
        let mut dec = self.dec.value() + d_dec;
        if dec > core::f64::consts::FRAC_PI_2 {
            dec = core::f64::consts::PI - dec;
            ra = (ra + core::f64::consts::PI).rem_euclid(TAU);
        } else if dec < -core::f64::consts::FRAC_PI_2 {
            dec = -core::f64::consts::PI - dec;
            ra = (ra + core::f64::consts::PI).rem_euclid(TAU);
        }
        (Radians::new(ra), Radians::new(dec))
    }
}

/// Constraints for a cone-search query.
#[derive(Debug, Clone, Copy, Default)]
pub struct CatalogFilter {
    /// If set, exclude records with `g_mag` greater than this.
    pub max_g_mag: Option<f64>,
    /// If set, only include records whose `quality_ok` matches.
    pub require_quality: Option<bool>,
}

pub(super) fn passes_filter(r: &CatalogRecord, f: &CatalogFilter) -> bool {
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

/// Haversine angular distance test on the unit sphere.
pub(super) fn inside_cone(
    ra: Radians,
    dec: Radians,
    c_ra: Radians,
    c_dec: Radians,
    radius: Radians,
) -> bool {
    let d_ra = ra.value() - c_ra.value();
    let s = (dec.value() - c_dec.value()).sin();
    let s2 = (d_ra / 2.0).sin();
    let a = (s / 2.0).powi(2) + dec.value().cos() * c_dec.value().cos() * s2 * s2;
    let c = 2.0 * a.sqrt().min(1.0).asin();
    c <= radius.value()
}
