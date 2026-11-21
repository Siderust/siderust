// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Catalog record type and filter
//!
//! [`CatalogRecord`] is the column-aligned row type for large stellar catalogs
//! (Gaia-style). Position is stored in [`crate::qtty::Radians`], and parallax
//! and radial velocity use typed quantities. Proper-motion components are
//! stored as scalars in milliarcseconds per Julian year.

use crate::qtty::{KmPerSeconds, MilliArcseconds, Radian, Radians};
#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Scalar, Transcendental};
use core::f64::consts::{FRAC_PI_2, TAU};

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
    /// Propagate the ICRS direction to a Julian-year epoch.
    ///
    /// Applies the local east (μᵅ⋆) and north (μᵟ) displacements linearly in
    /// the tangent plane and projects the result back onto the unit sphere.
    /// At a pole the stored RA defines the orientation of the tangent basis.
    /// This approximation omits parallax, radial velocity and perspective
    /// acceleration. The returned pair is right ascension and declination in
    /// radians. If either proper-motion component is missing, the input pair
    /// is returned unchanged. A non-finite displacement returns NaN angles.
    pub fn propagate_to(&self, target_epoch_jyr: f64) -> (Radians, Radians) {
        let (Some(mu_a), Some(mu_d)) = (self.pm_ra_cosdec, self.pm_dec) else {
            return (self.ra, self.dec);
        };
        let dt_yr = target_epoch_jyr - self.epoch_jyr;
        let east = MilliArcseconds::new(mu_a).to::<Radian>().value() * dt_yr;
        let north = MilliArcseconds::new(mu_d).to::<Radian>().value() * dt_yr;
        if !east.is_finite() || !north.is_finite() {
            return (Radians::new(f64::NAN), Radians::new(f64::NAN));
        }
        if east == 0.0 && north == 0.0 {
            return (self.ra, self.dec);
        }

        // Work in a Cartesian basis rotated by the reference RA about z.
        // r = (cos δ, 0, sin δ), east = (0, 1, 0), north = (-sin δ, 0, cos δ).
        let (sin_d, mut cos_d) = self.dec.sin_cos();
        if self.dec.value().abs() == FRAC_PI_2 {
            cos_d = 0.0;
        }
        let x = cos_d - north * sin_d;
        let z = sin_d + north * cos_d;
        let horizontal = x.hypot(east);
        if !horizontal.is_finite() {
            return (Radians::new(f64::NAN), Radians::new(f64::NAN));
        }
        let ra = (self.ra.value() + east.atan2(x)).rem_euclid(TAU);
        let dec = z.atan2(horizontal);
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

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    #[test]
    fn meridional_propagation_crosses_the_pole() {
        let offset = MilliArcseconds::new(2000.0).to::<Radian>().value().atan();
        let row = CatalogRecord {
            source_id: 1,
            ra: Radians::new(0.0),
            dec: Radians::new(FRAC_PI_2 - offset / 2.0),
            epoch_jyr: 2016.0,
            pm_ra_cosdec: Some(0.0),
            pm_dec: Some(2000.0),
            parallax: None,
            radial_velocity: None,
            g_mag: None,
            bp_mag: None,
            rp_mag: None,
            quality_ok: true,
        };
        let (ra, dec) = row.propagate_to(2017.0);
        assert_eq!(ra.value(), PI);
        assert!((dec.value() - row.dec.value()).abs() < 5e-16);

        let (ra, dec) = row.propagate_to(2015.0);
        assert_eq!(ra.value(), 0.0);
        assert!((dec.value() - (row.dec.value() - offset)).abs() < 5e-16);
    }
}
