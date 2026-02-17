// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Of-Date Ecliptic Coordinate Transformations
//!
//! This module provides traits for converting between equatorial and ecliptic
//! coordinates using the mean ecliptic plane of date (time-dependent obliquity).
//!
//! ## Approach
//!
//! The transformation uses:
//! - IAU 2006 precession to compute the mean equatorial frame of date
//! - Mean obliquity of date (without nutation)
//! - Rotation about the X-axis by the mean obliquity angle
//!
//! ## Time Requirements
//!
//! Ecliptic-of-date transformations require **TT** (Terrestrial Time) for:
//! - Precession calculations (IAU 2006 Fukushima-Williams)
//! - Mean obliquity evaluation
//!
//! ## Frame Compatibility
//!
//! These transformations connect:
//! - [`EquatorialMeanOfDate`] ↔ [`EclipticOfDate`]
//! - [`ICRS`]/[`GCRS`] ↔ [`EclipticOfDate`] (via precession matrix)
//!
//! Note: For J2000 ecliptic coordinates, use the time-independent [`Ecliptic`]
//! frame with [`TransformFrame`](crate::coordinates::transform::TransformFrame).
//!
//! ## Usage
//!
//! ```rust
//! use siderust::coordinates::cartesian::Direction;
//! use siderust::coordinates::frames::{EquatorialMeanOfDate, EclipticOfDate};
//! use siderust::coordinates::spherical;
//! use siderust::coordinates::transform::ecliptic_of_date::ToEclipticOfDate;
//! use siderust::time::JulianDate;
//! use qtty::*;
//!
//! let jd_tt = JulianDate::new(2_451_545.0);
//! let equatorial = spherical::Direction::<EquatorialMeanOfDate>::new(45.0 * DEG, 30.0 * DEG)
//!     .to_cartesian();
//!
//! let ecliptic: Direction<EclipticOfDate> = equatorial.to_ecliptic_of_date(&jd_tt);
//! ```

use crate::astro::precession;
use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::{EclipticOfDate, EquatorialMeanOfDate, GCRS, ICRS};
use crate::time::JulianDate;
use qtty::{Degrees, Radians};
use std::f64::consts::TAU;

// =============================================================================
// ToEclipticOfDate Trait
// =============================================================================

/// Convert coordinates to ecliptic-of-date frame using IAU 2006 precession.
///
/// This trait transforms directions to the mean ecliptic plane of date,
/// which is time-dependent due to precession.
pub trait ToEclipticOfDate {
    /// Convert this direction to ecliptic-of-date coordinates.
    ///
    /// # Arguments
    ///
    /// - `jd_tt`: Julian Date on the TT timescale (for precession/obliquity).
    ///
    /// # Returns
    ///
    /// A [`Direction<EclipticOfDate>`] with:
    /// - Ecliptic longitude normalized to `[0°, 360°)`
    /// - Ecliptic latitude in `[-90°, +90°]`
    fn to_ecliptic_of_date(&self, jd_tt: &JulianDate) -> Direction<EclipticOfDate>;
}

impl ToEclipticOfDate for Direction<EquatorialMeanOfDate> {
    #[inline]
    fn to_ecliptic_of_date(&self, jd_tt: &JulianDate) -> Direction<EclipticOfDate> {
        // Extract RA/Dec
        let spherical = self.to_spherical();
        let ra = Radians::from(spherical.azimuth);
        let dec = Radians::from(spherical.polar);

        // Convert to Cartesian
        let (sin_dec, cos_dec) = dec.sin_cos();
        let (sin_ra, cos_ra) = ra.sin_cos();
        let v_eq = [cos_dec * cos_ra, cos_dec * sin_ra, sin_dec];

        // Apply mean-of-date to ecliptic-of-date rotation
        let rot = precession::mean_equatorial_to_ecliptic_of_date_matrix(*jd_tt);
        let v_ecl = rot.apply_array(v_eq);

        // Convert back to spherical
        let lon = v_ecl[1].atan2(v_ecl[0]).rem_euclid(TAU);
        let lat = v_ecl[2].clamp(-1.0, 1.0).asin();

        let spherical_ecl = affn::spherical::Direction::<EclipticOfDate>::new_raw(
            Degrees::new(lat.to_degrees()),
            Degrees::new(lon.to_degrees()),
        );
        spherical_ecl.to_cartesian()
    }
}

impl ToEclipticOfDate for Direction<ICRS> {
    #[inline]
    fn to_ecliptic_of_date(&self, jd_tt: &JulianDate) -> Direction<EclipticOfDate> {
        // Extract RA/Dec
        let spherical = self.to_spherical();
        let ra = Radians::from(spherical.azimuth);
        let dec = Radians::from(spherical.polar);

        // Convert to Cartesian
        let (sin_dec, cos_dec) = dec.sin_cos();
        let (sin_ra, cos_ra) = ra.sin_cos();
        let v_eq = [cos_dec * cos_ra, cos_dec * sin_ra, sin_dec];

        // Apply ICRS/GCRS to ecliptic-of-date rotation
        let rot = precession::gcrs_to_ecliptic_of_date_matrix(*jd_tt);
        let v_ecl = rot.apply_array(v_eq);

        // Convert back to spherical
        let lon = v_ecl[1].atan2(v_ecl[0]).rem_euclid(TAU);
        let lat = v_ecl[2].clamp(-1.0, 1.0).asin();

        let spherical_ecl = affn::spherical::Direction::<EclipticOfDate>::new_raw(
            Degrees::new(lat.to_degrees()),
            Degrees::new(lon.to_degrees()),
        );
        spherical_ecl.to_cartesian()
    }
}

impl ToEclipticOfDate for Direction<GCRS> {
    #[inline]
    fn to_ecliptic_of_date(&self, jd_tt: &JulianDate) -> Direction<EclipticOfDate> {
        // Extract RA/Dec
        let spherical = self.to_spherical();
        let ra = Radians::from(spherical.azimuth);
        let dec = Radians::from(spherical.polar);

        // Convert to Cartesian
        let (sin_dec, cos_dec) = dec.sin_cos();
        let (sin_ra, cos_ra) = ra.sin_cos();
        let v_eq = [cos_dec * cos_ra, cos_dec * sin_ra, sin_dec];

        // Apply GCRS to ecliptic-of-date rotation
        let rot = precession::gcrs_to_ecliptic_of_date_matrix(*jd_tt);
        let v_ecl = rot.apply_array(v_eq);

        // Convert back to spherical
        let lon = v_ecl[1].atan2(v_ecl[0]).rem_euclid(TAU);
        let lat = v_ecl[2].clamp(-1.0, 1.0).asin();

        let spherical_ecl = affn::spherical::Direction::<EclipticOfDate>::new_raw(
            Degrees::new(lat.to_degrees()),
            Degrees::new(lon.to_degrees()),
        );
        spherical_ecl.to_cartesian()
    }
}

// =============================================================================
// FromEclipticOfDate Trait
// =============================================================================

/// Convert coordinates from ecliptic-of-date frame using IAU 2006 precession.
///
/// This trait transforms directions from the mean ecliptic plane of date
/// back to equatorial frames.
pub trait FromEclipticOfDate {
    /// Convert this ecliptic-of-date direction to equatorial mean-of-date coordinates.
    ///
    /// # Arguments
    ///
    /// - `jd_tt`: Julian Date on the TT timescale (for precession/obliquity).
    ///
    /// # Returns
    ///
    /// A [`Direction<EquatorialMeanOfDate>`] with:
    /// - Right ascension normalized to `[0°, 360°)`
    /// - Declination in `[-90°, +90°]`
    fn to_equatorial_mean_of_date(&self, jd_tt: &JulianDate) -> Direction<EquatorialMeanOfDate>;

    /// Convert this ecliptic-of-date direction to ICRS/GCRS coordinates.
    ///
    /// # Arguments
    ///
    /// - `jd_tt`: Julian Date on the TT timescale (for precession).
    ///
    /// # Returns
    ///
    /// A [`Direction<ICRS>`] with:
    /// - Right ascension normalized to `[0°, 360°)`
    /// - Declination in `[-90°, +90°]`
    fn to_icrs(&self, jd_tt: &JulianDate) -> Direction<ICRS>;
}

impl FromEclipticOfDate for Direction<EclipticOfDate> {
    #[inline]
    fn to_equatorial_mean_of_date(&self, jd_tt: &JulianDate) -> Direction<EquatorialMeanOfDate> {
        // Extract lon/lat
        let spherical = self.to_spherical();
        let lon = Radians::from(spherical.azimuth);
        let lat = Radians::from(spherical.polar);

        // Convert to Cartesian
        let (sin_lat, cos_lat) = lat.sin_cos();
        let (sin_lon, cos_lon) = lon.sin_cos();
        let v_ecl = [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat];

        // Apply ecliptic-of-date to mean-of-date rotation
        let rot = precession::ecliptic_of_date_to_mean_equatorial_matrix(*jd_tt);
        let v_eq = rot.apply_array(v_ecl);

        // Convert back to spherical
        let ra = v_eq[1].atan2(v_eq[0]).rem_euclid(TAU);
        let dec = v_eq[2].clamp(-1.0, 1.0).asin();

        let spherical_equ = affn::spherical::Direction::<EquatorialMeanOfDate>::new_raw(
            Degrees::new(dec.to_degrees()),
            Degrees::new(ra.to_degrees()),
        );
        spherical_equ.to_cartesian()
    }

    #[inline]
    fn to_icrs(&self, jd_tt: &JulianDate) -> Direction<ICRS> {
        // Extract lon/lat
        let spherical = self.to_spherical();
        let lon = Radians::from(spherical.azimuth);
        let lat = Radians::from(spherical.polar);

        // Convert to Cartesian
        let (sin_lat, cos_lat) = lat.sin_cos();
        let (sin_lon, cos_lon) = lon.sin_cos();
        let v_ecl = [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat];

        // Apply ecliptic-of-date to GCRS/ICRS rotation
        let rot = precession::ecliptic_of_date_to_gcrs_matrix(*jd_tt);
        let v_eq = rot.apply_array(v_ecl);

        // Convert back to spherical
        let ra = v_eq[1].atan2(v_eq[0]).rem_euclid(TAU);
        let dec = v_eq[2].clamp(-1.0, 1.0).asin();

        let spherical_equ = affn::spherical::Direction::<ICRS>::new_raw(
            Degrees::new(dec.to_degrees()),
            Degrees::new(ra.to_degrees()),
        );
        spherical_equ.to_cartesian()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::*;

    #[test]
    fn roundtrip_mean_of_date_is_stable() {
        let jd_tt = JulianDate::new(2_451_545.0);

        let ra = 1.0 * RAD;
        let dec = 0.5 * RAD;

        let spherical_equ = affn::spherical::Direction::<EquatorialMeanOfDate>::new_raw(
            Degrees::new(dec.value().to_degrees()),
            Degrees::new(ra.value().to_degrees())
        );
        let equatorial = spherical_equ.to_cartesian();
        let ecliptic = equatorial.to_ecliptic_of_date(&jd_tt);
        let back = ecliptic.to_equatorial_mean_of_date(&jd_tt);

        let orig_sph = equatorial.to_spherical();
        let back_sph = back.to_spherical();

        let dra = (Radians::from(back_sph.azimuth).value() - Radians::from(orig_sph.azimuth).value()).abs();
        let dra = dra.min(TAU - dra);
        let ddec = (Radians::from(back_sph.polar).value() - Radians::from(orig_sph.polar).value()).abs();

        assert!(dra < 1e-12, "RA roundtrip error too large: {dra}");
        assert!(ddec < 1e-12, "Dec roundtrip error too large: {ddec}");
    }

    #[test]
    fn roundtrip_icrs_is_stable() {
        let jd_tt = JulianDate::new(2_451_545.0);

        let ra = 1.0 * RAD;
        let dec = 0.5 * RAD;

        let spherical_icrs = affn::spherical::Direction::<ICRS>::new_raw(Degrees::new(dec.value().to_degrees()), Degrees::new(ra.value().to_degrees()));
        let icrs = spherical_icrs.to_cartesian();
        let ecliptic = icrs.to_ecliptic_of_date(&jd_tt);
        let back = ecliptic.to_icrs(&jd_tt);

        let orig_sph = icrs.to_spherical();
        let back_sph = back.to_spherical();

        let dra = (Radians::from(back_sph.azimuth).value() - Radians::from(orig_sph.azimuth).value()).abs();
        let dra = dra.min(TAU - dra);
        let ddec = (Radians::from(back_sph.polar).value() - Radians::from(orig_sph.polar).value()).abs();

        assert!(dra < 1e-12, "RA roundtrip error too large: {dra}");
        assert!(ddec < 1e-12, "Dec roundtrip error too large: {ddec}");
    }

    #[test]
    fn ecliptic_of_date_has_correct_obliquity() {
        // At J2000, the mean obliquity should be approximately 23.439279444444445°
        let jd_tt = JulianDate::J2000;

        // A point on the equatorial equator (RA = 0°, Dec = 0°) should be on the ecliptic
        let spherical_equ = affn::spherical::Direction::<EquatorialMeanOfDate>::new_raw(0.0 * DEG, 0.0 * DEG);
        let equatorial = spherical_equ.to_cartesian();
        let ecliptic = equatorial.to_ecliptic_of_date(&jd_tt);

        let sph = ecliptic.to_spherical();
        // Should be at ecliptic longitude ~0°, latitude ~0°
        assert!(Radians::from(sph.polar).abs() < 1e-10 * RAD, "Expected ecliptic latitude near 0");

        // A point at the north celestial pole (Dec = 90°) should be at the ecliptic pole
        let spherical_north = affn::spherical::Direction::<EquatorialMeanOfDate>::new_raw(90.0 * DEG, 0.0 * DEG);
        let north_pole = spherical_north.to_cartesian();
        let ecl_north = north_pole.to_ecliptic_of_date(&jd_tt);

        let sph_north = ecl_north.to_spherical();
        // Should be at ecliptic north pole (lat = 90° - obliquity)
        let obliquity_rad = 23.439279444444445_f64.to_radians();
        let expected_lat = std::f64::consts::FRAC_PI_2 - obliquity_rad;

        assert!(
            (Radians::from(sph_north.polar).value() - expected_lat).abs() < 1e-6,
            "Ecliptic latitude mismatch at north pole"
        );
    }
}
