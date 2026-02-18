// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # IAU 2006/GAST Horizontal Coordinate Transformations
//!
//! This module provides traits for converting between equatorial and horizontal
//! coordinates using IAU 2006 precession-nutation and GAST-based sidereal time.
//!
//! ## Approach
//!
//! The transformation follows ERFA conventions:
//! - GAST (Greenwich Apparent Sidereal Time) = `gast_iau2006(UT1, TT, dpsi, true_obliquity)`
//! - Hour Angle = GAST + observer_longitude - RA
//! - Azimuth uses astronomical convention: 0° = North, increasing clockwise through East
//!
//! ## Time Requirements
//!
//! Horizontal transformations require **both** UT1 and TT timescales:
//! - **UT1** for Earth rotation angle and sidereal time
//! - **TT** for precession/nutation models (IAU 2000B nutation)
//!
//! ## Observer Location
//!
//! The observer's geodetic location is provided via [`Geodetic<ECEF>`], which contains:
//! - Longitude (positive eastward)
//! - Latitude (positive northward)  
//! - Height above WGS84 ellipsoid
//!
//! ## Usage
//!
//! ```rust
//! use siderust::coordinates::cartesian::Direction;
//! use siderust::coordinates::centers::Geodetic;
//! use siderust::coordinates::frames::{ECEF, EquatorialTrueOfDate, Horizontal};
//! use siderust::coordinates::spherical;
//! use siderust::coordinates::transform::horizontal::ToHorizontal;
//! use siderust::time::JulianDate;
//! use qtty::*;
//!
//! let jd_ut1 = JulianDate::new(2_451_545.0);
//! let jd_tt = JulianDate::new(2_451_545.000_800_741);
//!
//! let equatorial = spherical::Direction::<EquatorialTrueOfDate>::new(45.0 * DEG, 30.0 * DEG)
//!     .to_cartesian();
//! let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.5 * DEG, 0.0 * M);
//!
//! let horizontal: Direction<Horizontal> = equatorial.to_horizontal(&jd_ut1, &jd_tt, &site);
//! ```

use crate::astro::{nutation, sidereal};
use crate::coordinates::cartesian::{Direction, Position};
use crate::coordinates::centers::{Geodetic, Topocentric};
use crate::coordinates::frames::{EquatorialTrueOfDate, Horizontal, ECEF};
use crate::time::JulianDate;
use qtty::{Degrees, LengthUnit, Radians};
use std::f64::consts::TAU;

// =============================================================================
// ToHorizontal Trait
// =============================================================================

/// Convert equatorial (true of date) coordinates to local horizontal coordinates.
///
/// This trait transforms directions from the [`EquatorialTrueOfDate`] frame to the
/// [`Horizontal`] frame using IAU 2006/GAST sidereal time.
///
/// The transformation is observer-dependent and requires geodetic location information.
pub trait ToHorizontal {
    /// Convert this equatorial direction to horizontal coordinates.
    ///
    /// # Arguments
    ///
    /// - `jd_ut1`: Julian Date on the UT1 timescale (for Earth rotation).
    /// - `jd_tt`: Julian Date on the TT timescale (for precession/nutation).
    /// - `site`: Observer's geodetic location.
    ///
    /// # Returns
    ///
    /// A [`Direction<Horizontal>`] with:
    /// - Azimuth normalized to `[0°, 360°)` (North = 0°, increasing clockwise)
    /// - Altitude in `[-90°, +90°]`
    fn to_horizontal(
        &self,
        jd_ut1: &JulianDate,
        jd_tt: &JulianDate,
        site: &Geodetic<ECEF>,
    ) -> Direction<Horizontal>;
}

impl ToHorizontal for Direction<EquatorialTrueOfDate> {
    #[inline]
    fn to_horizontal(
        &self,
        jd_ut1: &JulianDate,
        jd_tt: &JulianDate,
        site: &Geodetic<ECEF>,
    ) -> Direction<Horizontal> {
        // Extract RA/Dec from direction
        let spherical = self.to_spherical();
        let ra = Radians::from(spherical.azimuth);
        let dec = Radians::from(spherical.polar);

        // Convert observer location to radians
        let obs_lon = Radians::from(site.lon);
        let obs_lat = Radians::from(site.lat);

        // Compute GAST and hour angle
        let nut = nutation::nutation_iau2000b(*jd_tt);
        let gast = sidereal::gast_iau2006(*jd_ut1, *jd_tt, nut.dpsi, nut.true_obliquity());
        let ha = gast.value() + obs_lon.value() - ra.value();

        // Spherical trigonometry for equatorial → horizontal
        let (sh, ch) = ha.sin_cos();
        let (sd, cd) = dec.sin_cos();
        let (sp, cp) = obs_lat.sin_cos();

        // ERFA-compatible intermediate axes
        let x = -ch * cd * sp + sd * cp;
        let y = -sh * cd;
        let z = ch * cd * cp + sd * sp;
        let r = (x * x + y * y).sqrt();

        let az = if r != 0.0 { y.atan2(x) } else { 0.0 }.rem_euclid(TAU);
        let alt = z.atan2(r);

        // Construct horizontal direction from spherical coordinates
        let spherical_horiz = affn::spherical::Direction::<Horizontal>::new_raw(
            Degrees::new(alt.to_degrees()),
            Degrees::new(az.to_degrees()),
        );
        spherical_horiz.to_cartesian()
    }
}

// =============================================================================
// FromHorizontal Trait
// =============================================================================

/// Convert local horizontal coordinates to equatorial (true of date) coordinates.
///
/// This trait transforms directions from the [`Horizontal`] frame to the
/// [`EquatorialTrueOfDate`] frame using IAU 2006/GAST sidereal time.
///
/// The transformation is observer-dependent and requires geodetic location information.
pub trait FromHorizontal {
    /// Convert this horizontal direction to equatorial coordinates.
    ///
    /// # Arguments
    ///
    /// - `jd_ut1`: Julian Date on the UT1 timescale (for Earth rotation).
    /// - `jd_tt`: Julian Date on the TT timescale (for precession/nutation).
    /// - `site`: Observer's geodetic location.
    ///
    /// # Returns
    ///
    /// A [`Direction<EquatorialTrueOfDate>`] with:
    /// - Right ascension normalized to `[0°, 360°)`
    /// - Declination in `[-90°, +90°]`
    fn to_equatorial(
        &self,
        jd_ut1: &JulianDate,
        jd_tt: &JulianDate,
        site: &Geodetic<ECEF>,
    ) -> Direction<EquatorialTrueOfDate>;
}

impl FromHorizontal for Direction<Horizontal> {
    #[inline]
    fn to_equatorial(
        &self,
        jd_ut1: &JulianDate,
        jd_tt: &JulianDate,
        site: &Geodetic<ECEF>,
    ) -> Direction<EquatorialTrueOfDate> {
        // Extract Az/Alt from direction
        let spherical = self.to_spherical();
        let az = Radians::from(spherical.azimuth);
        let alt = Radians::from(spherical.polar);

        // Convert observer location to radians
        let obs_lon = Radians::from(site.lon);
        let obs_lat = Radians::from(site.lat);

        // Compute LAST (Local Apparent Sidereal Time)
        let nut = nutation::nutation_iau2000b(*jd_tt);
        let gast = sidereal::gast_iau2006(*jd_ut1, *jd_tt, nut.dpsi, nut.true_obliquity());
        let last = gast.value() + obs_lon.value();

        // Spherical trigonometry for horizontal → equatorial
        let (sp, cp) = obs_lat.sin_cos();
        let (sa, ca) = alt.sin_cos();
        let (sz, cz) = az.sin_cos();

        let dec = (sp * sa + cp * ca * cz).clamp(-1.0, 1.0).asin();
        let ha = (-sz * ca).atan2(sa * cp - ca * cz * sp);
        let ra = (last - ha).rem_euclid(TAU);

        // Construct equatorial direction from spherical coordinates
        let spherical_equ = affn::spherical::Direction::<EquatorialTrueOfDate>::new_raw(
            Degrees::new(dec.to_degrees()),
            Degrees::new(ra.to_degrees()),
        );
        spherical_equ.to_cartesian()
    }
}

// =============================================================================
// Topocentric Position Extensions
// =============================================================================

/// Extension trait for [`Position<Topocentric, EquatorialTrueOfDate, U>`].
///
/// This trait provides a convenience method for transforming topocentric equatorial
/// positions to horizontal coordinates, automatically extracting the observer site
/// from the position's center parameters.
pub trait TopocentricEquatorialExt<U: LengthUnit> {
    /// Convert this topocentric equatorial position to horizontal coordinates.
    ///
    /// The observer site is extracted from the position's [`Topocentric`] center parameters.
    ///
    /// # Arguments
    ///
    /// - `jd_ut1`: Julian Date on the UT1 timescale (for Earth rotation).
    /// - `jd_tt`: Julian Date on the TT timescale (for precession/nutation).
    ///
    /// # Returns
    ///
    /// A [`Position<Topocentric, Horizontal, U>`] with the same distance and observer site.
    fn to_horizontal_position(
        &self,
        jd_ut1: &JulianDate,
        jd_tt: &JulianDate,
    ) -> Position<Topocentric, Horizontal, U>;
}

impl<U: LengthUnit> TopocentricEquatorialExt<U> for Position<Topocentric, EquatorialTrueOfDate, U> {
    #[inline]
    fn to_horizontal_position(
        &self,
        jd_ut1: &JulianDate,
        jd_tt: &JulianDate,
    ) -> Position<Topocentric, Horizontal, U> {
        let site = *self.center_params();
        // Convert to spherical to get direction
        let sph = affn::Position::<Topocentric, EquatorialTrueOfDate, U>::to_spherical(self);
        let equ_dir = sph.direction().to_cartesian();
        let distance = sph.distance;

        // Transform direction
        let horiz_cart = equ_dir.to_horizontal(jd_ut1, jd_tt, &site);
        let horiz_sph = horiz_cart.to_spherical();

        // Create spherical position and convert back to cartesian
        let result_sph = horiz_sph.position_with_params::<Topocentric, U>(site, distance);
        affn::Position::<Topocentric, Horizontal, U>::from_spherical(&result_sph)
    }
}

/// Extension trait for [`Position<Topocentric, Horizontal, U>`].
///
/// This trait provides a convenience method for transforming topocentric horizontal
/// positions to equatorial coordinates, automatically extracting the observer site
/// from the position's center parameters.
pub trait TopocentricHorizontalExt<U: LengthUnit> {
    /// Convert this topocentric horizontal position to equatorial coordinates.
    ///
    /// The observer site is extracted from the position's [`Topocentric`] center parameters.
    ///
    /// # Arguments
    ///
    /// - `jd_ut1`: Julian Date on the UT1 timescale (for Earth rotation).
    /// - `jd_tt`: Julian Date on the TT timescale (for precession/nutation).
    ///
    /// # Returns
    ///
    /// A [`Position<Topocentric, EquatorialTrueOfDate, U>`] with the same distance and observer site.
    fn to_equatorial_position(
        &self,
        jd_ut1: &JulianDate,
        jd_tt: &JulianDate,
    ) -> Position<Topocentric, EquatorialTrueOfDate, U>;
}

impl<U: LengthUnit> TopocentricHorizontalExt<U> for Position<Topocentric, Horizontal, U> {
    #[inline]
    fn to_equatorial_position(
        &self,
        jd_ut1: &JulianDate,
        jd_tt: &JulianDate,
    ) -> Position<Topocentric, EquatorialTrueOfDate, U> {
        let site = *self.center_params();
        // Convert to spherical to get direction
        let sph = affn::Position::<Topocentric, Horizontal, U>::to_spherical(self);
        let horiz_dir = sph.direction().to_cartesian();
        let distance = sph.distance;

        // Transform direction
        let equ_cart = horiz_dir.to_equatorial(jd_ut1, jd_tt, &site);
        let equ_sph = equ_cart.to_spherical();

        // Create spherical position and convert back to cartesian
        let result_sph = equ_sph.position_with_params::<Topocentric, U>(site, distance);
        affn::Position::<Topocentric, EquatorialTrueOfDate, U>::from_spherical(&result_sph)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::*;

    #[test]
    fn roundtrip_equatorial_horizontal_is_stable() {
        let jd_ut1 = JulianDate::new(2_451_545.0);
        let jd_tt = JulianDate::new(2_451_545.000_800_741);

        let ra = 1.0 * RAD;
        let dec = 0.5 * RAD;
        let site = Geodetic::<ECEF>::new(
            0.0 * DEG,
            Degrees::new((0.6 * RAD).value().to_degrees()),
            0.0 * M,
        );

        let spherical_equ = affn::spherical::Direction::<EquatorialTrueOfDate>::new_raw(
            Degrees::new(dec.value().to_degrees()),
            Degrees::new(ra.value().to_degrees()),
        );
        let equatorial = spherical_equ.to_cartesian();
        let horizontal = equatorial.to_horizontal(&jd_ut1, &jd_tt, &site);
        let back = horizontal.to_equatorial(&jd_ut1, &jd_tt, &site);

        let spherical = equatorial.to_spherical();
        let back_spherical = back.to_spherical();

        let dra = (Radians::from(back_spherical.azimuth).value()
            - Radians::from(spherical.azimuth).value())
        .abs();
        let dra = dra.min(TAU - dra);
        let ddec = (Radians::from(back_spherical.polar).value()
            - Radians::from(spherical.polar).value())
        .abs();

        assert!(dra < 1e-12, "RA roundtrip error too large: {dra}");
        assert!(ddec < 1e-12, "Dec roundtrip error too large: {ddec}");
    }

    #[test]
    fn matches_erfa_reference_case() {
        // Reference generated from ERFA adapter in the lab benchmark suite
        let jd_ut1 = JulianDate::new(2_451_545.0);
        let jd_tt = JulianDate::new(2_451_545.000_800_741);

        let ra = 1.0 * RAD;
        let dec = 0.5 * RAD;
        let site = Geodetic::<ECEF>::new(
            0.0 * DEG,
            Degrees::new((0.6 * RAD).value().to_degrees()),
            0.0 * M,
        );

        let spherical_equ = affn::spherical::Direction::<EquatorialTrueOfDate>::new_raw(
            Degrees::new(dec.value().to_degrees()),
            Degrees::new(ra.value().to_degrees()),
        );
        let equatorial = spherical_equ.to_cartesian();
        let horizontal = equatorial.to_horizontal(&jd_ut1, &jd_tt, &site);

        let spherical = horizontal.to_spherical();
        let az = Radians::from(spherical.azimuth).value();
        let alt = Radians::from(spherical.polar).value();

        let az_ref = 0.670_382_099_942_421_8_f64;
        let alt_ref = -0.260_561_249_223_089_5_f64;

        assert!((az - az_ref).abs() < 1e-8, "Az mismatch: {az} vs {az_ref}");
        assert!(
            (alt - alt_ref).abs() < 1e-8,
            "Alt mismatch: {alt} vs {alt_ref}"
        );
    }

    #[test]
    fn topocentric_position_transformation() {
        use affn::spherical;
        use qtty::Kilometer;

        let jd_ut1 = JulianDate::new(2_451_545.0);
        let jd_tt = JulianDate::new(2_451_545.000_800_741);

        let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.5 * DEG, 100.0 * M);
        let ra = 45.0 * DEG;
        let dec = 30.0 * DEG;
        let distance = 1000.0 * KM;

        // Create equatorial position via spherical
        let equ_sph_dir = spherical::Direction::<EquatorialTrueOfDate>::new_raw(dec, ra);
        let equ_pos_sph =
            equ_sph_dir.position_with_params::<Topocentric, Kilometer>(site, distance);
        let equ_pos =
            affn::Position::<Topocentric, EquatorialTrueOfDate, Kilometer>::from_spherical(
                &equ_pos_sph,
            );

        // Transform to horizontal
        let horiz_pos = equ_pos.to_horizontal_position(&jd_ut1, &jd_tt);

        // Verify site preservation
        assert_eq!(horiz_pos.center_params(), &site);

        // Verify distance preservation (convert to comparable units)
        let horiz_dist = horiz_pos.distance();
        assert!(
            (horiz_dist.value() - distance.value()).abs() < 1e-10,
            "Distance mismatch"
        );

        // Transform back
        let back_pos = horiz_pos.to_equatorial_position(&jd_ut1, &jd_tt);

        // Verify roundtrip
        assert_eq!(back_pos.center_params(), &site);
        let back_dist = back_pos.distance();
        assert!(
            (back_dist.value() - distance.value()).abs() < 1e-10,
            "Roundtrip distance mismatch"
        );
    }
}
