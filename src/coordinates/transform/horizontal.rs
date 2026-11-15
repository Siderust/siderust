// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # IAU 2006/GAST Horizontal Coordinate Transformations
//!
//! This module provides traits for converting between equatorial and horizontal
//! coordinates using IAU 2006 precession-nutation and GAST-based sidereal time.
//!
//! ## Approach
//!
//! The transformation follows SOFA conventions:
//! - GAST (Greenwich Apparent Sidereal Time) = `gast_iau2006(UT1, TT, dpsi, mean_obliquity)`
//! - Hour Angle = GAST + observer_longitude - RA
//! - Azimuth uses astronomical convention: 0° = North, increasing clockwise through East
//!
//! ## Time Requirements
//!
//! Horizontal transformations require **both** UT1 and TT timescales:
//! - **UT1** for Earth rotation angle and sidereal time
//! - **TT** for precession/nutation models (IAU 2006A nutation)
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
//! use siderust::qtty::*;
//!
//! let jd_ut1 = tempoch::JulianDate::<tempoch::UT1>::new(2_451_545.0);
//! let jd_tt = siderust::JulianDate::new(2_451_545.000_800_741);
//!
//! let equatorial = spherical::Direction::<EquatorialTrueOfDate>::new(45.0 * DEG, 30.0 * DEG)
//!     .to_cartesian();
//! let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.5 * DEG, 0.0 * M);
//!
//! let horizontal: Direction<Horizontal> = equatorial.to_horizontal(&jd_ut1, &jd_tt, &site);
//! ```

use super::TransformContext;
use crate::astro::earth_rotation::jd_ut1_from_tt_eop;
use crate::astro::nutation::NutationModel;
use crate::astro::{nutation, sidereal};
use crate::coordinates::cartesian::{Direction, Position};
use crate::coordinates::centers::{Geodetic, Topocentric};
use crate::coordinates::frames::{EquatorialTrueOfDate, Horizontal, ECEF};
use crate::coordinates::spherical;
use crate::qtty::{Degrees, LengthUnit, Radians};
#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Scalar, Transcendental};
use crate::time::JulianDate;
use core::f64::consts::TAU;

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
        jd_ut1: &tempoch::JulianDate<tempoch::UT1>,
        jd_tt: &JulianDate,
        site: &Geodetic<ECEF>,
    ) -> Direction<Horizontal>;
}

impl ToHorizontal for Direction<EquatorialTrueOfDate> {
    #[inline]
    fn to_horizontal(
        &self,
        jd_ut1: &tempoch::JulianDate<tempoch::UT1>,
        jd_tt: &JulianDate,
        site: &Geodetic<ECEF>,
    ) -> Direction<Horizontal> {
        let nut = super::context::DefaultNutationModel::nutation(*jd_tt);
        let gast = sidereal::gast_iau2006(*jd_ut1, *jd_tt, nut.dpsi, nut.mean_obliquity);
        project(self.to_spherical(), site, gast).to_cartesian()
    }
}

/// Project a true-of-date equatorial direction onto the local horizon.
///
/// The context supplies Earth orientation for UT1 and the nutation model for
/// GAST. The direction must already describe the desired line of sight at the
/// TT epoch. This geometric projection applies no parallax, aberration or
/// refraction. It uses the context provider's infallible EOP lookup policy.
/// This equinox-based projection uses UT1 for GAST. For a body-fixed site
/// projection with polar motion, use
/// [`SurfaceHorizontalExt`](super::SurfaceHorizontalExt).
pub fn project_direction_with<Ctx: TransformContext>(
    direction: &spherical::Direction<EquatorialTrueOfDate>,
    site: &Geodetic<ECEF>,
    jd: JulianDate,
    ctx: &Ctx,
) -> spherical::Direction<Horizontal> {
    let eop = ctx.astro_context().eop_at_tt(jd);
    let ut1 = jd_ut1_from_tt_eop(jd, &eop);
    let nut = Ctx::Nut::nutation(jd);
    let gast = sidereal::gast_iau2006(ut1, jd, nut.dpsi, nut.mean_obliquity);
    project(*direction, site, gast)
}

/// Project a true-of-date topocentric position onto its observer's horizon.
///
/// The output preserves the site and distance. Use the same TT epoch and
/// context as for the preceding center and frame transformations. See
/// [`project_direction_with`] for the projection's scientific scope.
///
/// ```rust
/// use siderust::astro::nutation::Iau2000B;
/// use siderust::coordinates::{cartesian::Position, centers::{Geocentric, Geodetic, Topocentric},
///     frames::{EquatorialMeanJ2000, EquatorialTrueOfDate, ECEF},
///     transform::{AstroContext, PositionAstroExt, TransformCenter,
///         horizontal::project_position_with}};
/// use siderust::qtty::{Kilometer, DEG, M};
///
/// let jd = siderust::J2000;
/// let context: AstroContext = AstroContext::default();
/// let model = context.with_model::<Iau2000B>();
/// let site = Geodetic::<ECEF>::new(-17.9 * DEG, 28.8 * DEG, 2300.0 * M);
/// let target = Position::<Geocentric, EquatorialMeanJ2000, Kilometer>::new(384_400., 0., 0.);
/// let relative: Position<Topocentric, _, _> = target.to_center_with(site, jd, &model);
/// let equatorial = relative.to_frame_with::<EquatorialTrueOfDate, _>(&jd, &model);
/// let horizontal = project_position_with(&equatorial.to_spherical(), jd, &model);
/// assert_eq!(horizontal.center_params(), &site);
/// ```
///
/// A geocentric position must first be translated to the observer.
/// ```compile_fail
/// use siderust::coordinates::{centers::Geocentric, frames::EquatorialTrueOfDate,
///     spherical::Position, transform::{AstroContext, horizontal::project_position_with}};
/// use siderust::qtty::{Kilometer, DEG, KM};
///
/// let target: Position<Geocentric, EquatorialTrueOfDate, Kilometer> =
///     Position::new_unchecked(0.0 * DEG, 0.0 * DEG, 100.0 * KM);
/// let context: AstroContext = AstroContext::default();
/// project_position_with(&target, siderust::J2000, &context);
/// ```
pub fn project_position_with<U: LengthUnit, Ctx: TransformContext>(
    position: &spherical::Position<Topocentric, EquatorialTrueOfDate, U>,
    jd: JulianDate,
    ctx: &Ctx,
) -> spherical::Position<Topocentric, Horizontal, U> {
    let site = position.center_params();
    project_direction_with(&position.direction(), site, jd, ctx)
        .position_with_params::<Topocentric, U>(*site, position.distance)
}

pub(crate) fn project(
    direction: spherical::Direction<EquatorialTrueOfDate>,
    site: &Geodetic<ECEF>,
    gast: Radians,
) -> spherical::Direction<Horizontal> {
    let hour_angle = (gast + Radians::from(site.lon) - Radians::from(direction.azimuth)).value();
    let (sin_ha, cos_ha) = hour_angle.sin_cos();
    let (sin_dec, cos_dec) = Radians::from(direction.polar).sin_cos();
    let (sin_lat, cos_lat) = Radians::from(site.lat).sin_cos();

    let north = -cos_ha * cos_dec * sin_lat + sin_dec * cos_lat;
    let east = -sin_ha * cos_dec;
    let up = cos_ha * cos_dec * cos_lat + sin_dec * sin_lat;
    let planar = (north * north + east * east).sqrt();
    let az = if planar != 0.0 {
        east.atan2(north)
    } else {
        0.0
    }
    .rem_euclid(TAU);

    spherical::Direction::new_unchecked(
        Degrees::new(up.atan2(planar).to_degrees()),
        Degrees::new(az.to_degrees()),
    )
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
        jd_ut1: &tempoch::JulianDate<tempoch::UT1>,
        jd_tt: &JulianDate,
        site: &Geodetic<ECEF>,
    ) -> Direction<EquatorialTrueOfDate>;
}

impl FromHorizontal for Direction<Horizontal> {
    #[inline]
    fn to_equatorial(
        &self,
        jd_ut1: &tempoch::JulianDate<tempoch::UT1>,
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
        let nut =
            <super::context::DefaultNutationModel as nutation::NutationModel>::nutation(*jd_tt);
        let gast = sidereal::gast_iau2006(*jd_ut1, *jd_tt, nut.dpsi, nut.mean_obliquity);
        let last = (gast + obs_lon).value();

        // Spherical trigonometry for horizontal → equatorial
        let (sp, cp) = obs_lat.sin_cos();
        let (sa, ca) = alt.sin_cos();
        let (sz, cz) = az.sin_cos();

        let dec = (sp * sa + cp * ca * cz).clamp(-1.0, 1.0).asin();
        let ha = (-sz * ca).atan2(sa * cp - ca * cz * sp);
        let ra = (last - ha).rem_euclid(TAU);

        // Construct equatorial direction from spherical coordinates
        let spherical_equ = affn::spherical::Direction::<EquatorialTrueOfDate>::new_unchecked(
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
        jd_ut1: &tempoch::JulianDate<tempoch::UT1>,
        jd_tt: &JulianDate,
    ) -> Position<Topocentric, Horizontal, U>;

    /// Alias for [`to_horizontal_position`](Self::to_horizontal_position).
    #[inline]
    fn to_horizontal(
        &self,
        jd_ut1: &tempoch::JulianDate<tempoch::UT1>,
        jd_tt: &JulianDate,
    ) -> Position<Topocentric, Horizontal, U> {
        self.to_horizontal_position(jd_ut1, jd_tt)
    }
}

impl<U: LengthUnit> TopocentricEquatorialExt<U> for Position<Topocentric, EquatorialTrueOfDate, U> {
    #[inline]
    fn to_horizontal_position(
        &self,
        jd_ut1: &tempoch::JulianDate<tempoch::UT1>,
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
        jd_ut1: &tempoch::JulianDate<tempoch::UT1>,
        jd_tt: &JulianDate,
    ) -> Position<Topocentric, EquatorialTrueOfDate, U>;

    /// Alias for [`to_equatorial_position`](Self::to_equatorial_position).
    #[inline]
    fn to_equatorial(
        &self,
        jd_ut1: &tempoch::JulianDate<tempoch::UT1>,
        jd_tt: &JulianDate,
    ) -> Position<Topocentric, EquatorialTrueOfDate, U> {
        self.to_equatorial_position(jd_ut1, jd_tt)
    }
}

impl<U: LengthUnit> TopocentricHorizontalExt<U> for Position<Topocentric, Horizontal, U> {
    #[inline]
    fn to_equatorial_position(
        &self,
        jd_ut1: &tempoch::JulianDate<tempoch::UT1>,
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
    use crate::qtty::*;

    #[test]
    fn roundtrip_equatorial_horizontal_is_stable() {
        let jd_ut1 = tempoch::JulianDate::<tempoch::UT1>::new(2_451_545.0);
        let jd_tt = crate::time::JulianDate::new(2_451_545.000_800_741);

        let ra = 1.0 * RAD;
        let dec = 0.5 * RAD;
        let site = Geodetic::<ECEF>::new(
            0.0 * DEG,
            Degrees::new((0.6 * RAD).value().to_degrees()),
            0.0 * M,
        );

        let spherical_equ = affn::spherical::Direction::<EquatorialTrueOfDate>::new_unchecked(
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
    fn matches_sofa_reference_case() {
        // Reference generated from the SOFA-based lab benchmark suite
        let jd_ut1 = tempoch::JulianDate::<tempoch::UT1>::new(2_451_545.0);
        let jd_tt = crate::time::JulianDate::new(2_451_545.000_800_741);

        let ra = 1.0 * RAD;
        let dec = 0.5 * RAD;
        let site = Geodetic::<ECEF>::new(
            0.0 * DEG,
            Degrees::new((0.6 * RAD).value().to_degrees()),
            0.0 * M,
        );

        let spherical_equ = affn::spherical::Direction::<EquatorialTrueOfDate>::new_unchecked(
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
        use crate::qtty::Kilometer;
        use affn::spherical;

        let jd_ut1 = tempoch::JulianDate::<tempoch::UT1>::new(2_451_545.0);
        let jd_tt = crate::time::JulianDate::new(2_451_545.000_800_741);

        let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.5 * DEG, 100.0 * M);
        let ra = 45.0 * DEG;
        let dec = 30.0 * DEG;
        let distance = 1000.0 * KM;

        // Create equatorial position via spherical
        let equ_sph_dir = spherical::Direction::<EquatorialTrueOfDate>::new_unchecked(dec, ra);
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
