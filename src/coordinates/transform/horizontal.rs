// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Geometric Earth horizon transformations.
//!
//! The context-aware projection maps the selected equinox-of-date axes to
//! Earth-fixed axes, then resolves the sightline into local north, east and
//! up. Its inverse reverses the same basis change and Earth orientation.
//! Azimuth increases from north through east. At zenith and nadir, azimuth
//! is assigned zero when the horizontal component is numerically unresolved.
//!
//! The explicit-UT1 traits also support the equinox-based sidereal
//! transformation with a supplied UT1 and TT, using zero polar motion
//! and zero observed celestial pole offsets.

use super::TransformContext;
use super::{providers::frame_rotation_with, surface_horizontal::angles_from_horizontal_cartesian};
use crate::astro::nutation::NutationModel;
use crate::astro::sidereal;
use crate::coordinates::cartesian::{Direction, Position};
use crate::coordinates::centers::{Geodetic, Topocentric};
use crate::coordinates::frames::{EquatorialTrueOfDate, Horizontal, ECEF};
use crate::coordinates::{spherical, surface::EnuBasis};
use crate::qtty::{Degrees, LengthUnit};
#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Transcendental};
use crate::time::JulianDate;
use affn::Rotation3;

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
        let rotation = equinox_to_terrestrial(*jd_ut1, *jd_tt);
        let [n, e, u] = project_cartesian(self.as_array(), site, rotation);
        Direction::new_unchecked(n, e, u)
    }
}

/// Project a true-of-date equatorial direction onto the local horizon.
///
/// The frame provider uses the selected TT nutation model and EOP for
/// UT1, terrestrial polar motion and observed celestial pole offsets.
/// Supply the sightline at the reception epoch in the same model's
/// true equator and equinox of date. This projection only changes axes.
pub fn project_direction_with<Ctx: TransformContext>(
    direction: &spherical::Direction<EquatorialTrueOfDate>,
    site: &Geodetic<ECEF>,
    jd: JulianDate,
    ctx: &Ctx,
) -> spherical::Direction<Horizontal> {
    let rotation = frame_rotation_with::<EquatorialTrueOfDate, ECEF, _>(jd, ctx);
    project(direction, site, rotation)
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
/// use siderust::ephemeris::Ephemeris;
/// use siderust::coordinates::transform::context::DefaultEphemeris;
/// use siderust::qtty::{DEG, M};
///
/// let jd = siderust::JulianDate::new(2_460_310.25);
/// let context: AstroContext = AstroContext::default();
/// let model = context.with_model::<Iau2000B>();
/// let site = Geodetic::<ECEF>::new(-17.9 * DEG, 28.8 * DEG, 2300.0 * M);
/// let moon = DefaultEphemeris::moon_geocentric(jd);
/// let target: Position<Geocentric, EquatorialMeanJ2000, _> = moon.to_frame_with(&jd, &model);
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

/// Recover a sightline in the selected true equator and equinox of date.
///
/// Supply the original observer, reception epoch and context to reverse
/// [`project_direction_with`]. To express the result in other equatorial
/// axes, follow this operation with a typed frame rotation.
///
/// ```rust
/// use siderust::coordinates::{
///     centers::Geodetic,
///     frames::{Horizontal, ECEF},
///     spherical::Direction,
///     transform::{horizontal::unproject_direction_with, AstroContext},
/// };
/// use siderust::qtty::{DEG, M};
///
/// let site = Geodetic::<ECEF>::new(-17.9 * DEG, 28.8 * DEG, 2300.0 * M);
/// let context: AstroContext = AstroContext::default();
/// let horizontal = Direction::<Horizontal>::new(30.0 * DEG, 120.0 * DEG);
/// let equatorial = unproject_direction_with(&horizontal, &site, siderust::J2000, &context);
/// assert!(equatorial.ra().value().is_finite());
/// ```
pub fn unproject_direction_with<Ctx: TransformContext>(
    direction: &spherical::Direction<Horizontal>,
    site: &Geodetic<ECEF>,
    jd: JulianDate,
    ctx: &Ctx,
) -> spherical::Direction<EquatorialTrueOfDate> {
    let rotation = frame_rotation_with::<EquatorialTrueOfDate, ECEF, _>(jd, ctx);
    let [x, y, z] = unproject_cartesian(direction.to_cartesian().as_array(), site, rotation);
    // atan2 keeps the small transverse component near either celestial pole.
    let dec = z.atan2((x * x + y * y).sqrt());
    let ra = Degrees::new(y.atan2(x).to_degrees()).normalize();
    spherical::Direction::new_unchecked(Degrees::new(dec.to_degrees()), ra)
}

/// Reverse [`project_position_with`], preserving the stored observer and range.
pub fn unproject_position_with<U: LengthUnit, Ctx: TransformContext>(
    position: &spherical::Position<Topocentric, Horizontal, U>,
    jd: JulianDate,
    ctx: &Ctx,
) -> spherical::Position<Topocentric, EquatorialTrueOfDate, U> {
    let site = position.center_params();
    unproject_direction_with(&position.direction(), site, jd, ctx)
        .position_with_params::<Topocentric, U>(*site, position.distance)
}

pub(crate) fn project(
    direction: &spherical::Direction<EquatorialTrueOfDate>,
    site: &Geodetic<ECEF>,
    rotation: Rotation3,
) -> spherical::Direction<Horizontal> {
    let [north, east, up] = project_cartesian(direction.to_cartesian().as_array(), site, rotation);
    let (alt, az) = angles_from_horizontal_cartesian(north, east, up);
    spherical::Direction::new_unchecked(alt, az)
}

pub(super) fn project_cartesian(
    vector: [f64; 3],
    site: &Geodetic<ECEF>,
    rotation: Rotation3,
) -> [f64; 3] {
    let fixed = rotation.apply_array(vector);
    EnuBasis::geodetic(site.lon, site.lat).to_horizontal(fixed)
}

pub(super) fn unproject_cartesian(
    vector: [f64; 3],
    site: &Geodetic<ECEF>,
    rotation: Rotation3,
) -> [f64; 3] {
    let fixed = EnuBasis::geodetic(site.lon, site.lat).to_fixed(vector);
    rotation.inverse().apply_array(fixed)
}

fn equinox_to_terrestrial(ut1: tempoch::JulianDate<tempoch::UT1>, tt: JulianDate) -> Rotation3 {
    let nut = super::context::DefaultNutationModel::nutation(tt);
    Rotation3::rz(-sidereal::gast_iau2006(
        ut1,
        tt,
        nut.dpsi,
        nut.mean_obliquity,
    ))
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
        let rotation = equinox_to_terrestrial(*jd_ut1, *jd_tt);
        let [x, y, z] = unproject_cartesian(self.as_array(), site, rotation);
        Direction::new_unchecked(x, y, z)
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
        let rotation = equinox_to_terrestrial(*jd_ut1, *jd_tt);
        let [n, e, u] =
            project_cartesian(self.as_array().map(|axis| axis.value()), &site, rotation);
        Position::new_with_params(site, n, e, u)
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
        let rotation = equinox_to_terrestrial(*jd_ut1, *jd_tt);
        let [x, y, z] =
            unproject_cartesian(self.as_array().map(|axis| axis.value()), &site, rotation);
        Position::new_with_params(site, x, y, z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::*;
    use core::f64::consts::TAU;

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
