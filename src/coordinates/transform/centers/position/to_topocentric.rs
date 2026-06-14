// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Geocentric → Topocentric parallax correction.
//!
//! ## Scientific scope
//!
//! The **topocentric parallax** correction accounts for the difference between
//! a position measured from Earth's centre (geocentric) and the same position
//! as seen from a surface observer (topocentric). For the Moon the correction
//! reaches ~57 arcmin; for nearby asteroids it can be arcsec-level.
//!
//! The shift vector is the observer's geocentric position expressed in the
//! same Cartesian frame as the target, then subtracted from the target vector.
//!
//! ## Technical scope
//!
//! The observer's geodetic site (`Geodetic<ECEF>`) is first converted to an
//! ITRF Cartesian position, rotated into the target frame via the ITRS→Equatorial
//! rotation (GMST + polar motion), and then frame-transformed to match the
//! target coordinate frame before the subtraction.
//!
//! Two public entry points are provided:
//!
//! - [`to_topocentric_with_ctx`]: expert API accepting an explicit
//!   [`AstroContext`].
//! - [`to_topocentric_with`]: convenience wrapper accepting any
//!   [`TransformContext`].
//!
//! ## References
//!
//! - Urban, S. E. & Seidelmann, P. K. (2013). *Explanatory Supplement to the
//!   Astronomical Almanac*, 3rd ed. §7.2.
//! - IERS Conventions 2010, Chapter 4 (Earth rotation).

use crate::astro::earth_rotation_provider::itrs_to_equatorial_mean_j2000_rotation;
use crate::astro::eop::EopProvider;
use crate::astro::nutation::NutationModel;
use crate::coordinates::cartesian::Position;
use crate::coordinates::centers::{Geocentric, Geodetic, Topocentric};
use crate::coordinates::frames::{EquatorialMeanJ2000, MutableFrame, ECEF};
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::transform::context::{AstroContext, TransformContext};
use crate::ephemeris::Ephemeris;
use crate::qtty::{AstronomicalUnits, LengthUnit, Meter, Quantity};
use crate::time::JulianDate;

// =============================================================================
// Internal helper
// =============================================================================

#[inline]
fn observer_site_equatorial_mean_j2000_with_ctx<
    U: LengthUnit,
    Eph,
    Eop: EopProvider,
    Nut: NutationModel,
>(
    site: Geodetic<ECEF>,
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop>,
) -> Position<Geocentric, EquatorialMeanJ2000, U>
where
    Quantity<U>: From<Quantity<Meter>>,
{
    let site_itrf: Position<Geocentric, ECEF, U> = site.to_cartesian();
    let rot = itrs_to_equatorial_mean_j2000_rotation::<Eph, Eop, Nut>(jd, ctx);
    let [x_eq, y_eq, z_eq] = rot * [site_itrf.x(), site_itrf.y(), site_itrf.z()];
    Position::<Geocentric, EquatorialMeanJ2000, U>::new_with_params((), x_eq, y_eq, z_eq)
}

// =============================================================================
// Expert API: topocentric with custom ephemeris / EOP context
// =============================================================================

/// Transform a geocentric position to topocentric coordinates using a custom
/// [`AstroContext`].
///
/// This is the expert-level entry point when you need to override the Earth
/// Orientation Parameters (EOP), nutation model, or ephemeris backend.
///
/// For the default-precision path, call `pos.to_center(site, jd)` on any
/// `Position<Geocentric, F, U>`.
///
/// # Arguments
///
/// - `pos`: The geocentric position to transform.
/// - `site`: The observer's geodetic site coordinates in the ECEF frame.
/// - `jd`: The Julian Date (TT) of observation.
/// - `ctx`: The `AstroContext` carrying the EOP and ephemeris backends.
///
/// # Returns
///
/// The topocentric position in the same reference frame `F` and unit `U` as
/// the input, with the observer's geocentric position vector subtracted.
pub fn to_topocentric_with_ctx<F, U, Eph, Eop, Nut: NutationModel>(
    pos: &Position<Geocentric, F, U>,
    site: Geodetic<ECEF>,
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop>,
) -> Position<Topocentric, F, U>
where
    F: MutableFrame,
    U: LengthUnit,
    Eop: EopProvider,
    Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
    Position<Geocentric, EquatorialMeanJ2000, U>:
        crate::coordinates::transform::TransformFrame<Position<Geocentric, F, U>>,
{
    let site_equatorial =
        observer_site_equatorial_mean_j2000_with_ctx::<U, Eph, Eop, Nut>(site, jd, ctx);
    let site_in_frame: Position<Geocentric, F, U> =
        crate::coordinates::transform::TransformFrame::to_frame(&site_equatorial);
    let topo_vec = [
        pos.x() - site_in_frame.x(),
        pos.y() - site_in_frame.y(),
        pos.z() - site_in_frame.z(),
    ];
    Position::<Topocentric, F, U>::from_array(site, topo_vec)
}

/// Transform a geocentric position to topocentric coordinates using any
/// transform context wrapper.
///
/// # Arguments
///
/// - `pos`: The geocentric position to transform.
/// - `site`: The observer's geodetic site coordinates in the ECEF frame.
/// - `jd`: The Julian Date (TT) of observation.
/// - `ctx`: Any type implementing [`TransformContext`].
///
/// # Returns
///
/// The topocentric position in the same reference frame `F` and unit `U` as
/// the input.
pub fn to_topocentric_with<F, U, Ctx>(
    pos: &Position<Geocentric, F, U>,
    site: Geodetic<ECEF>,
    jd: JulianDate,
    ctx: &Ctx,
) -> Position<Topocentric, F, U>
where
    F: MutableFrame,
    U: LengthUnit,
    Ctx: TransformContext,
    Ctx::Eph: Ephemeris,
    Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
    Position<Geocentric, EquatorialMeanJ2000, U>:
        crate::coordinates::transform::TransformFrame<Position<Geocentric, F, U>>,
{
    to_topocentric_with_ctx::<F, U, Ctx::Eph, Ctx::Eop, Ctx::Nut>(
        pos,
        site,
        jd,
        ctx.astro_context(),
    )
}

// =============================================================================
// Geocentric → Topocentric  (TransformCenter impl)
// =============================================================================

// Accept custom model/eop context types for expert callers.
impl<F: MutableFrame, U: LengthUnit> TransformCenter<Topocentric, F, U>
    for Position<Geocentric, F, U>
where
    Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
    Position<Geocentric, EquatorialMeanJ2000, U>:
        crate::coordinates::transform::TransformFrame<Position<Geocentric, F, U>>,
{
    fn to_center_as<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        &self,
        params: Geodetic<ECEF>,
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop>,
    ) -> Position<Topocentric, F, U> {
        to_topocentric_with_ctx::<F, U, Eph, Eop, Nut>(self, params, jd, ctx)
    }
}

// =============================================================================
// Topocentric → Geocentric  (TransformCenter impl)
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Geocentric, F, U>
    for Position<Topocentric, F, U>
where
    Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
    Position<Geocentric, EquatorialMeanJ2000, U>:
        crate::coordinates::transform::TransformFrame<Position<Geocentric, F, U>>,
{
    /// Inverse parallax: recovers the geocentric position from topocentric.
    ///
    /// The observer site is read from the stored center params.
    fn to_center_as<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        &self,
        _params: (),
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop>,
    ) -> Position<Geocentric, F, U> {
        let site = self.center_params();
        let site_equatorial =
            observer_site_equatorial_mean_j2000_with_ctx::<U, Eph, Eop, Nut>(*site, jd, _ctx);
        let site_in_frame: Position<Geocentric, F, U> =
            crate::coordinates::transform::TransformFrame::to_frame(&site_equatorial);
        let geo_vec = [
            self.x() + site_in_frame.x(),
            self.y() + site_in_frame.y(),
            self.z() + site_in_frame.z(),
        ];
        Position::<Geocentric, F, U>::from_array_origin(geo_vec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::cartesian::position;
    use crate::coordinates::transform::TransformCenter;
    use crate::qtty::*;

    #[test]
    fn test_observer_site_geocentric_position() {
        let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
        let pos: Position<Geocentric, ECEF, Kilometer> = site.to_cartesian();

        assert!(
            (pos.x() - Kilometers::new(3980.0)).abs() < Kilometers::new(100.0),
            "x={}",
            pos.x()
        );
        assert!(pos.y().abs().value() < 10.0, "y={}", pos.y());
        assert!(
            (pos.z() - Kilometers::new(4970.0)).abs() < Kilometers::new(100.0),
            "z={}",
            pos.z()
        );

        let r = pos.distance();
        assert!(
            (r - Kilometers::new(6371.0)).abs() < Kilometers::new(50.0),
            "distance={}",
            r
        );
    }

    #[test]
    fn test_topocentric_parallax_moon_like() {
        let moon_geo =
            position::EquatorialMeanJ2000::<Kilometer, Geocentric>::new(384_400.0, 0.0, 0.0);

        let site = Geodetic::<ECEF>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M);
        let jd = crate::J2000;

        let moon_topo: position::EquatorialMeanJ2000<Kilometer, Topocentric> =
            moon_geo.to_center((site, jd));

        let diff = (moon_geo.x() - moon_topo.x()).abs();
        assert!(
            diff.value() > 1000.0 && diff.value() < 10000.0,
            "Parallax shift = {} km, expected ~6371 km",
            diff
        );
    }

    #[test]
    fn test_topocentric_roundtrip() {
        let geo = position::EquatorialMeanJ2000::<Kilometer, Geocentric>::new(
            100_000.0, 50_000.0, 25_000.0,
        );

        let site = Geodetic::<ECEF>::new(10.0 * DEG, 45.0 * DEG, 100.0 * M);
        let jd = crate::J2000;

        let topo: position::EquatorialMeanJ2000<Kilometer, Topocentric> = geo.to_center((site, jd));
        let geo_recovered: position::EquatorialMeanJ2000<Kilometer, Geocentric> =
            topo.to_center(jd);

        assert!(
            (geo.x() - geo_recovered.x()).abs().value() < 1e-6,
            "x: {} vs {}",
            geo.x(),
            geo_recovered.x()
        );
        assert!(
            (geo.y() - geo_recovered.y()).abs().value() < 1e-6,
            "y: {} vs {}",
            geo.y(),
            geo_recovered.y()
        );
        assert!(
            (geo.z() - geo_recovered.z()).abs().value() < 1e-6,
            "z: {} vs {}",
            geo.z(),
            geo_recovered.z()
        );
    }

    #[test]
    fn test_distant_object_small_parallax() {
        let star_geo = position::EquatorialMeanJ2000::<Au, Geocentric>::new(206265.0, 0.0, 0.0);

        let site = Geodetic::<ECEF>::new(0.0 * DEG, 45.0 * DEG, 0.0 * M);
        let jd = crate::J2000;

        let star_topo: position::EquatorialMeanJ2000<Au, Topocentric> =
            star_geo.to_center((site, jd));

        let rel_diff = (star_geo.x() - star_topo.x()).abs() / star_geo.x();
        assert!(rel_diff < 1e-6, "Relative parallax too large: {}", rel_diff);
    }
}
