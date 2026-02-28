// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::astro::earth_rotation_provider::itrs_to_equatorial_mean_j2000_rotation;
use crate::astro::eop::EopProvider;
use crate::coordinates::cartesian::Position;
use crate::coordinates::centers::{Geocentric, Geodetic, Topocentric};
use crate::coordinates::frames::{EquatorialMeanJ2000, MutableFrame, ECEF};
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::transform::context::AstroContext;
use crate::time::JulianDate;
use qtty::{AstronomicalUnits, LengthUnit, Meter, Quantity};

// =============================================================================
// Internal helper
// =============================================================================

#[inline]
fn observer_site_equatorial_mean_j2000_with_ctx<U: LengthUnit, Eph, Eop: EopProvider, Nut>(
    site: Geodetic<ECEF>,
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop, Nut>,
) -> Position<Geocentric, EquatorialMeanJ2000, U>
where
    Quantity<U>: From<Quantity<Meter>>,
{
    let site_itrf: Position<Geocentric, ECEF, U> = site.to_cartesian();
    let rot = itrs_to_equatorial_mean_j2000_rotation(jd, ctx);
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
pub fn to_topocentric_with_ctx<F, U, Eph, Eop, Nut>(
    pos: &Position<Geocentric, F, U>,
    site: Geodetic<ECEF>,
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop, Nut>,
) -> Position<Topocentric, F, U>
where
    F: MutableFrame,
    U: LengthUnit,
    Eop: EopProvider,
    Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
    Position<Geocentric, EquatorialMeanJ2000, U>:
        crate::coordinates::transform::TransformFrame<Position<Geocentric, F, U>>,
{
    let site_equatorial = observer_site_equatorial_mean_j2000_with_ctx(site, jd, ctx);
    let site_in_frame: Position<Geocentric, F, U> =
        crate::coordinates::transform::TransformFrame::to_frame(&site_equatorial);
    let topo_vec = nalgebra::Vector3::new(
        pos.x() - site_in_frame.x(),
        pos.y() - site_in_frame.y(),
        pos.z() - site_in_frame.z(),
    );
    Position::<Topocentric, F, U>::from_vec3(site, topo_vec)
}

// =============================================================================
// Geocentric → Topocentric  (TransformCenter impl)
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Topocentric, F, U>
    for Position<Geocentric, F, U>
where
    Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
    Position<Geocentric, EquatorialMeanJ2000, U>:
        crate::coordinates::transform::TransformFrame<Position<Geocentric, F, U>>,
{
    /// Transform to topocentric coordinates.
    ///
    /// `params` is the observer's geographic location ([`Geodetic<ECEF>`]).
    ///
    /// # Example
    ///
    /// ```rust
    /// use siderust::coordinates::cartesian::Position;
    /// use siderust::coordinates::centers::{Geocentric, Topocentric, Geodetic};
    /// use siderust::coordinates::frames::{ECEF, EquatorialMeanJ2000};
    /// use siderust::coordinates::transform::TransformCenter;
    /// use siderust::time::JulianDate;
    /// use qtty::*;
    ///
    /// let moon_geo = Position::<Geocentric, EquatorialMeanJ2000, Kilometer>::new(
    ///     384_400.0, 0.0, 0.0
    /// );
    /// let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
    /// let moon_topo = moon_geo.to_center((site, JulianDate::J2000));
    /// ```
    fn to_center_with(
        &self,
        params: Geodetic<ECEF>,
        jd: JulianDate,
        ctx: &AstroContext,
    ) -> Position<Topocentric, F, U> {
        to_topocentric_with_ctx(self, params, jd, ctx)
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
    fn to_center_with(
        &self,
        _params: (),
        jd: JulianDate,
        _ctx: &AstroContext,
    ) -> Position<Geocentric, F, U> {
        let site = self.center_params();
        let ctx: AstroContext = AstroContext::default();
        let site_equatorial = observer_site_equatorial_mean_j2000_with_ctx(*site, jd, &ctx);
        let site_in_frame: Position<Geocentric, F, U> =
            crate::coordinates::transform::TransformFrame::to_frame(&site_equatorial);
        let geo_vec = nalgebra::Vector3::new(
            self.x() + site_in_frame.x(),
            self.y() + site_in_frame.y(),
            self.z() + site_in_frame.z(),
        );
        Position::<Geocentric, F, U>::from_vec3_origin(geo_vec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::cartesian::position;
    use crate::coordinates::transform::TransformCenter;
    use qtty::*;

    #[test]
    fn test_observer_site_geocentric_position() {
        let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
        let pos: Position<Geocentric, ECEF, Kilometer> = site.to_cartesian();

        assert!(
            (pos.x() - Kilometers::new(3980.0)).abs() < Kilometers::new(100.0),
            "x={}",
            pos.x()
        );
        assert!(pos.y().abs() < 10.0, "y={}", pos.y());
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
        let jd = JulianDate::J2000;

        let moon_topo: position::EquatorialMeanJ2000<Kilometer, Topocentric> =
            moon_geo.to_center((site, jd));

        let diff = (moon_geo.x() - moon_topo.x()).abs();
        assert!(
            diff > 1000.0 && diff < 10000.0,
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
        let jd = JulianDate::J2000;

        let topo: position::EquatorialMeanJ2000<Kilometer, Topocentric> = geo.to_center((site, jd));
        let geo_recovered: position::EquatorialMeanJ2000<Kilometer, Geocentric> =
            topo.to_center(jd);

        assert!(
            (geo.x() - geo_recovered.x()).abs() < 1e-6,
            "x: {} vs {}",
            geo.x(),
            geo_recovered.x()
        );
        assert!(
            (geo.y() - geo_recovered.y()).abs() < 1e-6,
            "y: {} vs {}",
            geo.y(),
            geo_recovered.y()
        );
        assert!(
            (geo.z() - geo_recovered.z()).abs() < 1e-6,
            "z: {} vs {}",
            geo.z(),
            geo_recovered.z()
        );
    }

    #[test]
    fn test_distant_object_small_parallax() {
        let star_geo = position::EquatorialMeanJ2000::<Au, Geocentric>::new(206265.0, 0.0, 0.0);

        let site = Geodetic::<ECEF>::new(0.0 * DEG, 45.0 * DEG, 0.0 * M);
        let jd = JulianDate::J2000;

        let star_topo: position::EquatorialMeanJ2000<Au, Topocentric> =
            star_geo.to_center((site, jd));

        let rel_diff = ((star_geo.x() - star_topo.x()).abs() / star_geo.x()).simplify();
        assert!(rel_diff < 1e-6, "Relative parallax too large: {}", rel_diff);
    }
}
