// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::astro::eop::EopProvider;
use crate::astro::earth_rotation_provider::itrs_to_equatorial_mean_j2000_rotation;
use crate::coordinates::cartesian::Position;
use crate::coordinates::centers::{Geodetic, Geocentric, Topocentric};
use crate::coordinates::frames::{EquatorialMeanJ2000, MutableFrame, ECEF};
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::transform::context::AstroContext;
use crate::time::JulianDate;
use qtty::{AstronomicalUnits, LengthUnit, Meter, Quantity};

// =============================================================================
// Extension Trait for Topocentric Transforms
// =============================================================================

/// Extension trait for transforming geocentric positions to topocentric coordinates.
///
/// Topocentric coordinates are measured from the observer's location on Earth's
/// surface. This is a real translation:
///
///   r_topo = r_geo - r_site
///
/// where r_site is the observer's geocentric position at the observation time.
///
/// For nearby objects (Moon, satellites, planets), the parallax is significant.
/// For distant stars, it's negligible but the math is still correct.
pub trait ToTopocentricExt<F: MutableFrame, U: LengthUnit> {
    /// Transform to topocentric coordinates with a specific observer site and context.
    ///
    /// Uses the full IAU 2006 terrestrial→celestial chain:
    /// polar motion `W(xp, yp, s')`, Earth rotation `ERA(UT1)`, and CIO/CIP
    /// orientation `Q(X, Y, s)` with optional EOP `dX/dY` corrections.
    fn to_topocentric_with_ctx<Eph, Eop: EopProvider, Nut>(
        &self,
        site: Geodetic<ECEF>,
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Position<Topocentric, F, U>;

    /// Transform to topocentric coordinates with a specific observer site.
    ///
    /// This applies a real parallax correction by translating the geocentric
    /// position by the observer's location. The observer's ITRF position is
    /// rotated to the celestial frame using the default EOP-backed context.
    ///
    /// # Arguments
    ///
    /// * `site` - The observer's geographic location
    /// * `jd` - The Julian Date of observation (for Earth rotation)
    ///
    /// # Returns
    ///
    /// The position as seen from the observer's location (topocentric).
    fn to_topocentric(&self, site: Geodetic<ECEF>, jd: JulianDate) -> Position<Topocentric, F, U> {
        let ctx: AstroContext = AstroContext::default();
        self.to_topocentric_with_ctx(site, jd, &ctx)
    }
}

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
// Geocentric → Topocentric (real parallax translation)
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> ToTopocentricExt<F, U> for Position<Geocentric, F, U>
where
    Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
    Position<Geocentric, EquatorialMeanJ2000, U>:
        crate::coordinates::transform::TransformFrame<Position<Geocentric, F, U>>,
{
    /// Transform to topocentric coordinates with a specific observer site.
    ///
    /// # Example
    ///
    /// ```rust
    /// use siderust::coordinates::cartesian::Position;
    /// use siderust::coordinates::centers::{Geocentric, Topocentric, Geodetic};
    /// use siderust::coordinates::frames::{ECEF, EquatorialMeanJ2000};
    /// use siderust::coordinates::transform::centers::ToTopocentricExt;
    /// use siderust::time::JulianDate;
    /// use qtty::*;
    ///
    /// // Moon at roughly 384,400 km from Earth's center
    /// let moon_geo = Position::<Geocentric, EquatorialMeanJ2000, Kilometer>::new(
    ///     384_400.0, 0.0, 0.0
    /// );
    ///
    /// // Observer in Greenwich
    /// let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
    ///
    /// // Get topocentric position (will differ by observer offset)
    /// let moon_topo = moon_geo.to_topocentric(site, JulianDate::J2000);
    /// ```
    fn to_topocentric_with_ctx<Eph, Eop: EopProvider, Nut>(
        &self,
        site: Geodetic<ECEF>,
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Position<Topocentric, F, U> {
        let site_equatorial = observer_site_equatorial_mean_j2000_with_ctx(site, jd, ctx);

        // Transform observer position to target frame
        let site_in_frame: Position<Geocentric, F, U> =
            crate::coordinates::transform::TransformFrame::to_frame(&site_equatorial);

        // Apply parallax: r_topo = r_geo - r_site
        let topo_vec = nalgebra::Vector3::new(
            self.x() - site_in_frame.x(),
            self.y() - site_in_frame.y(),
            self.z() - site_in_frame.z(),
        );

        Position::<Topocentric, F, U>::from_vec3(site, topo_vec)
    }
}

// =============================================================================
// Topocentric → Geocentric (inverse parallax translation)
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Geocentric, F, U>>
    for Position<Topocentric, F, U>
where
    Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
    Position<Geocentric, EquatorialMeanJ2000, U>:
        crate::coordinates::transform::TransformFrame<Position<Geocentric, F, U>>,
{
    /// Transform back to geocentric coordinates.
    ///
    /// This is the inverse of `to_topocentric`: it adds back the observer's
    /// geocentric position to recover the geocentric position.
    ///
    /// # Arguments
    ///
    /// * `jd` - The Julian Date of observation (for Earth rotation)
    fn to_center(&self, jd: JulianDate) -> Position<Geocentric, F, U> {
        // Get the observer site from the stored parameters
        let site = self.center_params();
        let ctx: AstroContext = AstroContext::default();
        let site_equatorial = observer_site_equatorial_mean_j2000_with_ctx(*site, jd, &ctx);

        // Transform observer position to target frame
        let site_in_frame: Position<Geocentric, F, U> =
            crate::coordinates::transform::TransformFrame::to_frame(&site_equatorial);

        // Inverse parallax: r_geo = r_topo + r_site
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
    use qtty::*;

    #[test]
    fn test_observer_site_geocentric_position() {
        // Greenwich at sea level
        let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
        let pos: Position<Geocentric, ECEF, Kilometer> = site.to_cartesian();

        // Greenwich is at ~3980 km from center in x, ~0 in y, ~4970 km in z
        // (roughly, for latitude ~51.5°)
        assert!(
            (pos.x() - Kilometers::new(3980.0)).abs() < Kilometers::new(100.0),
            "x={}",
            pos.x()
        );
        assert!(pos.y().abs() < 10.0, "y={}", pos.y()); // Near prime meridian
        assert!(
            (pos.z() - Kilometers::new(4970.0)).abs() < Kilometers::new(100.0),
            "z={}",
            pos.z()
        );

        // Distance should be roughly Earth's radius
        let r = pos.distance();
        assert!(
            (r - Kilometers::new(6371.0)).abs() < Kilometers::new(50.0),
            "distance={}",
            r
        );
    }

    #[test]
    fn test_topocentric_parallax_moon_like() {
        // Simulate Moon-like object at 384,400 km along x-axis
        let moon_geo =
            position::EquatorialMeanJ2000::<Kilometer, Geocentric>::new(384_400.0, 0.0, 0.0);

        // Observer at equator, prime meridian
        let site = Geodetic::<ECEF>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M);
        let jd = JulianDate::J2000;

        let moon_topo = moon_geo.to_topocentric(site, jd);

        // Topocentric and geocentric should differ by about Earth's radius (~6371 km)
        let diff = (moon_geo.x() - moon_topo.x()).abs();

        // The difference should be on the order of Earth's radius
        // (exact value depends on Earth-rotation angle at J2000)
        assert!(
            diff > 1000.0 && diff < 10000.0,
            "Parallax shift = {} km, expected ~6371 km",
            diff
        );
    }

    #[test]
    fn test_topocentric_roundtrip() {
        // Create a geocentric position
        let geo = position::EquatorialMeanJ2000::<Kilometer, Geocentric>::new(
            100_000.0, 50_000.0, 25_000.0,
        );

        let site = Geodetic::<ECEF>::new(10.0 * DEG, 45.0 * DEG, 100.0 * M);
        let jd = JulianDate::J2000;

        // Convert to topocentric and back
        let topo = geo.to_topocentric(site, jd);
        let geo_recovered: position::EquatorialMeanJ2000<Kilometer, Geocentric> =
            topo.to_center(jd);

        // Should recover original position
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
        // Very distant object (star-like, 100 pc = ~3e15 km)
        let star_geo = position::EquatorialMeanJ2000::<Au, Geocentric>::new(
            206265.0, 0.0, 0.0, // ~1 parsec
        );

        let site = Geodetic::<ECEF>::new(0.0 * DEG, 45.0 * DEG, 0.0 * M);
        let jd = JulianDate::J2000;

        let star_topo = star_geo.to_topocentric(site, jd);

        // Relative difference should be tiny (~Earth_radius / 1_pc)
        let rel_diff = ((star_geo.x() - star_topo.x()).abs() / star_geo.x()).simplify();

        // Earth radius / 1 pc ≈ 6371 km / (3.086e13 km) ≈ 2e-10
        assert!(rel_diff < 1e-6, "Relative parallax too large: {}", rel_diff);
    }
}
