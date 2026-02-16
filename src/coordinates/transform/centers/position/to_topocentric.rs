// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::astro::sidereal::gmst_iau2006;
use crate::coordinates::cartesian::Position;
use crate::coordinates::centers::{Geocentric, ObserverSite, Topocentric};
use crate::coordinates::frames::{EquatorialMeanJ2000, MutableFrame, ECEF};
use crate::coordinates::transform::centers::TransformCenter;
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
    /// Transform to topocentric coordinates with a specific observer site.
    ///
    /// This applies a real parallax correction by translating the geocentric
    /// position by the observer's location. The observer's ITRF position is
    /// rotated to the celestial frame using Greenwich Mean Sidereal Time.
    ///
    /// # Arguments
    ///
    /// * `site` - The observer's geographic location
    /// * `jd` - The Julian Date of observation (for Earth rotation)
    ///
    /// # Returns
    ///
    /// The position as seen from the observer's location (topocentric).
    fn to_topocentric(&self, site: ObserverSite, jd: JulianDate) -> Position<Topocentric, F, U>;
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
    /// use siderust::coordinates::centers::{Geocentric, Topocentric, ObserverSite};
    /// use siderust::coordinates::frames::EquatorialMeanJ2000;
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
    /// let site = ObserverSite::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
    ///
    /// // Get topocentric position (will differ by observer offset)
    /// let moon_topo = moon_geo.to_topocentric(site, JulianDate::J2000);
    /// ```
    fn to_topocentric(&self, site: ObserverSite, jd: JulianDate) -> Position<Topocentric, F, U> {
        // Get observer's ITRF position
        let site_itrf: Position<Geocentric, ECEF, U> = site.geocentric_itrf();

        // Rotate ITRF to celestial frame using GMST (IAU 2006 ERA-based)
        // For now, jd_ut1 ≈ jd_tt (NullEop equivalent); Phase 4 threads real EOP.
        let gmst_rad = gmst_iau2006(jd, jd).value();

        // Rotate from ECEF (x toward Greenwich, z toward pole) to equatorial
        // R_z(-GMST) transforms ECEF to equatorial
        let (sin_g, cos_g) = gmst_rad.sin_cos();

        let x_eq = site_itrf.x() * cos_g - site_itrf.y() * sin_g;
        let y_eq = site_itrf.x() * sin_g + site_itrf.y() * cos_g;
        let z_eq = site_itrf.z();

        let site_equatorial =
            Position::<Geocentric, EquatorialMeanJ2000, U>::new_with_params((), x_eq, y_eq, z_eq);

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

        // Get observer's ITRF position
        let site_itrf: Position<Geocentric, ECEF, U> = site.geocentric_itrf();

        // Rotate ITRF to celestial frame using GMST (IAU 2006 ERA-based)
        let gmst_rad = gmst_iau2006(jd, jd).value();

        let (sin_g, cos_g) = gmst_rad.sin_cos();

        let x_eq = site_itrf.x() * cos_g - site_itrf.y() * sin_g;
        let y_eq = site_itrf.x() * sin_g + site_itrf.y() * cos_g;
        let z_eq = site_itrf.z();

        let site_equatorial =
            Position::<Geocentric, EquatorialMeanJ2000, U>::new_with_params((), x_eq, y_eq, z_eq);

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
        let site = ObserverSite::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
        let pos: Position<Geocentric, ECEF, Kilometer> = site.geocentric_itrf();

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
        let site = ObserverSite::new(0.0 * DEG, 0.0 * DEG, 0.0 * M);
        let jd = JulianDate::J2000;

        let moon_topo = moon_geo.to_topocentric(site, jd);

        // Topocentric and geocentric should differ by about Earth's radius (~6371 km)
        let diff = (moon_geo.x() - moon_topo.x()).abs();

        // The difference should be on the order of Earth's radius
        // (exact value depends on GMST at J2000)
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

        let site = ObserverSite::new(10.0 * DEG, 45.0 * DEG, 100.0 * M);
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

        let site = ObserverSite::new(0.0 * DEG, 45.0 * DEG, 0.0 * M);
        let jd = JulianDate::J2000;

        let star_topo = star_geo.to_topocentric(site, jd);

        // Relative difference should be tiny (~Earth_radius / 1_pc)
        let rel_diff = ((star_geo.x() - star_topo.x()).abs() / star_geo.x()).simplify();

        // Earth radius / 1 pc ≈ 6371 km / (3.086e13 km) ≈ 2e-10
        assert!(rel_diff < 1e-6, "Relative parallax too large: {}", rel_diff);
    }
}
