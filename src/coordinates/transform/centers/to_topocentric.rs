//! Center transformation: Geocentric → Topocentric
//!
//! This module implements transformations between Geocentric and Topocentric
//! reference centers. The Topocentric center represents coordinates as seen
//! from a specific observer location on Earth's surface.

use crate::astro::JulianDate;
use crate::coordinates::cartesian::Vector;
use crate::coordinates::centers::{Geocentric, ObserverSite, Topocentric};
use crate::coordinates::frames::MutableFrame;
use crate::coordinates::spherical::direction::DirectionUnit;
use crate::coordinates::transform::centers::TransformCenter;
use qtty::LengthUnit;

// =============================================================================
// Geocentric → Topocentric (Direction only)
// =============================================================================

/// For directions (unit vectors), the transformation from geocentric to topocentric
/// is essentially an identity transformation since we're not accounting for parallax.
/// The main change is attaching the observer site information.
///
/// For precise work with nearby objects (Moon, artificial satellites), topocentric
/// parallax corrections should be applied. This implementation is suitable for
/// stars and distant objects.
impl<F: MutableFrame> TransformCenter<Vector<Topocentric, F, DirectionUnit>>
    for Vector<Geocentric, F, DirectionUnit>
{
    fn to_center(&self, _jd: JulianDate) -> Vector<Topocentric, F, DirectionUnit> {
        // For directions of distant objects, geocentric ≈ topocentric
        // The site must be provided externally or use default
        Vector::<Topocentric, F, DirectionUnit>::from_vec3(ObserverSite::default(), self.as_vec3())
    }
}

/// Transform with a specific observer site.
pub trait TransformToTopocentric<Coord> {
    fn to_topocentric(&self, site: ObserverSite, jd: JulianDate) -> Coord;
}

impl<F: MutableFrame> TransformToTopocentric<Vector<Topocentric, F, DirectionUnit>>
    for Vector<Geocentric, F, DirectionUnit>
{
    fn to_topocentric(
        &self,
        site: ObserverSite,
        _jd: JulianDate,
    ) -> Vector<Topocentric, F, DirectionUnit> {
        // For distant objects, just attach the site info
        Vector::<Topocentric, F, DirectionUnit>::from_vec3(site, self.as_vec3())
    }
}

// =============================================================================
// Geocentric → Topocentric (Position with distance)
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> TransformToTopocentric<Vector<Topocentric, F, U>>
    for Vector<Geocentric, F, U>
{
    fn to_topocentric(&self, site: ObserverSite, _jd: JulianDate) -> Vector<Topocentric, F, U> {
        // For positions, we should ideally subtract the observer's position vector
        // from the geocentric position. For distant objects, this is negligible.
        // TODO: Implement proper topocentric parallax for nearby objects (Moon, etc.)
        Vector::<Topocentric, F, U>::from_vec3(site, self.as_vec3())
    }
}

// =============================================================================
// Topocentric → Geocentric
// =============================================================================

impl<F: MutableFrame> TransformCenter<Vector<Geocentric, F, DirectionUnit>>
    for Vector<Topocentric, F, DirectionUnit>
{
    fn to_center(&self, _jd: JulianDate) -> Vector<Geocentric, F, DirectionUnit> {
        // For directions of distant objects, topocentric ≈ geocentric
        Vector::<Geocentric, F, DirectionUnit>::from_vec3_origin(self.as_vec3())
    }
}

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Vector<Geocentric, F, U>>
    for Vector<Topocentric, F, U>
{
    fn to_center(&self, _jd: JulianDate) -> Vector<Geocentric, F, U> {
        // For positions, we should add back the observer's position
        // TODO: Implement proper inverse topocentric parallax
        Vector::<Geocentric, F, U>::from_vec3_origin(self.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::frames::Equatorial;
    use qtty::{AstronomicalUnit, DEG, M};

    #[test]
    fn test_geocentric_to_topocentric_direction() {
        let geo_vec = Vector::<Geocentric, Equatorial, DirectionUnit>::new_with_params(
            (),
            1.0,
            0.0,
            0.0,
        );

        let site = ObserverSite::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
        let jd = JulianDate::J2000;

        let topo_vec: Vector<Topocentric, Equatorial, DirectionUnit> =
            geo_vec.to_topocentric(site.clone(), jd);

        assert_eq!(topo_vec.x().value(), geo_vec.x().value());
        assert_eq!(topo_vec.y().value(), geo_vec.y().value());
        assert_eq!(topo_vec.z().value(), geo_vec.z().value());
        assert_eq!(*topo_vec.center_params(), site);
    }

    #[test]
    fn test_topocentric_to_geocentric_direction() {
        let site = ObserverSite::new(-17.8925 * DEG, 28.7543 * DEG, 2396.0 * M);
        let topo_vec = Vector::<Topocentric, Equatorial, DirectionUnit>::new_with_params(
            site,
            0.5,
            0.5,
            0.707,
        );

        let jd = JulianDate::J2000;
        let geo_vec: Vector<Geocentric, Equatorial, DirectionUnit> = topo_vec.to_center(jd);

        assert_eq!(topo_vec.x().value(), geo_vec.x().value());
        assert_eq!(topo_vec.y().value(), geo_vec.y().value());
        assert_eq!(topo_vec.z().value(), geo_vec.z().value());
    }

    #[test]
    fn test_roundtrip_geocentric_topocentric() {
        let site = ObserverSite::new(10.0 * DEG, 45.0 * DEG, 100.0 * M);
        let jd = JulianDate::J2000;

        let original =
            Vector::<Geocentric, Equatorial, AstronomicalUnit>::new(1.0, 2.0, 3.0);

        let topo: Vector<Topocentric, Equatorial, AstronomicalUnit> =
            original.to_topocentric(site, jd);
        let back: Vector<Geocentric, Equatorial, AstronomicalUnit> = topo.to_center(jd);

        assert!((original.x().value() - back.x().value()).abs() < 1e-10);
        assert!((original.y().value() - back.y().value()).abs() < 1e-10);
        assert!((original.z().value() - back.z().value()).abs() < 1e-10);
    }
}
