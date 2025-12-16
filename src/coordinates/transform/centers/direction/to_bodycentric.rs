//! # Bodycentric Direction Transformations
//!
//! This module provides transformations for directions to body-centered coordinate systems.
//!
//! ## Overview
//!
//! For distant objects (stars, galaxies), the direction from different centers is essentially
//! the same since the angular parallax is negligible. For nearby objects, the direction would
//! need to account for the body's position.
//!
//! This implementation provides methods to:
//! 1. Attach body-centric parameters to directions (for distant objects)
//! 2. Properly compute directions for nearby objects using position-based calculations

use crate::astro::JulianDate;
use crate::coordinates::cartesian::Vector;
use crate::coordinates::centers::{
    Barycentric, Bodycentric, BodycentricParams, Geocentric, Heliocentric,
};
use crate::coordinates::frames::MutableFrame;
use crate::coordinates::spherical::direction::DirectionUnit;
use crate::coordinates::transform::centers::TransformCenter;

// =============================================================================
// Geocentric → Bodycentric (Direction)
// =============================================================================

impl<F: MutableFrame> Vector<Geocentric, F, DirectionUnit> {
    /// Transform direction to body-centric coordinates.
    ///
    /// For distant objects (stars), the direction is essentially the same from any
    /// position within the solar system, so this just attaches the body parameters.
    ///
    /// For nearby objects, consider using position-based calculations and then
    /// extracting the direction.
    ///
    /// # Arguments
    ///
    /// - `body_params`: The orbital parameters of the body to use as the new center.
    /// - `_jd`: The Julian date (currently unused for distant objects).
    ///
    /// # Example
    ///
    /// ```rust
    /// use siderust::coordinates::centers::{BodycentricParams, Geocentric};
    /// use siderust::coordinates::cartesian::Direction;
    /// use siderust::coordinates::frames;
    /// use siderust::astro::orbit::Orbit;
    /// use siderust::astro::JulianDate;
    /// use qtty::*;
    ///
    /// let satellite_orbit = Orbit::new(
    ///     0.0000426 * AU,
    ///     0.0,
    ///     Degrees::new(51.6),
    ///     Degrees::new(0.0),
    ///     Degrees::new(0.0),
    ///     Degrees::new(0.0),
    ///     JulianDate::J2000,
    /// );
    ///
    /// // Create a cartesian direction (unit vector)
    /// let star_dir: Direction<Geocentric, frames::Equatorial> = Direction::new(0.5, 0.5, 0.707);
    /// let sat_params = BodycentricParams::geocentric(satellite_orbit);
    /// let star_from_sat = star_dir.to_bodycentric(sat_params, JulianDate::J2000);
    /// ```
    pub fn to_bodycentric(
        &self,
        body_params: BodycentricParams,
        _jd: JulianDate,
    ) -> Vector<Bodycentric, F, DirectionUnit> {
        // For distant objects, direction is the same
        Vector::<Bodycentric, F, DirectionUnit>::from_vec3(body_params, self.as_vec3())
    }
}

// =============================================================================
// Heliocentric → Bodycentric (Direction)
// =============================================================================

impl<F: MutableFrame> Vector<Heliocentric, F, DirectionUnit> {
    /// Transform direction to body-centric coordinates.
    ///
    /// For distant objects (stars), the direction is essentially the same from any
    /// position within the solar system.
    pub fn to_bodycentric(
        &self,
        body_params: BodycentricParams,
        _jd: JulianDate,
    ) -> Vector<Bodycentric, F, DirectionUnit> {
        Vector::<Bodycentric, F, DirectionUnit>::from_vec3(body_params, self.as_vec3())
    }
}

// =============================================================================
// Barycentric → Bodycentric (Direction)
// =============================================================================

impl<F: MutableFrame> Vector<Barycentric, F, DirectionUnit> {
    /// Transform direction to body-centric coordinates.
    ///
    /// For distant objects (stars), the direction is essentially the same from any
    /// position within the solar system.
    pub fn to_bodycentric(
        &self,
        body_params: BodycentricParams,
        _jd: JulianDate,
    ) -> Vector<Bodycentric, F, DirectionUnit> {
        Vector::<Bodycentric, F, DirectionUnit>::from_vec3(body_params, self.as_vec3())
    }
}

// =============================================================================
// Bodycentric → Other Centers (Direction)
// =============================================================================

impl<F: MutableFrame> TransformCenter<Vector<Geocentric, F, DirectionUnit>>
    for Vector<Bodycentric, F, DirectionUnit>
{
    /// Transform from body-centric to geocentric direction.
    ///
    /// For distant objects, the direction is essentially the same.
    fn to_center(&self, _jd: JulianDate) -> Vector<Geocentric, F, DirectionUnit> {
        Vector::<Geocentric, F, DirectionUnit>::from_vec3_origin(self.as_vec3())
    }
}

impl<F: MutableFrame> TransformCenter<Vector<Heliocentric, F, DirectionUnit>>
    for Vector<Bodycentric, F, DirectionUnit>
{
    /// Transform from body-centric to heliocentric direction.
    fn to_center(&self, _jd: JulianDate) -> Vector<Heliocentric, F, DirectionUnit> {
        Vector::<Heliocentric, F, DirectionUnit>::from_vec3_origin(self.as_vec3())
    }
}

impl<F: MutableFrame> TransformCenter<Vector<Barycentric, F, DirectionUnit>>
    for Vector<Bodycentric, F, DirectionUnit>
{
    /// Transform from body-centric to barycentric direction.
    fn to_center(&self, _jd: JulianDate) -> Vector<Barycentric, F, DirectionUnit> {
        Vector::<Barycentric, F, DirectionUnit>::from_vec3_origin(self.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::orbit::Orbit;
    use crate::coordinates::cartesian::Direction;
    use crate::coordinates::frames;
    use qtty::{Degrees, AU};

    #[test]
    fn test_direction_geocentric_to_bodycentric() {
        // Create a geocentric direction using cartesian type
        let star: Direction<Geocentric, frames::Equatorial> = Direction::new(0.5, 0.5, 0.707);

        let satellite_orbit = Orbit::new(
            0.0001 * AU,
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );

        let sat_params = BodycentricParams::geocentric(satellite_orbit);
        let star_from_sat = star.to_bodycentric(sat_params, JulianDate::J2000);

        // For distant stars, direction should be essentially the same
        assert!((star.x().value() - star_from_sat.x().value()).abs() < 1e-10);
        assert!((star.y().value() - star_from_sat.y().value()).abs() < 1e-10);
        assert!((star.z().value() - star_from_sat.z().value()).abs() < 1e-10);
    }

    #[test]
    fn test_direction_roundtrip() {
        let original: Direction<Geocentric, frames::Equatorial> = Direction::new(0.6, 0.0, 0.8);

        let orbit = Orbit::new(
            1.0 * AU,
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );

        let params = BodycentricParams::heliocentric(orbit);

        // Transform to bodycentric and back
        let bodycentric = original.to_bodycentric(params, JulianDate::J2000);
        let recovered: Direction<Geocentric, frames::Equatorial> =
            bodycentric.to_center(JulianDate::J2000);

        // Should recover original direction
        assert!((original.x().value() - recovered.x().value()).abs() < 1e-10);
        assert!((original.y().value() - recovered.y().value()).abs() < 1e-10);
        assert!((original.z().value() - recovered.z().value()).abs() < 1e-10);
    }
}
