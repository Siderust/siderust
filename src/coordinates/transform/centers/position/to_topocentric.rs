use crate::astro::JulianDate;
use crate::coordinates::cartesian::Vector;
use crate::coordinates::centers::{Geocentric, ObserverSite, Topocentric};
use crate::coordinates::frames::MutableFrame;
use crate::coordinates::transform::centers::TransformCenter;
use qtty::LengthUnit;

// =============================================================================
// Geocentric → Topocentric (Position with distance)
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> Vector<Geocentric, F, U> {
    /// Transform to topocentric coordinates with a specific observer site.
    ///
    /// For positions, this attaches the observer site information. For distant objects,
    /// the parallax effect is negligible.
    ///
    /// # TODO
    /// Implement proper topocentric parallax for nearby objects (Moon, satellites, etc.)
    pub fn to_topocentric(&self, site: ObserverSite, _jd: JulianDate) -> Vector<Topocentric, F, U> {
        Vector::<Topocentric, F, U>::from_vec3(site, self.as_vec3())
    }
}

// =============================================================================
// Topocentric → Geocentric (Position with distance)
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Vector<Geocentric, F, U>>
    for Vector<Topocentric, F, U>
{
    fn to_center(&self, _jd: JulianDate) -> Vector<Geocentric, F, U> {
        // For positions, we should add back the observer's position
        // TODO: Implement proper inverse topocentric parallax
        Vector::<Geocentric, F, U>::from_vec3_origin(self.as_vec3())
    }
}
