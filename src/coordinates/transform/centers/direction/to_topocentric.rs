use crate::astro::JulianDate;
use crate::coordinates::cartesian::Vector;
use crate::coordinates::centers::{Geocentric, ObserverSite, Topocentric};
use crate::coordinates::frames::MutableFrame;
use crate::coordinates::spherical::direction::DirectionUnit;
use crate::coordinates::transform::centers::TransformCenter;

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

impl<F: MutableFrame> Vector<Geocentric, F, DirectionUnit> {
    /// Transform to topocentric coordinates with a specific observer site.
    ///
    /// For directions of distant objects, this attaches the observer site information
    /// without applying parallax corrections.
    pub fn to_topocentric(
        &self,
        site: ObserverSite,
        _jd: JulianDate,
    ) -> Vector<Topocentric, F, DirectionUnit> {
        Vector::<Topocentric, F, DirectionUnit>::from_vec3(site, self.as_vec3())
    }
}

// =============================================================================
// Topocentric → Geocentric (Direction only)
// =============================================================================

impl<F: MutableFrame> TransformCenter<Vector<Geocentric, F, DirectionUnit>>
    for Vector<Topocentric, F, DirectionUnit>
{
    fn to_center(&self, _jd: JulianDate) -> Vector<Geocentric, F, DirectionUnit> {
        // For directions of distant objects, topocentric ≈ geocentric
        Vector::<Geocentric, F, DirectionUnit>::from_vec3_origin(self.as_vec3())
    }
}
