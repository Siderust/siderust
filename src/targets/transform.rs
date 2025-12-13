use super::Target;
use crate::coordinates::{cartesian, centers::*, frames::*, spherical, transform::Transform};
use qtty::Unit;

/// Blanket implementation to allow chaining two consecutive `Transform` operations.
///
/// This implementation allows converting a `Target` in Cartesian coordinates from one
/// reference center and frame (`C1`, `F1`) to another (`C2`, `F2`) by applying two
/// transformations:
/// 1. Frame transformation (within the same center)
/// 2. Center transformation (within the new frame)
impl<C1, F1, C2, F2, U> From<&Target<cartesian::Vector<C1, F1, U>>>
    for Target<cartesian::Vector<C2, F2, U>>
where
    cartesian::Vector<C1, F1, U>: Transform<cartesian::Vector<C1, F2, U>>, // transform frame
    cartesian::Vector<C1, F2, U>: Transform<cartesian::Vector<C2, F2, U>>, // transform center
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    U: Unit,
{
    fn from(orig: &Target<cartesian::Vector<C1, F1, U>>) -> Self {
        // Step 1: Transform to new frame, keeping the original center.
        // Step 2: Transform to new center, now using the new frame.
        Self::new_raw(
            orig.position.transform(orig.time).transform(orig.time),
            orig.time,
            orig.proper_motion.clone(),
        )
    }
}

/// Blanket implementation for transforming `Target` in spherical coordinates,
/// involving frame and center changes. Internally uses Cartesian conversions.
///
/// The transformation follows these steps:
/// 1. Convert spherical coordinates to Cartesian.
/// 2. Apply frame transformation.
/// 3. Apply center transformation.
/// 4. Convert back to spherical coordinates.
impl<C1, F1, C2, F2, U> From<&Target<spherical::SphericalCoord<C1, F1, U>>>
    for Target<spherical::SphericalCoord<C2, F2, U>>
where
    cartesian::Vector<C1, F1, U>: Transform<cartesian::Vector<C1, F2, U>>, // transform frame
    cartesian::Vector<C1, F2, U>: Transform<cartesian::Vector<C2, F2, U>>, // transform center
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    U: Unit,
{
    fn from(orig: &Target<spherical::SphericalCoord<C1, F1, U>>) -> Self {
        // Step 1: Convert spherical to Cartesian
        // Step 2: Transform to new frame
        // Step 3: Transform to new center
        // Step 4: Convert back to spherical
        Self::new_raw(
            orig.position
                .to_cartesian()
                .transform(orig.time)
                .transform(orig.time)
                .to_spherical(),
            orig.time,
            orig.proper_motion.clone(),
        )
    }
}
