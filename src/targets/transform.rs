use super::Target;
use crate::coordinates::{
    cartesian,
    spherical,
    centers::*,
    frames::*,
    transform::Transform
};
use crate::units::LengthUnit;

/// Blanket implementation to allow chaining two consecutive `Transform` operations.
///
/// This implementation allows converting a `Target` in Cartesian coordinates from one
/// reference center and frame (`C1`, `F1`) to another (`C2`, `F2`) by applying two 
/// transformations:
/// 1. Frame transformation (within the same center)
/// 2. Center transformation (within the new frame)
impl<C1, F1, C2, F2, U> From<&Target<cartesian::Position<C1, F1, U>>> for Target<cartesian::Position<C2, F2, U>>
where
    cartesian::Position<C1, F1, U>: Transform<cartesian::Position<C1, F2, U>>, // transform frame
    cartesian::Position<C1, F2, U>: Transform<cartesian::Position<C2, F2, U>>, // transform center
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    U: LengthUnit,
{
    fn from(orig: &Target<cartesian::Position<C1, F1, U>>) -> Self {
        // Step 1: Transform to new frame, keeping the original center.
        // Step 2: Transform to new center, now using the new frame.
        Self::new_raw(
            orig.position.transform(orig.time)
                         .transform(orig.time),
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
impl<C1, F1, C2, F2, U> From<&Target<spherical::Position<C1, F1, U>>> for Target<spherical::Position<C2, F2, U>>
where
    cartesian::Position<C1, F1, U>: Transform<cartesian::Position<C1, F2, U>>, // transform frame
    cartesian::Position<C1, F2, U>: Transform<cartesian::Position<C2, F2, U>>, // transform center
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    U: LengthUnit,
{
    fn from(orig: &Target<spherical::Position<C1, F1, U>>) -> Self {
        // Step 1: Convert spherical to Cartesian
        // Step 2: Transform to new frame
        // Step 3: Transform to new center
        // Step 4: Convert back to spherical
        Self::new_raw(
            orig.position.to_cartesian()
                         .transform(orig.time)
                         .transform(orig.time)
                         .to_spherical(),
            orig.time,
            orig.proper_motion.clone(),
        )
    }
}



/// Blanket implementation to allow chaining two consecutive `Transform` operations.
///
/// This implementation allows converting a `Target` in Cartesian coordinates from one
/// reference center and frame (`C1`, `F1`) to another (`C2`, `F2`) by applying two 
/// transformations:
/// 1. Frame transformation (within the same center)
/// 2. Center transformation (within the new frame)
impl<C1, F1, C2, F2> From<&Target<cartesian::Direction<C1, F1>>> for Target<cartesian::Direction<C2, F2>>
where
    cartesian::Direction<C1, F1>: Transform<cartesian::Direction<C1, F2>>, // transform frame
    cartesian::Direction<C1, F2>: Transform<cartesian::Direction<C2, F2>>, // transform center
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
{
    fn from(orig: &Target<cartesian::Direction<C1, F1>>) -> Self {
        // Step 1: Transform to new frame, keeping the original center.
        // Step 2: Transform to new center, now using the new frame.
        Self::new_raw(
            orig.position.transform(orig.time)
                         .transform(orig.time),
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
impl<C1, F1, C2, F2> From<&Target<spherical::Direction<C1, F1>>> for Target<spherical::Direction<C2, F2>>
where
    cartesian::Direction<C1, F1>: Transform<cartesian::Direction<C1, F2>>, // transform frame
    cartesian::Direction<C1, F2>: Transform<cartesian::Direction<C2, F2>>, // transform center
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
{
    fn from(orig: &Target<spherical::Direction<C1, F1>>) -> Self {
        // Step 1: Convert spherical to Cartesian
        // Step 2: Transform to new frame
        // Step 3: Transform to new center
        // Step 4: Convert back to spherical
        Self::new_raw(
            orig.position.to_cartesian()
                         .transform(orig.time)
                         .transform(orig.time)
                         .to_spherical(),
            orig.time,
            orig.proper_motion.clone(),
        )
    }
}
