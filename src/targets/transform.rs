use super::Target;
use crate::coordinates::{
    cartesian::Vector,
    spherical::SphericalCoord,
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
impl<C1, F1, C2, F2, U> From<&Target<Vector<C1, F1, U>>> for Target<Vector<C2, F2, U>>
where
    Vector<C1, F1, U>: Transform<Vector<C1, F2, U>>, // transform frame
    Vector<C1, F2, U>: Transform<Vector<C2, F2, U>>, // transform center
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    U: LengthUnit,
{
    fn from(orig: &Target<Vector<C1, F1, U>>) -> Self {
        // Step 1: Transform to new frame, keeping the original center.
        let mid: Vector<C1, F2, U> = orig.position.transform(orig.time);
        // Step 2: Transform to new center, now using the new frame.
        Self::new_raw(
            mid.transform(orig.time),
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
impl<C1, F1, C2, F2, U> From<&Target<SphericalCoord<C1, F1, U>>> for Target<SphericalCoord<C2, F2, U>>
where
    Vector<C1, F1, U>: Transform<Vector<C1, F2, U>>, // transform frame
    Vector<C1, F2, U>: Transform<Vector<C2, F2, U>>, // transform center
    Vector<C1, F1, U>: for<'a> From<&'a SphericalCoord<C1, F1, U>>, // to_cartesian
    SphericalCoord<C2, F2, U>: for<'a> From<&'a Vector<C2, F2, U>>, // to_spherical
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    U: LengthUnit,
{
    fn from(orig: &Target<SphericalCoord<C1, F1, U>>) -> Self {
        // Step 1: Convert spherical to Cartesian
        let cart: Vector<C1, F1, U> = orig.position.to_cartesian();
        // Step 2: Transform to new frame
        let cart_mid: Vector<C1, F2, U> = cart.transform(orig.time);
        // Step 3: Transform to new center
        let cart_dest: Vector<C2, F2, U> = cart_mid.transform(orig.time);
        // Step 4: Convert back to spherical
        let mid: SphericalCoord<C2, F2, U> = cart_dest.to_spherical();
        // Construct the new Target
        Self::new_raw(
            mid,
            orig.time,
            orig.proper_motion.clone(),
        )
    }
}
