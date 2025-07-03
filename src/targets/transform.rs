use super::Target;
use crate::coordinates::{
    cartesian::Vector,
    spherical::SphericalCoord,
    centers::*,
    frames::*,
    kinds::Kind,
    transform::Transform
};
use crate::units::Distance;

/// Blanket implementation to allow chaining two consecutive `Transform` operations.
///
/// This implementation allows converting a `Target` in Cartesian coordinates from one
/// reference center and frame (`C1`, `F1`) to another (`C2`, `F2`) by applying two 
/// transformations:
/// 1. Frame transformation (within the same center)
/// 2. Center transformation (within the new frame)
impl<C1, F1, C2, F2, U, K> From<&Target<Vector<C1, F1, U, K>>> for Target<Vector<C2, F2, U, K>>
where
    Vector<C1, F1, U, K>: Transform<Vector<C1, F2, U, K>>, // transform frame
    Vector<C1, F2, U, K>: Transform<Vector<C2, F2, U, K>>, // transform center
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    U: Distance,
    K: Kind,
{
    fn from(orig: &Target<Vector<C1, F1, U, K>>) -> Self {
        // Step 1: Transform to new frame, keeping the original center.
        let mid: Vector<C1, F2, U, K> = orig.position.transform(orig.time);
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
impl<C1, F1, C2, F2, U, K> From<&Target<SphericalCoord<C1, F1, U, K>>> for Target<SphericalCoord<C2, F2, U, K>>
where
    Vector<C1, F1, U, K>: Transform<Vector<C1, F2, U, K>>, // transform frame
    Vector<C1, F2, U, K>: Transform<Vector<C2, F2, U, K>>, // transform center
    Vector<C1, F1, U, K>: for<'a> From<&'a SphericalCoord<C1, F1, U, K>>, // to_cartesian
    SphericalCoord<C2, F2, U, K>: for<'a> From<&'a Vector<C2, F2, U, K>>, // to_spherical
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    U: Distance,
    K: Kind,
{
    fn from(orig: &Target<SphericalCoord<C1, F1, U, K>>) -> Self {
        // Step 1: Convert spherical to Cartesian
        let cart: Vector<C1, F1, U, K> = orig.position.to_cartesian();
        // Step 2: Transform to new frame
        let cart_mid: Vector<C1, F2, U, K> = cart.transform(orig.time);
        // Step 3: Transform to new center
        let cart_dest: Vector<C2, F2, U, K> = cart_mid.transform(orig.time);
        // Step 4: Convert back to spherical
        let mid: SphericalCoord<C2, F2, U, K> = cart_dest.to_spherical();
        // Construct the new Target
        Self::new_raw(
            mid,
            orig.time,
            orig.proper_motion.clone(),
        )
    }
}
