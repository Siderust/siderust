use super::TransformFrame;
use crate::coordinates::math::rotations;
use crate::coordinates::{cartesian::Vector, centers::ReferenceCenter, frames};
use qtty::Unit;

/// Rotate an ecliptic‐J2000 Cartesian vector into the mean equatorial‐J2000 frame.
///
/// The transformation is a right‐hand rotation about +X by the obliquity ε.
impl<C: ReferenceCenter, U: Unit> TransformFrame<Vector<C, frames::Equatorial, U>>
    for Vector<C, frames::Ecliptic, U>
{
    fn to_frame(&self) -> Vector<C, frames::Equatorial, U> {
        let r = rotations::rotate_ecliptic_to_equatorial(
            self.x().value(),
            self.y().value(),
            self.z().value(),
        );
        Vector::from_vec3(
            self.center_params().clone(),
            nalgebra::Vector3::new(
                qtty::Quantity::new(r.x),
                qtty::Quantity::new(r.y),
                qtty::Quantity::new(r.z),
            ),
        )
    }
}

// Implement Transform trait for ICRS -> Equatorial (identity)
impl<C: ReferenceCenter, U: Unit> TransformFrame<Vector<C, frames::Equatorial, U>>
    for Vector<C, frames::ICRS, U>
{
    fn to_frame(&self) -> Vector<C, frames::Equatorial, U> {
        Vector::from_vec3(self.center_params().clone(), self.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use crate::coordinates::{centers, frames, spherical::Position};
    use crate::macros::assert_spherical_eq;
    use qtty::{AstronomicalUnit, Degrees};

    const EPS: f64 = 1.0e-12;

    #[test]
    fn round_trip_ecliptic_equatorial() {
        let ecliptic_orig =
            Position::<centers::Barycentric, frames::Ecliptic, AstronomicalUnit>::new(
                Degrees::new(123.4),
                Degrees::new(-21.0),
                2.7,
            );
        let equatorial =
            Position::<centers::Barycentric, frames::Equatorial, AstronomicalUnit>::from(
                &ecliptic_orig,
            );
        let ecliptic_rec =
            Position::<centers::Barycentric, frames::Ecliptic, AstronomicalUnit>::from(&equatorial);

        assert_spherical_eq!(ecliptic_orig, ecliptic_rec, EPS);
    }
}
