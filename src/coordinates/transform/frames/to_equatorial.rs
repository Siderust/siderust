use super::TransformFrame;
use crate::coordinates::{cartesian::Position, centers::ReferenceCenter, frames};
use qtty::LengthUnit;

/// Rotate an ecliptic‐J2000 Cartesian vector into the mean equatorial‐J2000 frame.
///
/// The transformation is a right‐hand rotation about +X by the obliquity ε.
impl<C: ReferenceCenter, U: LengthUnit> TransformFrame<Position<C, frames::Equatorial, U>>
    for Position<C, frames::Ecliptic, U>
{
    fn to_frame(&self) -> Position<C, frames::Equatorial, U> {
        let eps = 23.439281_f64.to_radians(); // obliquity in radians
        let (sin_e, cos_e) = (eps.sin(), eps.cos());

        let y = self.y();
        let z = self.z();
        Position::from_vec3(
            self.center_params().clone(),
            nalgebra::Vector3::new(self.x(), cos_e * y - sin_e * z, sin_e * y + cos_e * z),
        )
    }
}

// Implement Transform trait for ICRS -> Equatorial (identity)
impl<C: ReferenceCenter, U: LengthUnit> TransformFrame<Position<C, frames::Equatorial, U>>
    for Position<C, frames::ICRS, U>
{
    fn to_frame(&self) -> Position<C, frames::Equatorial, U> {
        Position::from_vec3(self.center_params().clone(), *self.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use crate::astro::JulianDate;
    use crate::coordinates::transform::Transform;
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
        let equatorial: Position<centers::Barycentric, frames::Equatorial, AstronomicalUnit> =
            ecliptic_orig.transform(JulianDate::J2000);
        let ecliptic_rec: Position<centers::Barycentric, frames::Ecliptic, AstronomicalUnit> =
            equatorial.transform(JulianDate::J2000);

        assert_spherical_eq!(ecliptic_orig, ecliptic_rec, EPS);
    }
}
