use crate::coordinates::{
    cartesian::Vector,
    centers::ReferenceCenter,
    kinds::Kind,
    frames
};
use crate::coordinates::transform::Transform;
use crate::units::LengthUnit;

/// Rotate an ecliptic‐J2000 Cartesian vector into the mean equatorial‐J2000 frame.
///
/// The transformation is a right‐hand rotation about +X by the obliquity ε.
impl<C: ReferenceCenter, K: Kind, U: LengthUnit> Transform<Vector<C, frames::Equatorial, U, K>> for Vector<C, frames::Ecliptic, U, K> {
    fn transform(&self, _jd: crate::units::JulianDay) -> Vector<C, frames::Equatorial, U, K> {
        let eps = 23.439281_f64.to_radians(); // obliquity in radians
        let (sin_e, cos_e) = (eps.sin(), eps.cos());

        let y = self.y();
        let z = self.z();
        Vector::new(
            self.x(),
            cos_e * y - sin_e * z,
            sin_e * y + cos_e * z
        )
    }
}

// Implement Transform trait for ICRS -> Equatorial (identity)
impl<C: ReferenceCenter, K: Kind, U: LengthUnit> Transform<Vector<C, frames::Equatorial, U, K>> for Vector<C, frames::ICRS, U, K> {
    fn transform(&self, _jd: crate::units::JulianDay) -> Vector<C, frames::Equatorial, U, K> {
        Vector::new(self.x(), self.y(), self.z())
    }
}


#[cfg(test)]
mod tests {
    use crate::coordinates::{
        spherical::Position,
        centers, frames
    };
    use crate::units::Degrees;
    use crate::macros::assert_spherical_eq;

    const EPS: f64 = 1.0e-12;

    #[test]
    fn round_trip_ecliptic_equatorial() {
        let ecliptic_orig = Position::<centers::Barycentric, frames::Ecliptic>::new(
            Degrees::new(123.4),
            Degrees::new(-21.0),
            2.7,
        );
        let equatorial  = Position::<centers::Barycentric, frames::Equatorial>::from(&ecliptic_orig);
        let ecliptic_rec = Position::<centers::Barycentric, frames::Ecliptic>::from(&equatorial);

        assert_spherical_eq!(ecliptic_orig, ecliptic_rec, EPS);
    }

}
