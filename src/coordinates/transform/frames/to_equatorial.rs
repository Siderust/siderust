use crate::coordinates::{
    cartesian::Vector,
    centers::ReferenceCenter,
    kinds::Kind,
    frames
};
use crate::coordinates::transform::Transform;
use crate::units::Unit;

/// Rotate an ecliptic‐J2000 Cartesian vector into the mean equatorial‐J2000 frame.
///
/// The transformation is a right‐hand rotation about +X by the obliquity ε.
impl<C: ReferenceCenter, K: Kind, U: Unit> Transform<Vector<C, frames::Equatorial, U, K>> for Vector<C, frames::Ecliptic, U, K> {
    fn transform(&self, _jd: crate::units::JulianDay) -> Vector<C, frames::Equatorial, U, K> {
        let x:f64 = self.x().into();
        let y:f64 = self.y().into();
        let z:f64 = self.z().into();
        let eps = 23.439281_f64.to_radians(); // obliquity in radians
        let (sin_e, cos_e) = (eps.sin(), eps.cos());

        let x_eq = x;
        let y_eq = cos_e * y - sin_e * z;
        let z_eq = sin_e * y + cos_e * z;

        Vector::new(x_eq.into(), y_eq.into(), z_eq.into())
    }
}

// Implement Transform trait for ICRS -> Equatorial (identity)
impl<C: ReferenceCenter, K: Kind, U: Unit> Transform<Vector<C, frames::Equatorial, U, K>> for Vector<C, frames::ICRS, U, K> {
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
