use crate::coordinates::{
    CartesianCoord,
    centers::ReferenceCenter,
    frames
};
use crate::coordinates::transform::Transform;

/// Rotate an ecliptic‐J2000 Cartesian vector into the mean equatorial‐J2000 frame.
///
/// The transformation is a right‐hand rotation about +X by the obliquity ε.
impl<C: ReferenceCenter> Transform<CartesianCoord<C, frames::Equatorial>> for CartesianCoord<C, frames::Ecliptic> {
    fn transform(&self, _jd: crate::units::JulianDay) -> CartesianCoord<C, frames::Equatorial> {
        let eps = 23.439281_f64.to_radians(); // obliquity in radians
        let (sin_e, cos_e) = (eps.sin(), eps.cos());

        let x_eq = self.x();
        let y_eq = cos_e * self.y() - sin_e * self.z();
        let z_eq = sin_e * self.y() + cos_e * self.z();

        CartesianCoord::new(x_eq, y_eq, z_eq)
    }
}

// Implement Transform trait for ICRS -> Equatorial (identity)
impl<C: ReferenceCenter> Transform<CartesianCoord<C, frames::Equatorial>> for CartesianCoord<C, frames::ICRS> {
    fn transform(&self, _jd: crate::units::JulianDay) -> CartesianCoord<C, frames::Equatorial> {
        CartesianCoord::new(self.x(), self.y(), self.z())
    }
}


#[cfg(test)]
mod tests {
    use crate::units::Degrees;
    use crate::coordinates::{
        SphericalCoord,
        centers, frames
    };
    use crate::macros::assert_spherical_eq;

    const EPS: f64 = 1.0e-12;

    #[test]
    fn round_trip_ecliptic_equatorial() {
        let ecliptic_orig = SphericalCoord::<centers::Barycentric, frames::Ecliptic>::new(
            Degrees::new(123.4),
            Degrees::new(-21.0),
            2.7,
        );
        let equatorial  = SphericalCoord::<centers::Barycentric, frames::Equatorial>::from(&ecliptic_orig);
        let ecliptic_rec = SphericalCoord::<centers::Barycentric, frames::Ecliptic>::from(&equatorial);

        assert_spherical_eq!(ecliptic_orig, ecliptic_rec, EPS);
    }

}
