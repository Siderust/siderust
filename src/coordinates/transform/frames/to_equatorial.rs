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

/*
#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::Degrees;
    use crate::coordinates::{CartesianCoord, SphericalCoord};
    use crate::coordinates::centers;
    use crate::coordinates::frames;
    use approx::assert_abs_diff_eq;

    // Un valor de tolerancia único para todos los tests
    const EPS: f64 = 1.0e-12;

    /// Comprueba que Spherical -> Cartesian -> Spherical es un casi-identidad
    #[test]
    fn round_trip_ecliptic_equatorial() {
        let sph0 = SphericalCoord::<centers::Barycentric, frames::Ecliptic>::new(
            Degrees::new(123.4),
            Degrees::new(-21.0),
            2.7,
        );
        let cart  = CartesianCoord::from(&sph0);
        let sph1  = SphericalCoord::from(&cart);

        // Radio
        assert_abs_diff_eq!(sph0.radial_distance, sph1.radial_distance, epsilon = EPS);
        // Inclinación / polar (φ)
        assert_abs_diff_eq!(sph0.polar.0, sph1.polar.0, epsilon = EPS);
        // Azimut (θ): normaliza diferencia en [0, 360)
        let d_theta = (sph0.azimuth.0 - sph1.azimuth.0).rem_euclid(360.0);
        assert!(d_theta.min(360.0 - d_theta) < 1.0e-10);
    }

    /// Comprueba que Cartesian -> Spherical -> Cartesian es un casi-identidad
    #[test]
    fn round_trip_cartesian_spherical() {
        let cart0 = CartesianCoord::<centers::Barycentric, frames::Ecliptic>::new(
            0.8, -1.2, 2.4,
        );
        let sph   = SphericalCoord::from(&cart0);
        let cart1 = CartesianCoord::from(&sph);

        assert_abs_diff_eq!(cart0.x(), cart1.x(), epsilon = EPS);
        assert_abs_diff_eq!(cart0.y(), cart1.y(), epsilon = EPS);
        assert_abs_diff_eq!(cart0.z(), cart1.z(), epsilon = EPS);
    }

    /// Comprueba que la rotación Eclíptica -> Ecuatorial es la esperada
    #[test]
    fn ecliptic_to_equatorial_rotation() {
        // Vector puro +Y eclíptico debe rotar a (0, cos ε, sin ε) en ecuatorial
        let v_ecl = CartesianCoord::<centers::Barycentric, frames::Ecliptic>::new(0.0, 1.0, 0.0);
        let v_eq  : CartesianCoord<_, frames::Equatorial> = v_ecl.transform(crate::units::JD_J2000);

        let eps = 23.439_281_f64.to_radians();
        assert_abs_diff_eq!(v_eq.x(), 0.0, epsilon = EPS);
        assert_abs_diff_eq!(v_eq.y(),  eps.cos(), epsilon = EPS);
        assert_abs_diff_eq!(v_eq.z(),  eps.sin(), epsilon = EPS);

        // Conserva la norma (rotación ortogonal)
        assert_abs_diff_eq!(v_eq.distance_from_origin(),
                            v_ecl.distance_from_origin(), epsilon = EPS);
    }
}
*/