use crate::coordinates::{
    CartesianCoord,
    centers, frames
};
use crate::coordinates::transform::Transform;

// Implement Transform trait for ICRS -> Ecliptic
impl<C: centers::ReferenceCenter> Transform<CartesianCoord<C, frames::Ecliptic>> for CartesianCoord<C, frames::ICRS> {
    fn transform(&self, _jd: crate::units::JulianDay) -> CartesianCoord<C, frames::Ecliptic> {
        let eps = 23.439281_f64.to_radians(); // obliquity in radians
        let (sin_e, cos_e) = (eps.sin(), eps.cos());

        let x_ecl = self.x();
        let y_ecl = cos_e * self.y() + sin_e * self.z();
        let z_ecl = -sin_e * self.y() + cos_e * self.z();

        CartesianCoord::new(x_ecl, y_ecl, z_ecl)
    }
}

// Implement Transform trait for Equatorial -> Ecliptic
impl<C: centers::ReferenceCenter> Transform<CartesianCoord<C, frames::Ecliptic>> for CartesianCoord<C, frames::Equatorial> {
    fn transform(&self, _jd: crate::units::JulianDay) -> CartesianCoord<C, frames::Ecliptic> {
        let eps = 23.439281_f64.to_radians(); // obliquity in radians
        let (sin_e, cos_e) = (eps.sin(), eps.cos());

        let x_ecl = self.x();
        let y_ecl = cos_e * self.y() + sin_e * self.z();
        let z_ecl = -sin_e * self.y() + cos_e * self.z();

        CartesianCoord::new(x_ecl, y_ecl, z_ecl)
    }
}

#[cfg(test)]
mod tests {
    use crate::coordinates::*;
    use crate::coordinates::frames::*;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    /// Helper function to compare floating-point values within a small tolerance
    fn approx_eq(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    fn serialize(ecl: &CartesianCoord::<centers::Heliocentric, Ecliptic>) -> CartesianCoord::<centers::Heliocentric, Ecliptic> {
        let hcrs: HCRS = ecl.into(); // Convert to ICRS
        let ecl_back: CartesianCoord::<centers::Heliocentric, Ecliptic> = (&hcrs).into(); // Convert back to Ecliptic
        ecl_back
    }

    /// Check if two coordinates are approximately equal
    fn coords_approx_eq(a: &CartesianCoord<impl centers::ReferenceCenter, impl ReferenceFrame>, 
                        b: &CartesianCoord<impl centers::ReferenceCenter, impl ReferenceFrame>, 
                        epsilon: f64) -> bool {
        approx_eq(a.x(), b.x(), epsilon) &&
        approx_eq(a.y(), b.y(), epsilon) &&
        approx_eq(a.z(), b.z(), epsilon)
    }

    /// **Test 1: Identity transformation (Zero vector)**
    #[test]
    fn test_zero_vector_transformation() {
        let zero_ecl = CartesianCoord::<centers::Heliocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        let zero_ecl_back = serialize(&zero_ecl);

        assert!(coords_approx_eq(&zero_ecl, &zero_ecl_back, EPSILON), 
                "Zero vector transformation should be reversible.");
    }

    /// **Test 3: Edge case - Aligned along X-axis (Should not change)**
    #[test]
    fn test_x_axis_aligned() {
        let coord_ecl = CartesianCoord::<centers::Heliocentric, Ecliptic>::new(1.0, 0.0, 0.0);
        let coord_ecl_back = serialize(&coord_ecl);

        assert!(coords_approx_eq(&coord_ecl, &coord_ecl_back, EPSILON), 
                "X-aligned vector should remain unchanged after transformation.");
    }

    /// **Test 4: Edge case - Aligned along Y-axis**
    #[test]
    fn test_y_axis_aligned() {
        let coord_ecl = CartesianCoord::<centers::Heliocentric, Ecliptic>::new(0.0, 1.0, 0.0);
        let coord_ecl_back = serialize(&coord_ecl);

        assert!(coords_approx_eq(&coord_ecl, &coord_ecl_back, EPSILON), 
                "Y-aligned vector should recover after round-trip transformation.");
    }

    /// **Test 5: Edge case - Aligned along Z-axis**
    #[test]
    fn test_z_axis_aligned() {
        let coord_ecl = CartesianCoord::<centers::Heliocentric, Ecliptic>::new(0.0, 0.0, 1.0);
        let coord_ecl_back = serialize(&coord_ecl);

        assert!(coords_approx_eq(&coord_ecl, &coord_ecl_back, EPSILON), 
                "Z-aligned vector should recover after round-trip transformation.");
    }

    /// **Test 6: Transformation with extreme values**
    #[test]
    fn test_large_values() {
        let coord_ecl = CartesianCoord::<centers::Heliocentric, Ecliptic>::new(1e10, -1e10, 5e9);
        let coord_ecl_back = serialize(&coord_ecl);

        assert!(coords_approx_eq(&coord_ecl, &coord_ecl_back, EPSILON), 
                "Large values should not cause precision errors.");
    }

    /// **Test 7: Transformation with small values (Precision test)**
    #[test]
    fn test_small_values() {
        let coord_ecl = CartesianCoord::<centers::Heliocentric, Ecliptic>::new(1e-10, -1e-10, 5e-11);
        let coord_ecl_back = serialize(&coord_ecl);

        assert!(coords_approx_eq(&coord_ecl, &coord_ecl_back, EPSILON), 
                "Small values should not cause precision errors.");
    }
}
