use crate::coordinates::{
    cartesian::Vector,
    centers::ReferenceCenter,
    kinds::Kind,
    frames
};
use crate::coordinates::transform::Transform;
use crate::units::LengthUnit;

// Implement Transform trait for ICRS -> Ecliptic

impl<C: ReferenceCenter, K: Kind, U> Transform<Vector<C, frames::Ecliptic, U, K>>
    for Vector<C, frames::ICRS, U, K>
where
    U: LengthUnit,
{
    fn transform(&self, _jd: crate::astro::JulianDate) -> Vector<C, frames::Ecliptic, U, K> {
        let eps = 23.439281_f64.to_radians();
        let (sin_e, cos_e) = (eps.sin(), eps.cos());

        let y = self.y();
        let z = self.z();
        Vector::new(
            self.x(),
             cos_e * y + sin_e * z,
            -sin_e * y + cos_e * z,
        )
    }
}

// Implement Transform trait for Equatorial -> Ecliptic
impl<C: ReferenceCenter, K: Kind, U: LengthUnit> Transform<Vector<C, frames::Ecliptic, U, K>> for Vector<C, frames::Equatorial, U, K> {
    fn transform(&self, _jd: crate::astro::JulianDate) -> Vector<C, frames::Ecliptic, U, K> {
        let eps = 23.439281_f64.to_radians(); // obliquity in radians
        let (sin_e, cos_e) = (eps.sin(), eps.cos());

        let y = self.y();
        let z = self.z();
        Vector::new(
            self.x(),
             cos_e * y + sin_e * z,
            -sin_e * y + cos_e * z,
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::coordinates::*;
    use crate::coordinates::cartesian::position::*;
    use crate::units::Degrees;
    use crate::macros::assert_cartesian_eq;
    use crate::macros::assert_spherical_eq;
    use crate::units::{Quantity, LengthUnit, AstronomicalUnit};

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    fn serialize<U: LengthUnit>(ecl: &Ecliptic<U>) -> Ecliptic<U>
    where
        Quantity<U>: From<Quantity<AstronomicalUnit>>
    {
        let hcrs: ICRS<U> = ecl.into(); // Convert to ICRS
        let ecl_back: Ecliptic<U> = (&hcrs).into(); // Convert back to Ecliptic
        ecl_back
    }

    /// **Test 1: Identity transformation (Zero vector)**
    #[test]
    fn test_zero_vector_transformation() {
        let zero_ecl = Ecliptic::<AstronomicalUnit>::CENTER;
        let zero_ecl_back = serialize(&zero_ecl);

        assert_cartesian_eq!(zero_ecl, zero_ecl_back, EPSILON, "Zero vector transformation should be reversible.");
    }

    /// **Test 3: Edge case - Aligned along X-axis (Should not change)**
    #[test]
    fn test_x_axis_aligned() {
        let coord_ecl = Ecliptic::<AstronomicalUnit>::new(1.0, 0.0, 0.0);
        let coord_ecl_back = serialize(&coord_ecl);

        assert_cartesian_eq!(coord_ecl, coord_ecl_back, EPSILON, "X-aligned vector should remain unchanged after transformation.");
    }

    /// **Test 4: Edge case - Aligned along Y-axis**
    #[test]
    fn test_y_axis_aligned() {
        let coord_ecl = Ecliptic::<AstronomicalUnit>::new(0.0, 1.0, 0.0);
        let coord_ecl_back = serialize(&coord_ecl);

        assert_cartesian_eq!(coord_ecl, coord_ecl_back, EPSILON, "Y-aligned vector should recover after round-trip transformation.");
    }

    /// **Test 5: Edge case - Aligned along Z-axis**
    #[test]
    fn test_z_axis_aligned() {
        let coord_ecl = Ecliptic::<AstronomicalUnit>::new(0.0, 0.0, 1.0);
        let coord_ecl_back = serialize(&coord_ecl);

        assert_cartesian_eq!(coord_ecl, coord_ecl_back, EPSILON, "Z-aligned vector should recover after round-trip transformation.");
    }

    /// **Test 6: Transformation with extreme values**
    #[test]
    fn test_large_values() {
        let coord_ecl = Ecliptic::<AstronomicalUnit>::new(1e10, -1e10, 5e9);
        let coord_ecl_back = serialize(&coord_ecl);

        assert_cartesian_eq!(coord_ecl, coord_ecl_back, EPSILON, "Large values should not cause precision errors.");
    }

    /// **Test 7: Transformation with small values (Precision test)**
    #[test]
    fn test_small_values() {
        let coord_ecl = Ecliptic::<AstronomicalUnit>::new(1e-10, -1e-10, 5e-11);
        let coord_ecl_back = serialize(&coord_ecl);

        assert_cartesian_eq!(coord_ecl, coord_ecl_back, EPSILON, "Small values should not cause precision errors.");
    }

    #[test]
    fn round_trip_equatorial_ecliptic() {
        let equatorial_orig = spherical::Position::<centers::Barycentric, frames::Equatorial, AstronomicalUnit>::new(
            Degrees::new(123.4),
            Degrees::new(-21.0),
            2.7,
        );
        let ecliptic  = spherical::Position::<centers::Barycentric, frames::Ecliptic, AstronomicalUnit>::from(&equatorial_orig);
        let equatorial_rec = spherical::Position::<centers::Barycentric, frames::Equatorial, AstronomicalUnit>::from(&ecliptic);

        assert_spherical_eq!(equatorial_orig, equatorial_rec, 1e-10);
    }
}
