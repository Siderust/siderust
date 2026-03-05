// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use super::bias::frame_bias_icrs_to_j2000;
use super::TransformFrame;
use crate::coordinates::{cartesian::Position, centers::ReferenceCenter, frames};
use affn::Rotation3;
use qtty::LengthUnit;

// Implement Transform trait for ICRS -> EclipticMeanJ2000
impl<C: ReferenceCenter, U> TransformFrame<Position<C, frames::EclipticMeanJ2000, U>>
    for Position<C, frames::ICRS, U>
where
    U: LengthUnit,
{
    fn to_frame(&self) -> Position<C, frames::EclipticMeanJ2000, U> {
        // J2000 mean obliquity ε₀ (IAU 2006): 84381.406″
        let eps = (84381.406_f64 / 3600.0).to_radians();
        let (sin_e, cos_e) = (eps.sin(), eps.cos());

        let rot: Rotation3 = frame_bias_icrs_to_j2000();
        let [bx, by, bz] = rot * [self.x(), self.y(), self.z()];
        Position::from_vec3(
            self.center_params().clone(),
            nalgebra::Vector3::new(bx, cos_e * by + sin_e * bz, -sin_e * by + cos_e * bz),
        )
    }
}

// Implement Transform trait for EquatorialMeanJ2000 -> EclipticMeanJ2000
impl<C: ReferenceCenter, U: LengthUnit> TransformFrame<Position<C, frames::EclipticMeanJ2000, U>>
    for Position<C, frames::EquatorialMeanJ2000, U>
{
    fn to_frame(&self) -> Position<C, frames::EclipticMeanJ2000, U> {
        // J2000 mean obliquity ε₀ (IAU 2006): 84381.406″
        let eps = (84381.406_f64 / 3600.0).to_radians();
        let (sin_e, cos_e) = (eps.sin(), eps.cos());

        let y = self.y();
        let z = self.z();
        Position::from_vec3(
            self.center_params().clone(),
            nalgebra::Vector3::new(self.x(), cos_e * y + sin_e * z, -sin_e * y + cos_e * z),
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::coordinates::cartesian::position::*;
    use crate::coordinates::transform::Transform;
    use crate::coordinates::*;
    use crate::macros::assert_cartesian_eq;
    use crate::macros::assert_spherical_eq;
    use crate::time::JulianDate;
    use qtty::Degrees;
    use qtty::*;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    fn serialize<U: LengthUnit>(ecl: &EclipticMeanJ2000<U>) -> EclipticMeanJ2000<U> {
        use crate::coordinates::transform::TransformFrame;
        // Use to_frame() (pure frame rotation, no center shift) so the
        // Heliocentric center is preserved across the round-trip.
        // Using transform() with the ICRS<U> default (Barycentric) would shift
        // the center and introduce the Sun's barycentric offset as residual.
        let hcrs: ICRS<U, centers::Heliocentric> = ecl.to_frame();
        hcrs.to_frame()
    }

    /// **Test 1: Identity transformation (Zero vector)**
    #[test]
    fn test_zero_vector_transformation() {
        let zero_ecl = EclipticMeanJ2000::<AstronomicalUnit>::CENTER;
        let zero_ecl_back = serialize(&zero_ecl);

        assert_cartesian_eq!(
            zero_ecl,
            zero_ecl_back,
            EPSILON,
            "Zero vector transformation should be reversible."
        );
    }

    /// **Test 3: Edge case - Aligned along X-axis (Should not change)**
    #[test]
    fn test_x_axis_aligned() {
        let coord_ecl = EclipticMeanJ2000::<AstronomicalUnit>::new(1.0, 0.0, 0.0);
        let coord_ecl_back = serialize(&coord_ecl);

        assert_cartesian_eq!(
            coord_ecl,
            coord_ecl_back,
            EPSILON,
            "X-aligned vector should remain unchanged after transformation."
        );
    }

    /// **Test 4: Edge case - Aligned along Y-axis**
    #[test]
    fn test_y_axis_aligned() {
        let coord_ecl = EclipticMeanJ2000::<AstronomicalUnit>::new(0.0, 1.0, 0.0);
        let coord_ecl_back = serialize(&coord_ecl);

        assert_cartesian_eq!(
            coord_ecl,
            coord_ecl_back,
            EPSILON,
            "Y-aligned vector should recover after round-trip transformation."
        );
    }

    /// **Test 5: Edge case - Aligned along Z-axis**
    #[test]
    fn test_z_axis_aligned() {
        let coord_ecl = EclipticMeanJ2000::<AstronomicalUnit>::new(0.0, 0.0, 1.0);
        let coord_ecl_back = serialize(&coord_ecl);

        assert_cartesian_eq!(
            coord_ecl,
            coord_ecl_back,
            EPSILON,
            "Z-aligned vector should recover after round-trip transformation."
        );
    }

    /// **Test 6: Transformation with extreme values**
    #[test]
    fn test_large_values() {
        let coord_ecl = EclipticMeanJ2000::<AstronomicalUnit>::new(1e10, -1e10, 5e9);
        let coord_ecl_back = serialize(&coord_ecl);

        assert_cartesian_eq!(
            coord_ecl,
            coord_ecl_back,
            5e-5,
            "Large values should not cause precision errors."
        );
    }

    /// **Test 7: Transformation with small values (Precision test)**
    #[test]
    fn test_small_values() {
        let coord_ecl = EclipticMeanJ2000::<AstronomicalUnit>::new(1e-10, -1e-10, 5e-11);
        let coord_ecl_back = serialize(&coord_ecl);

        assert_cartesian_eq!(
            coord_ecl,
            coord_ecl_back,
            EPSILON,
            "Small values should not cause precision errors."
        );
    }

    #[test]
    fn round_trip_equatorial_ecliptic() {
        let equatorial_orig = spherical::Position::<
            centers::Barycentric,
            frames::EquatorialMeanJ2000,
            AstronomicalUnit,
        >::new(Degrees::new(123.4), Degrees::new(-21.0), 2.7);
        let ecliptic: spherical::Position<centers::Barycentric, frames::EclipticMeanJ2000, Au> =
            equatorial_orig.transform(JulianDate::J2000);
        let equatorial_rec: spherical::Position<
            centers::Barycentric,
            frames::EquatorialMeanJ2000,
            Au,
        > = ecliptic.transform(JulianDate::J2000);

        assert_spherical_eq!(equatorial_orig, equatorial_rec, 1e-10);
    }
}
