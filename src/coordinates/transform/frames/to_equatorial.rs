// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use super::bias::frame_bias_icrs_to_j2000;
use super::TransformFrame;
use crate::coordinates::{cartesian::Position, centers::ReferenceCenter, frames};
use affn::Rotation3;
use qtty::LengthUnit;

/// Rotate an ecliptic‐J2000 Cartesian vector into the mean equatorial‐J2000 frame.
///
/// The transformation is a right‐hand rotation about +X by the obliquity ε.
impl<C: ReferenceCenter, U: LengthUnit> TransformFrame<Position<C, frames::EquatorialMeanJ2000, U>>
    for Position<C, frames::EclipticMeanJ2000, U>
{
    fn to_frame(&self) -> Position<C, frames::EquatorialMeanJ2000, U> {
        // J2000 mean obliquity ε₀ (IAU 2006): 84381.406″
        let eps = (84381.406_f64 / 3600.0).to_radians();
        let (sin_e, cos_e) = (eps.sin(), eps.cos());

        let y = self.y();
        let z = self.z();
        Position::from_vec3(
            self.center_params().clone(),
            nalgebra::Vector3::new(self.x(), cos_e * y - sin_e * z, sin_e * y + cos_e * z),
        )
    }
}

// Implement Transform trait for ICRS -> EquatorialMeanJ2000 (frame bias)
impl<C: ReferenceCenter, U: LengthUnit> TransformFrame<Position<C, frames::EquatorialMeanJ2000, U>>
    for Position<C, frames::ICRS, U>
{
    fn to_frame(&self) -> Position<C, frames::EquatorialMeanJ2000, U> {
        let rot: Rotation3 = frame_bias_icrs_to_j2000();
        let [x, y, z] = rot * [self.x(), self.y(), self.z()];
        Position::from_vec3(
            self.center_params().clone(),
            nalgebra::Vector3::new(x, y, z),
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::coordinates::transform::Transform;
    use crate::coordinates::{centers, frames, spherical::Position};
    use crate::macros::assert_spherical_eq;
    use crate::time::JulianDate;
    use qtty::{AstronomicalUnit, Degrees};

    const EPS: f64 = 1.0e-12;

    #[test]
    fn round_trip_ecliptic_equatorial() {
        let ecliptic_orig =
            Position::<centers::Barycentric, frames::EclipticMeanJ2000, AstronomicalUnit>::new(
                Degrees::new(123.4),
                Degrees::new(-21.0),
                2.7,
            );
        let equatorial: Position<
            centers::Barycentric,
            frames::EquatorialMeanJ2000,
            AstronomicalUnit,
        > = ecliptic_orig.transform(JulianDate::J2000);
        let ecliptic_rec: Position<centers::Barycentric, frames::EclipticMeanJ2000, AstronomicalUnit> =
            equatorial.transform(JulianDate::J2000);

        assert_spherical_eq!(ecliptic_orig, ecliptic_rec, EPS);
    }
}
