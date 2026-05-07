// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame rotation implementations targeting [`EquatorialMeanJ2000`](crate::coordinates::frames::EquatorialMeanJ2000).
//!
//! ## Scientific scope
//!
//! The **J2000.0 mean equatorial frame** has its X-axis pointing toward the
//! mean vernal equinox at J2000.0 and its Z-axis aligned with Earth's mean
//! rotation pole at J2000.0. It is the traditional epoch-based equatorial
//! frame used for star catalogues and astrometric reductions.
//!
//! This module implements two rotation chains targeting this frame:
//!
//! - **EclipticMeanJ2000 → EquatorialMeanJ2000**: inverse obliquity rotation.
//! - **ICRS → EquatorialMeanJ2000**: IAU 2006 frame-bias matrix `rb`.
//!
//! ## Technical scope
//!
//! Both impls apply a fixed 3×3 rotation from [`bias`]:
//!
//! - `obliquity_ecl_to_eq()` = `Rx(−ε₀)`, tilting the ecliptic plane to the
//!   equatorial plane.
//! - `frame_bias_icrs_to_j2000()` = the `rb` matrix from `eraBp06(J2000)`.
//!
//! Spherical types are handled by the blanket impl in the parent module via
//! a Cartesian round-trip.
//!
//! ## References
//!
//! - IAU 2006 resolution B1 (frame-bias definition).
//! - SOFA routines `iauBp06`, `iauObl06`.

use super::bias;
use super::TransformFrame;
use crate::coordinates::{cartesian::Position, centers::ReferenceCenter, frames};
use crate::qtty::LengthUnit;

/// Rotate an ecliptic‐J2000 Cartesian vector into the mean equatorial‐J2000 frame.
impl<C: ReferenceCenter, U: LengthUnit> TransformFrame<Position<C, frames::EquatorialMeanJ2000, U>>
    for Position<C, frames::EclipticMeanJ2000, U>
{
    fn to_frame(&self) -> Position<C, frames::EquatorialMeanJ2000, U> {
        let rot = bias::obliquity_ecl_to_eq();
        let [x, y, z] = rot * [self.x(), self.y(), self.z()];
        Position::from_array(
            self.center_params().clone(),
            [x, y, z],
        )
    }
}

// Implement Transform trait for ICRS -> EquatorialMeanJ2000 (frame bias)
impl<C: ReferenceCenter, U: LengthUnit> TransformFrame<Position<C, frames::EquatorialMeanJ2000, U>>
    for Position<C, frames::ICRS, U>
{
    fn to_frame(&self) -> Position<C, frames::EquatorialMeanJ2000, U> {
        let rot = bias::frame_bias_icrs_to_j2000();
        let [x, y, z] = rot * [self.x(), self.y(), self.z()];
        Position::from_array(
            self.center_params().clone(),
            [x, y, z],
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::coordinates::transform::Transform;
    use crate::coordinates::{centers, frames, spherical::Position};
    use crate::macros::assert_spherical_eq;
    use crate::qtty::{AstronomicalUnit, Degrees};
    use crate::time::JulianDate;

    const EPS: f64 = 1.0e-12;

    #[test]
    fn round_trip_ecliptic_equatorial() {
        let ecliptic_orig = Position::<
            centers::Barycentric,
            frames::EclipticMeanJ2000,
            AstronomicalUnit,
        >::new(Degrees::new(123.4), Degrees::new(-21.0), 2.7);
        let equatorial: Position<
            centers::Barycentric,
            frames::EquatorialMeanJ2000,
            AstronomicalUnit,
        > = ecliptic_orig.transform(JulianDate::J2000);
        let ecliptic_rec: Position<
            centers::Barycentric,
            frames::EclipticMeanJ2000,
            AstronomicalUnit,
        > = equatorial.transform(JulianDate::J2000);

        assert_spherical_eq!(ecliptic_orig, ecliptic_rec, EPS);
    }
}
