// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame transformations for Velocity types.
//!
//! Velocities are free vectors, so frame transformations are pure rotations
//! (same as directions, but with velocity units instead of dimensionless units).

use crate::coordinates::cartesian::Velocity;
use crate::coordinates::frames::{self, MutableFrame};
use crate::coordinates::transform::frames::bias;
use crate::coordinates::transform::TransformFrame;
use nalgebra::Vector3;
use crate::qtty::Unit;

/// Identity frame transform for velocities.
impl<F, U> TransformFrame<Velocity<F, U>> for Velocity<F, U>
where
    F: MutableFrame,
    U: Unit,
{
    fn to_frame(&self) -> Velocity<F, U> {
        Self::from_vec3(*self.as_vec3())
    }
}

/// Frame transform from EclipticMeanJ2000 to EquatorialMeanJ2000 for velocities.
impl<U: Unit> TransformFrame<Velocity<frames::EquatorialMeanJ2000, U>>
    for Velocity<frames::EclipticMeanJ2000, U>
{
    fn to_frame(&self) -> Velocity<frames::EquatorialMeanJ2000, U> {
        let rot = bias::obliquity_ecl_to_eq();
        let [x, y, z] = rot * [self.x(), self.y(), self.z()];
        Velocity::<frames::EquatorialMeanJ2000, U>::from_vec3(Vector3::new(x, y, z))
    }
}

/// Frame transform from EquatorialMeanJ2000 to EclipticMeanJ2000 for velocities.
impl<U: Unit> TransformFrame<Velocity<frames::EclipticMeanJ2000, U>>
    for Velocity<frames::EquatorialMeanJ2000, U>
{
    fn to_frame(&self) -> Velocity<frames::EclipticMeanJ2000, U> {
        let rot = bias::obliquity_eq_to_ecl();
        let [x, y, z] = rot * [self.x(), self.y(), self.z()];
        Velocity::<frames::EclipticMeanJ2000, U>::from_vec3(Vector3::new(x, y, z))
    }
}

/// Frame transform from ICRS to EquatorialMeanJ2000 for velocities (frame bias).
impl<U: Unit> TransformFrame<Velocity<frames::EquatorialMeanJ2000, U>>
    for Velocity<frames::ICRS, U>
{
    fn to_frame(&self) -> Velocity<frames::EquatorialMeanJ2000, U> {
        let rot = bias::frame_bias_icrs_to_j2000();
        let [x, y, z] = rot * [self.x(), self.y(), self.z()];
        Velocity::<frames::EquatorialMeanJ2000, U>::from_vec3(Vector3::new(x, y, z))
    }
}

/// Frame transform from EquatorialMeanJ2000 to ICRS for velocities (frame bias).
impl<U: Unit> TransformFrame<Velocity<frames::ICRS, U>>
    for Velocity<frames::EquatorialMeanJ2000, U>
{
    fn to_frame(&self) -> Velocity<frames::ICRS, U> {
        let rot = bias::frame_bias_j2000_to_icrs();
        let [x, y, z] = rot * [self.x(), self.y(), self.z()];
        Velocity::<frames::ICRS, U>::from_vec3(Vector3::new(x, y, z))
    }
}

/// Frame transform from EclipticMeanJ2000 to ICRS for velocities.
impl<U: Unit> TransformFrame<Velocity<frames::ICRS, U>> for Velocity<frames::EclipticMeanJ2000, U> {
    fn to_frame(&self) -> Velocity<frames::ICRS, U> {
        // EclipticMeanJ2000 -> EquatorialMeanJ2000 -> ICRS
        let eq: Velocity<frames::EquatorialMeanJ2000, U> = self.to_frame();
        eq.to_frame()
    }
}

/// Frame transform from ICRS to EclipticMeanJ2000 for velocities.
impl<U: Unit> TransformFrame<Velocity<frames::EclipticMeanJ2000, U>> for Velocity<frames::ICRS, U> {
    fn to_frame(&self) -> Velocity<frames::EclipticMeanJ2000, U> {
        // ICRS -> EquatorialMeanJ2000 -> EclipticMeanJ2000
        let eq: Velocity<frames::EquatorialMeanJ2000, U> = self.to_frame();
        eq.to_frame()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::{AstronomicalUnit, Day, Per, Quantity};

    type AuPerDay = Per<AstronomicalUnit, Day>;

    fn vel_ecl(x: f64, y: f64, z: f64) -> Velocity<frames::EclipticMeanJ2000, AuPerDay> {
        Velocity::from_vec3(Vector3::new(
            Quantity::new(x),
            Quantity::new(y),
            Quantity::new(z),
        ))
    }

    fn assert_velocity_close<F, U>(a: &Velocity<F, U>, b: &Velocity<F, U>, eps: f64)
    where
        F: frames::MutableFrame,
        U: Unit,
    {
        let tol = Quantity::<U>::new(eps);
        assert!((a.x() - b.x()).abs() < tol);
        assert!((a.y() - b.y()).abs() < tol);
        assert!((a.z() - b.z()).abs() < tol);
    }

    #[test]
    fn ecliptic_to_equatorial_and_back_preserves_velocity() {
        let v_ecl = vel_ecl(0.1, -0.2, 0.3);
        let v_eq: Velocity<frames::EquatorialMeanJ2000, AuPerDay> = v_ecl.to_frame();
        let back: Velocity<frames::EclipticMeanJ2000, AuPerDay> = v_eq.to_frame();

        assert_velocity_close(&v_ecl, &back, 1e-12);
    }

    #[test]
    fn icrs_bias_roundtrip_for_velocity_is_identity() {
        let icrs = Velocity::<frames::ICRS, AuPerDay>::from_vec3(Vector3::new(
            Quantity::new(-1.0),
            Quantity::new(0.5),
            Quantity::new(0.25),
        ));

        let eq: Velocity<frames::EquatorialMeanJ2000, AuPerDay> = icrs.to_frame();
        let back: Velocity<frames::ICRS, AuPerDay> = eq.to_frame();

        assert_velocity_close(&icrs, &back, 1e-14);
    }

    #[test]
    fn ecliptic_to_icrs_composition_matches_direct_roundtrip() {
        let v_ecl = vel_ecl(0.01, 0.02, -0.05);
        let icrs: Velocity<frames::ICRS, AuPerDay> = v_ecl.to_frame();
        let back: Velocity<frames::EclipticMeanJ2000, AuPerDay> = icrs.to_frame();

        assert_velocity_close(&v_ecl, &back, 1e-12);
    }
}
