//! Frame transformations for Velocity types.
//!
//! Velocities are free vectors, so frame transformations are pure rotations
//! (same as directions, but with velocity units instead of dimensionless units).

use crate::coordinates::cartesian::Velocity;
use crate::coordinates::frames::{self, MutableFrame};
use crate::coordinates::transform::frames::bias::{frame_bias_icrs_to_j2000, frame_bias_j2000_to_icrs};
use crate::coordinates::transform::TransformFrame;
use affn::Rotation3;
use nalgebra::Vector3;
use qtty::{Quantity, Unit};

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

/// Frame transform from Ecliptic to EquatorialMeanJ2000 for velocities.
/// Rotation about +X by the obliquity ε ≈ 23.439281°.
impl<U: Unit> TransformFrame<Velocity<frames::EquatorialMeanJ2000, U>>
    for Velocity<frames::Ecliptic, U>
{
    fn to_frame(&self) -> Velocity<frames::EquatorialMeanJ2000, U> {
        let eps = 23.439281_f64.to_radians();
        let (sin_eps, cos_eps) = (eps.sin(), eps.cos());

        let x = self.x().value();
        let y = self.y().value();
        let z = self.z().value();

        // Same rotation as Direction: about +X axis
        let new_x = x;
        let new_y = y * cos_eps - z * sin_eps;
        let new_z = y * sin_eps + z * cos_eps;

        Velocity::<frames::EquatorialMeanJ2000, U>::from_vec3(Vector3::new(
            Quantity::new(new_x),
            Quantity::new(new_y),
            Quantity::new(new_z),
        ))
    }
}

/// Frame transform from EquatorialMeanJ2000 to Ecliptic for velocities.
/// Inverse rotation about +X by the obliquity ε.
impl<U: Unit> TransformFrame<Velocity<frames::Ecliptic, U>>
    for Velocity<frames::EquatorialMeanJ2000, U>
{
    fn to_frame(&self) -> Velocity<frames::Ecliptic, U> {
        let eps = 23.439281_f64.to_radians();
        let (sin_eps, cos_eps) = (eps.sin(), eps.cos());

        let x = self.x().value();
        let y = self.y().value();
        let z = self.z().value();

        // Inverse rotation (transpose)
        let new_x = x;
        let new_y = y * cos_eps + z * sin_eps;
        let new_z = -y * sin_eps + z * cos_eps;

        Velocity::<frames::Ecliptic, U>::from_vec3(Vector3::new(
            Quantity::new(new_x),
            Quantity::new(new_y),
            Quantity::new(new_z),
        ))
    }
}

/// Frame transform from ICRS to EquatorialMeanJ2000 for velocities (frame bias).
impl<U: Unit> TransformFrame<Velocity<frames::EquatorialMeanJ2000, U>>
    for Velocity<frames::ICRS, U>
{
    fn to_frame(&self) -> Velocity<frames::EquatorialMeanJ2000, U> {
        let rot: Rotation3 = frame_bias_icrs_to_j2000();
        let [x, y, z] = rot.apply_array([self.x().value(), self.y().value(), self.z().value()]);
        Velocity::<frames::EquatorialMeanJ2000, U>::from_vec3(Vector3::new(
            Quantity::new(x),
            Quantity::new(y),
            Quantity::new(z),
        ))
    }
}

/// Frame transform from EquatorialMeanJ2000 to ICRS for velocities (frame bias).
impl<U: Unit> TransformFrame<Velocity<frames::ICRS, U>>
    for Velocity<frames::EquatorialMeanJ2000, U>
{
    fn to_frame(&self) -> Velocity<frames::ICRS, U> {
        let rot: Rotation3 = frame_bias_j2000_to_icrs();
        let [x, y, z] = rot.apply_array([self.x().value(), self.y().value(), self.z().value()]);
        Velocity::<frames::ICRS, U>::from_vec3(Vector3::new(
            Quantity::new(x),
            Quantity::new(y),
            Quantity::new(z),
        ))
    }
}

/// Frame transform from Ecliptic to ICRS for velocities.
impl<U: Unit> TransformFrame<Velocity<frames::ICRS, U>> for Velocity<frames::Ecliptic, U> {
    fn to_frame(&self) -> Velocity<frames::ICRS, U> {
        // Ecliptic -> EquatorialMeanJ2000 -> ICRS
        let eq: Velocity<frames::EquatorialMeanJ2000, U> = self.to_frame();
        eq.to_frame()
    }
}

/// Frame transform from ICRS to Ecliptic for velocities.
impl<U: Unit> TransformFrame<Velocity<frames::Ecliptic, U>> for Velocity<frames::ICRS, U> {
    fn to_frame(&self) -> Velocity<frames::Ecliptic, U> {
        // ICRS -> EquatorialMeanJ2000 -> Ecliptic
        let eq: Velocity<frames::EquatorialMeanJ2000, U> = self.to_frame();
        eq.to_frame()
    }
}
