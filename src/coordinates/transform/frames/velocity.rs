//! Frame transformations for Velocity types.
//!
//! Velocities are free vectors, so frame transformations are pure rotations
//! (same as directions, but with velocity units instead of dimensionless units).

use crate::coordinates::cartesian::Velocity;
use crate::coordinates::frames::{self, MutableFrame};
use crate::coordinates::math::rotations;
use crate::coordinates::transform::TransformFrame;
use nalgebra::Vector3;
use qtty::{Quantity, Unit};

/// Identity frame transform for velocities.
impl<F, U> TransformFrame<Velocity<F, U>> for Velocity<F, U>
where
    F: MutableFrame,
    U: Unit,
{
    fn to_frame(&self) -> Velocity<F, U> {
        Self::from_vec3(self.as_vec3())
    }
}

/// Frame transform from Ecliptic to Equatorial for velocities.
/// Rotation about +X by the obliquity ε ≈ 23.439281°.
impl<U: Unit> TransformFrame<Velocity<frames::Equatorial, U>> for Velocity<frames::Ecliptic, U> {
    fn to_frame(&self) -> Velocity<frames::Equatorial, U> {
        let r = rotations::rotate_ecliptic_to_equatorial(
            self.x().value(),
            self.y().value(),
            self.z().value(),
        );
        Velocity::<frames::Equatorial, U>::from_vec3(Vector3::new(
            Quantity::new(r.x),
            Quantity::new(r.y),
            Quantity::new(r.z),
        ))
    }
}

/// Frame transform from Equatorial to Ecliptic for velocities.
/// Inverse rotation about +X by the obliquity ε.
impl<U: Unit> TransformFrame<Velocity<frames::Ecliptic, U>> for Velocity<frames::Equatorial, U> {
    fn to_frame(&self) -> Velocity<frames::Ecliptic, U> {
        let r = rotations::rotate_equatorial_to_ecliptic(
            self.x().value(),
            self.y().value(),
            self.z().value(),
        );
        Velocity::<frames::Ecliptic, U>::from_vec3(Vector3::new(
            Quantity::new(r.x),
            Quantity::new(r.y),
            Quantity::new(r.z),
        ))
    }
}

/// Frame transform from ICRS to Equatorial for velocities (identity for now).
impl<U: Unit> TransformFrame<Velocity<frames::Equatorial, U>> for Velocity<frames::ICRS, U> {
    fn to_frame(&self) -> Velocity<frames::Equatorial, U> {
        Velocity::<frames::Equatorial, U>::from_vec3(self.as_vec3())
    }
}

/// Frame transform from Equatorial to ICRS for velocities (identity for now).
impl<U: Unit> TransformFrame<Velocity<frames::ICRS, U>> for Velocity<frames::Equatorial, U> {
    fn to_frame(&self) -> Velocity<frames::ICRS, U> {
        Velocity::<frames::ICRS, U>::from_vec3(self.as_vec3())
    }
}

/// Frame transform from Ecliptic to ICRS for velocities.
impl<U: Unit> TransformFrame<Velocity<frames::ICRS, U>> for Velocity<frames::Ecliptic, U> {
    fn to_frame(&self) -> Velocity<frames::ICRS, U> {
        // Ecliptic -> Equatorial -> ICRS
        let eq: Velocity<frames::Equatorial, U> = self.to_frame();
        eq.to_frame()
    }
}

/// Frame transform from ICRS to Ecliptic for velocities.
impl<U: Unit> TransformFrame<Velocity<frames::Ecliptic, U>> for Velocity<frames::ICRS, U> {
    fn to_frame(&self) -> Velocity<frames::Ecliptic, U> {
        // ICRS -> Equatorial -> Ecliptic
        let eq: Velocity<frames::Equatorial, U> = self.to_frame();
        eq.to_frame()
    }
}
