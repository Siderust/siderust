//! Frame transformations for Direction types.
//!
//! Directions are free unit vectors, so frame transformations are pure rotations
//! without any translation or center dependence.

use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::{self, MutableFrame};
use crate::coordinates::transform::frames::bias::{frame_bias_icrs_to_j2000, frame_bias_j2000_to_icrs};
use crate::coordinates::transform::TransformFrame;
use affn::Rotation3;
use nalgebra::Vector3;

/// Identity frame transform for directions.
impl<F: MutableFrame> TransformFrame<Direction<F>> for Direction<F> {
    fn to_frame(&self) -> Direction<F> {
        Direction::<F>::from_vec3(self.as_vec3())
    }
}

/// Frame transform from Ecliptic to EquatorialMeanJ2000 for directions.
/// Rotation about +X by the J2000 mean obliquity ε₀ (IAU 2006): 84381.406″.
impl TransformFrame<Direction<frames::EquatorialMeanJ2000>> for Direction<frames::Ecliptic> {
    fn to_frame(&self) -> Direction<frames::EquatorialMeanJ2000> {
        let eps = (84381.406_f64 / 3600.0).to_radians();
        let (sin_eps, cos_eps) = (eps.sin(), eps.cos());

        let x = self.x();
        let y = self.y();
        let z = self.z();

        // Rotation matrix about +X axis:
        // [ 1    0        0     ] [ x ]
        // [ 0  cos(ε)  -sin(ε) ] [ y ]
        // [ 0  sin(ε)   cos(ε) ] [ z ]
        let new_x = x;
        let new_y = y * cos_eps - z * sin_eps;
        let new_z = y * sin_eps + z * cos_eps;

        Direction::<frames::EquatorialMeanJ2000>::from_vec3(Vector3::new(new_x, new_y, new_z))
    }
}

/// Frame transform from EquatorialMeanJ2000 to Ecliptic for directions.
/// Inverse rotation about +X by the obliquity ε.
impl TransformFrame<Direction<frames::Ecliptic>> for Direction<frames::EquatorialMeanJ2000> {
    fn to_frame(&self) -> Direction<frames::Ecliptic> {
        let eps = (84381.406_f64 / 3600.0).to_radians();
        let (sin_eps, cos_eps) = (eps.sin(), eps.cos());

        let x = self.x();
        let y = self.y();
        let z = self.z();

        // Inverse rotation matrix (transpose of the forward rotation):
        // [ 1    0        0     ] [ x ]
        // [ 0  cos(ε)   sin(ε) ] [ y ]
        // [ 0 -sin(ε)   cos(ε) ] [ z ]
        let new_x = x;
        let new_y = y * cos_eps + z * sin_eps;
        let new_z = -y * sin_eps + z * cos_eps;

        Direction::<frames::Ecliptic>::from_vec3(Vector3::new(new_x, new_y, new_z))
    }
}

/// Frame transform from ICRS to EquatorialMeanJ2000 for directions (frame bias).
impl TransformFrame<Direction<frames::EquatorialMeanJ2000>> for Direction<frames::ICRS> {
    fn to_frame(&self) -> Direction<frames::EquatorialMeanJ2000> {
        let rot: Rotation3 = frame_bias_icrs_to_j2000();
        let [x, y, z] = rot.apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::EquatorialMeanJ2000>::from_vec3(Vector3::new(x, y, z))
    }
}

/// Frame transform from EquatorialMeanJ2000 to ICRS for directions (frame bias).
impl TransformFrame<Direction<frames::ICRS>> for Direction<frames::EquatorialMeanJ2000> {
    fn to_frame(&self) -> Direction<frames::ICRS> {
        let rot: Rotation3 = frame_bias_j2000_to_icrs();
        let [x, y, z] = rot.apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::ICRS>::from_vec3(Vector3::new(x, y, z))
    }
}

/// Frame transform from Ecliptic to ICRS for directions.
/// Composed via EquatorialMeanJ2000 intermediate frame.
impl TransformFrame<Direction<frames::ICRS>> for Direction<frames::Ecliptic> {
    fn to_frame(&self) -> Direction<frames::ICRS> {
        // Ecliptic -> EquatorialMeanJ2000 -> ICRS
        let eq: Direction<frames::EquatorialMeanJ2000> = self.to_frame();
        eq.to_frame()
    }
}

/// Frame transform from ICRS to Ecliptic for directions.
/// Composed via EquatorialMeanJ2000 intermediate frame.
impl TransformFrame<Direction<frames::Ecliptic>> for Direction<frames::ICRS> {
    fn to_frame(&self) -> Direction<frames::Ecliptic> {
        // ICRS -> EquatorialMeanJ2000 -> Ecliptic
        let eq: Direction<frames::EquatorialMeanJ2000> = self.to_frame();
        eq.to_frame()
    }
}
