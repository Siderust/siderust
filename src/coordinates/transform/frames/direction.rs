//! Frame transformations for Direction types.
//!
//! Directions are free unit vectors, so frame transformations are pure rotations
//! without any translation or center dependence.

use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::{self, MutableFrame};
use crate::coordinates::transform::TransformFrame;
use nalgebra::Vector3;

/// Identity frame transform for directions.
impl<F: MutableFrame> TransformFrame<Direction<F>> for Direction<F> {
    fn to_frame(&self) -> Direction<F> {
        Direction::<F>::from_vec3(self.as_vec3())
    }
}

/// Frame transform from Ecliptic to Equatorial for directions.
/// Rotation about +X by the obliquity ε ≈ 23.439281°.
impl TransformFrame<Direction<frames::Equatorial>> for Direction<frames::Ecliptic> {
    fn to_frame(&self) -> Direction<frames::Equatorial> {
        let eps = 23.439281_f64.to_radians();
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

        Direction::<frames::Equatorial>::from_vec3(Vector3::new(new_x, new_y, new_z))
    }
}

/// Frame transform from Equatorial to Ecliptic for directions.
/// Inverse rotation about +X by the obliquity ε.
impl TransformFrame<Direction<frames::Ecliptic>> for Direction<frames::Equatorial> {
    fn to_frame(&self) -> Direction<frames::Ecliptic> {
        let eps = 23.439281_f64.to_radians();
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

/// Frame transform from ICRS to Equatorial for directions (identity for now).
/// In a full implementation, this would include precession/nutation corrections.
impl TransformFrame<Direction<frames::Equatorial>> for Direction<frames::ICRS> {
    fn to_frame(&self) -> Direction<frames::Equatorial> {
        Direction::<frames::Equatorial>::from_vec3(self.as_vec3())
    }
}

/// Frame transform from Equatorial to ICRS for directions (identity for now).
impl TransformFrame<Direction<frames::ICRS>> for Direction<frames::Equatorial> {
    fn to_frame(&self) -> Direction<frames::ICRS> {
        Direction::<frames::ICRS>::from_vec3(self.as_vec3())
    }
}

/// Frame transform from Ecliptic to ICRS for directions.
/// Composed via Equatorial intermediate frame.
impl TransformFrame<Direction<frames::ICRS>> for Direction<frames::Ecliptic> {
    fn to_frame(&self) -> Direction<frames::ICRS> {
        // Ecliptic -> Equatorial -> ICRS
        let eq: Direction<frames::Equatorial> = self.to_frame();
        eq.to_frame()
    }
}

/// Frame transform from ICRS to Ecliptic for directions.
/// Composed via Equatorial intermediate frame.
impl TransformFrame<Direction<frames::Ecliptic>> for Direction<frames::ICRS> {
    fn to_frame(&self) -> Direction<frames::Ecliptic> {
        // ICRS -> Equatorial -> Ecliptic
        let eq: Direction<frames::Equatorial> = self.to_frame();
        eq.to_frame()
    }
}
