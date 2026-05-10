// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Local orbital frames built from a Cartesian inertial state.
//!
//! Unlike the global astronomical frames in [`crate::coordinates::frames`]
//! (GCRS, ICRS, TEME, …), local orbital frames are *state-dependent*: their
//! orientation only makes sense relative to a particular spacecraft state.
//!
//! This module provides:
//!
//! - The marker frame types [`RTN`], [`VNC`], and [`LVLH`]; each is a
//!   zero-sized type implementing [`affn::frames::ReferenceFrame`].  The
//!   type tag prevents accidentally mixing displacements expressed in
//!   different local frames or with raw inertial vectors.
//! - The typed object [`LocalOrbitalFrame<M>`], which captures the rotation
//!   from the inertial parent frame ([`GCRS`]) into the local frame `M` for
//!   a specific [`OrbitState`].  Build it via
//!   [`LocalOrbitalFrame::<RTN>::from_state`] /
//!   [`LocalOrbitalFrame::<LVLH>::from_state`] /
//!   [`LocalOrbitalFrame::<VNC>::from_state`] and then rotate many vectors
//!   per state without rebuilding the basis.
//!
//! ## Conventions
//!
//! - **RTN** (also called RIC): R = along position vector (radial),
//!   N = orbit normal (`r×v`), T = N×R (transverse, along-track-ish).  The
//!   frame is right-handed and R is exactly along position.
//! - **VNC**: V = along velocity, N = orbit normal (`r×v`), C = V×N
//!   (co-normal).
//! - **LVLH**: Z = radial inward (`−R̂`), X = velocity projected onto the
//!   local horizontal, Y = Z×X.
//!
//! ## Usage
//!
//! ```rust
//! use siderust::astro::dynamics::{
//!     OrbitState, Position, Velocity,
//!     frames::{LocalOrbitalFrame, RTN},
//! };
//! use siderust::coordinates::frames::GCRS;
//! use siderust::qtty::unit::Kilometer;
//! use siderust::time::JulianDate;
//! use affn::cartesian::Displacement;
//!
//! let s = OrbitState::new(
//!     JulianDate::new(2_451_545.0),
//!     Position::<GCRS>::new(7000.0, 0.0, 0.0),
//!     Velocity::<GCRS>::new(0.0, 7.5, 0.0),
//! );
//!
//! let rtn = LocalOrbitalFrame::<RTN>::from_state(&s);
//! let d_gcrs = Displacement::<GCRS, Kilometer>::new(1.0, 0.0, 0.0);
//! let d_rtn = rtn.to_local(d_gcrs);
//! assert!((d_rtn.x().value() - 1.0).abs() < 1e-12);
//! ```

use core::marker::PhantomData;

use affn::cartesian::{Direction, Displacement};
use affn::DeriveReferenceFrame;
use affn::Rotation3;

use crate::coordinates::frames::GCRS;
use crate::qtty::unit::Kilometer;

use super::state::OrbitState;

// =============================================================================
// Local orbital frame marker types
// =============================================================================

/// Radial / Transverse / Normal local orbital frame.
///
/// - **R**: along the position vector (radial, geocentre → satellite)
/// - **T**: transverse, = N×R (along-track-ish; *not* the velocity direction)
/// - **N**: orbit normal, = normalize(r×v)
///
/// The frame is right-handed.  Use
/// [`LocalOrbitalFrame::<RTN>::from_state`] to materialise it from an
/// inertial state.
#[derive(Debug, Copy, Clone, DeriveReferenceFrame)]
pub struct RTN;

/// Velocity / Normal / Co-normal local orbital frame.
///
/// - **V**: along the velocity vector
/// - **N**: orbit normal, = normalize(r×v)
/// - **C**: co-normal, = V×N
#[derive(Debug, Copy, Clone, DeriveReferenceFrame)]
pub struct VNC;

/// Local-Vertical / Local-Horizontal frame.
///
/// - **Z**: radial inward (= −R̂)
/// - **X**: velocity projected onto the local horizontal
/// - **Y**: = Z×X
#[derive(Debug, Copy, Clone, DeriveReferenceFrame)]
pub struct LVLH;

// =============================================================================
// Typed local orbital frame object
// =============================================================================

/// Materialised local orbital frame for a specific [`OrbitState`].
///
/// `M` is the local frame marker (one of [`RTN`], [`VNC`], [`LVLH`]).
/// The captured [`Rotation3`] takes vectors from the inertial parent frame
/// ([`GCRS`]) into the local frame.
///
/// Construct with [`Self::from_state`] (specialised on `M`).  Use
/// [`Self::to_local`] / [`Self::to_inertial`] to rotate displacements while
/// keeping the frame tag in the type system.
#[derive(Debug, Clone, Copy)]
pub struct LocalOrbitalFrame<M> {
    rotation: Rotation3,
    _marker: PhantomData<M>,
}

impl<M> LocalOrbitalFrame<M> {
    /// Wrap an existing GCRS→`M` rotation matrix.
    ///
    /// Prefer the specialised `from_state` constructors below; this is for
    /// callers that already have a basis matrix from another source.
    #[inline]
    pub fn from_rotation(rotation: Rotation3) -> Self {
        Self { rotation, _marker: PhantomData }
    }

    /// The captured GCRS→`M` rotation matrix.
    #[inline]
    pub fn rotation(&self) -> Rotation3 {
        self.rotation
    }

    /// The transposed (`M`→GCRS) rotation matrix.
    #[inline]
    pub fn rotation_inverse(&self) -> Rotation3 {
        self.rotation.transpose()
    }
}

/// Velocity unit direction extracted from an orbit state.
///
/// The velocity quantity has units `Per<Km, Second>`, which is not a
/// `LengthUnit`, so `normalize()` is not available on the vector directly.
/// We extract the raw scalar components and construct a `Direction`.
#[inline]
fn velocity_direction(s: &OrbitState) -> Direction<GCRS> {
    Direction::new(
        s.velocity.x().value(),
        s.velocity.y().value(),
        s.velocity.z().value(),
    )
}

impl LocalOrbitalFrame<RTN> {
    /// Build the GCRS→RTN rotation for the given state.
    pub fn from_state(s: &OrbitState) -> Self {
        let r_hat: Direction<GCRS> = s.position.direction_unchecked();
        let v_hat: Direction<GCRS> = velocity_direction(s);
        let n_hat = r_hat
            .cross(&v_hat)
            .expect("position and velocity must not be parallel");
        let t_hat = n_hat
            .cross(&r_hat)
            .expect("n and r are orthogonal by construction");
        let rotation = Rotation3::from_matrix_unchecked([
            r_hat.as_array(),
            t_hat.as_array(),
            n_hat.as_array(),
        ]);
        Self::from_rotation(rotation)
    }
}

impl LocalOrbitalFrame<VNC> {
    /// Build the GCRS→VNC rotation for the given state.
    pub fn from_state(s: &OrbitState) -> Self {
        let v_hat: Direction<GCRS> = velocity_direction(s);
        let r_hat: Direction<GCRS> = s.position.direction_unchecked();
        let n_hat = r_hat
            .cross(&v_hat)
            .expect("position and velocity must not be parallel");
        let c_hat = v_hat
            .cross(&n_hat)
            .expect("v and n are orthogonal by construction");
        let rotation = Rotation3::from_matrix_unchecked([
            v_hat.as_array(),
            n_hat.as_array(),
            c_hat.as_array(),
        ]);
        Self::from_rotation(rotation)
    }
}

impl LocalOrbitalFrame<LVLH> {
    /// Build the GCRS→LVLH rotation for the given state.
    pub fn from_state(s: &OrbitState) -> Self {
        let r_hat: Direction<GCRS> = s.position.direction_unchecked();
        let v_hat: Direction<GCRS> = velocity_direction(s);
        let z_hat = r_hat.negate();
        let x_hat = v_hat
            .cross(&r_hat)
            .and_then(|n| n.cross(&z_hat))
            .expect("velocity and position must not be parallel");
        let y_hat = z_hat
            .cross(&x_hat)
            .expect("z and x are orthogonal by construction");
        let rotation = Rotation3::from_matrix_unchecked([
            x_hat.as_array(),
            y_hat.as_array(),
            z_hat.as_array(),
        ]);
        Self::from_rotation(rotation)
    }
}

// =============================================================================
// Typed transform helpers
// =============================================================================

impl<M: affn::frames::ReferenceFrame> LocalOrbitalFrame<M> {
    /// Rotate a GCRS displacement into the local frame `M`.
    ///
    /// The result carries the local frame tag in its type, so it cannot be
    /// accidentally mixed with inertial displacements.
    #[inline]
    pub fn to_local(&self, v: Displacement<GCRS, Kilometer>) -> Displacement<M, Kilometer> {
        (self.rotation * v).reinterpret_frame()
    }

    /// Rotate a local-frame displacement back into GCRS.
    #[inline]
    pub fn to_inertial(&self, v: Displacement<M, Kilometer>) -> Displacement<GCRS, Kilometer> {
        (self.rotation.transpose() * v).reinterpret_frame()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::time::JulianDate;

    fn circular_orbit() -> OrbitState {
        OrbitState::new(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn rtn_basis_orthonormal_for_circular_orbit() {
        let s = circular_orbit();
        let r = LocalOrbitalFrame::<RTN>::from_state(&s).rotation();
        let m = r.as_matrix();
        assert!((m[0][0] - 1.0).abs() < 1e-12);
        assert!((m[1][1] - 1.0).abs() < 1e-12);
        assert!((m[2][2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn vnc_basis_aligned_with_velocity_for_circular_orbit() {
        let s = circular_orbit();
        let r = LocalOrbitalFrame::<VNC>::from_state(&s).rotation();
        let m = r.as_matrix();
        assert!((m[0][1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn rotate_gcrs_to_rtn_typed_result() {
        let s = circular_orbit();
        let f = LocalOrbitalFrame::<RTN>::from_state(&s);
        let d = Displacement::<GCRS, Kilometer>::new(1.0, 0.0, 0.0);
        let d_rtn = f.to_local(d);
        assert!((d_rtn.x().value() - 1.0).abs() < 1e-12);
        assert!(d_rtn.y().value().abs() < 1e-12);
        assert!(d_rtn.z().value().abs() < 1e-12);
    }

    #[test]
    fn round_trip_through_local_frame() {
        let s = circular_orbit();
        let f = LocalOrbitalFrame::<RTN>::from_state(&s);
        let d = Displacement::<GCRS, Kilometer>::new(1.5, -0.3, 0.7);
        let d_back = f.to_inertial(f.to_local(d));
        assert!((d_back.x().value() - 1.5).abs() < 1e-12);
        assert!((d_back.y().value() - (-0.3)).abs() < 1e-12);
        assert!((d_back.z().value() - 0.7).abs() < 1e-12);
    }
}
