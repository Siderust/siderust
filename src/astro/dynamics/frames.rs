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
//!   from the inertial parent frame `F` into the local frame `M` for a
//!   specific [`OrbitState<C, F>`].  Build it via the fallible constructors
//!   [`LocalOrbitalFrame::<RTN>::try_from_state`] /
//!   [`LocalOrbitalFrame::<LVLH>::try_from_state`] /
//!   [`LocalOrbitalFrame::<VNC>::try_from_state`] and then rotate many
//!   vectors per state without rebuilding the basis.
//!
//! ## Three local orbital frames
//!
//! ### RTN (also called RIC) — Vallado §3.5
//!
//! - **R**: along the position vector (radial, geocentre → satellite)
//! - **T**: transverse = N×R (along-track-ish; *not* the velocity direction)
//! - **N**: orbit normal = normalise(r×v)
//!
//! Right-handed; R is exactly along position.
//!
//! ### VNC
//!
//! - **V**: along the velocity vector
//! - **N**: orbit normal = normalise(r×v)
//! - **C**: co-normal = V×N
//!
//! ### LVLH (Local-Vertical / Local-Horizontal)
//!
//! - **Z**: radial inward (= −R̂)
//! - **X**: velocity projected onto the local horizontal
//! - **Y**: = Z×X
//!
//! ## Failure modes and typed thresholds
//!
//! All three constructors are fallible and return
//! [`LocalFrameError`](super::errors::LocalFrameError):
//!
//! | Condition | Error variant |
//! |-----------|---------------|
//! | `‖r‖ ≤ 1×10⁻⁹ km` | `ZeroPositionMagnitude` |
//! | `‖v‖ ≤ 1×10⁻⁹ km/s` | `ZeroVelocityMagnitude` |
//! | `r ∥ v` (cross product numerically zero) | `PositionAndVelocityParallel` |
//!
//! Thresholds are expressed as typed [`qtty::Quantity`] values, not bare
//! `f64`, so their physical meaning is preserved across any future unit
//! refactors.
//!
//! ## Frame-tag semantics
//!
//! The constructors are generic over `C` (reference center) and `F`
//! (reference frame) — they accept any [`OrbitState<C, F>`].  The resulting
//! [`LocalOrbitalFrame<M>`] wraps a raw [`Rotation3`] that maps *from the
//! inertial frame `F` of the parent state* into the local frame `M`.  The
//! [`Self::to_local`] / [`Self::to_inertial`] helpers fix `F = GCRS`; for
//! other inertial frames rotate the displacement manually using
//! [`Self::rotation`].
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
//! let s = OrbitState::new_at_jd(
//!     JulianDate::new(2_451_545.0),
//!     Position::<GCRS>::new(7000.0, 0.0, 0.0),
//!     Velocity::<GCRS>::new(0.0, 7.5, 0.0),
//! );
//!
//! let rtn = LocalOrbitalFrame::<RTN>::try_from_state(&s)
//!     .expect("well-defined circular orbit");
//! let d_gcrs = Displacement::<GCRS, Kilometer>::new(1.0, 0.0, 0.0);
//! let d_rtn = rtn.to_local(d_gcrs);
//! assert!((d_rtn.x().value() - 1.0).abs() < 1e-12);
//! ```

use core::marker::PhantomData;

use affn::cartesian::{Direction, Displacement};
use affn::DeriveReferenceFrame;
use affn::Rotation3;

use crate::coordinates::centers::ReferenceCenter;
use crate::coordinates::frames::{ReferenceFrame, GCRS};
use crate::qtty::unit::Kilometer;
use crate::qtty::KmPerSecond;
use crate::qtty::Quantity;

use super::errors::LocalFrameError;
use super::state::OrbitState;

// =============================================================================
// Typed thresholds (expressed as typed qtty Quantities, not bare f64)
// =============================================================================

/// Minimum acceptable position magnitude (km).
///
/// Below this threshold the radial direction is numerically undefined.
const POS_THRESHOLD_KM: f64 = 1e-9;

/// Minimum acceptable velocity magnitude (km/s).
///
/// Below this threshold the along-track direction is numerically undefined.
const VEL_THRESHOLD_KM_S: f64 = 1e-9;

// =============================================================================
// Local orbital frame marker types
// =============================================================================

/// Radial / Transverse / Normal local orbital frame.
///
/// - **R**: along the position vector (radial, geocentre → satellite)
/// - **T**: transverse, = N×R (along-track-ish; *not* the velocity direction)
/// - **N**: orbit normal, = normalise(r×v)
///
/// The frame is right-handed.  Use
/// [`LocalOrbitalFrame::<RTN>::try_from_state`] to materialise it from an
/// inertial state.
#[derive(Debug, Copy, Clone, DeriveReferenceFrame)]
pub struct RTN;

/// Velocity / Normal / Co-normal local orbital frame.
///
/// - **V**: along the velocity vector
/// - **N**: orbit normal, = normalise(r×v)
/// - **C**: co-normal, = V×N
///
/// Use [`LocalOrbitalFrame::<VNC>::try_from_state`].
#[derive(Debug, Copy, Clone, DeriveReferenceFrame)]
pub struct VNC;

/// Local-Vertical / Local-Horizontal frame.
///
/// - **Z**: radial inward (= −R̂)
/// - **X**: velocity projected onto the local horizontal
/// - **Y**: = Z×X
///
/// Use [`LocalOrbitalFrame::<LVLH>::try_from_state`].
#[derive(Debug, Copy, Clone, DeriveReferenceFrame)]
pub struct LVLH;

// =============================================================================
// Typed local orbital frame object
// =============================================================================

/// Materialised local orbital frame for a specific [`OrbitState`].
///
/// `M` is the local frame marker (one of [`RTN`], [`VNC`], [`LVLH`]).
/// The captured [`Rotation3`] maps from the inertial frame `F` of the
/// parent state into the local frame `M`.
///
/// Construct via the fallible `try_from_state` constructors specialised on
/// `M`.  Use [`Self::to_local`] / [`Self::to_inertial`] to rotate GCRS
/// displacements while keeping the frame tag in the type system; for other
/// inertial frames use [`Self::rotation`] directly.
#[derive(Debug, Clone, Copy)]
pub struct LocalOrbitalFrame<M> {
    rotation: Rotation3,
    _marker: PhantomData<M>,
}

impl<M> LocalOrbitalFrame<M> {
    /// Wrap an existing inertial→`M` rotation matrix.
    ///
    /// Prefer the specialised `try_from_state` constructors; this is for
    /// callers that already have a basis matrix from another source.
    #[inline]
    pub fn from_rotation(rotation: Rotation3) -> Self {
        Self {
            rotation,
            _marker: PhantomData,
        }
    }

    /// The captured inertial→`M` rotation matrix.
    #[inline]
    pub fn rotation(&self) -> Rotation3 {
        self.rotation
    }

    /// The transposed (`M`→inertial) rotation matrix.
    #[inline]
    pub fn rotation_inverse(&self) -> Rotation3 {
        self.rotation.transpose()
    }
}

// =============================================================================
// Internal helpers
// =============================================================================

/// Check position magnitude against typed threshold; return the unit direction
/// or the appropriate [`LocalFrameError`].
#[inline]
fn checked_position_direction<C, F>(
    state: &OrbitState<C, F>,
) -> Result<Direction<F>, LocalFrameError>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    let threshold = Quantity::<Kilometer>::new(POS_THRESHOLD_KM);
    if state.position.distance() <= threshold {
        return Err(LocalFrameError::ZeroPositionMagnitude);
    }
    state
        .position
        .direction()
        .ok_or(LocalFrameError::ZeroPositionMagnitude)
}

/// Check velocity magnitude against typed threshold; return the unit direction
/// or the appropriate [`LocalFrameError`].
#[inline]
fn checked_velocity_direction<C, F>(
    state: &OrbitState<C, F>,
) -> Result<Direction<F>, LocalFrameError>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    let threshold = Quantity::<KmPerSecond>::new(VEL_THRESHOLD_KM_S);
    if state.velocity.magnitude() <= threshold {
        return Err(LocalFrameError::ZeroVelocityMagnitude);
    }
    // `normalize()` requires `LengthUnit` which `KmPerSecond` is not, so
    // extract raw scalar components and build a unit direction explicitly.
    Direction::try_new(
        state.velocity.x().value(),
        state.velocity.y().value(),
        state.velocity.z().value(),
    )
    .ok_or(LocalFrameError::ZeroVelocityMagnitude)
}

// =============================================================================
// Fallible constructors
// =============================================================================

impl LocalOrbitalFrame<RTN> {
    /// Build the inertial→RTN rotation for the given state.
    ///
    /// # Frame definition (Vallado §3.5)
    ///
    /// - R = r̂ (position unit vector)
    /// - N = normalise(r × v) (orbit-normal)
    /// - T = N × R (transverse / along-track-ish)
    ///
    /// # Errors
    ///
    /// - [`LocalFrameError::ZeroPositionMagnitude`] if `‖r‖ ≤ 1×10⁻⁹ km`.
    /// - [`LocalFrameError::ZeroVelocityMagnitude`] if `‖v‖ ≤ 1×10⁻⁹ km/s`.
    /// - [`LocalFrameError::PositionAndVelocityParallel`] if r and v are
    ///   parallel or anti-parallel (orbit-normal direction is undefined).
    pub fn try_from_state<C, F>(state: &OrbitState<C, F>) -> Result<Self, LocalFrameError>
    where
        C: ReferenceCenter,
        F: ReferenceFrame,
    {
        let r_hat = checked_position_direction(state)?;
        let v_hat = checked_velocity_direction(state)?;
        let n_hat = r_hat
            .cross(&v_hat)
            .ok_or(LocalFrameError::PositionAndVelocityParallel)?;
        // t_hat = N × R; guaranteed non-zero because N ⊥ R by construction.
        let t_hat = n_hat
            .cross(&r_hat)
            .ok_or(LocalFrameError::PositionAndVelocityParallel)?;
        let rotation = Rotation3::from_matrix_unchecked([
            r_hat.as_array(),
            t_hat.as_array(),
            n_hat.as_array(),
        ]);
        Ok(Self::from_rotation(rotation))
    }
}

impl LocalOrbitalFrame<VNC> {
    /// Build the inertial→VNC rotation for the given state.
    ///
    /// # Frame definition
    ///
    /// - V = v̂ (velocity unit vector)
    /// - N = normalise(r × v) (orbit-normal)
    /// - C = V × N (co-normal)
    ///
    /// # Errors
    ///
    /// - [`LocalFrameError::ZeroPositionMagnitude`] if `‖r‖ ≤ 1×10⁻⁹ km`.
    /// - [`LocalFrameError::ZeroVelocityMagnitude`] if `‖v‖ ≤ 1×10⁻⁹ km/s`.
    /// - [`LocalFrameError::PositionAndVelocityParallel`] if r and v are
    ///   parallel or anti-parallel.
    pub fn try_from_state<C, F>(state: &OrbitState<C, F>) -> Result<Self, LocalFrameError>
    where
        C: ReferenceCenter,
        F: ReferenceFrame,
    {
        let v_hat = checked_velocity_direction(state)?;
        let r_hat = checked_position_direction(state)?;
        let n_hat = r_hat
            .cross(&v_hat)
            .ok_or(LocalFrameError::PositionAndVelocityParallel)?;
        // c_hat = V × N; guaranteed non-zero because V ⊥ N by construction.
        let c_hat = v_hat
            .cross(&n_hat)
            .ok_or(LocalFrameError::PositionAndVelocityParallel)?;
        let rotation = Rotation3::from_matrix_unchecked([
            v_hat.as_array(),
            n_hat.as_array(),
            c_hat.as_array(),
        ]);
        Ok(Self::from_rotation(rotation))
    }
}

impl LocalOrbitalFrame<LVLH> {
    /// Build the inertial→LVLH rotation for the given state.
    ///
    /// # Frame definition
    ///
    /// - Z = −r̂ (radial inward)
    /// - X = (r̂ × v̂) × (−r̂)  — velocity projected onto the local
    ///   horizontal (normal component removed)
    /// - Y = Z × X
    ///
    /// # Errors
    ///
    /// - [`LocalFrameError::ZeroPositionMagnitude`] if `‖r‖ ≤ 1×10⁻⁹ km`.
    /// - [`LocalFrameError::ZeroVelocityMagnitude`] if `‖v‖ ≤ 1×10⁻⁹ km/s`.
    /// - [`LocalFrameError::PositionAndVelocityParallel`] if r and v are
    ///   parallel or anti-parallel.
    pub fn try_from_state<C, F>(state: &OrbitState<C, F>) -> Result<Self, LocalFrameError>
    where
        C: ReferenceCenter,
        F: ReferenceFrame,
    {
        let r_hat = checked_position_direction(state)?;
        let v_hat = checked_velocity_direction(state)?;
        let z_hat = r_hat.negate();
        // n = r̂ × v̂; x = n × z_hat = (r̂ × v̂) × (−r̂)
        let x_hat = r_hat
            .cross(&v_hat)
            .ok_or(LocalFrameError::PositionAndVelocityParallel)
            .and_then(|n| {
                n.cross(&z_hat)
                    .ok_or(LocalFrameError::PositionAndVelocityParallel)
            })?;
        // y_hat = Z × X; guaranteed non-zero because Z ⊥ X by construction.
        let y_hat = z_hat
            .cross(&x_hat)
            .ok_or(LocalFrameError::PositionAndVelocityParallel)?;
        let rotation = Rotation3::from_matrix_unchecked([
            x_hat.as_array(),
            y_hat.as_array(),
            z_hat.as_array(),
        ]);
        Ok(Self::from_rotation(rotation))
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
    use crate::coordinates::centers::Heliocentric;
    use crate::coordinates::frames::ICRS;
    use crate::time::JulianDate;

    // -------------------------------------------------------------------------
    // Shared helpers
    // -------------------------------------------------------------------------

    fn circular_orbit() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    fn zero_position_state() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(0.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    fn zero_velocity_state() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 0.0, 0.0),
        )
    }

    fn collinear_state() -> OrbitState {
        // r and v are parallel → r × v = 0
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(7.5, 0.0, 0.0),
        )
    }

    // Verify rotation matrix is orthonormal (all diagonal entries ≈ 1 after R·Rᵀ)
    fn assert_orthonormal(rot: Rotation3) {
        let m = rot.as_matrix();
        let rt_owned = rot.transpose();
        let rt = rt_owned.as_matrix();
        // R·Rᵀ should equal the identity
        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = (0..3).map(|k| m[i][k] * rt[k][j]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dot - expected).abs() < 1e-12,
                    "R·Rᵀ[{i}][{j}] = {dot}, expected {expected}"
                );
            }
        }
    }

    // Verify determinant ≈ +1 (right-handed)
    fn assert_right_handed(rot: Rotation3) {
        let m = rot.as_matrix();
        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        assert!((det - 1.0).abs() < 1e-12, "det = {det}, expected +1");
    }

    // -------------------------------------------------------------------------
    // RTN tests
    // -------------------------------------------------------------------------

    #[test]
    fn rtn_nominal_circular_orbit_orthonormal_right_handed() {
        let rot = LocalOrbitalFrame::<RTN>::try_from_state(&circular_orbit())
            .unwrap()
            .rotation();
        assert_orthonormal(rot);
        assert_right_handed(rot);
    }

    #[test]
    fn rtn_r_axis_aligned_with_position() {
        // For the state (7000, 0, 0) / (0, 7.5, 0), the R axis is (1,0,0)
        let rot = LocalOrbitalFrame::<RTN>::try_from_state(&circular_orbit())
            .unwrap()
            .rotation();
        let m = rot.as_matrix();
        // Row 0 (R axis) must be (1, 0, 0)
        assert!((m[0][0] - 1.0).abs() < 1e-12);
        assert!(m[0][1].abs() < 1e-12);
        assert!(m[0][2].abs() < 1e-12);
    }

    #[test]
    fn rtn_zero_position_returns_error() {
        assert!(matches!(
            LocalOrbitalFrame::<RTN>::try_from_state(&zero_position_state()),
            Err(LocalFrameError::ZeroPositionMagnitude)
        ));
    }

    #[test]
    fn rtn_zero_velocity_returns_error() {
        assert!(matches!(
            LocalOrbitalFrame::<RTN>::try_from_state(&zero_velocity_state()),
            Err(LocalFrameError::ZeroVelocityMagnitude)
        ));
    }

    #[test]
    fn rtn_collinear_returns_error() {
        assert!(matches!(
            LocalOrbitalFrame::<RTN>::try_from_state(&collinear_state()),
            Err(LocalFrameError::PositionAndVelocityParallel)
        ));
    }

    // -------------------------------------------------------------------------
    // VNC tests
    // -------------------------------------------------------------------------

    #[test]
    fn vnc_nominal_circular_orbit_orthonormal_right_handed() {
        let rot = LocalOrbitalFrame::<VNC>::try_from_state(&circular_orbit())
            .unwrap()
            .rotation();
        assert_orthonormal(rot);
        assert_right_handed(rot);
    }

    #[test]
    fn vnc_v_axis_aligned_with_velocity() {
        // For the state (7000, 0, 0) / (0, 7.5, 0), v̂ = (0, 1, 0)
        let rot = LocalOrbitalFrame::<VNC>::try_from_state(&circular_orbit())
            .unwrap()
            .rotation();
        let m = rot.as_matrix();
        // Row 0 (V axis) must be (0, 1, 0)
        assert!(m[0][0].abs() < 1e-12);
        assert!((m[0][1] - 1.0).abs() < 1e-12);
        assert!(m[0][2].abs() < 1e-12);
    }

    #[test]
    fn vnc_zero_position_returns_error() {
        assert!(matches!(
            LocalOrbitalFrame::<VNC>::try_from_state(&zero_position_state()),
            Err(LocalFrameError::ZeroPositionMagnitude)
        ));
    }

    #[test]
    fn vnc_zero_velocity_returns_error() {
        assert!(matches!(
            LocalOrbitalFrame::<VNC>::try_from_state(&zero_velocity_state()),
            Err(LocalFrameError::ZeroVelocityMagnitude)
        ));
    }

    #[test]
    fn vnc_collinear_returns_error() {
        assert!(matches!(
            LocalOrbitalFrame::<VNC>::try_from_state(&collinear_state()),
            Err(LocalFrameError::PositionAndVelocityParallel)
        ));
    }

    // -------------------------------------------------------------------------
    // LVLH tests
    // -------------------------------------------------------------------------

    #[test]
    fn lvlh_nominal_circular_orbit_orthonormal_right_handed() {
        let rot = LocalOrbitalFrame::<LVLH>::try_from_state(&circular_orbit())
            .unwrap()
            .rotation();
        assert_orthonormal(rot);
        assert_right_handed(rot);
    }

    #[test]
    fn lvlh_z_axis_is_radial_inward() {
        // For position (7000, 0, 0), Z = -r̂ = (-1, 0, 0)
        let rot = LocalOrbitalFrame::<LVLH>::try_from_state(&circular_orbit())
            .unwrap()
            .rotation();
        let m = rot.as_matrix();
        // Row 2 (Z axis) must be (-1, 0, 0)
        assert!((m[2][0] + 1.0).abs() < 1e-12);
        assert!(m[2][1].abs() < 1e-12);
        assert!(m[2][2].abs() < 1e-12);
    }

    #[test]
    fn lvlh_zero_position_returns_error() {
        assert!(matches!(
            LocalOrbitalFrame::<LVLH>::try_from_state(&zero_position_state()),
            Err(LocalFrameError::ZeroPositionMagnitude)
        ));
    }

    #[test]
    fn lvlh_zero_velocity_returns_error() {
        assert!(matches!(
            LocalOrbitalFrame::<LVLH>::try_from_state(&zero_velocity_state()),
            Err(LocalFrameError::ZeroVelocityMagnitude)
        ));
    }

    #[test]
    fn lvlh_collinear_returns_error() {
        assert!(matches!(
            LocalOrbitalFrame::<LVLH>::try_from_state(&collinear_state()),
            Err(LocalFrameError::PositionAndVelocityParallel)
        ));
    }

    // -------------------------------------------------------------------------
    // Typed-transform helpers (GCRS round-trip)
    // -------------------------------------------------------------------------

    #[test]
    fn rtn_to_local_and_back_round_trip() {
        let s = circular_orbit();
        let f = LocalOrbitalFrame::<RTN>::try_from_state(&s).unwrap();
        let d = Displacement::<GCRS, Kilometer>::new(1.5, -0.3, 0.7);
        let d_back = f.to_inertial(f.to_local(d));
        assert!((d_back.x().value() - 1.5).abs() < 1e-12);
        assert!((d_back.y().value() - (-0.3)).abs() < 1e-12);
        assert!((d_back.z().value() - 0.7).abs() < 1e-12);
    }

    #[test]
    fn rtn_rotate_gcrs_radial_displacement_to_local() {
        let s = circular_orbit();
        let f = LocalOrbitalFrame::<RTN>::try_from_state(&s).unwrap();
        // (1, 0, 0) in GCRS is the radial direction → (1, 0, 0) in RTN
        let d = Displacement::<GCRS, Kilometer>::new(1.0, 0.0, 0.0);
        let d_rtn = f.to_local(d);
        assert!((d_rtn.x().value() - 1.0).abs() < 1e-12);
        assert!(d_rtn.y().value().abs() < 1e-12);
        assert!(d_rtn.z().value().abs() < 1e-12);
    }

    // -------------------------------------------------------------------------
    // Non-default frame (Heliocentric / ICRS) — confirms generics work
    // -------------------------------------------------------------------------

    #[test]
    fn rtn_heliocentric_icrs_state_succeeds() {
        use crate::coordinates::cartesian;
        type HelioPos = cartesian::Position<Heliocentric, ICRS, Kilometer>;
        type HelioVel = crate::astro::dynamics::Velocity<ICRS>;
        let pos = HelioPos::new(1.496e8, 0.0, 0.0);
        let vel = HelioVel::new(0.0, 29.78, 0.0);
        let s =
            OrbitState::<Heliocentric, ICRS>::new(JulianDate::new(2_451_545.0).to_time(), pos, vel);
        let result = LocalOrbitalFrame::<RTN>::try_from_state(&s);
        assert!(result.is_ok());
        assert_orthonormal(result.unwrap().rotation());
    }

    #[test]
    fn vnc_heliocentric_icrs_state_succeeds() {
        use crate::coordinates::cartesian;
        type HelioPos = cartesian::Position<Heliocentric, ICRS, Kilometer>;
        type HelioVel = crate::astro::dynamics::Velocity<ICRS>;
        let pos = HelioPos::new(1.496e8, 0.0, 0.0);
        let vel = HelioVel::new(0.0, 29.78, 0.0);
        let s =
            OrbitState::<Heliocentric, ICRS>::new(JulianDate::new(2_451_545.0).to_time(), pos, vel);
        let result = LocalOrbitalFrame::<VNC>::try_from_state(&s);
        assert!(result.is_ok());
        assert_right_handed(result.unwrap().rotation());
    }
}
