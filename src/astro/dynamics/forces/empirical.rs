// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Constant empirical acceleration in the RTN local orbital frame.
//!
//! ## Scope
//!
//! Provides [`EmpiricalAcceleration`], a force model that injects a constant
//! user-supplied acceleration vector expressed in the Radial/Transverse/Normal
//! (RTN) local orbital frame, rotated into GCRS at each call.
//!
//! ## Equations
//!
//! The RTN frame is defined as:
//! - **R**: along the position vector (radial direction)
//! - **T**: transverse = N×R (along-track component)
//! - **N**: orbit normal = normalise(r×v)
//!
//! The force model rotates the constant RTN acceleration vector into GCRS:
//!
//! ```text
//! a_GCRS = R_{GCRS←RTN} · a_RTN
//! ```
//!
//! ## Units & frames
//!
//! RTN components in km/s².  Rotation matrix is dimensionless and frame-preserving.
//! Output acceleration in km/s² (GCRS).
//!
//! ## Failure modes
//!
//! Returns [`DynamicsError::DegenerateGeometry`] when the RTN frame cannot be
//! constructed:
//! - Zero position magnitude (`|r| = 0`)
//! - Zero velocity magnitude (`|v| = 0`)
//! - Position parallel to velocity (`r ∥ v`)
//!
//! ## References
//!
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, §3.5.
//! * Montenbruck & Gill, *Satellite Orbits*, §2.2.

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::frames::{LocalOrbitalFrame, RTN};
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::frames::GCRS;
use crate::qtty::KmPerSecondsSquared;

use super::traits::ForceModel;

/// Constant empirical acceleration expressed in the RTN local orbital frame.
///
/// The (R, T, N) components are rotated into GCRS at each evaluation using
/// the instantaneous state.  Useful for fitting unmodelled accelerations in
/// orbit determination.
#[derive(Debug, Clone, Copy)]
pub struct EmpiricalAcceleration {
    /// Radial component (km/s²).
    pub radial: KmPerSecondsSquared,
    /// Transverse component (km/s²).
    pub transverse: KmPerSecondsSquared,
    /// Normal (out-of-plane) component (km/s²).
    pub normal: KmPerSecondsSquared,
}

impl EmpiricalAcceleration {
    /// Build an empirical acceleration from typed RTN components.
    pub fn rtn(
        radial: KmPerSecondsSquared,
        transverse: KmPerSecondsSquared,
        normal: KmPerSecondsSquared,
    ) -> Self {
        Self {
            radial,
            transverse,
            normal,
        }
    }
}

impl ForceModel for EmpiricalAcceleration {
    fn acceleration(
        &self,
        s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        // Build the inertial→RTN rotation.
        let rtn_frame = LocalOrbitalFrame::<RTN>::try_from_state(s).map_err(|e| {
            DynamicsError::DegenerateGeometry {
                reason: match e {
                    crate::astro::dynamics::errors::LocalFrameError::ZeroPositionMagnitude => {
                        "zero position magnitude — RTN frame undefined"
                    }
                    crate::astro::dynamics::errors::LocalFrameError::ZeroVelocityMagnitude => {
                        "zero velocity magnitude — RTN frame undefined"
                    }
                    crate::astro::dynamics::errors::LocalFrameError::PositionAndVelocityParallel => {
                        "position and velocity are parallel — RTN frame undefined"
                    }
                },
            }
        })?;

        // Acceleration in RTN frame as raw array [R, T, N].
        let a_rtn = [
            self.radial.value(),
            self.transverse.value(),
            self.normal.value(),
        ];

        // Rotate from RTN → inertial (GCRS): use the transpose (RTN→inertial).
        let rot_to_inertial = rtn_frame.rotation_inverse();
        let a_gcrs = rot_to_inertial.apply_array(a_rtn);

        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            a_gcrs[0], a_gcrs[1], a_gcrs[2],
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::state::OrbitState;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::time::JulianDate;

    fn leo() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::from_raw_unchecked(qtty::Day::new(2_451_545.0)),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    fn zero_velocity_state() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::from_raw_unchecked(qtty::Day::new(2_451_545.0)),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 0.0, 0.0),
        )
    }

    #[test]
    fn empirical_magnitude_preserved() {
        let r = 1e-9_f64;
        let t = 2e-9_f64;
        let n = 3e-9_f64;
        let expected_mag = (r * r + t * t + n * n).sqrt();

        let f = EmpiricalAcceleration::rtn(
            KmPerSecondsSquared::new(r),
            KmPerSecondsSquared::new(t),
            KmPerSecondsSquared::new(n),
        );
        let ctx = DynamicsContext::empty();
        let a = f.acceleration(&leo(), &ctx).unwrap();

        let ax = a.x().value();
        let ay = a.y().value();
        let az = a.z().value();
        let got_mag = (ax * ax + ay * ay + az * az).sqrt();

        assert!(
            (got_mag - expected_mag).abs() < 1e-15,
            "magnitude should be preserved under rotation: expected {expected_mag:.3e}, got {got_mag:.3e}"
        );
    }

    #[test]
    fn pure_radial_empirical_is_collinear_with_position() {
        let f = EmpiricalAcceleration::rtn(
            KmPerSecondsSquared::new(1e-9),
            KmPerSecondsSquared::new(0.0),
            KmPerSecondsSquared::new(0.0),
        );
        let ctx = DynamicsContext::empty();
        let state = leo();
        let a = f.acceleration(&state, &ctx).unwrap();

        let rx = state.position.x().value();
        let ry = state.position.y().value();
        let rz = state.position.z().value();

        let ax = a.x().value();
        let ay = a.y().value();
        let az = a.z().value();

        // cross product r × a must be zero for collinear vectors
        let cross_x = ry * az - rz * ay;
        let cross_y = rz * ax - rx * az;
        let cross_z = rx * ay - ry * ax;
        let cross_mag = (cross_x * cross_x + cross_y * cross_y + cross_z * cross_z).sqrt();

        assert!(
            cross_mag < 1e-15,
            "pure radial empirical acceleration should be parallel to position, cross={cross_mag:.3e}"
        );
    }

    #[test]
    fn degenerate_zero_velocity_returns_error() {
        let f = EmpiricalAcceleration::rtn(
            KmPerSecondsSquared::new(1e-9),
            KmPerSecondsSquared::new(0.0),
            KmPerSecondsSquared::new(0.0),
        );
        let ctx = DynamicsContext::empty();
        let result = f.acceleration(&zero_velocity_state(), &ctx);
        assert!(
            matches!(result, Err(DynamicsError::DegenerateGeometry { .. })),
            "expected DegenerateGeometry error for zero-velocity state"
        );
    }
}
