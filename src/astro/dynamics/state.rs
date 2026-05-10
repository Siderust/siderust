// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Spacecraft and orbit state primitives.
//!
//! These types carry only the *minimum* information force models and
//! integrators need; richer state (including spacecraft mass, parameter
//! blocks, and provenance) lives in [`SpacecraftState`].
//!
//! ## Typed fields
//!
//! Position, velocity, and acceleration are stored directly as [`affn`]
//! typed values so that frame and unit constraints are enforced at compile
//! time:
//!
//! ```rust
//! use siderust::astro::dynamics::{OrbitState, Position, Velocity};
//! use siderust::time::JulianDate;
//! use siderust::coordinates::frames::GCRS;
//!
//! let pos = Position::<GCRS>::new(7000.0, 0.0, 0.0);
//! let vel = Velocity::<GCRS>::new(0.0, 7.5, 0.0);
//! let s = OrbitState::new(JulianDate::new(2_451_545.0), pos, vel);
//! assert!((s.position.x().value() - 7000.0).abs() < 1e-12);
//! ```

use crate::coordinates::cartesian;
use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::GCRS;
use crate::qtty::unit::{Kilometer, Per, Second as SecondUnit};
use crate::qtty::{Kilograms, Second, SquareMeters};
use crate::time::JulianDate;
use affn::cartesian::Vector;

// =============================================================================
// Type aliases for the propagated state
// =============================================================================

/// Geocentric inertial position, default `Kilometer`.
///
/// Frame defaults to [`GCRS`] for use as `Position::<GCRS>` /
/// `Position::<GCRS, Kilometer>`.
pub type Position<S = GCRS, U = Kilometer> = cartesian::Position<Geocentric, S, U>;

/// Inertial velocity vector, default `km/s` in [`GCRS`].
pub type Velocity<S = GCRS, U = Per<Kilometer, SecondUnit>> = cartesian::Velocity<S, U>;

/// Inertial acceleration vector, default `km/s²` in [`GCRS`].
///
/// Stored as the canonical `Per<velocity, time>` unit so that it composes
/// dimensionally with [`Velocity`] under [`qtty`] arithmetic.
pub type Acceleration<S = GCRS, U = Per<Per<Kilometer, SecondUnit>, SecondUnit>> = Vector<S, U>;

/// Default unit of [`Velocity`] used by the propagator.
pub type VelocityUnit = Per<Kilometer, SecondUnit>;

/// Default unit of [`Acceleration`] used by the propagator.
pub type AccelerationUnit = Per<Per<Kilometer, SecondUnit>, SecondUnit>;

// =============================================================================
// OrbitState
// =============================================================================

/// Cartesian inertial position + velocity in km / (km/s) in [`GCRS`].
///
/// Position and velocity are stored as typed [`affn`] vectors, giving
/// compile-time frame and unit guarantees.
#[derive(Debug, Clone, Copy)]
pub struct OrbitState {
    /// Epoch (TT scale, Julian Date).
    pub epoch_tt: JulianDate,
    /// Position in GCRS, km.
    pub position: Position<GCRS, Kilometer>,
    /// Velocity in GCRS, km/s.
    pub velocity: Velocity<GCRS, VelocityUnit>,
}

impl PartialEq for OrbitState {
    fn eq(&self, other: &Self) -> bool {
        self.epoch_tt == other.epoch_tt
            && self.position.x() == other.position.x()
            && self.position.y() == other.position.y()
            && self.position.z() == other.position.z()
            && self.velocity.x() == other.velocity.x()
            && self.velocity.y() == other.velocity.y()
            && self.velocity.z() == other.velocity.z()
    }
}

impl OrbitState {
    /// Construct from typed GCRS position and velocity.
    #[inline]
    pub fn new(
        epoch_tt: JulianDate,
        position: Position<GCRS>,
        velocity: Velocity<GCRS>,
    ) -> Self {
        Self {
            epoch_tt,
            position,
            velocity,
        }
    }

    /// Advance position and velocity by `dt` along `deriv`.
    ///
    /// The epoch is **not** updated — the caller is responsible for advancing
    /// `epoch_tt` to the new time. This mirrors the mathematical step
    /// `x(t + h) ≈ x(t) + h · ẋ(t)`.
    #[inline]
    pub fn advance(&self, deriv: &StateDerivative, dt: Second) -> Self {
        let dt_s = dt.value();
        let new_pos = Position::<GCRS, Kilometer>::new(
            self.position.x().value() + dt_s * deriv.vel.x().value(),
            self.position.y().value() + dt_s * deriv.vel.y().value(),
            self.position.z().value() + dt_s * deriv.vel.z().value(),
        );
        let new_vel = Velocity::<GCRS, VelocityUnit>::new(
            self.velocity.x().value() + dt_s * deriv.acc.x().value(),
            self.velocity.y().value() + dt_s * deriv.acc.y().value(),
            self.velocity.z().value() + dt_s * deriv.acc.z().value(),
        );
        Self {
            epoch_tt: self.epoch_tt,
            position: new_pos,
            velocity: new_vel,
        }
    }

    /// Build the typed velocity `(km/s)` from the orbit state.
    ///
    /// Convenience accessor: returns a clone of [`OrbitState::velocity`] with
    /// the canonical default unit, useful when constructing
    /// [`StateDerivative::position_rate_from_state`].
    #[inline]
    pub fn velocity_typed(&self) -> Velocity<GCRS, VelocityUnit> {
        self.velocity
    }
}

// =============================================================================
// StateDerivative
// =============================================================================

/// Time derivative of an [`OrbitState`]: the 6-vector `[dr/dt, dv/dt]`.
///
/// `vel` is the position rate (= the orbit-state velocity, km/s) and `acc`
/// is the velocity rate (= inertial acceleration in km/s²). Both are typed
/// [`affn`] vectors in [`GCRS`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StateDerivative {
    /// Position rate (= velocity), km/s.
    pub vel: Velocity<GCRS, VelocityUnit>,
    /// Velocity rate (= acceleration), km/s².
    pub acc: Acceleration<GCRS, AccelerationUnit>,
}

impl StateDerivative {
    /// Construct from a typed velocity and a typed acceleration.
    #[inline]
    pub fn new(
        vel: Velocity<GCRS, VelocityUnit>,
        acc: Acceleration<GCRS, AccelerationUnit>,
    ) -> Self {
        Self { vel, acc }
    }

    /// Convenience constructor from raw `f64` km/s and km/s² components.
    ///
    /// Provided as a transitional shim while [`crate::astro::dynamics`]
    /// force models and integrators are migrated to the typed API. New
    /// code should call [`StateDerivative::new`] with typed vectors.
    #[inline]
    pub fn from_components(vel_kms: [f64; 3], acc_km_s2: [f64; 3]) -> Self {
        let vel = Velocity::<GCRS, VelocityUnit>::new(vel_kms[0], vel_kms[1], vel_kms[2]);
        let acc =
            Acceleration::<GCRS, AccelerationUnit>::new(acc_km_s2[0], acc_km_s2[1], acc_km_s2[2]);
        Self { vel, acc }
    }

    /// Return the velocity (position rate) as a typed vector.
    #[inline]
    pub fn velocity(&self) -> Velocity<GCRS, VelocityUnit> {
        self.vel
    }

    /// Return the acceleration (velocity rate) as a typed vector.
    #[inline]
    pub fn acceleration(&self) -> Acceleration<GCRS, AccelerationUnit> {
        self.acc
    }

    /// Weighted RK4 combination of four stages: `(k1 + 2 k2 + 2 k3 + k4) / 6`.
    #[inline]
    pub fn rk4_combine(k1: &Self, k2: &Self, k3: &Self, k4: &Self) -> Self {
        Self {
            vel: (k1.vel + k2.vel.scale(2.0) + k3.vel.scale(2.0) + k4.vel).scale(1.0 / 6.0),
            acc: (k1.acc + k2.acc.scale(2.0) + k3.acc.scale(2.0) + k4.acc).scale(1.0 / 6.0),
        }
    }

    /// Scale this derivative by a dimensionless factor.
    #[inline]
    pub fn scaled(&self, factor: f64) -> Self {
        Self {
            vel: self.vel.scale(factor),
            acc: self.acc.scale(factor),
        }
    }

    /// Element-wise addition.
    #[inline]
    pub fn add(&self, other: &Self) -> Self {
        Self {
            vel: self.vel + other.vel,
            acc: self.acc + other.acc,
        }
    }
}

// =============================================================================
// SpacecraftProperties + SpacecraftState
// =============================================================================

/// Spacecraft properties carried alongside the orbit state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpacecraftProperties {
    /// Total mass.
    pub mass: Kilograms,
    /// Cross-section for atmospheric drag.
    pub drag_area: SquareMeters,
    /// Drag coefficient (dimensionless).
    pub cd: f64,
    /// Cross-section for solar radiation pressure.
    pub srp_area: SquareMeters,
    /// SRP coefficient (dimensionless).
    pub cr: f64,
}

impl SpacecraftProperties {
    /// Reasonable demo defaults for a small LEO platform.
    pub fn demo_leo() -> Self {
        Self {
            mass: Kilograms::new(500.0),
            drag_area: SquareMeters::new(2.0),
            cd: 2.2,
            srp_area: SquareMeters::new(2.0),
            cr: 1.3,
        }
    }
}

/// Combined spacecraft state used by propagators and estimators.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpacecraftState {
    /// Inertial Cartesian orbit state.
    pub orbit: OrbitState,
    /// Spacecraft physical properties.
    pub properties: SpacecraftProperties,
}

// =============================================================================
// (No inherent impls on the type aliases — methods live on the underlying
// `affn` types via the existing `cartesian::Position` / `affn::Vector` impls.)
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn typed_roundtrip_preserves_values() {
        let epoch = JulianDate::new(2_451_545.0);
        let pos = Position::<GCRS>::new(7000.0, 100.0, -200.0);
        let vel = Velocity::<GCRS>::new(0.5, 7.4, -0.1);

        let s = OrbitState::new(epoch, pos, vel);

        assert!((s.position.x().value() - 7000.0).abs() < f64::EPSILON);
        assert!((s.position.y().value() - 100.0).abs() < f64::EPSILON);
        assert!((s.position.z().value() - (-200.0)).abs() < f64::EPSILON);
        assert!((s.velocity.x().value() - 0.5).abs() < f64::EPSILON);
        assert!((s.velocity.y().value() - 7.4).abs() < f64::EPSILON);
        assert!((s.velocity.z().value() - (-0.1)).abs() < f64::EPSILON);
    }

    #[test]
    fn advance_applies_derivative_correctly() {
        let epoch = JulianDate::new(2_451_545.0);
        let pos = Position::<GCRS>::new(7000.0, 0.0, 0.0);
        let vel = Velocity::<GCRS>::new(0.0, 7.5, 0.0);
        let s = OrbitState::new(epoch, pos, vel);

        let deriv = StateDerivative::from_components([0.0, 7.5, 0.0], [0.0, 0.0, -9.8e-3]);
        let dt = Second::new(10.0);
        let s2 = s.advance(&deriv, dt);

        assert!((s2.position.x().value() - 7000.0).abs() < 1e-10);
        assert!((s2.position.y().value() - 75.0).abs() < 1e-10);
        assert!((s2.velocity.z().value() - (-0.098)).abs() < 1e-10);
        // Epoch is unchanged by advance.
        assert_eq!(s2.epoch_tt, epoch);
    }

    #[test]
    fn rk4_combine_is_weighted_mean() {
        let k = StateDerivative::from_components([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let combined = StateDerivative::rk4_combine(&k, &k, &k, &k);
        // (1 + 2 + 2 + 1) / 6 = 1
        assert!((combined.vel.x().value() - 1.0).abs() < 1e-12);
        assert!((combined.acc.y().value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn spacecraft_properties_demo_leo() {
        let p = SpacecraftProperties::demo_leo();
        assert!((p.mass.value() - 500.0).abs() < f64::EPSILON);
        assert!((p.drag_area.value() - 2.0).abs() < f64::EPSILON);
        assert!((p.srp_area.value() - 2.0).abs() < f64::EPSILON);
    }
}
