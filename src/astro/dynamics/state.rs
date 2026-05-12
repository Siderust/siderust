// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Spacecraft and orbit state primitives.
//!
//! These types carry only the *minimum* information force models and
//! integrators need; richer state (including spacecraft mass, parameter
//! blocks, and provenance) lives in [`SpacecraftState`].
//!
//! ## Scope
//!
//! Provides [`OrbitState<C, F>`] — a frame- and center-tagged Cartesian state
//! `(r, v, epoch)` — and [`StateDerivative<F>`] — the velocity and acceleration
//! time derivatives used by all integrators.
//!
//! ## Equations
//!
//! The Cartesian equations of motion are:
//!
//! ```text
//! dr/dt = v
//! dv/dt = a(r, v, t) = [force model output]
//! ```
//!
//! These are propagated by RK4, DOPRI5, or DOP853 integrators.
//!
//! ## Generic over center and frame
//!
//! [`OrbitState<C, F>`] is parameterised over the reference center `C` and
//! the reference frame `F`.  The defaults (`Geocentric`, `GCRS`) give the
//! conventional geocentric inertial state used by most force models.
//!
//! ## Units & frames
//!
//! - Position: **km**, GCRS frame (Geocentric, inertial)
//! - Velocity: **km/s**, GCRS frame
//! - Acceleration: **km/s²**, GCRS frame
//! - Epoch: **seconds TT** (Terrestrial Time) as a [`Time<TT>`](crate::time::Time)
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
//! let s = OrbitState::new_at_jd(JulianDate::new(2_451_545.0), pos, vel);
//! assert!((s.position.x().value() - 7000.0).abs() < 1e-12);
//! ```
//!
//! ## Failure modes
//!
//! Construction is infallible for valid positions and velocities.
//! Failures arise only in integrators when force models or providers fail.
//!
//! ## Removed API
//!
//! `StateDerivative::from_components` was removed in this release because it
//! bypassed the typed vector API.  Use [`StateDerivative::new`] with typed
//! [`Velocity`] and [`Acceleration`] vectors instead.
//!
//! ```compile_fail
//! use siderust::astro::dynamics::state::StateDerivative;
//! // `from_components` no longer exists — this must not compile:
//! let _ = StateDerivative::from_components([0.0; 3], [0.0; 3]);
//! ```

use crate::coordinates::cartesian;
use crate::coordinates::centers::{Geocentric, ReferenceCenter};
use crate::coordinates::frames::{ReferenceFrame, GCRS};
use crate::qtty::unit::Kilometer;
use crate::qtty::{
    AreaToMass, DragCoefficient, Kilograms, KmPerSecond, KmPerSecondSquared, Second, SquareMeters,
    SrpCoefficient,
};
use crate::time::{JulianDate, Time, JD, TT};
use affn::cartesian::Vector;

// =============================================================================
// Type aliases for the propagated state
// =============================================================================

/// Inertial position in center `C` and frame `F`, default `Geocentric`/`GCRS`/`km`.
///
/// For the usual geocentric-inertial case, use `Position::<GCRS>` or
/// `Position::<GCRS, Kilometer>`.
pub type Position<F = GCRS, U = Kilometer> = cartesian::Position<Geocentric, F, U>;

/// Inertial velocity vector, default `km/s` in [`GCRS`].
///
/// The unit `KmPerSecond` = `Per<Kilometer, Second>` is the canonical
/// astrodynamics velocity unit.
pub type Velocity<F = GCRS, U = KmPerSecond> = cartesian::Velocity<F, U>;

/// Inertial acceleration vector, default `km/s²` in [`GCRS`].
///
/// The unit `KmPerSecondSquared` = `Per<Per<Kilometer, Second>, Second>`.
pub type Acceleration<F = GCRS, U = KmPerSecondSquared> = Vector<F, U>;

/// Default unit of [`Velocity`] used by the propagator (`km/s`).
pub type VelocityUnit = KmPerSecond;

/// Default unit of [`Acceleration`] used by the propagator (`km/s²`).
pub type AccelerationUnit = KmPerSecondSquared;

// =============================================================================
// OrbitState
// =============================================================================

/// Generic Cartesian position + velocity state parameterised over center `C`
/// and frame `F`.
///
/// Defaults to `Geocentric`/[`GCRS`], the standard geocentric inertial frame
/// used by most LEO/MEO force models.  Use the explicit type parameters to
/// express heliocentric, barycentric, or other states:
///
/// ```rust
/// use siderust::astro::dynamics::state::OrbitState;
/// use siderust::coordinates::centers::Heliocentric;
/// use siderust::coordinates::frames::ICRS;
/// use siderust::astro::dynamics::state::{Position, Velocity};
/// use siderust::time::JulianDate;
///
/// // Heliocentric/ICRS state (for inner-planet propagation).
/// type HelioPos = siderust::coordinates::cartesian::Position<
///     Heliocentric,
///     ICRS,
///     siderust::qtty::unit::Kilometer,
/// >;
/// type HelioVel = siderust::coordinates::cartesian::Velocity<
///     ICRS,
///     siderust::qtty::KmPerSecond,
/// >;
/// let pos = HelioPos::new(1.496e8, 0.0, 0.0);
/// let vel = HelioVel::new(0.0, 29.78, 0.0);
/// let s = OrbitState::<Heliocentric, ICRS>::new(
///     JulianDate::new(2_451_545.0).to_time(),
///     pos,
///     vel,
/// );
/// assert!((s.position.x().value() - 1.496e8).abs() < 1.0);
/// ```
///
/// # Units
///
/// - `epoch`: continuous TT instant ([`Time<TT>`]).
/// - `position`: km in frame `F` centred on `C`.
/// - `velocity`: km/s in frame `F`.
#[derive(Debug, Clone)]
pub struct OrbitState<C = Geocentric, F = GCRS>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    /// Epoch as a continuous TT instant.
    ///
    /// Use [`OrbitState::epoch_jd`] to get the Julian Date encoding, and
    /// [`OrbitState::new_at_jd`] to construct from a [`JulianDate`].
    pub epoch: Time<TT>,
    /// Position in frame `F` centred on `C`, km.
    pub position: cartesian::Position<C, F, Kilometer>,
    /// Velocity in frame `F`, km/s.
    pub velocity: Velocity<F, VelocityUnit>,
}

// Manual Copy impl with explicit Params bound (the derive macro cannot express this).
impl<C, F> Copy for OrbitState<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    C::Params: Copy,
{
}

impl<C, F> PartialEq for OrbitState<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn eq(&self, other: &Self) -> bool {
        self.epoch == other.epoch
            && self.position.x() == other.position.x()
            && self.position.y() == other.position.y()
            && self.position.z() == other.position.z()
            && self.velocity.x() == other.velocity.x()
            && self.velocity.y() == other.velocity.y()
            && self.velocity.z() == other.velocity.z()
    }
}

impl<C, F> OrbitState<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    /// Construct from a continuous [`Time<TT>`] epoch, typed position, and
    /// typed velocity.
    ///
    /// # Example
    ///
    /// ```rust
    /// use siderust::astro::dynamics::state::{OrbitState, Position, Velocity};
    /// use siderust::time::{JulianDate, Time, TT};
    ///
    /// let epoch = JulianDate::new(2_451_545.0).to_time();
    /// let pos = Position::new(7000.0, 0.0, 0.0);
    /// let vel = Velocity::new(0.0, 7.5, 0.0);
    /// let s = OrbitState::new(epoch, pos, vel);
    /// ```
    #[inline]
    pub fn new(
        epoch: Time<TT>,
        position: cartesian::Position<C, F, Kilometer>,
        velocity: Velocity<F, VelocityUnit>,
    ) -> Self {
        Self {
            epoch,
            position,
            velocity,
        }
    }

    /// Boundary constructor: build from a [`JulianDate`] (TT) by converting
    /// it to the internal [`Time<TT>`] representation.
    ///
    /// This is the preferred constructor when the epoch comes from an
    /// ephemeris or other source that expresses time as a JD.
    ///
    /// # Example
    ///
    /// ```rust
    /// use siderust::astro::dynamics::state::{OrbitState, Position, Velocity};
    /// use siderust::time::JulianDate;
    ///
    /// let s = OrbitState::new_at_jd(
    ///     JulianDate::new(2_451_545.0),
    ///     Position::new(7000.0, 0.0, 0.0),
    ///     Velocity::new(0.0, 7.5, 0.0),
    /// );
    /// assert!((s.epoch_jd().jd_value() - 2_451_545.0).abs() < 1e-9);
    /// ```
    #[inline]
    pub fn new_at_jd(
        epoch_jd: JulianDate,
        position: cartesian::Position<C, F, Kilometer>,
        velocity: Velocity<F, VelocityUnit>,
    ) -> Self {
        Self::new(epoch_jd.to_time(), position, velocity)
    }

    /// Return the epoch encoded as a TT Julian Date.
    ///
    /// This is the boundary helper used by ephemeris callers that require a
    /// [`JulianDate`] scalar.
    #[inline]
    pub fn epoch_jd(&self) -> JulianDate {
        self.epoch.to::<JD>()
    }

    /// Advance position and velocity by `dt` along `deriv`.
    ///
    /// The epoch is **not** updated — the caller is responsible for advancing
    /// `epoch` to the new time. This mirrors the mathematical step:
    /// `x(t + h) ≈ x(t) + h · ẋ(t)`.
    #[inline]
    pub fn advance(&self, deriv: &StateDerivative<F>, dt: Second) -> Self {
        let dt_s = dt.value();
        let new_pos = cartesian::Position::<C, F, Kilometer>::new_with_params(
            self.position.center_params().clone(),
            self.position.x().value() + dt_s * deriv.vel.x().value(),
            self.position.y().value() + dt_s * deriv.vel.y().value(),
            self.position.z().value() + dt_s * deriv.vel.z().value(),
        );
        let new_vel = Velocity::<F, VelocityUnit>::new(
            self.velocity.x().value() + dt_s * deriv.acc.x().value(),
            self.velocity.y().value() + dt_s * deriv.acc.y().value(),
            self.velocity.z().value() + dt_s * deriv.acc.z().value(),
        );
        Self {
            epoch: self.epoch,
            position: new_pos,
            velocity: new_vel,
        }
    }

    /// Return a copy of the velocity as a typed `km/s` vector.
    ///
    /// Convenience accessor useful when constructing a [`StateDerivative`].
    #[inline]
    pub fn velocity_typed(&self) -> Velocity<F, VelocityUnit> {
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
/// [`affn`] vectors in frame `F`.
///
/// Generic over the reference frame `F` so that derivatives can be formed in
/// any inertial frame without losing the frame tag.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StateDerivative<F = GCRS>
where
    F: ReferenceFrame,
{
    /// Position rate (= velocity), km/s in frame `F`.
    pub vel: Velocity<F, VelocityUnit>,
    /// Velocity rate (= acceleration), km/s² in frame `F`.
    pub acc: Acceleration<F, AccelerationUnit>,
}

impl<F: ReferenceFrame> StateDerivative<F> {
    /// Construct from a typed velocity (km/s) and a typed acceleration (km/s²).
    ///
    /// Both vectors must be expressed in the same frame `F`.
    #[inline]
    pub fn new(vel: Velocity<F, VelocityUnit>, acc: Acceleration<F, AccelerationUnit>) -> Self {
        Self { vel, acc }
    }

    /// Return the velocity (position rate) as a typed vector.
    #[inline]
    pub fn velocity(&self) -> Velocity<F, VelocityUnit> {
        self.vel
    }

    /// Return the acceleration (velocity rate) as a typed vector.
    #[inline]
    pub fn acceleration(&self) -> Acceleration<F, AccelerationUnit> {
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

/// Physical properties of a spacecraft, used by force models.
///
/// Raw area/mass fields are kept for legibility alongside the precomputed
/// area-to-mass ratios consumed by drag and SRP force models.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpacecraftProperties {
    /// Total wet mass (kg).
    pub mass: Kilograms,
    /// Drag cross-section (m²).
    pub drag_area: SquareMeters,
    /// Drag coefficient C_D (dimensionless).
    pub cd: DragCoefficient,
    /// SRP cross-section (m²).
    pub srp_area: SquareMeters,
    /// SRP coefficient C_R (dimensionless).
    pub cr: SrpCoefficient,
    /// Precomputed drag area-to-mass ratio (m²/kg).
    pub area_to_mass_drag: AreaToMass,
    /// Precomputed SRP area-to-mass ratio (m²/kg).
    pub area_to_mass_srp: AreaToMass,
}

impl SpacecraftProperties {
    /// Construct from all physical parameters.
    ///
    /// The area-to-mass ratios are computed from the supplied values so they
    /// are always consistent with `drag_area`, `srp_area`, and `mass`.
    #[inline]
    pub fn new(
        mass: Kilograms,
        drag_area: SquareMeters,
        cd: DragCoefficient,
        srp_area: SquareMeters,
        cr: SrpCoefficient,
    ) -> Self {
        let area_to_mass_drag = AreaToMass::new(drag_area.value() / mass.value());
        let area_to_mass_srp = AreaToMass::new(srp_area.value() / mass.value());
        Self {
            mass,
            drag_area,
            cd,
            srp_area,
            cr,
            area_to_mass_drag,
            area_to_mass_srp,
        }
    }

    /// Drag area-to-mass ratio: `drag_area / mass` (m²/kg).
    #[inline]
    pub fn drag_area_to_mass(&self) -> AreaToMass {
        self.area_to_mass_drag
    }

    /// SRP area-to-mass ratio: `srp_area / mass` (m²/kg).
    #[inline]
    pub fn srp_area_to_mass(&self) -> AreaToMass {
        self.area_to_mass_srp
    }

    /// Reasonable demo defaults for a small LEO platform.
    ///
    /// - mass = 500 kg
    /// - drag area = 2 m², C_D = 2.2
    /// - SRP area  = 2 m², C_R = 1.3
    pub fn demo_leo() -> Self {
        Self::new(
            Kilograms::new(500.0),
            SquareMeters::new(2.0),
            DragCoefficient::new(2.2),
            SquareMeters::new(2.0),
            SrpCoefficient::new(1.3),
        )
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
    use crate::coordinates::centers::Heliocentric;
    use crate::coordinates::frames::ICRS;

    // ── OrbitState round-trip ─────────────────────────────────────────────────

    #[test]
    fn typed_roundtrip_preserves_values() {
        let epoch = JulianDate::new(2_451_545.0);
        let pos = Position::<GCRS>::new(7000.0, 100.0, -200.0);
        let vel = Velocity::<GCRS>::new(0.5, 7.4, -0.1);

        let s = OrbitState::new_at_jd(epoch, pos, vel);

        assert!((s.position.x().value() - 7000.0).abs() < f64::EPSILON);
        assert!((s.position.y().value() - 100.0).abs() < f64::EPSILON);
        assert!((s.position.z().value() - (-200.0)).abs() < f64::EPSILON);
        assert!((s.velocity.x().value() - 0.5).abs() < f64::EPSILON);
        assert!((s.velocity.y().value() - 7.4).abs() < f64::EPSILON);
        assert!((s.velocity.z().value() - (-0.1)).abs() < f64::EPSILON);
    }

    #[test]
    fn new_at_jd_epoch_roundtrip() {
        let jd_in = JulianDate::new(2_451_545.0);
        let s = OrbitState::new_at_jd(
            jd_in,
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        );
        let jd_out = s.epoch_jd();
        assert!(
            (jd_out.jd_value() - 2_451_545.0).abs() < 1e-9,
            "JD round-trip error: {}",
            jd_out.jd_value() - 2_451_545.0
        );
    }

    // ── advance ───────────────────────────────────────────────────────────────

    #[test]
    fn advance_applies_derivative_correctly() {
        let epoch = JulianDate::new(2_451_545.0).to_time();
        let pos = Position::<GCRS>::new(7000.0, 0.0, 0.0);
        let vel = Velocity::<GCRS>::new(0.0, 7.5, 0.0);
        let s = OrbitState::new(epoch, pos, vel);

        let deriv = StateDerivative::new(
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
            Acceleration::<GCRS>::new(0.0, 0.0, -9.8e-3),
        );
        let dt = Second::new(10.0);
        let s2 = s.advance(&deriv, dt);

        assert!((s2.position.x().value() - 7000.0).abs() < 1e-10);
        assert!((s2.position.y().value() - 75.0).abs() < 1e-10);
        assert!((s2.velocity.z().value() - (-0.098)).abs() < 1e-10);
        // Epoch is unchanged by advance.
        assert_eq!(s2.epoch, epoch);
    }

    #[test]
    fn advance_leaves_epoch_unchanged() {
        let epoch = JulianDate::new(2_451_545.0).to_time();
        let s = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        );
        let deriv = StateDerivative::new(
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
            Acceleration::<GCRS>::new(0.0, 0.0, 0.0),
        );
        let s2 = s.advance(&deriv, Second::new(60.0));
        assert_eq!(s2.epoch, epoch, "advance must not mutate the epoch");
    }

    // ── StateDerivative ───────────────────────────────────────────────────────

    #[test]
    fn rk4_combine_is_weighted_mean() {
        let k = StateDerivative::new(
            Velocity::<GCRS>::new(1.0, 0.0, 0.0),
            Acceleration::<GCRS>::new(0.0, 1.0, 0.0),
        );
        let combined = StateDerivative::rk4_combine(&k, &k, &k, &k);
        // (1 + 2 + 2 + 1) / 6 = 1
        assert!((combined.vel.x().value() - 1.0).abs() < 1e-12);
        assert!((combined.acc.y().value() - 1.0).abs() < 1e-12);
    }

    // ── SpacecraftProperties ──────────────────────────────────────────────────

    #[test]
    fn spacecraft_properties_demo_leo() {
        let p = SpacecraftProperties::demo_leo();
        assert!((p.mass.value() - 500.0).abs() < f64::EPSILON);
        assert!((p.drag_area.value() - 2.0).abs() < f64::EPSILON);
        assert!((p.srp_area.value() - 2.0).abs() < f64::EPSILON);
        assert!((p.cd.value() - 2.2).abs() < f64::EPSILON);
        assert!((p.cr.value() - 1.3).abs() < f64::EPSILON);
        // area-to-mass: 2 m² / 500 kg = 0.004 m²/kg
        assert!((p.drag_area_to_mass().value() - 0.004).abs() < 1e-12);
        assert!((p.srp_area_to_mass().value() - 0.004).abs() < 1e-12);
    }

    // ── Generic instantiation ─────────────────────────────────────────────────

    #[test]
    fn generic_heliocentric_icrs_state_compiles() {
        use crate::coordinates::cartesian;
        use crate::qtty::unit::Kilometer;

        type HelioPos = cartesian::Position<Heliocentric, ICRS, Kilometer>;
        type HelioVel = cartesian::Velocity<ICRS, KmPerSecond>;

        let pos = HelioPos::new(1.496e8, 0.0, 0.0);
        let vel = HelioVel::new(0.0, 29.78, 0.0);
        let s = OrbitState::<Heliocentric, ICRS>::new_at_jd(JulianDate::new(2_451_545.0), pos, vel);
        assert!((s.position.x().value() - 1.496e8).abs() < 1.0);
        assert!((s.velocity.y().value() - 29.78).abs() < 1e-10);
    }

    #[test]
    fn state_derivative_velocity_and_acceleration_accessors() {
        let k = StateDerivative::new(
            Velocity::<GCRS>::new(1.0, 2.0, 3.0),
            Acceleration::<GCRS>::new(4.0, 5.0, 6.0),
        );
        assert!((k.velocity().x().value() - 1.0).abs() < 1e-12);
        assert!((k.velocity().y().value() - 2.0).abs() < 1e-12);
        assert!((k.acceleration().x().value() - 4.0).abs() < 1e-12);
        assert!((k.acceleration().z().value() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn state_derivative_add_and_scaled() {
        let k = StateDerivative::new(
            Velocity::<GCRS>::new(1.0, 0.0, 0.0),
            Acceleration::<GCRS>::new(0.0, 1.0, 0.0),
        );
        let sum = k.add(&k);
        assert!((sum.velocity().x().value() - 2.0).abs() < 1e-12);
        let scaled = k.scaled(0.5);
        assert!((scaled.velocity().x().value() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn velocity_typed_returns_velocity_field() {
        let s = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(1.0, 2.0, 3.0),
        );
        let v = s.velocity_typed();
        assert!((v.x().value() - 1.0).abs() < 1e-12);
        assert!((v.z().value() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn orbit_state_partial_eq_distinguishes_states() {
        let s1 = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        );
        let s2 = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7001.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        );
        assert_ne!(s1, s2);
        assert_eq!(s1, s1);
    }

    #[test]
    fn spacecraft_state_constructs_and_compares() {
        let orbit = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        );
        let props = SpacecraftProperties::demo_leo();
        let sc = SpacecraftState { orbit, properties: props };
        let sc2 = SpacecraftState { orbit, properties: props };
        assert_eq!(sc, sc2);
        assert!((sc.properties.mass.value() - 500.0).abs() < f64::EPSILON);
    }
}
