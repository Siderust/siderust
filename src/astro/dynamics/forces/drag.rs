// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Cannonball atmospheric drag — [`DragForce`] force model.
//!
//! ## Scope
//!
//! Provides a cannonball aerodynamic drag model that reads atmospheric mass
//! density from the [`DynamicsContext`]'s atmosphere provider.  Any
//! [`super::super::atmosphere::DensityProvider`] implementation can be injected at context-build time.
//!
//! ## Equations
//!
//! ```text
//! a_drag = − ½ · Cd · (A/m) · ρ(h) · |v_rel|_SI · v_rel
//! ```
//!
//! where:
//!
//! * `Cd` — drag coefficient (dimensionless);
//! * `A/m` — area-to-mass ratio (m²/kg);
//! * `ρ(h)` — atmospheric mass density at geodetic altitude `h` (kg/m³),
//!   read from the context's [`super::super::atmosphere::DensityProvider`];
//! * `v_rel = v − ω_⊕ × r` — velocity relative to the co-rotating atmosphere
//!   (km/s), subtracted element-wise in raw f64 after extracting typed values;
//! * `|v_rel|_SI` — magnitude of `v_rel` converted to m/s for unit consistency.
//!
//! ## Unit arithmetic
//!
//! ```text
//! a [km/s²] = −½ · Cd · (A/m)[m²/kg] · ρ[kg/m³] · |v_rel|[m/s] · v_rel[km/s]
//!           = [m²/kg · kg/m³ · m/s · km/s]
//!           = [1/m · m/s · km/s]
//!           = [km/s²]  ✓
//! ```
//!
//! ## Frame/center assumptions
//!
//! All vectors are in `<Geocentric, GCRS>`.  The Earth rotation is modelled
//! as a fixed sidereal rate (default [`OMEGA_EARTH_RAD_S`]).  Geodetic
//! altitude uses a placeholder identity GCRS↔ITRF rotation (error ≤ 20 km at
//! the poles — well within the model uncertainty).
//!
//! ## References
//!
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, §8.
//! * Montenbruck & Gill, *Satellite Orbits*, §3.4.

use affn::cartesian::Vector;

use crate::astro::dynamics::atmosphere::geodetic_altitude;
use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::GCRS;
use crate::ext_qtty::dynamics::InverseSecond;
use crate::qtty::unit::Kilometer;
use crate::qtty::{AreaToMass, DragCoefficient, InverseSeconds, Kilometers, Ratios};

use super::traits::{ForceModel, OMEGA_EARTH_RAD_S};

/// WGS-84 equatorial radius in km.
const WGS84_A_KM: f64 = 6_378.137;
/// WGS-84 flattening.
const WGS84_F: f64 = 1.0 / 298.257_223_563;

// =============================================================================
// DragForce
// =============================================================================

/// Cannonball atmospheric drag acceleration.
///
/// ```text
/// a_drag = − ½ · Cd · (A/m) · ρ(h) · |v_rel|_SI · v_rel
/// ```
///
/// The atmosphere density model is **not embedded** in this struct; it is read
/// at propagation time from [`DynamicsContext::require_atmosphere()`].  Inject
/// a provider via [`DynamicsContext::builder()`] + `with_atmosphere(…)`.
///
/// ## Constructor
///
/// ```rust
/// use siderust::astro::dynamics::forces::DragForce;
/// use siderust::qtty::{AreaToMass, DragCoefficient, InverseSeconds};
///
/// // Default Earth rotation (7.292115e-5 rad/s about +Z):
/// let drag = DragForce::new(DragCoefficient::new(2.2), AreaToMass::new(0.01));
///
/// // Custom omega_earth (e.g. for tests or other bodies):
/// let drag_custom = DragForce::with_omega(
///     DragCoefficient::new(2.2),
///     AreaToMass::new(0.01),
///     InverseSeconds::new(7.292115e-5),
/// );
/// ```
#[derive(Debug, Clone)]
pub struct DragForce {
    /// Drag coefficient `C_D` (dimensionless, typical LEO ≈ 2.2).
    pub cd: DragCoefficient,
    /// Effective area-to-mass ratio `A/m` (m²/kg).
    pub area_to_mass: AreaToMass,
    /// Earth angular velocity vector `ω_⊕` (rad/s) in GCRS.
    ///
    /// Default: `[0, 0, 7.292115e-5]` rad/s (IAU sidereal rate).
    pub omega_earth: Vector<GCRS, InverseSecond>,
}

impl DragForce {
    /// Construct a [`DragForce`] using the IAU Earth rotation rate about +Z.
    ///
    /// # Arguments
    ///
    /// * `cd`           — Drag coefficient (e.g. `DragCoefficient::new(2.2)`).
    /// * `area_to_mass` — Area-to-mass ratio (m²/kg, e.g. `AreaToMass::new(0.01)`).
    pub fn new(cd: DragCoefficient, area_to_mass: AreaToMass) -> Self {
        Self::with_omega(cd, area_to_mass, OMEGA_EARTH_RAD_S)
    }

    /// Construct a [`DragForce`] with a custom Earth angular-velocity magnitude.
    ///
    /// The rotation axis is assumed to be aligned with +Z (GCRS).
    ///
    /// # Arguments
    ///
    /// * `cd` — Drag coefficient.
    /// * `area_to_mass` — Area-to-mass ratio (m²/kg).
    /// * `omega_z` — Angular velocity magnitude about +Z (rad/s), typed as
    ///   [`InverseSeconds`].  Use `InverseSeconds::new(0.0)` to disable the
    ///   rotating-atmosphere correction.
    pub fn with_omega(
        cd: DragCoefficient,
        area_to_mass: AreaToMass,
        omega_z: InverseSeconds,
    ) -> Self {
        Self {
            cd,
            area_to_mass,
            omega_earth: Vector::<GCRS, InverseSecond>::new(0.0_f64, 0.0_f64, omega_z.value()),
        }
    }
}

impl ForceModel for DragForce {
    /// Compute the aerodynamic drag acceleration (km/s²).
    ///
    /// # Errors
    ///
    /// * [`DynamicsError::Provider`] — context has no atmosphere provider.
    /// * [`DynamicsError::AltitudeBelowSurface`] — geocentric altitude < 0.
    /// * [`DynamicsError::AtmosphereProviderError`] — density model error.
    #[inline]
    fn acceleration(
        &self,
        s: &OrbitState,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        // --- get density provider from context ---
        let provider = ctx.require_atmosphere()?;

        // --- geodetic altitude ---
        let pos_typed: &affn::cartesian::Position<Geocentric, GCRS, Kilometer> = &s.position;
        let h = geodetic_altitude(pos_typed, Kilometers::new(WGS84_A_KM), Ratios::new(WGS84_F));

        // --- guard: below surface ---
        if h.value() < 0.0 {
            return Err(DynamicsError::AltitudeBelowSurface {
                altitude_km: h.value(),
            });
        }

        // --- density ---
        let rho = provider.density(h)?.value();

        // --- rotating-atmosphere correction: ω × r ---
        let rx = s.position.x().value();
        let ry = s.position.y().value();
        let oz = self.omega_earth.z().value(); // ω_z (rad/s)

        // ω × r = (ω_z · ry directed -x, ω_z · rx directed +y, 0)  [km/s]
        //   (since ω = [0, 0, ω_z])
        let omega_cross_r_x = -oz * ry;
        let omega_cross_r_y = oz * rx;

        // --- relative velocity (km/s) ---
        let vx = s.velocity.x().value() - omega_cross_r_x;
        let vy = s.velocity.y().value() - omega_cross_r_y;
        let vz = s.velocity.z().value(); // ω × r has no z component

        // --- drag acceleration ---
        // |v_rel| in m/s for (A/m)[m²/kg]·ρ[kg/m³]·|v_rel|[m/s]·v_rel[km/s] → [km/s²]
        let v_mag_m_s = (vx * vx + vy * vy + vz * vz).sqrt() * 1_000.0;
        let pre = -0.5 * self.cd.value() * self.area_to_mass.value() * rho * v_mag_m_s;

        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            pre * vx,
            pre * vy,
            pre * vz,
        ))
    }
}

/// Type alias: drag model with context-injected atmosphere (same as [`DragForce`]).
///
/// Retained for API compatibility; previously referred to
/// `DragForce<ExponentialAtmosphere>`.
pub type ExponentialDrag = DragForce;

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::astro::dynamics::atmosphere::{
        ConstantDensity, DensityProvider, ExponentialAtmosphere, Nrlmsise00LiteApprox,
    };
    use crate::astro::dynamics::context::{DynamicsContext, DynamicsContextBuilder};
    use crate::astro::dynamics::integrators::rk4_propagate;
    use crate::astro::dynamics::state::OrbitState;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::qtty::{
        AreaToMass, DragCoefficient, InverseSeconds, KilogramsPerCubicMeter, Second,
    };
    use crate::time::JulianDate;

    use super::super::composite::CompositeForce;
    use super::super::traits::R_EARTH;
    use super::super::two_body::TwoBody;

    fn ctx_with_exponential() -> DynamicsContext {
        let atm: Arc<dyn DensityProvider + Send + Sync> =
            Arc::new(ExponentialAtmosphere::LEO_500KM);
        DynamicsContextBuilder::new().with_atmosphere(atm).build()
    }

    fn ctx_with_constant(rho: f64) -> DynamicsContext {
        let atm: Arc<dyn DensityProvider + Send + Sync> = Arc::new(ConstantDensity {
            rho: KilogramsPerCubicMeter::new(rho),
        });
        DynamicsContextBuilder::new().with_atmosphere(atm).build()
    }

    // ── no provider → error ───────────────────────────────────────────────────

    #[test]
    fn drag_error_when_no_provider() {
        let drag = DragForce::new(DragCoefficient::new(2.2), AreaToMass::new(0.01));
        let ctx = DynamicsContext::empty();
        let r0 = R_EARTH.value() + 400.0;
        let s = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.7, 0.0),
        );
        let result = drag.acceleration(&s, &ctx);
        assert!(
            result.is_err(),
            "expected error with no atmosphere provider"
        );
    }

    // ── drag direction opposes velocity ───────────────────────────────────────

    #[test]
    fn drag_direction_opposes_velocity() {
        let drag = DragForce::new(DragCoefficient::new(2.2), AreaToMass::new(0.01));
        let ctx = ctx_with_exponential();
        let r0 = R_EARTH.value() + 400.0;
        // prograde velocity along +Y
        let s = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.7, 0.0),
        );
        let acc = drag.acceleration(&s, &ctx).unwrap();
        // drag opposes prograde motion → acceleration along -Y
        assert!(
            acc.y().value() < 0.0,
            "drag should oppose +Y velocity, got ay = {:.3e}",
            acc.y().value()
        );
        // x and z components should be tiny for near-prograde orbit
        assert!(acc.x().value().abs() < acc.y().value().abs().abs());
    }

    // ── below surface → error ─────────────────────────────────────────────────

    #[test]
    fn drag_below_surface_returns_error() {
        let drag = DragForce::new(DragCoefficient::new(2.2), AreaToMass::new(0.01));
        let ctx = ctx_with_constant(1.0e-9);
        // Position well inside Earth
        let r_underground = R_EARTH.value() - 100.0;
        let s = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r_underground, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.0, 0.0),
        );
        let result = drag.acceleration(&s, &ctx);
        assert!(
            matches!(result, Err(DynamicsError::AltitudeBelowSurface { .. })),
            "expected AltitudeBelowSurface, got {result:?}"
        );
    }

    // ── rotating atmosphere correction reduces drag slightly ──────────────────

    #[test]
    fn rotating_atmosphere_changes_drag_magnitude() {
        let r0 = R_EARTH.value() + 400.0;
        let v_circ = 7.7_f64; // km/s prograde +Y

        // No rotation
        let drag_norot = DragForce::with_omega(
            DragCoefficient::new(2.2),
            AreaToMass::new(0.01),
            InverseSeconds::new(0.0),
        );
        // With Earth rotation
        let drag_rot = DragForce::new(DragCoefficient::new(2.2), AreaToMass::new(0.01));

        let ctx = ctx_with_constant(1.0e-11);
        let s = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v_circ, 0.0),
        );

        let acc_norot = drag_norot.acceleration(&s, &ctx).unwrap();
        let acc_rot = drag_rot.acceleration(&s, &ctx).unwrap();

        let mag_norot = (acc_norot.x().value().powi(2)
            + acc_norot.y().value().powi(2)
            + acc_norot.z().value().powi(2))
        .sqrt();
        let mag_rot = (acc_rot.x().value().powi(2)
            + acc_rot.y().value().powi(2)
            + acc_rot.z().value().powi(2))
        .sqrt();

        // At 400 km LEO, ω × r ≈ 0.46 km/s eastward.
        // For prograde orbit (+Y), r on +X: ω × r = +Y component → v_rel smaller.
        // So rotating atmosphere reduces the drag magnitude.
        assert!(
            mag_rot < mag_norot,
            "mag_rot={mag_rot:.3e} should be less than mag_norot={mag_norot:.3e}"
        );

        let frac_diff = (mag_norot - mag_rot) / mag_norot;
        assert!(
            frac_diff > 0.01 && frac_diff < 0.20,
            "fractional drag difference {frac_diff:.3} should be ~6 % at 400 km LEO"
        );
    }

    // ── orbit decay under drag ────────────────────────────────────────────────

    #[test]
    fn drag_decays_orbit_altitude() {
        let mu: f64 = 398_600.441_8;
        let r0 = R_EARTH.value() + 350.0;
        let v0 = (mu / r0).sqrt();
        let s0 = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v0, 0.0),
        );

        // Inject a dense constant atmosphere to ensure measurable decay
        let atm: Arc<dyn DensityProvider + Send + Sync> = Arc::new(ConstantDensity {
            rho: KilogramsPerCubicMeter::new(1.0e-11),
        });
        let ctx = DynamicsContextBuilder::new().with_atmosphere(atm).build();

        let force = CompositeForce::empty()
            .push(Box::new(TwoBody::earth()))
            .push(Box::new(DragForce::with_omega(
                DragCoefficient::new(2.2),
                AreaToMass::new(5.0),
                OMEGA_EARTH_RAD_S,
            )));

        let s_end = rk4_propagate(&force, s0, Second::new(30.0), 360, &ctx).unwrap();
        let r_end = s_end.position.distance().value();
        assert!(
            r_end < r0,
            "expected drag-driven decay; r0={r0:.3}, r_end={r_end:.3}"
        );
    }

    // ── NRLMSISE via context ──────────────────────────────────────────────────

    #[test]
    fn drag_with_nrlmsise_provider() {
        let atm: Arc<dyn DensityProvider + Send + Sync> = Arc::new(Nrlmsise00LiteApprox);
        let ctx = DynamicsContextBuilder::new().with_atmosphere(atm).build();

        let drag = DragForce::new(DragCoefficient::new(2.2), AreaToMass::new(0.01));
        let r0 = R_EARTH.value() + 400.0;
        let s = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.7, 0.0),
        );

        let acc = drag.acceleration(&s, &ctx).unwrap();
        // drag should be non-trivial negative Y acceleration
        assert!(
            acc.y().value() < -1e-9,
            "expected meaningful drag at 400 km, got ay = {:.3e}",
            acc.y().value()
        );
    }
}
