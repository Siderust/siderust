// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Core trait definitions and shared Jacobian building blocks for force models.
//!
//! ## Scope
//!
//! This file owns:
//!
//! - [`ForcePartials<F>`] — the frame-tagged `(A_r, A_v)` Jacobian block.
//! - [`ForceModel<C, F>`] — the central force-model trait.
//! - Constants shared by multiple force models (Earth radius, rotation rate, etc.).
//!
//! ## Equations
//!
//! For a conservative central force `a = a(r)` the variational equation is
//!
//! ```text
//! δ̈r = A_r · δr,   A_r = ∂a/∂r,   A_v = 0
//! ```
//!
//! ## Units
//!
//! All accelerations are in km/s².  `A_r` has units s⁻², `A_v` has units s⁻¹.
//!
//! ## Frame/center assumptions
//!
//! The default frame is [`GCRS`] and the default center is [`Geocentric`].
//!
//! ## References
//!
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, §1, §8.
//! * Montenbruck & Gill, *Satellite Orbits*, §3.2, §3.4.

use affn::matrix3::FrameMatrix3;

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::centers::{Geocentric, ReferenceCenter};
use crate::coordinates::frames::{ReferenceFrame, GCRS};
use crate::qtty::{AstronomicalUnit, InverseSeconds, Kilometer, Kilometers, Pascals, Unit};

// Re-export so force-model subfiles can use them without extra imports.
pub use crate::astro::dynamics::units::{GravitationalParameter, GM_EARTH, GM_MOON, GM_SUN};

// =============================================================================
// Earth constants (used by multiple force models)
// =============================================================================

/// Earth mean equatorial radius (GRS-80 / WGS-84), km.
pub const R_EARTH: Kilometers = Kilometers::new(6_378.137);

/// Earth mean rotation rate (sidereal), rad/s (IAU 2000).
///
/// Typed as an angular rate (`InverseSeconds` ≡ `Quantity<Per<Ratio, Second>>`).
/// Use `.value()` when you need the raw `f64`.
pub const OMEGA_EARTH_RAD_S: InverseSeconds = InverseSeconds::new(7.292_115_146_706_979e-5);

/// Minimum geocentric radius below which force-model computations are
/// considered degenerate (satellite inside Earth).  100 km is well below any
/// physical orbit but avoids false positives for very low-perigee trajectories.
pub const DEGENERATE_RADIUS_KM: f64 = 100.0;

/// Solar radiation pressure at 1 AU, N/m².
pub const P0: Pascals = Pascals::new(4.560e-6);

/// Standard gravitational parameter of the Sun, km³/s².
pub const MU_SUN: GravitationalParameter = GM_SUN;
/// Standard gravitational parameter of the Moon, km³/s².
pub const MU_MOON: GravitationalParameter = GM_MOON;

/// Astronomical unit in km (derived from unit ratios, not a magic number).
pub(super) const AU_IN_KM: f64 = AstronomicalUnit::RATIO / Kilometer::RATIO;

// =============================================================================
// ForcePartials
// =============================================================================

/// Frame-tagged Jacobian blocks of the acceleration: `∂a/∂[r, v]`.
///
/// These are the `A_r` and `A_v` matrices that form the lower half of the
/// variational dynamics matrix `F(t)`:
///
/// ```text
/// F(t) = [ 0   I  ]
///        [ A_r A_v ]
/// ```
///
/// where:
/// * `A_r = ∂a/∂r` (units: km/s² per km = s⁻²)
/// * `A_v = ∂a/∂v` (units: km/s² per km/s = s⁻¹)
///
/// Both matrices are stored as raw `f64` values in the tagged frame `F`.
/// For conservative forces such as gravity, `A_v = 0`.
#[derive(Debug, Clone, Copy)]
pub struct ForcePartials<F = GCRS> {
    /// `∂a/∂r` in frame `F` (units: s⁻²).
    pub d_acc_d_pos: FrameMatrix3<F>,
    /// `∂a/∂v` in frame `F` (units: s⁻¹). Zero for conservative forces.
    pub d_acc_d_vel: FrameMatrix3<F>,
}

impl<F> ForcePartials<F> {
    /// Zero partials — both Jacobian blocks are zero matrices.
    ///
    /// Used as a neutral element when compositing forces whose partials are
    /// not implemented.
    pub fn zero() -> Self {
        Self {
            d_acc_d_pos: FrameMatrix3::zero(),
            d_acc_d_vel: FrameMatrix3::zero(),
        }
    }

    /// Analytic `∂a/∂r` for the Newtonian two-body acceleration `a = -μ r / |r|³`.
    ///
    /// The formula is:
    ///
    /// ```text
    /// A_r = -μ/r³ I + 3μ/r⁵ r rᵀ
    /// ```
    ///
    /// `A_v = 0` because central gravity has no velocity dependence.
    ///
    /// # Arguments
    ///
    /// * `mu`  — standard gravitational parameter μ = GM (km³/s²).
    /// * `r`   — Cartesian position vector `[x, y, z]` (km) in frame `F`.
    pub fn two_body(mu: GravitationalParameter, r: [f64; 3]) -> Self {
        let r_norm = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
        if r_norm == 0.0 {
            return Self::zero();
        }
        let r3 = r_norm.powi(3);
        let r5 = r_norm.powi(5);
        let mu_val = mu.value();
        let diag = -mu_val / r3;
        let scale = 3.0 * mu_val / r5;
        let mut d_acc_d_pos = FrameMatrix3::from_diagonal([diag; 3]);
        d_acc_d_pos.add_outer_product_in_place(r, [scale * r[0], scale * r[1], scale * r[2]]);
        Self {
            d_acc_d_pos,
            d_acc_d_vel: FrameMatrix3::zero(),
        }
    }

    /// Element-wise sum of `self` and `other`, returning a new [`ForcePartials`].
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        let mut d_acc_d_pos = FrameMatrix3::from_array(*self.d_acc_d_pos.as_array());
        d_acc_d_pos.add_in_place(&other.d_acc_d_pos);
        let mut d_acc_d_vel = FrameMatrix3::from_array(*self.d_acc_d_vel.as_array());
        d_acc_d_vel.add_in_place(&other.d_acc_d_vel);
        Self {
            d_acc_d_pos,
            d_acc_d_vel,
        }
    }

    /// Element-wise add `other` into `self` in place.
    pub fn add_in_place(&mut self, other: &Self) {
        self.d_acc_d_pos.add_in_place(&other.d_acc_d_pos);
        self.d_acc_d_vel.add_in_place(&other.d_acc_d_vel);
    }
}

// =============================================================================
// ForceModel trait
// =============================================================================

/// A force model evaluated on an inertial [`OrbitState<C, F>`].
///
/// All providers (ephemeris, atmosphere, gravity field) are accessed via
/// `ctx`.  Implementors store only their tunable physical parameters.
///
/// # Type parameters
///
/// * `C` — reference center (default [`Geocentric`]).
/// * `F` — reference frame (default [`GCRS`]).
///
/// # Units
///
/// The returned acceleration is in km/s² in frame `F`.
///
/// # Errors
///
/// Returns a [`DynamicsError`] when:
/// * a required provider is absent from `ctx` (e.g. ephemeris for third-body),
/// * the spacecraft is below the surface,
/// * geometry is degenerate.
///
/// # Partial derivatives
///
/// The default `partials` implementation returns
/// `Err(DynamicsError::Provider(_))` ("analytic partials not implemented").
/// Override it in force models that have a closed-form Jacobian.
///
/// The `partials` convention is the lower half of `F(t) = [[0, I], [A_r, A_v]]`:
///
/// * `A_r = ∂a/∂r` (3×3, units: s⁻²)
/// * `A_v = ∂a/∂v` (3×3, units: s⁻¹; zero for conservative forces)
pub trait ForceModel<C = Geocentric, F = GCRS>: Send + Sync
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    /// Compute the inertial acceleration acting on `state`.
    ///
    /// # Frame
    ///
    /// The result is expressed in frame `F` with units km/s².
    fn acceleration(
        &self,
        state: &OrbitState<C, F>,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<F, AccelerationUnit>, DynamicsError>;

    /// Analytic partial derivatives of the acceleration with respect to
    /// position and velocity: `(∂a/∂r, ∂a/∂v)` in frame `F`.
    ///
    /// The default implementation returns an error indicating that analytic
    /// partials are not available for this model.  Override this method in
    /// force models that have a closed-form Jacobian.
    fn partials(
        &self,
        _state: &OrbitState<C, F>,
        _ctx: &DynamicsContext,
    ) -> Result<ForcePartials<F>, DynamicsError> {
        Err(DynamicsError::Provider(Box::new(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "analytic partials not implemented for this force model",
        ))))
    }

    /// Human-readable name identifying this force model.
    ///
    /// The default returns the fully-qualified Rust type name.
    fn name(&self) -> &'static str {
        std::any::type_name::<Self>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::time::JulianDate;

    fn sample_state() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn force_partials_zero_has_zero_matrices() {
        let z = ForcePartials::<GCRS>::zero();
        let a = z.d_acc_d_pos.as_array();
        let v = z.d_acc_d_vel.as_array();
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(a[i][j], 0.0);
                assert_eq!(v[i][j], 0.0);
            }
        }
    }

    #[test]
    fn force_partials_two_body_diagonal_sign() {
        let mu = GM_EARTH;
        let r = [7000.0_f64, 0.0, 0.0];
        let p = ForcePartials::<GCRS>::two_body(mu, r);
        let a = p.d_acc_d_pos.as_array();
        assert!(
            (a[0][1]).abs() < 1e-30,
            "off-diagonal [0][1] should be zero"
        );
        let v = p.d_acc_d_vel.as_array();
        for row in v {
            for &val in row {
                assert_eq!(val, 0.0, "d_acc_d_vel must be zero for two-body");
            }
        }
    }

    #[test]
    fn force_partials_two_body_zero_position_returns_zero() {
        let p = ForcePartials::<GCRS>::two_body(GM_EARTH, [0.0, 0.0, 0.0]);
        let a = p.d_acc_d_pos.as_array();
        for row in a {
            for &val in row {
                assert_eq!(val, 0.0);
            }
        }
    }

    #[test]
    fn force_partials_add_sums_elements() {
        let a = ForcePartials::<GCRS>::two_body(GM_EARTH, [7000.0, 0.0, 0.0]);
        let b = ForcePartials::<GCRS>::two_body(GM_EARTH, [7000.0, 0.0, 0.0]);
        let c = a.add(&b);
        let ca = c.d_acc_d_pos.as_array();
        let aa = a.d_acc_d_pos.as_array();
        for i in 0..3 {
            for j in 0..3 {
                let rel = (ca[i][j] - 2.0 * aa[i][j]).abs();
                assert!(
                    rel < 1e-20,
                    "add must sum element-wise: [{i}][{j}] rel={rel}"
                );
            }
        }
    }

    #[test]
    fn force_partials_add_in_place_sums_elements() {
        let mut a = ForcePartials::<GCRS>::two_body(GM_EARTH, [7000.0, 0.0, 0.0]);
        let b = ForcePartials::<GCRS>::two_body(GM_EARTH, [7000.0, 0.0, 0.0]);
        let orig = *a.d_acc_d_pos.as_array();
        a.add_in_place(&b);
        let ca = a.d_acc_d_pos.as_array();
        for i in 0..3 {
            for j in 0..3 {
                let rel = (ca[i][j] - 2.0 * orig[i][j]).abs();
                assert!(
                    rel < 1e-20,
                    "add_in_place must sum in place: [{i}][{j}] rel={rel}"
                );
            }
        }
    }

    struct ConstantForce;
    impl ForceModel for ConstantForce {
        fn acceleration(
            &self,
            _state: &OrbitState,
            _ctx: &DynamicsContext,
        ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
            Ok(Acceleration::<GCRS, AccelerationUnit>::new(1e-6, 0.0, 0.0))
        }
    }

    #[test]
    fn default_partials_returns_error() {
        let f = ConstantForce;
        let ctx = DynamicsContext::empty();
        let result = f.partials(&sample_state(), &ctx);
        assert!(result.is_err(), "default partials must return an error");
    }

    #[test]
    fn default_name_returns_type_name() {
        let f = ConstantForce;
        let name = f.name();
        assert!(!name.is_empty(), "type name must not be empty");
    }
}
