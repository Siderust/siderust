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
use crate::qtty::{AstronomicalUnit, Kilometer, Kilometers, Pascals, Unit};

// Re-export so force-model subfiles can use them without extra imports.
pub use crate::astro::dynamics::units::{GravitationalParameter, GM_EARTH, GM_MOON, GM_SUN};

// =============================================================================
// Earth constants (used by multiple force models)
// =============================================================================

/// Earth mean equatorial radius (GRS-80 / WGS-84), km.
pub const R_EARTH: Kilometers = Kilometers::new(6_378.137);

/// Earth mean rotation rate (sidereal), rad/s (IAU 2000).
pub const OMEGA_EARTH_RAD_S: f64 = 7.292_115_146_706_979e-5;

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
        let mut data = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                data[i][j] = scale * r[i] * r[j] + if i == j { diag } else { 0.0 };
            }
        }
        Self {
            d_acc_d_pos: FrameMatrix3::from_array(data),
            d_acc_d_vel: FrameMatrix3::zero(),
        }
    }

    /// Element-wise sum of `self` and `other`, returning a new [`ForcePartials`].
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        let a_pos = self.d_acc_d_pos.as_array();
        let b_pos = other.d_acc_d_pos.as_array();
        let mut out_pos = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                out_pos[i][j] = a_pos[i][j] + b_pos[i][j];
            }
        }
        let a_vel = self.d_acc_d_vel.as_array();
        let b_vel = other.d_acc_d_vel.as_array();
        let mut out_vel = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                out_vel[i][j] = a_vel[i][j] + b_vel[i][j];
            }
        }
        Self {
            d_acc_d_pos: FrameMatrix3::from_array(out_pos),
            d_acc_d_vel: FrameMatrix3::from_array(out_vel),
        }
    }

    /// Element-wise add `other` into `self` in place.
    pub fn add_in_place(&mut self, other: &Self) {
        let b_pos = *other.d_acc_d_pos.as_array();
        let cur_pos = *self.d_acc_d_pos.as_array();
        let mut out_pos = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                out_pos[i][j] = cur_pos[i][j] + b_pos[i][j];
            }
        }
        self.d_acc_d_pos = FrameMatrix3::from_array(out_pos);

        let b_vel = *other.d_acc_d_vel.as_array();
        let cur_vel = *self.d_acc_d_vel.as_array();
        let mut out_vel = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                out_vel[i][j] = cur_vel[i][j] + b_vel[i][j];
            }
        }
        self.d_acc_d_vel = FrameMatrix3::from_array(out_vel);
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
