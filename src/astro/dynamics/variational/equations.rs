// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Variational equations for Cartesian orbital mechanics.
//!
//! ## Mathematical formulation
//!
//! The full propagated state is `(y, Φ)` where:
//!
//! - `y ∈ ℝ⁶` — the Cartesian orbit state `[r; v]`
//! - `Φ ∈ ℝ^{6×6}` — the state-transition matrix (STM), `Φ(t₀) = I`
//!
//! The equations of motion are:
//!
//! ```text
//! dy/dt = f(y)         =  [v ; a(r, v)]
//! dΦ/dt = A(y) · Φ,   Φ(t₀) = I
//! ```
//!
//! where `A = ∂f/∂y` is the 6×6 dynamics Jacobian:
//!
//! ```text
//! A = [  0   I  ]
//!     [ A_r A_v ]
//! ```
//!
//! with:
//! - `A_r = ∂a/∂r` (3×3, s⁻²) — from [`ForcePartials::d_acc_d_pos`]
//! - `A_v = ∂a/∂v` (3×3, s⁻¹) — from [`ForcePartials::d_acc_d_vel`]
//!
//! ## Velocity-partial assumption
//!
//! For conservative forces (two-body, J2, third-body gravity) `A_v = 0`.
//! Atmospheric drag has a non-zero `A_v`; however, the [`ForceModel::partials`]
//! default implementation returns an error for drag, so drag is effectively
//! excluded from the analytic variational system in this release.
//! Any force model that *does* expose analytic velocity partials via
//! [`ForcePartials::d_acc_d_vel`] will have them included automatically.
//!
//! ## Internal representation
//!
//! All numerics in this module use raw `f64` arrays for efficiency.
//! Typed boundaries are provided by the propagator in `super::propagator`.

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::forces::traits::{ForceModel, ForcePartials};
use crate::astro::dynamics::state::OrbitState;
use crate::coordinates::frames::GCRS;

// =============================================================================
// A-matrix construction
// =============================================================================

/// Build the 6×6 dynamics Jacobian `A = ∂f/∂y` from the force partials.
///
/// Layout (row-major, `[r, v]` ordering):
///
/// ```text
/// A = [ 0_{3×3}   I_{3×3} ]
///     [   A_r       A_v   ]
/// ```
///
/// - Rows/cols 0–2: position subspace
/// - Rows/cols 3–5: velocity subspace
pub(super) fn build_a_matrix(partials: &ForcePartials<GCRS>) -> [[f64; 6]; 6] {
    let ar = partials.d_acc_d_pos.as_array();
    let av = partials.d_acc_d_vel.as_array();
    let mut a = [[0.0_f64; 6]; 6];

    // Top-right block: identity (∂v/∂v = I)
    a[0][3] = 1.0;
    a[1][4] = 1.0;
    a[2][5] = 1.0;

    // Bottom-left block: ∂a/∂r
    for i in 0..3 {
        for j in 0..3 {
            a[3 + i][j] = ar[i][j];
        }
    }

    // Bottom-right block: ∂a/∂v
    for i in 0..3 {
        for j in 0..3 {
            a[3 + i][3 + j] = av[i][j];
        }
    }

    a
}

// =============================================================================
// 6×6 matrix operations
// =============================================================================

/// Row-major 6×6 matrix multiplication: `C = A · B`.
#[inline]
pub(super) fn mat6_mul(a: &[[f64; 6]; 6], b: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut c = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for k in 0..6 {
            let aik = a[i][k];
            if aik == 0.0 {
                continue;
            }
            for j in 0..6 {
                c[i][j] += aik * b[k][j];
            }
        }
    }
    c
}

/// Element-wise scale a 6×6 matrix: `B = s · A`.
#[inline]
pub(super) fn mat6_scale(m: &[[f64; 6]; 6], s: f64) -> [[f64; 6]; 6] {
    let mut out = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            out[i][j] = m[i][j] * s;
        }
    }
    out
}

/// Element-wise addition of two 6×6 matrices: `C = A + B`.
#[inline]
pub(super) fn mat6_add(a: &[[f64; 6]; 6], b: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut c = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}

/// 6×6 identity matrix (raw array).
#[inline]
pub(super) fn identity_6x6() -> [[f64; 6]; 6] {
    let mut m = [[0.0_f64; 6]; 6];
    for (i, row) in m.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    m
}

// =============================================================================
// Variational derivative
// =============================================================================

/// Compute the time derivative of the combined `(y, Φ)` variational state.
///
/// Returns `(dy/dt, dΦ/dt)` where:
/// - `dy/dt = [v; a(r,v)]` (raw 6-array, km and km/s)
/// - `dΦ/dt = A(y) · Φ` (raw 6×6 array)
///
/// Errors if the force model does not provide analytic partials.
pub(super) fn variational_derivative<FM>(
    force: &FM,
    state: &OrbitState,
    phi: &[[f64; 6]; 6],
    ctx: &DynamicsContext,
) -> Result<([f64; 6], [[f64; 6]; 6]), DynamicsError>
where
    FM: ForceModel,
{
    let acc = force.acceleration(state, ctx)?;
    let dy = [
        state.velocity.x().value(),
        state.velocity.y().value(),
        state.velocity.z().value(),
        acc.x().value(),
        acc.y().value(),
        acc.z().value(),
    ];

    let partials = force.partials(state, ctx)?;
    let a_mat = build_a_matrix(&partials);
    let dphi = mat6_mul(&a_mat, phi);

    Ok((dy, dphi))
}
