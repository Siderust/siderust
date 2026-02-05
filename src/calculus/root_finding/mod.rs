// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Root-Finding Algorithms
//!
//! This module provides numerical methods for finding roots of scalar functions,
//! primarily used for refining threshold crossings in astronomical calculations.
//!
//! ## Available Methods
//!
//! | Method     | Type              | Convergence | Function Evals/Iter | Best For |
//! |------------|-------------------|-------------|---------------------|----------|
//! | Newton     | Derivative-based  | Quadratic   | ~3 (finite-diff)    | Smooth functions with good initial guess |
//! | Bisection  | Bracketing        | Linear      | 1                   | Guaranteed convergence, discontinuous functions |
//! | Brent      | Bracketing+Hybrid | Superlinear | 1-2                 | **Best general-purpose**: robust + fast |
//!
//! ## Recommended Usage
//!
//! - **[`brent::refine_root`]**: Preferred for most astronomical applications.
//!   Combines the robustness of bisection with superlinear convergence. Requires
//!   a valid bracket (sign change) but typically converges in fewer function
//!   evaluations than Newton with finite differences.
//!
//! - **[`newton::refine_root`]**: Use when you have a good initial guess and the
//!   function is smooth. Note: uses 3 function evaluations per iteration due to
//!   finite-difference derivative estimation.
//!
//! - **[`bisection::refine_root`]**: Fallback for pathological cases. Guaranteed
//!   to converge but slower than Brent.
//!
//! ## Module Structure
//!
//! Each method is implemented in its own submodule with dedicated constants and tests:
//! - [`newton`] — Newton-Raphson with finite-difference derivatives
//! - [`bisection`] — Classic bisection method
//! - [`brent`] — Brent's method (hybrid: bisection + secant + inverse quadratic)

use crate::astro::JulianDate;
use qtty::{Days, Seconds};

pub mod bisection;
pub mod brent;
pub mod newton;

// =============================================================================
// Public API Constants
// =============================================================================

/// Standard time step for finite-difference derivative estimation (1 second).
/// 
/// This is exposed for external modules that need to sample functions at a
/// consistent time scale for derivative-like operations.
pub const DERIVATIVE_STEP: Days = Seconds::new(1.0).to::<qtty::Day>();

// =============================================================================
// Re-exports for Convenience
// =============================================================================

pub use bisection::refine_root as refine_root_bisection;
pub use brent::refine_root as refine_root_brent;
pub use newton::refine_root as refine_root_newton;

// =============================================================================
// High-Level Crossing Finders
// =============================================================================

/// Finds a threshold crossing between two JD values.
///
/// Tries Newton-Raphson first, falls back to bisection if Newton fails.
pub fn find_crossing<F>(
    jd_a: JulianDate,
    jd_b: JulianDate,
    scalar_fn: &F,
    threshold: f64,
) -> Option<JulianDate>
where
    F: Fn(JulianDate) -> f64,
{
    let guess = JulianDate::new((jd_a.value() + jd_b.value()) * 0.5);

    newton::refine_root(guess, scalar_fn, threshold)
        .or_else(|| bisection::refine_root(jd_a, jd_b, scalar_fn, threshold))
}

/// Finds a threshold crossing using Brent's method (derivative-free).
///
/// This is the **recommended method** when you already have a bracketed interval
/// (i.e., from scanning for sign changes). It provides:
/// - Robust convergence (guaranteed when bracket is valid)
/// - ~1-2 function evaluations per iteration (vs ~3 for Newton with finite-diff)
/// - Superlinear convergence in practice
///
/// # Arguments
/// - `jd_a`, `jd_b`: Bracket endpoints (must have sign change in `scalar_fn - threshold`)
/// - `scalar_fn`: Function to find root of
/// - `threshold`: Target value (finds where `scalar_fn(x) == threshold`)
///
/// # Performance
/// For expensive functions like VSOP87 ephemeris calculations, Brent typically
/// reduces total evaluations by 50-70% compared to Newton+finite-difference.
pub fn find_crossing_brent<F>(
    jd_a: JulianDate,
    jd_b: JulianDate,
    scalar_fn: &F,
    threshold: f64,
) -> Option<JulianDate>
where
    F: Fn(JulianDate) -> f64,
{
    brent::refine_root(jd_a, jd_b, scalar_fn, threshold)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn find_crossing_prefers_newton_and_returns_root() {
        let threshold = 200.0;
        let lo = JulianDate::new(threshold - 2.0);
        let hi = JulianDate::new(threshold + 2.0);
        let scalar = |jd: JulianDate| jd.value();

        let found = find_crossing(lo, hi, &scalar, threshold);
        assert!((found.unwrap().value() - threshold).abs() < 1e-9);
    }

    #[test]
    fn find_crossing_falls_back_to_bisection_on_spiky_function() {
        let lo = JulianDate::new(-1.0);
        let hi = JulianDate::new(1.0);
        let scalar = |jd: JulianDate| {
            if jd.value() < 0.0 {
                -1.0
            } else {
                1.0
            }
        };

        let found = find_crossing(lo, hi, &scalar, 0.0).expect("expected fallback to succeed");
        assert!(found.value().abs() < 1e-6);
    }

    #[test]
    fn find_crossing_brent_finds_sine_root() {
        let lo = JulianDate::new(3.0);
        let hi = JulianDate::new(4.0);
        let scalar = |jd: JulianDate| jd.value().sin();

        let result = find_crossing_brent(lo, hi, &scalar, 0.0).expect("should find pi");
        assert!((result.value() - std::f64::consts::PI).abs() < 1e-9);
    }
}
