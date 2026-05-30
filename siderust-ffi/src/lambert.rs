// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! C FFI for Lambert's two-point boundary-value problem.
//!
//! Wraps [`keplerian::lambert::lambert`] (Izzo's algorithm) with a flat C ABI.
//! All inputs/outputs use raw `f64` arrays; no Rust generics cross the boundary.
//!
//! ## Lifecycle
//!
//! This module exposes a single stateless function — no allocation or free.
//!
//! 1. Fill `r1_km[3]` and `r2_km[3]` with departure/arrival positions in km.
//! 2. Call [`siderust_lambert_solve`].
//! 3. Read `out_v1_kms[3]` and `out_v2_kms[3]`.
//! 4. Optionally inspect [`SiderustLambertDiagnostics`].
//!
//! ## Status codes
//!
//! | Value | Meaning |
//! |-------|---------|
//! | 0 | Success |
//! | 1 | A required pointer was null |
//! | 9 | Solver did not converge or geometry is degenerate |

use affn::cartesian::Position;
use affn::frames::ICRS;
use keplerian::lambert::{lambert, LambertBranch};
use qtty::length::Kilometer;
use qtty::Second;

use crate::error::SiderustStatus;

// ─── Gravitational parameter ──────────────────────────────────────────────────

use qtty::dynamics::GravitationalParameter;

// ─── Public types ─────────────────────────────────────────────────────────────

/// Branch selector for the Lambert solver.
///
/// | Value | Meaning |
/// |-------|---------|
/// | 0 | Prograde (direct) transfer |
/// | 1 | Retrograde transfer |
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustLambertBranch {
    /// Prograde (co-rotating) single-revolution transfer.
    Prograde = 0,
    /// Retrograde single-revolution transfer.
    Retrograde = 1,
}

/// Householder-iteration diagnostics accompanying a Lambert solution.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustLambertDiagnostics {
    /// Number of Householder iterations performed.
    pub iterations: u32,
    /// Final residual from the Householder solver (dimensionless).
    pub residual: f64,
    /// Number of complete revolutions (0 for single-revolution transfers).
    pub revolutions: u32,
}

// ─── Entry point ──────────────────────────────────────────────────────────────

/// Solve Lambert's single-revolution two-point boundary-value problem.
///
/// All position/velocity arrays follow the standard ordering `[x, y, z]`.
///
/// # Parameters
///
/// * `r1_km`       — departure position, km (3 elements).
/// * `r2_km`       — arrival position, km (3 elements).
/// * `tof_s`       — time of flight, seconds.
/// * `mu_km3_s2`   — gravitational parameter of the central body, km³·s⁻².
/// * `branch`      — 0 = prograde, 1 = retrograde.
/// * `out_v1_kms`  — receives departure velocity, km/s (3 elements, **not** null).
/// * `out_v2_kms`  — receives arrival velocity, km/s (3 elements, **not** null).
/// * `out_diag`    — receives solver diagnostics; may be null.
///
/// # Returns
///
/// [`SiderustStatus::Ok`] on success.
/// [`SiderustStatus::NullPointer`] if any required pointer is null.
/// [`SiderustStatus::InvalidArgument`] if the solver fails to converge or
/// the input geometry is degenerate (e.g. anti-parallel positions, zero ToF).
#[no_mangle]
pub extern "C" fn siderust_lambert_solve(
    r1_km: *const f64,
    r2_km: *const f64,
    tof_s: f64,
    mu_km3_s2: f64,
    branch: i32,
    out_v1_kms: *mut f64,
    out_v2_kms: *mut f64,
    out_diag: *mut SiderustLambertDiagnostics,
) -> SiderustStatus {
    crate::ffi_guard!({
        if r1_km.is_null() || r2_km.is_null() {
            return SiderustStatus::NullPointer;
        }
        crate::check_out!(out_v1_kms);
        crate::check_out!(out_v2_kms);

        // SAFETY: caller guarantees r1_km/r2_km point to at least 3 f64 values.
        let r1_arr = unsafe { std::slice::from_raw_parts(r1_km, 3) };
        let r2_arr = unsafe { std::slice::from_raw_parts(r2_km, 3) };

        let r1 = Position::<(), ICRS, Kilometer>::new(r1_arr[0], r1_arr[1], r1_arr[2]);
        let r2 = Position::<(), ICRS, Kilometer>::new(r2_arr[0], r2_arr[1], r2_arr[2]);

        let tof = Second::new(tof_s);
        let mu = GravitationalParameter::new(mu_km3_s2);
        let lb = if branch == 1 {
            LambertBranch::Retrograde
        } else {
            LambertBranch::Prograde
        };

        match lambert(r1, r2, tof, mu, lb) {
            Ok(sol) => {
                // SAFETY: out_v1_kms / out_v2_kms checked non-null above.
                unsafe {
                    *out_v1_kms.add(0) = sol.v1.x().value();
                    *out_v1_kms.add(1) = sol.v1.y().value();
                    *out_v1_kms.add(2) = sol.v1.z().value();
                    *out_v2_kms.add(0) = sol.v2.x().value();
                    *out_v2_kms.add(1) = sol.v2.y().value();
                    *out_v2_kms.add(2) = sol.v2.z().value();
                    if !out_diag.is_null() {
                        *out_diag = SiderustLambertDiagnostics {
                            iterations: sol.diagnostics.iterations,
                            residual: sol.diagnostics.residual,
                            revolutions: sol.diagnostics.revolutions,
                        };
                    }
                }
                SiderustStatus::Ok
            }
            Err(_) => SiderustStatus::InvalidArgument,
        }
    })
}
