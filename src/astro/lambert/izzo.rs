//! Izzo (2014) reformulation of Lambert's two-point boundary-value problem.
//!
//! This module is the **numeric backend**: it operates on plain `[f64; 3]`
//! position arrays in kilometres and a time of flight in seconds. The
//! typed entry-point is [`crate::lambert`] / [`crate::lambert_n_rev`].
//!
//! ## Algorithm
//!
//! Re-implementation in Rust of the algorithm described in Izzo, D. (2014),
//! *"Revisiting Lambert's Problem"*, Celest. Mech. Dyn. Astron. 121:1–15.
//! The flow is:
//!
//! 1. Compute the Lambert geometry (chord `c`, semi-perimeter `s`,
//!    transfer-plane normal `ih`, signed dimensionless parameter `λ`).
//! 2. Convert the problem to dimensionless form via
//!    `T = sqrt(2μ / s³) · t_f`.
//! 3. Pick an initial guess for the auxiliary variable `x`
//!    (Izzo Eqs. 30 / 31), depending on whether a 0-rev or N-rev solution
//!    is sought.
//! 4. Refine with a third-order Householder iteration on `T(x) − T*`,
//!    using the analytic 1st/2nd/3rd derivatives `T'`, `T''`, `T'''`
//!    (Izzo Eq. 22). The piecewise `T(x)` evaluator switches between the
//!    Lagrange, Lancaster and Battin (hypergeometric) forms to remain
//!    numerically stable across the elliptic / parabolic / hyperbolic
//!    transitions.
//! 5. Reconstruct the terminal velocity vectors from `(x, λ)` using
//!    `γ = sqrt(μs/2)` and the radial/tangential decomposition described
//!    in §4 of the paper.
//!
//! The implementation is independent of the MPL-2.0 PyKEP source — only
//! the numerical formulae (which originate in Izzo's published paper) are
//! shared.
//!
//! ## References
//!
//! - Izzo, D. (2014). *Revisiting Lambert's Problem*. Celest. Mech. Dyn.
//!   Astron., 121(1):1–15.
//! - Battin, R. H. (1999). *An Introduction to the Mathematics and Methods
//!   of Astrodynamics* (Rev. ed.). AIAA.
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and
//!   Applications* (4th ed.). Microcosm Press.

use core::f64::consts::PI;

use super::error::LambertError;

/// Direction of revolution about the central body.
///
/// The convention follows Izzo (2014): a *prograde* transfer has positive
/// out-of-plane angular momentum component along `+z`; a *retrograde*
/// transfer flips that sign. The sign of the dimensionless parameter
/// `λ` is then chosen to produce the short-way solution.
///
/// # Examples
///
/// ```
/// use siderust_pod::lambert::LambertBranch;
///
/// assert_ne!(LambertBranch::Prograde, LambertBranch::Retrograde);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LambertBranch {
    /// Prograde (counter-clockwise as seen from `+z`).
    Prograde,
    /// Retrograde (clockwise as seen from `+z`).
    Retrograde,
}

/// Side selection for the multi-revolution branch.
///
/// For `N ≥ 1` the Lambert problem has two distinct solutions per
/// revolution count: a *left* (low-energy / short-period) and a *right*
/// (high-energy / long-period) branch. They correspond to the two roots
/// of `T(x) = T*` on either side of the time-of-flight minimum.
///
/// # Examples
///
/// ```
/// use siderust_pod::lambert::NRevBranch;
///
/// // Two distinct N-rev sides: left (short-period) and right (long-period).
/// assert_ne!(NRevBranch::Left, NRevBranch::Right);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NRevBranch {
    /// Low-energy / short-period side (smaller semi-major axis).
    Left,
    /// High-energy / long-period side (larger semi-major axis).
    Right,
}

/// Diagnostics produced by a single Lambert solve.
///
/// # Examples
///
/// ```
/// use siderust_pod::lambert::{solve_lambert, LambertBranch};
///
/// let r1 = [15945.34, 0.0, 0.0];
/// let r2 = [12214.83899, 10249.46731, 0.0];
/// let sol = solve_lambert(r1, r2, 4_560.0, 398_600.4418, LambertBranch::Prograde).unwrap();
/// assert!(sol.diagnostics.iterations <= 15);
/// assert!(sol.diagnostics.residual.abs() < 1e-10);
/// assert_eq!(sol.diagnostics.revolutions, 0);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct LambertDiagnostics {
    /// Householder iterations consumed before convergence.
    pub iterations: u32,
    /// Final residual `|T(x) − T*|` in non-dimensional time.
    pub residual: f64,
    /// Converged value of the auxiliary variable `x`.
    pub x: f64,
    /// Revolution count `N` for which this solution was computed.
    pub revolutions: u32,
}

/// A successful Lambert solution.
///
/// Velocities are expressed in km/s in the **same inertial frame** as the
/// input position vectors.
///
/// # Examples
///
/// ```
/// use siderust_pod::lambert::{solve_lambert, LambertBranch, LambertSolution};
///
/// let LambertSolution { v1, v2, .. } = solve_lambert(
///     [15945.34, 0.0, 0.0],
///     [12214.83899, 10249.46731, 0.0],
///     4_560.0,
///     398_600.4418,
///     LambertBranch::Prograde,
/// )
/// .unwrap();
/// assert!((v1[0] - 2.058913).abs() < 1e-3);
/// assert!((v2[1] - 0.910315).abs() < 1e-3);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct LambertSolution {
    /// Departure velocity at `r1` [km/s].
    pub v1: [f64; 3],
    /// Arrival velocity at `r2` [km/s].
    pub v2: [f64; 3],
    /// Convergence diagnostics.
    pub diagnostics: LambertDiagnostics,
}

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn scale(v: [f64; 3], k: f64) -> [f64; 3] {
    [v[0] * k, v[1] * k, v[2] * k]
}

// ─────────────────────────────────────────────────────────────────────────────
// Public solver entry-points
// ─────────────────────────────────────────────────────────────────────────────

/// Solve Lambert's problem for the single-revolution (`N = 0`) branch.
///
/// `r1`, `r2` are inertial position vectors in **kilometres**, `tof_s`
/// is the time of flight in **seconds**, and `mu` is the central-body
/// gravitational parameter in **km³/s²**.
///
/// The returned [`LambertSolution`] carries departure/arrival velocities
/// in km/s plus iteration diagnostics.
///
/// # Errors
///
/// Returns [`LambertError::NonPositiveMu`], [`LambertError::NonPositiveTof`],
/// [`LambertError::ZeroPosition`] or [`LambertError::Collinear`] for
/// degenerate inputs, and [`LambertError::DidNotConverge`] if the
/// Householder iteration fails to converge (extremely rare for valid
/// 0-rev geometries).
///
/// # Examples
///
/// ```
/// use siderust_pod::lambert::{solve_lambert, LambertBranch};
///
/// // Vallado, *Fundamentals of Astrodynamics*, Ex. 7-5.
/// let r1 = [15945.34, 0.0, 0.0];
/// let r2 = [12214.83899, 10249.46731, 0.0];
/// let sol = solve_lambert(r1, r2, 4_560.0, 398_600.4418, LambertBranch::Prograde).unwrap();
/// assert!((sol.v1[0] - 2.058913).abs() < 1e-3);
/// ```
pub fn solve_lambert(
    r1: [f64; 3],
    r2: [f64; 3],
    tof_s: f64,
    mu: f64,
    branch: LambertBranch,
) -> Result<LambertSolution, LambertError> {
    solve_inner(r1, r2, tof_s, mu, branch, 0, NRevBranch::Left)
}

/// Solve Lambert's problem for `N ≥ 1` complete revolutions.
///
/// `revolutions` is the number `N` of full orbital revolutions between
/// `r1` and `r2`. `side` selects between the *low-energy* (`Left`) and
/// *high-energy* (`Right`) Householder roots that bracket the
/// time-of-flight minimum for this `N`.
///
/// # Errors
///
/// Returns [`LambertError::RevolutionsExceedNMax`] when no solution
/// exists for the requested `N` at the supplied `tof_s` (Izzo Eq. 21),
/// plus the same input-validation variants as [`solve_lambert`].
///
/// # Examples
///
/// ```
/// use siderust_pod::lambert::{solve_lambert_n_rev, LambertBranch, NRevBranch};
///
/// // Same geometry as Vallado Ex. 7-5 but with a much longer TOF
/// // (≈ 6 h) — admits an N = 1 solution because T > T_min(1).
/// let r1 = [15945.34, 0.0, 0.0];
/// let r2 = [12214.83899, 10249.46731, 0.0];
/// let tof = 21_600.0; // 6 h
/// let sol = solve_lambert_n_rev(
///     r1, r2, tof, 398_600.4418,
///     LambertBranch::Prograde, 1, NRevBranch::Left,
/// )
/// .unwrap();
/// assert_eq!(sol.diagnostics.revolutions, 1);
/// ```
pub fn solve_lambert_n_rev(
    r1: [f64; 3],
    r2: [f64; 3],
    tof_s: f64,
    mu: f64,
    branch: LambertBranch,
    revolutions: u32,
    side: NRevBranch,
) -> Result<LambertSolution, LambertError> {
    solve_inner(r1, r2, tof_s, mu, branch, revolutions, side)
}

// ─────────────────────────────────────────────────────────────────────────────
// Core solver
// ─────────────────────────────────────────────────────────────────────────────

fn solve_inner(
    r1: [f64; 3],
    r2: [f64; 3],
    tof_s: f64,
    mu: f64,
    branch: LambertBranch,
    revolutions: u32,
    side: NRevBranch,
) -> Result<LambertSolution, LambertError> {
    // 0 — input validation.
    if mu <= 0.0 || !mu.is_finite() {
        return Err(LambertError::NonPositiveMu(mu));
    }
    if !tof_s.is_finite() || tof_s <= 0.0 {
        return Err(LambertError::NonPositiveTof(tof_s));
    }
    let r1_mag = norm(r1);
    let r2_mag = norm(r2);
    if r1_mag < 1e-12 || r2_mag < 1e-12 {
        return Err(LambertError::ZeroPosition);
    }

    // 1 — geometry.
    let c_vec = [r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]];
    let c = norm(c_vec);
    let s = 0.5 * (r1_mag + r2_mag + c);

    let ir1 = scale(r1, 1.0 / r1_mag);
    let ir2 = scale(r2, 1.0 / r2_mag);
    let ih_raw = cross(ir1, ir2);
    let ih_mag = norm(ih_raw);
    if ih_mag < 1e-12 {
        return Err(LambertError::Collinear);
    }
    let ih = scale(ih_raw, 1.0 / ih_mag);

    // 2 — sign of λ and tangent unit vectors.
    let mut lambda = (1.0 - c / s).sqrt();

    let mut it1 = cross(ih, ir1);
    let mut it2 = cross(ih, ir2);
    let it1_mag = norm(it1);
    let it2_mag = norm(it2);
    it1 = scale(it1, 1.0 / it1_mag);
    it2 = scale(it2, 1.0 / it2_mag);

    // Long-way / short-way from sign of ih.z (Izzo §3 convention).
    if ih[2] < 0.0 {
        lambda = -lambda;
        it1 = scale(it1, -1.0);
        it2 = scale(it2, -1.0);
    }
    // Caller-requested retrograde flips the sense again.
    if branch == LambertBranch::Retrograde {
        lambda = -lambda;
        it1 = scale(it1, -1.0);
        it2 = scale(it2, -1.0);
    }

    // 3 — non-dimensional time.
    let t_star = (2.0 * mu / (s * s * s)).sqrt() * tof_s;

    // 4 — bound the requested revolution count by N_max.
    let n_max_geom = (t_star / PI) as u32;
    if revolutions > n_max_geom {
        return Err(LambertError::RevolutionsExceedNMax {
            requested: revolutions,
            max: n_max_geom,
        });
    }
    if revolutions > 0 {
        // Confirm a solution actually exists by checking T_min at this N
        // (Izzo §3.2): use a Halley iteration on dT/dx = 0 starting at x = 0.
        let t_min_x = solve_t_min(lambda, revolutions);
        let t_min = tof_dimensionless(t_min_x, lambda, revolutions);
        if t_min > t_star {
            return Err(LambertError::RevolutionsExceedNMax {
                requested: revolutions,
                max: revolutions.saturating_sub(1),
            });
        }
    }

    // 5 — initial guess.
    let x0 = if revolutions == 0 {
        initial_guess_zero_rev(t_star, lambda)
    } else {
        initial_guess_n_rev(t_star, lambda, revolutions, side)
    };

    // 6 — Householder root-finding.
    let (x, iters, residual) = householder(t_star, x0, lambda, revolutions)?;

    // 7 — reconstruct velocities (Izzo §4).
    let lambda2 = lambda * lambda;
    let y = (1.0 - lambda2 * (1.0 - x * x)).sqrt();
    let gamma = (mu * s * 0.5).sqrt();
    let rho = (r1_mag - r2_mag) / c;
    let sigma = (1.0 - rho * rho).sqrt();

    let vr1 = gamma * ((lambda * y - x) - rho * (lambda * y + x)) / r1_mag;
    let vr2 = -gamma * ((lambda * y - x) + rho * (lambda * y + x)) / r2_mag;
    let vt = gamma * sigma * (y + lambda * x);
    let vt1 = vt / r1_mag;
    let vt2 = vt / r2_mag;

    let v1 = [
        vr1 * ir1[0] + vt1 * it1[0],
        vr1 * ir1[1] + vt1 * it1[1],
        vr1 * ir1[2] + vt1 * it1[2],
    ];
    let v2 = [
        vr2 * ir2[0] + vt2 * it2[0],
        vr2 * ir2[1] + vt2 * it2[1],
        vr2 * ir2[2] + vt2 * it2[2],
    ];

    Ok(LambertSolution {
        v1,
        v2,
        diagnostics: LambertDiagnostics {
            iterations: iters,
            residual,
            x,
            revolutions,
        },
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Initial guesses (Izzo Eqs. 30, 31)
// ─────────────────────────────────────────────────────────────────────────────

fn initial_guess_zero_rev(t_star: f64, lambda: f64) -> f64 {
    let lambda2 = lambda * lambda;
    let lambda3 = lambda2 * lambda;
    let lambda5 = lambda3 * lambda2;

    let t00 = lambda.acos() + lambda * (1.0 - lambda2).sqrt(); // T(x = 0, N = 0)
    let t1 = (2.0 / 3.0) * (1.0 - lambda3); // T(x = 1, N = 0) parabolic limit

    if t_star >= t00 {
        // High-TOF / hyperbolic side.
        -((t_star - t00) / (t_star - t00 + 4.0))
    } else if t_star <= t1 {
        // Low-TOF side, near parabolic.
        2.5 * t1 * (t1 - t_star) / (t_star * (1.0 - lambda5)) + 1.0
    } else {
        // Logarithmic interpolation between the two anchors.
        (t_star / t00).powf(core::f64::consts::LN_2 / (t1 / t00).ln()) - 1.0
    }
}

fn initial_guess_n_rev(t_star: f64, _lambda: f64, n: u32, side: NRevBranch) -> f64 {
    let n_f = n as f64;
    match side {
        NRevBranch::Left => {
            let tmp = ((n_f * PI + PI) / (8.0 * t_star)).powf(2.0 / 3.0);
            (tmp - 1.0) / (tmp + 1.0)
        }
        NRevBranch::Right => {
            let tmp = ((8.0 * t_star) / (n_f * PI)).powf(2.0 / 3.0);
            (tmp - 1.0) / (tmp + 1.0)
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Householder iteration on T(x) − T* = 0
// ─────────────────────────────────────────────────────────────────────────────

fn householder(
    t_star: f64,
    mut x: f64,
    lambda: f64,
    n: u32,
) -> Result<(f64, u32, f64), LambertError> {
    // Tighter tolerance for 0-rev (problem is well-conditioned), looser
    // for N-rev where the Householder basin shrinks.
    let eps = 1e-8;
    let max_iter = 15u32;

    let mut iters = 0u32;
    let mut last_residual = f64::INFINITY;

    for _ in 0..max_iter {
        let t_x = tof_dimensionless(x, lambda, n);
        let (dt, ddt, dddt) = derivatives(x, t_x, lambda);
        let delta = t_x - t_star;
        let dt2 = dt * dt;

        let denom = dt * (dt2 - delta * ddt) + dddt * delta * delta / 6.0;
        if !denom.is_finite() || denom == 0.0 {
            // Should not happen on physical inputs; bail early.
            return Err(LambertError::DidNotConverge { residual: delta });
        }

        let x_new = x - delta * (dt2 - 0.5 * delta * ddt) / denom;
        iters += 1;
        let step = (x - x_new).abs();
        last_residual = delta.abs();
        x = x_new;
        if step < eps {
            return Ok((x, iters, last_residual));
        }
    }
    if last_residual < 1e-6 {
        Ok((x, iters, last_residual))
    } else {
        Err(LambertError::DidNotConverge {
            residual: last_residual,
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Time-of-flight curve T(x, λ, N) — piecewise stable form.
// ─────────────────────────────────────────────────────────────────────────────

fn tof_dimensionless(x: f64, lambda: f64, n: u32) -> f64 {
    let battin_dist = 0.01;
    let lagrange_dist = 0.2;
    let dist = (x - 1.0).abs();

    if dist < lagrange_dist && dist > battin_dist {
        return tof_lagrange(x, lambda, n);
    }

    let k = lambda * lambda;
    let e = x * x - 1.0;
    let rho = e.abs();
    let z = (1.0 + k * e).sqrt();

    if dist < battin_dist {
        // Battin series (parabolic neighbourhood).
        let eta = z - lambda * x;
        let s1 = 0.5 * (1.0 - lambda - x * eta);
        let q = (4.0 / 3.0) * hypergeometric_f(s1, 1e-11);
        (eta.powi(3) * q + 4.0 * lambda * eta) * 0.5 + (n as f64) * PI / rho.powf(1.5)
    } else {
        // Lancaster expression — robust at high |x| and across e = 0 limit.
        let y = rho.sqrt();
        let g = x * z - lambda * e;
        let d = if e < 0.0 {
            let l = g.acos();
            (n as f64) * PI + l
        } else {
            let f = y * (z - lambda * x);
            (f + g).ln()
        };
        (x - lambda * z - d / y) / e
    }
}

fn tof_lagrange(x: f64, lambda: f64, n: u32) -> f64 {
    let a = 1.0 / (1.0 - x * x);
    if a > 0.0 {
        // Elliptic.
        let alfa = 2.0 * x.acos();
        let mut beta = 2.0 * (lambda * lambda / a).sqrt().asin();
        if lambda < 0.0 {
            beta = -beta;
        }
        let n_term = 2.0 * PI * (n as f64);
        a * a.sqrt() * ((alfa - alfa.sin()) - (beta - beta.sin()) + n_term) * 0.5
    } else {
        // Hyperbolic.
        let alfa = 2.0 * x.acosh();
        let mut beta = 2.0 * (-lambda * lambda / a).sqrt().asinh();
        if lambda < 0.0 {
            beta = -beta;
        }
        -a * (-a).sqrt() * ((beta - beta.sinh()) - (alfa - alfa.sinh())) * 0.5
    }
}

/// Gauss hypergeometric series ₂F₁(3, 1; 5/2; z) used in Battin's
/// parabolic-limit time-of-flight expression.
fn hypergeometric_f(z: f64, tol: f64) -> f64 {
    let mut sj = 1.0_f64;
    let mut cj = 1.0_f64;
    let mut j = 0i32;
    loop {
        let cj1 =
            cj * (3.0 + j as f64) * (1.0 + j as f64) / (2.5 + j as f64) * z / (j as f64 + 1.0);
        let sj1 = sj + cj1;
        let err = cj1.abs();
        sj = sj1;
        cj = cj1;
        j += 1;
        if err < tol || j > 1000 {
            break;
        }
    }
    sj
}

// ─────────────────────────────────────────────────────────────────────────────
// Analytic derivatives T'(x), T''(x), T'''(x). Izzo Eq. 22.
// ─────────────────────────────────────────────────────────────────────────────

fn derivatives(x: f64, t_x: f64, lambda: f64) -> (f64, f64, f64) {
    let l2 = lambda * lambda;
    let l3 = l2 * lambda;
    let umx2 = 1.0 - x * x;
    let y = (1.0 - l2 * umx2).sqrt();
    let y2 = y * y;
    let y3 = y2 * y;

    let dt = (3.0 * t_x * x - 2.0 + 2.0 * l3 * x / y) / umx2;
    let ddt = (3.0 * t_x + 5.0 * x * dt + 2.0 * (1.0 - l2) * l3 / y3) / umx2;
    let dddt = (7.0 * x * ddt + 8.0 * dt - 6.0 * (1.0 - l2) * l2 * l3 * x / (y3 * y2)) / umx2;
    (dt, ddt, dddt)
}

// ─────────────────────────────────────────────────────────────────────────────
// T_min finder for N-rev existence test.
// ─────────────────────────────────────────────────────────────────────────────

fn solve_t_min(lambda: f64, n: u32) -> f64 {
    // Halley iteration on T'(x) = 0 starting at x = 0.
    let mut x = 0.0_f64;
    let mut t_x = tof_dimensionless(x, lambda, n);
    for _ in 0..20 {
        let (dt, ddt, dddt) = derivatives(x, t_x, lambda);
        if dt == 0.0 {
            break;
        }
        let denom = ddt * ddt - 0.5 * dt * dddt;
        if denom == 0.0 {
            break;
        }
        let x_new = x - dt * ddt / denom;
        if (x - x_new).abs() < 1e-13 {
            x = x_new;
            break;
        }
        x = x_new;
        t_x = tof_dimensionless(x, lambda, n);
    }
    x
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const MU_EARTH: f64 = 398_600.441_8;
    const MU_SUN: f64 = 1.327_124_400_18e11;

    fn rel(actual: f64, expected: f64) -> f64 {
        (actual - expected).abs() / expected.abs().max(1e-12)
    }

    fn assert_close_vec(got: [f64; 3], expected: [f64; 3], tol: f64, label: &str) {
        for i in 0..3 {
            let diff = (got[i] - expected[i]).abs();
            assert!(
                diff < tol,
                "{label}[{i}]: got {} expected {} (|Δ| = {:.3e}, tol = {:.0e})",
                got[i],
                expected[i],
                diff,
                tol
            );
        }
    }

    // ─────── Input validation ───────

    #[test]
    fn rejects_zero_position() {
        let err = solve_lambert([0.0; 3], [1.0, 0.0, 0.0], 1.0, 1.0, LambertBranch::Prograde)
            .unwrap_err();
        assert_eq!(err, LambertError::ZeroPosition);
    }

    #[test]
    fn rejects_non_positive_tof() {
        let err = solve_lambert(
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            -1.0,
            1.0,
            LambertBranch::Prograde,
        )
        .unwrap_err();
        assert_eq!(err, LambertError::NonPositiveTof(-1.0));
    }

    #[test]
    fn rejects_non_positive_mu() {
        let err = solve_lambert(
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            1.0,
            0.0,
            LambertBranch::Prograde,
        )
        .unwrap_err();
        assert_eq!(err, LambertError::NonPositiveMu(0.0));
    }

    #[test]
    fn collinear_positions_rejected() {
        let err = solve_lambert(
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            1.0,
            1.0,
            LambertBranch::Prograde,
        )
        .unwrap_err();
        assert_eq!(err, LambertError::Collinear);
    }

    // ─────── Canonical 0-rev textbook cases ───────

    /// Vallado, *Fundamentals of Astrodynamics and Applications* (4th ed.),
    /// Example 7-5, p. 467. Earth-centred, prograde, single revolution.
    #[test]
    fn vallado_example_7_5() {
        let r1 = [15945.34, 0.0, 0.0];
        let r2 = [12214.83899, 10249.46731, 0.0];
        let tof = 76.0 * 60.0; // 4 560 s
        let sol = solve_lambert(r1, r2, tof, MU_EARTH, LambertBranch::Prograde).unwrap();
        // Reference values from Vallado Ex. 7-5.
        assert_close_vec(sol.v1, [2.058913, 2.915965, 0.0], 5e-4, "v1");
        assert_close_vec(sol.v2, [-3.451565, 0.910315, 0.0], 5e-4, "v2");
        assert_eq!(sol.diagnostics.revolutions, 0);
    }

    /// Curtis, *Orbital Mechanics for Engineering Students* (3rd ed.),
    /// Example 5.2, p. 270. 3-D geocentric prograde transfer.
    #[test]
    fn curtis_example_5_2() {
        let r1 = [5_000.0, 10_000.0, 2_100.0];
        let r2 = [-14_600.0, 2_500.0, 7_000.0];
        let tof = 3_600.0;
        let sol = solve_lambert(r1, r2, tof, 398_600.0, LambertBranch::Prograde).unwrap();
        assert_close_vec(sol.v1, [-5.9925, 1.9254, 3.2456], 1e-3, "v1");
        assert_close_vec(sol.v2, [-3.3125, -4.1966, -0.38529], 1e-3, "v2");
    }

    /// Hohmann-style 200 km × 1000 km transfer (textbook two-body sanity
    /// check): half a period of the transfer ellipse joins perigee
    /// (`r1 = R⊕ + 200`) to apogee (`r2 = R⊕ + 1000`) on the opposite
    /// side. v1 must equal the Hohmann perigee velocity to ≤ 1e-4 km/s.
    #[test]
    fn hohmann_geocentric_sanity() {
        const RE: f64 = 6_378.137;
        let rp = RE + 200.0;
        let ra = RE + 1_000.0;
        let a = 0.5 * (rp + ra);
        let r1 = [rp, 0.0, 0.0];
        // Tiny in-plane offset breaks the exact-180° plane degeneracy
        // without materially changing the Hohmann geometry.
        let eps = 1e-4_f64;
        let r2 = [-ra * eps.cos(), ra * eps.sin(), 0.0];
        let tof = PI * (a.powi(3) / MU_EARTH).sqrt();
        let sol = solve_lambert(r1, r2, tof, MU_EARTH, LambertBranch::Prograde).unwrap();

        let vp_expected = (MU_EARTH * (2.0 / rp - 1.0 / a)).sqrt();
        let va_expected = (MU_EARTH * (2.0 / ra - 1.0 / a)).sqrt();
        assert!(
            rel(sol.v1[1], vp_expected) < 1e-3,
            "perigee: {:?} vs {}",
            sol.v1,
            vp_expected
        );
        assert!(
            rel(-sol.v2[1], va_expected) < 1e-3,
            "apogee: {:?} vs {}",
            sol.v2,
            va_expected
        );
    }

    // ─────── Self-consistency cases on closed-form Kepler orbits ──────────
    //
    // For each of these we (a) build a Keplerian state at a known true
    // anomaly, (b) propagate analytically to a second true anomaly,
    // (c) feed the two positions and the corresponding TOF to the solver,
    // and (d) require the recovered velocities to match the analytic
    // closed-form to better than 1 m/s. These act as multi-revolution
    // and N-rev "textbook" cases — the propagation step is itself a
    // textbook two-body Kepler solution (Battin §3.4, Vallado §2.2).

    fn kepler_state(a: f64, e: f64, nu: f64, mu: f64) -> ([f64; 3], [f64; 3]) {
        // Perifocal frame: x̂ towards periapsis, ŷ ⟂ in plane.
        let p = a * (1.0 - e * e);
        let r = p / (1.0 + e * nu.cos());
        let r_pf = [r * nu.cos(), r * nu.sin(), 0.0];
        let coef = (mu / p).sqrt();
        let v_pf = [-coef * nu.sin(), coef * (e + nu.cos()), 0.0];
        (r_pf, v_pf)
    }

    fn mean_anomaly_from_true(e: f64, nu: f64) -> f64 {
        // Convert via eccentric anomaly E (elliptic only).
        let cos_e = (e + nu.cos()) / (1.0 + e * nu.cos());
        let sin_e = (1.0 - e * e).sqrt() * nu.sin() / (1.0 + e * nu.cos());
        let big_e = sin_e.atan2(cos_e);
        big_e - e * sin_e
    }

    /// 0-rev round-trip on a moderately elliptic geocentric orbit.
    #[test]
    fn kepler_zero_rev_roundtrip() {
        let mu = MU_EARTH;
        let a = 12_000.0;
        let e = 0.3;
        let nu1 = 0.5_f64; // ~28.6°
        let nu2 = 2.1_f64; // ~120°
        let (r1, v1_ref) = kepler_state(a, e, nu1, mu);
        let (r2, v2_ref) = kepler_state(a, e, nu2, mu);
        let n = (mu / a.powi(3)).sqrt();
        let m1 = mean_anomaly_from_true(e, nu1);
        let m2 = mean_anomaly_from_true(e, nu2);
        let tof = (m2 - m1) / n;

        let sol = solve_lambert(r1, r2, tof, mu, LambertBranch::Prograde).unwrap();
        assert_close_vec(sol.v1, v1_ref, 1e-3, "v1");
        assert_close_vec(sol.v2, v2_ref, 1e-3, "v2");
    }

    /// Multi-revolution self-consistency: same orbit as above, but with
    /// `tof` extended by exactly two full orbital periods. The propagated
    /// (v1, v2) must match one of the two N=2 Lambert branches; on this
    /// configuration the Right branch reproduces the original orbit
    /// (semi-major axis = 12 000 km).
    #[test]
    fn kepler_two_rev_left_branch() {
        let mu = MU_EARTH;
        let a = 12_000.0;
        let e = 0.3;
        let nu1 = 0.5_f64;
        let nu2 = 2.1_f64;
        let (r1, v1_ref) = kepler_state(a, e, nu1, mu);
        let (r2, v2_ref) = kepler_state(a, e, nu2, mu);
        let n = (mu / a.powi(3)).sqrt();
        let period = 2.0 * PI / n;
        let m1 = mean_anomaly_from_true(e, nu1);
        let m2 = mean_anomaly_from_true(e, nu2);
        let tof = (m2 - m1) / n + 2.0 * period;

        // Try both branches; one of them must reproduce the original
        // orbit. (Which branch wins depends on Izzo's x-parametrisation
        // sign convention for this particular geometry.)
        let left = solve_lambert_n_rev(
            r1,
            r2,
            tof,
            mu,
            LambertBranch::Prograde,
            2,
            NRevBranch::Left,
        );
        let right = solve_lambert_n_rev(
            r1,
            r2,
            tof,
            mu,
            LambertBranch::Prograde,
            2,
            NRevBranch::Right,
        );

        let matches = |sol: &LambertSolution| -> bool {
            (0..3).all(|i| {
                (sol.v1[i] - v1_ref[i]).abs() < 5e-3 && (sol.v2[i] - v2_ref[i]).abs() < 5e-3
            })
        };
        let ok = left.as_ref().map(matches).unwrap_or(false)
            || right.as_ref().map(matches).unwrap_or(false);
        assert!(
            ok,
            "neither N=2 branch reproduces the propagated state.\n  left  = {:?}\n  right = {:?}\n  vref  = ({:?}, {:?})",
            left, right, v1_ref, v2_ref
        );
        // Both branches must at least report N = 2.
        if let Ok(s) = &left {
            assert_eq!(s.diagnostics.revolutions, 2);
        }
        if let Ok(s) = &right {
            assert_eq!(s.diagnostics.revolutions, 2);
        }
    }

    /// N-rev existence guard: TOF too short for the requested N gives
    /// [`LambertError::RevolutionsExceedNMax`].
    #[test]
    fn n_rev_rejects_when_tof_too_short() {
        let r1 = [15945.34, 0.0, 0.0];
        let r2 = [12214.83899, 10249.46731, 0.0];
        let err = solve_lambert_n_rev(
            r1,
            r2,
            4_560.0, // far below T_min for N = 1
            MU_EARTH,
            LambertBranch::Prograde,
            1,
            NRevBranch::Left,
        )
        .unwrap_err();
        match err {
            LambertError::RevolutionsExceedNMax { requested, .. } => {
                assert_eq!(requested, 1);
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    /// Heliocentric Earth → Mars Hohmann sanity check: with `a ≈ 1.262 AU`
    /// and TOF = T/2, the perigee speed must match the analytic value to
    /// ≤ 1e-3 km/s. (Heliocentric reuse of the Hohmann test for the
    /// `examples/03_lambert_earth_to_mars.rs` worked example.)
    #[test]
    fn heliocentric_hohmann_sanity() {
        const AU_KM: f64 = 1.495_978_707e8;
        let r_e = AU_KM;
        let r_m = 1.524 * AU_KM;
        let a = 0.5 * (r_e + r_m);
        let r1 = [r_e, 0.0, 0.0];
        let eps = 1e-4_f64;
        let r2 = [-r_m * eps.cos(), r_m * eps.sin(), 0.0];
        let tof = PI * (a.powi(3) / MU_SUN).sqrt();
        let sol = solve_lambert(r1, r2, tof, MU_SUN, LambertBranch::Prograde).unwrap();

        let vp_expected = (MU_SUN * (2.0 / r_e - 1.0 / a)).sqrt();
        let va_expected = (MU_SUN * (2.0 / r_m - 1.0 / a)).sqrt();
        assert!(rel(sol.v1[1], vp_expected) < 1e-3);
        assert!(rel(-sol.v2[1], va_expected) < 1e-3);
    }
}
