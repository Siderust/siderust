//! # Gauss-Newton nonlinear estimation
//!
//! ## Scientific scope
//!
//! This module implements a classic Gauss-Newton outer loop for mildly
//! nonlinear POD problems, where the design matrix must be rebuilt around
//! the current best estimate at each iteration. It is appropriate for
//! short-arc orbit-state corrections and similar parameter updates when the
//! residual surface is locally well behaved.
//!
//! The algorithm assumes the caller can re-linearize the problem after each
//! update. It does not perform trust-region control, robust editing, or
//! global convergence safeguards beyond simple stopping criteria on step
//! size and reduced chi-square change.
//!
//! ## Technical scope
//!
//! The entry points are `NonlinearOptions`, `NonlinearReport`,
//! `NonlinearError`, and `gauss_newton`. Callers provide a closure that
//! assembles fresh `NormalEquations` from the current parameter vector, and
//! the solver returns an updated parameter vector plus iteration
//! diagnostics in raw solver coordinates.
//!
//! This module does not own measurement units or force-model semantics.
//! Those stay with the upstream model code that produces the normal
//! equations.
//!
//! ## References
//!
//! - Tapley, B. D., Schutz, B. E., & Born, G. H. (2004). Statistical Orbit
//!   Determination. Elsevier Academic Press.
//! - Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications
//!   (4th ed.). Microcosm Press.
use super::wls::{NormalEquations, WlsResult, WlsSolverError};
use thiserror::Error;

/// Convergence options for [`gauss_newton`].
#[derive(Debug, Clone, Copy)]
pub struct NonlinearOptions {
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Relative-update convergence tolerance (parameter L2-norm).
    pub tol_rel: f64,
    /// Relative-χ² convergence tolerance.
    pub tol_chi2_rel: f64,
}

impl Default for NonlinearOptions {
    fn default() -> Self {
        Self {
            max_iter: 12,
            tol_rel: 1e-9,
            tol_chi2_rel: 1e-6,
        }
    }
}

/// Errors during nonlinear iteration.
#[derive(Debug, Error)]
pub enum NonlinearError {
    /// WLS solve failed.
    #[error(transparent)]
    Solver(#[from] WlsSolverError),
    /// Reached `max_iter` without convergence.
    #[error("did not converge in {0} iterations")]
    DidNotConverge(usize),
}

/// Result of a nonlinear iteration loop.
#[derive(Debug, Clone)]
pub struct NonlinearReport {
    /// Final parameter vector.
    pub parameters: Vec<f64>,
    /// Final WLS solve.
    pub last: WlsResult,
    /// Number of iterations performed.
    pub iterations: usize,
}

/// Run Gauss-Newton iterations.
///
/// `assemble(params) -> NormalEquations` is invoked at every iteration
/// with the *current* best parameter vector and is responsible for
/// producing prefit residuals and partials around that point.
///
/// # Errors
///
/// Returns [`NonlinearError::Solver`] wrapping a [`WlsSolverError::Other`]
/// if `opts.max_iter == 0`, any tolerance is non-finite or negative, or
/// any element of `initial` is non-finite.  Returns
/// [`NonlinearError::DidNotConverge`] if the loop completes without
/// satisfying the convergence criteria.
pub fn gauss_newton<F>(
    initial: Vec<f64>,
    opts: NonlinearOptions,
    mut assemble: F,
) -> Result<NonlinearReport, NonlinearError>
where
    F: FnMut(&[f64]) -> Result<NormalEquations, NonlinearError>,
{
    use super::wls::WlsSolverError;

    if opts.max_iter == 0 {
        return Err(NonlinearError::Solver(WlsSolverError::other(
            "max_iter must be > 0",
        )));
    }
    if !opts.tol_rel.is_finite() || opts.tol_rel < 0.0 {
        return Err(NonlinearError::Solver(WlsSolverError::other(format!(
            "tol_rel must be finite and ≥ 0 (got {})",
            opts.tol_rel
        ))));
    }
    if !opts.tol_chi2_rel.is_finite() || opts.tol_chi2_rel < 0.0 {
        return Err(NonlinearError::Solver(WlsSolverError::other(format!(
            "tol_chi2_rel must be finite and ≥ 0 (got {})",
            opts.tol_chi2_rel
        ))));
    }
    if let Some(bad) = initial.iter().copied().find(|v| !v.is_finite()) {
        return Err(NonlinearError::Solver(WlsSolverError::other(format!(
            "initial vector must be finite (got {bad})"
        ))));
    }
    let mut params = initial;
    let mut last_chi2 = f64::INFINITY;
    #[allow(unused_assignments)]
    let mut last_result: Option<WlsResult> = None;
    for it in 0..opts.max_iter {
        let ne = assemble(&params)?;
        if ne.n_params() != params.len() {
            return Err(NonlinearError::Solver(WlsSolverError::other(format!(
                "assembler returned {} parameters but initial vector has {} entries",
                ne.n_params(),
                params.len()
            ))));
        }
        let result = ne.solve()?;
        let mut max_rel = 0.0_f64;
        for (i, dp) in result.update.iter().enumerate() {
            let scale = params[i].abs().max(1.0);
            let rel = dp.abs() / scale;
            if rel > max_rel {
                max_rel = rel;
            }
            params[i] += dp;
        }
        let chi2 = result.reduced_chi2();
        let chi2_rel = if last_chi2.is_finite() && last_chi2 > 0.0 {
            (chi2 - last_chi2).abs() / last_chi2
        } else {
            f64::INFINITY
        };
        last_chi2 = chi2;
        last_result = Some(result);
        if max_rel < opts.tol_rel && chi2_rel < opts.tol_chi2_rel {
            return Ok(NonlinearReport {
                parameters: params,
                last: last_result.unwrap(),
                iterations: it + 1,
            });
        }
    }
    Err(NonlinearError::DidNotConverge(opts.max_iter))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pod::estimation::wls::NormalEquations;

    fn default_opts() -> NonlinearOptions {
        NonlinearOptions {
            max_iter: 20,
            tol_rel: 1e-9,
            tol_chi2_rel: 1e-6,
        }
    }

    #[test]
    fn converges_for_1d_linear_problem() {
        // y = x, measurement y = 3.0, sigma = 0.1
        // After 1 Gauss-Newton step params[0]=3.0; use tol_chi2_rel=2.0 so iter 2 terminates.
        let opts = NonlinearOptions {
            max_iter: 20,
            tol_rel: 1e-9,
            tol_chi2_rel: 2.0,
        };
        let report = gauss_newton(vec![0.0], opts, |params| {
            let mut ne = NormalEquations::new(1);
            ne.add_row(&[(0, 1.0)], 3.0 - params[0], 0.1).unwrap();
            Ok(ne)
        })
        .unwrap();
        assert!((report.parameters[0] - 3.0).abs() < 1e-6);
        assert!(report.iterations >= 1);
    }

    #[test]
    fn max_iter_zero_is_error() {
        let opts = NonlinearOptions {
            max_iter: 0,
            ..default_opts()
        };
        let err = gauss_newton(vec![0.0], opts, |_| Ok(NormalEquations::new(1))).unwrap_err();
        assert!(matches!(err, NonlinearError::Solver(_)));
    }

    #[test]
    fn non_finite_tol_rel_is_error() {
        let opts = NonlinearOptions {
            tol_rel: f64::NAN,
            ..default_opts()
        };
        let err = gauss_newton(vec![0.0], opts, |_| Ok(NormalEquations::new(1))).unwrap_err();
        assert!(matches!(err, NonlinearError::Solver(_)));
    }

    #[test]
    fn negative_tol_chi2_is_error() {
        let opts = NonlinearOptions {
            tol_chi2_rel: -1.0,
            ..default_opts()
        };
        let err = gauss_newton(vec![0.0], opts, |_| Ok(NormalEquations::new(1))).unwrap_err();
        assert!(matches!(err, NonlinearError::Solver(_)));
    }

    #[test]
    fn non_finite_initial_is_error() {
        let err = gauss_newton(vec![f64::NAN], default_opts(), |_| {
            Ok(NormalEquations::new(1))
        })
        .unwrap_err();
        assert!(matches!(err, NonlinearError::Solver(_)));
    }

    #[test]
    fn dimension_mismatch_is_error() {
        // Initial has 2 params but assembler returns NormalEquations::new(1)
        let err = gauss_newton(vec![0.0, 0.0], default_opts(), |_| {
            Ok(NormalEquations::new(1))
        })
        .unwrap_err();
        assert!(matches!(err, NonlinearError::Solver(_)));
    }

    #[test]
    fn did_not_converge_after_max_iter() {
        // Residual 100.0 – x starting from x=0 with max_iter=1 won't converge
        let opts = NonlinearOptions {
            max_iter: 1,
            tol_rel: 0.0,
            tol_chi2_rel: 0.0,
        };
        let err = gauss_newton(vec![0.0], opts, |params| {
            let mut ne = NormalEquations::new(1);
            ne.add_row(&[(0, 1.0)], 100.0 - params[0], 0.1).unwrap();
            Ok(ne)
        })
        .unwrap_err();
        assert!(matches!(err, NonlinearError::DidNotConverge(1)));
    }
}
