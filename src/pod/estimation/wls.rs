//! # Weighted least-squares normal equations
//!
//! ## Scientific scope
//!
//! Batch POD commonly reduces to solving a symmetric positive-definite
//! normal system assembled from many scalar observation equations. This
//! module implements that algebraic core for short deterministic estimation
//! arcs where the caller has already formed residuals and design-matrix
//! coefficients.
//!
//! The numerical method is dense Cholesky factorization, which is
//! appropriate for the small-to-medium problems exercised by the current
//! workspace. It is not intended as a sparse or square-root information
//! solver for very large operational networks.
//!
//! ## Technical scope
//!
//! The public surface centers on `NormalEquations`, `WlsResult`, and
//! `WlsSolverError`. Callers stream rows into the accumulator, then solve
//! for the parameter update and inverse normal matrix in raw `f64` form.
//!
//! Observation modelling, parameter semantics, and any typed quantity
//! handling remain outside this module. It only owns the linear algebra and
//! bookkeeping of the solve.
//!
//! ## References
//!
//! - Tapley, B. D., Schutz, B. E., & Born, G. H. (2004). Statistical Orbit
//!   Determination. Elsevier Academic Press.
//! - Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications
//!   (4th ed.). Microcosm Press.
use crate::coordinates::frames::GCRS;
use affn::matrix3::{FrameMatrix3, SymmetricFrameMatrix3};
use faer::linalg::solvers::Solve;
use faer::{Mat, Side};
use principia::StateCovariance;
use thiserror::Error;

/// Errors emerging from WLS assembly or solve.
#[derive(Debug, Error)]
pub enum WlsSolverError {
    /// Cholesky failed (matrix indefinite / rank deficient).
    #[error("normal equations not positive definite: {0}")]
    NotPositiveDefinite(String),
    /// Index out of bounds for the configured number of parameters.
    #[error("invalid parameter index {index} (n_params = {n_params})")]
    InvalidIndex {
        /// Offending index.
        index: usize,
        /// Total parameter count.
        n_params: usize,
    },
    /// Wraps an underlying propagation/STM/IO failure that prevents row assembly.
    #[error("solver row assembly failed: {0}")]
    Other(String),
}

impl WlsSolverError {
    /// Convenience constructor for [`WlsSolverError::Other`].
    pub fn other<S: Into<String>>(msg: S) -> Self {
        Self::Other(msg.into())
    }
}

/// Accumulator for the normal equations.
#[derive(Debug, Clone)]
pub struct NormalEquations {
    n: usize,
    /// Upper triangle is filled; the matrix is symmetrised on solve.
    pub n_matrix: Mat<f64>,
    /// Right-hand side, shape (n_params, 1).
    pub b: Mat<f64>,
    /// Weighted χ²: Σ (rᵢ / σᵢ)².
    pub chi2: f64,
    /// Number of measurement rows accumulated.
    pub n_obs: usize,
}

impl NormalEquations {
    /// Allocate for `n_params` total parameters.
    pub fn new(n_params: usize) -> Self {
        Self {
            n: n_params,
            n_matrix: Mat::zeros(n_params, n_params),
            b: Mat::zeros(n_params, 1),
            chi2: 0.0,
            n_obs: 0,
        }
    }

    /// Total parameter count.
    pub fn n_params(&self) -> usize {
        self.n
    }

    /// Add one weighted measurement row.
    ///
    /// `partials` is a sparse vector of `(index, value)` pairs, `residual`
    /// is the measurement minus prediction, and `sigma` is its standard
    /// deviation (must be > 0).
    pub fn add_row(
        &mut self,
        partials: &[(usize, f64)],
        residual: f64,
        sigma: f64,
    ) -> Result<(), WlsSolverError> {
        if sigma <= 0.0 {
            return Err(WlsSolverError::NotPositiveDefinite(format!(
                "sigma must be > 0 (got {sigma})"
            )));
        }
        let w = 1.0 / (sigma * sigma);
        for &(i, hi) in partials {
            if i >= self.n {
                return Err(WlsSolverError::InvalidIndex {
                    index: i,
                    n_params: self.n,
                });
            }
            // Accumulate into the upper triangle; symmetrise at solve time.
            self.b[(i, 0)] += w * hi * residual;
            for &(j, hj) in partials {
                if j < i {
                    continue;
                }
                self.n_matrix[(i, j)] += w * hi * hj;
            }
        }
        self.chi2 += w * residual * residual;
        self.n_obs += 1;
        Ok(())
    }

    /// Solve the normal equations and return parameter update + covariance.
    pub fn solve(self) -> Result<WlsResult, WlsSolverError> {
        // Symmetrise.
        let mut a = self.n_matrix.clone();
        for i in 0..self.n {
            for j in i + 1..self.n {
                let v = a[(i, j)];
                a[(j, i)] = v;
            }
        }

        let llt = a
            .as_ref()
            .llt(Side::Lower)
            .map_err(|e| WlsSolverError::NotPositiveDefinite(format!("{e:?}")))?;

        // Δp = N⁻¹ b.
        let dp = llt.solve(&self.b);

        // Covariance = N⁻¹ via solving against identity.
        let mut id = Mat::<f64>::zeros(self.n, self.n);
        for i in 0..self.n {
            id[(i, i)] = 1.0;
        }
        let cov = llt.solve(&id);

        let mut update = vec![0.0; self.n];
        for i in 0..self.n {
            update[i] = dp[(i, 0)];
        }
        let mut cov_arr = vec![vec![0.0f64; self.n]; self.n];
        for i in 0..self.n {
            for j in 0..self.n {
                cov_arr[i][j] = cov[(i, j)];
            }
        }
        Ok(WlsResult {
            update,
            covariance: cov_arr,
            chi2: self.chi2,
            n_obs: self.n_obs,
            n_params: self.n,
        })
    }
}

/// Output of one WLS iteration.
#[derive(Debug, Clone)]
pub struct WlsResult {
    /// Parameter update vector, length `n_params`.
    pub update: Vec<f64>,
    /// Posterior covariance, `n_params × n_params` (row-major).
    pub covariance: Vec<Vec<f64>>,
    /// Weighted χ² of the prefit residuals used in this iteration.
    pub chi2: f64,
    /// Number of observation rows accumulated.
    pub n_obs: usize,
    /// Total parameter count.
    pub n_params: usize,
}

impl WlsResult {
    /// Reduced χ² = χ² / max(1, n_obs − n_params).
    pub fn reduced_chi2(&self) -> f64 {
        let dof = self.n_obs.saturating_sub(self.n_params).max(1) as f64;
        self.chi2 / dof
    }

    /// Return the posterior covariance as a typed [`StateCovariance<GCRS>`]
    /// when `n_params == 6` (position + velocity solution).
    ///
    /// Returns `None` for any other parameter count.
    pub fn to_state_covariance(&self) -> Option<StateCovariance<GCRS>> {
        if self.n_params != 6 {
            return None;
        }
        let c = &self.covariance;
        let rr = SymmetricFrameMatrix3::<GCRS>::from_upper([
            [c[0][0], c[0][1], c[0][2]],
            [c[1][0], c[1][1], c[1][2]],
            [c[2][0], c[2][1], c[2][2]],
        ]);
        let rv = FrameMatrix3::<GCRS>::from_array([
            [c[0][3], c[0][4], c[0][5]],
            [c[1][3], c[1][4], c[1][5]],
            [c[2][3], c[2][4], c[2][5]],
        ]);
        let vv = SymmetricFrameMatrix3::<GCRS>::from_upper([
            [c[3][3], c[3][4], c[3][5]],
            [c[4][3], c[4][4], c[4][5]],
            [c[5][3], c[5][4], c[5][5]],
        ]);
        Some(StateCovariance::<GCRS>::from_blocks(rr, rv, vv))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solves_simple_2d_system() {
        // Fit y = a + b*x to (0,1), (1,2), (2,3) with σ=1.
        // Truth: a=1, b=1.
        let mut ne = NormalEquations::new(2);
        for (x, y) in [(0.0, 1.0), (1.0, 2.0), (2.0, 3.0)] {
            ne.add_row(&[(0, 1.0), (1, x)], y, 1.0).unwrap();
        }
        let r = ne.solve().unwrap();
        assert!((r.update[0] - 1.0).abs() < 1e-12);
        assert!((r.update[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn flags_singular_system() {
        // Two identical rows ⇒ rank deficient.
        let mut ne = NormalEquations::new(2);
        ne.add_row(&[(0, 1.0)], 1.0, 1.0).unwrap();
        ne.add_row(&[(0, 1.0)], 1.0, 1.0).unwrap();
        assert!(matches!(
            ne.solve(),
            Err(WlsSolverError::NotPositiveDefinite(_))
        ));
    }
}
