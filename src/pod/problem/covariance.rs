//! Parameter-vector covariance.
//!
//! State-only, frame/center-aware covariance lives in `siderust` as
//! `crate::astro::dynamics::covariance::StateCovariance` and is re-used
//! directly. This module adds [`ParameterCovariance`] for the *full*
//! estimator parameter vector, which mixes state, receiver clocks, biases,
//! scales, and ambiguities (see [`crate::pod::problem::parameter::ParameterOrdering`]).

use crate::pod::problem::parameter::ParameterOrdering;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Symmetric covariance matrix over the full estimator parameter vector.
///
/// Stored row-major in `data` with `params.len()` rows and columns. Symmetry
/// is *not* enforced at construction; callers should validate via
/// [`ParameterCovariance::validate_symmetric`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ParameterCovariance {
    /// Parameter ordering defining the row/column correspondence.
    pub params: ParameterOrdering,
    /// Row-major `n x n` data, where `n = params.len()`.
    pub data: Vec<f64>,
}

/// Errors when validating or transforming a parameter covariance.
#[derive(Debug, thiserror::Error)]
pub enum CovarianceError {
    /// Row/column count does not match the expected square dimension.
    #[error("dimension mismatch: expected {expected}x{expected}, got {actual} entries")]
    Dimension {
        /// Expected dimension (will be squared for matrix size).
        expected: usize,
        /// Actual number of elements provided.
        actual: usize,
    },
    /// Matrix is not symmetric within tolerance.
    #[error("not symmetric: |C[{i},{j}] - C[{j},{i}]| = {delta} > tol {tol}")]
    NotSymmetric {
        /// Row index.
        i: usize,
        /// Column index.
        j: usize,
        /// Absolute difference between symmetric pairs.
        delta: f64,
        /// Tolerance threshold that was exceeded.
        tol: f64,
    },
    /// Diagonal element is negative, violating positive semi-definiteness.
    #[error("not positive semi-definite: diagonal entry {i} = {value} < 0")]
    NotPsdDiagonal {
        /// Diagonal index.
        i: usize,
        /// Non-negative value found.
        value: f64,
    },
}

impl ParameterCovariance {
    /// Dimension of the parameter covariance matrix (number of parameters).
    pub fn n(&self) -> usize {
        self.params.len()
    }

    /// Verify dimensions match `n x n`.
    pub fn validate_dim(&self) -> Result<(), CovarianceError> {
        let n = self.n();
        let expected = n * n;
        if self.data.len() != expected {
            return Err(CovarianceError::Dimension {
                expected: n,
                actual: self.data.len(),
            });
        }
        Ok(())
    }

    /// Verify symmetry within `tol` element-wise.
    pub fn validate_symmetric(&self, tol: f64) -> Result<(), CovarianceError> {
        self.validate_dim()?;
        let n = self.n();
        for i in 0..n {
            for j in (i + 1)..n {
                let a = self.data[i * n + j];
                let b = self.data[j * n + i];
                let delta = (a - b).abs();
                if delta > tol {
                    return Err(CovarianceError::NotSymmetric { i, j, delta, tol });
                }
            }
        }
        Ok(())
    }

    /// Cheap PSD necessary-condition check: all diagonal entries non-negative.
    /// A full PSD test (eigendecomposition) belongs in the estimator crate.
    pub fn validate_diagonal_nonneg(&self) -> Result<(), CovarianceError> {
        self.validate_dim()?;
        let n = self.n();
        for i in 0..n {
            let v = self.data[i * n + i];
            if v < 0.0 {
                return Err(CovarianceError::NotPsdDiagonal { i, value: v });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pod::problem::parameter::{Parameter, ParameterKind};

    fn ord(n: usize) -> ParameterOrdering {
        ParameterOrdering {
            params: (0..n)
                .map(|_| Parameter {
                    kind: ParameterKind::DragScale,
                    initial_value: 1.0,
                    apriori_sigma: None,
                })
                .collect(),
        }
    }

    #[test]
    fn validates_symmetry() {
        let cov = ParameterCovariance {
            params: ord(2),
            data: vec![1.0, 0.5, 0.5, 2.0],
        };
        cov.validate_symmetric(1e-12).unwrap();
        cov.validate_diagonal_nonneg().unwrap();
    }

    #[test]
    fn detects_asymmetry() {
        let cov = ParameterCovariance {
            params: ord(2),
            data: vec![1.0, 0.5, 0.6, 2.0],
        };
        assert!(cov.validate_symmetric(1e-12).is_err());
    }

    #[test]
    fn detects_negative_diagonal() {
        let cov = ParameterCovariance {
            params: ord(1),
            data: vec![-1.0],
        };
        assert!(cov.validate_diagonal_nonneg().is_err());
    }
}
