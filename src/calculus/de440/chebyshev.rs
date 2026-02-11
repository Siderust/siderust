// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Chebyshev polynomial evaluation using the Clenshaw recurrence.
//!
//! Given coefficients `c[0..n]` and normalized argument `tau ∈ [-1, 1]`,
//! evaluates:
//!
//! ```text
//! f(tau) = c[0]*T_0(tau) + c[1]*T_1(tau) + ... + c[n-1]*T_{n-1}(tau)
//! ```
//!
//! where `T_k` are Chebyshev polynomials of the first kind.

/// Evaluate a Chebyshev polynomial using the Clenshaw algorithm.
///
/// - `coeffs`: Chebyshev coefficients `c[0..n]`.
/// - `tau`: normalized time in `[-1, 1]`.
///
/// Returns the polynomial value.
#[inline]
pub fn evaluate(coeffs: &[f64], tau: f64) -> f64 {
    let n = coeffs.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return coeffs[0];
    }

    // Clenshaw recurrence (backwards):
    //   b_{n+1} = 0
    //   b_n     = 0
    //   b_k     = 2*tau*b_{k+1} - b_{k+2} + c[k]   for k = n-1, n-2, ..., 1
    //   result  = tau*b_1 - b_2 + c[0]
    //
    // Note: The Chebyshev convention uses c[0]*T_0 + c[1]*T_1 + ...
    // where T_0(tau) = 1, so the final step is:
    //   result = c[0] + tau*b_1 - b_2

    let two_tau = 2.0 * tau;
    let mut b_kp1 = 0.0; // b_{k+1}
    let mut b_kp2 = 0.0; // b_{k+2}

    for k in (1..n).rev() {
        let b_k = two_tau * b_kp1 - b_kp2 + coeffs[k];
        b_kp2 = b_kp1;
        b_kp1 = b_k;
    }

    coeffs[0] + tau * b_kp1 - b_kp2
}

/// Evaluate the derivative of a Chebyshev polynomial.
///
/// Returns `df/dtau` — the derivative with respect to the normalized
/// argument `tau`. To get `df/dt` in physical time units, multiply by
/// `dtau/dt = 1/RADIUS`.
///
/// Uses the recurrence for Chebyshev derivatives:
/// ```text
/// T_0'(tau) = 0
/// T_1'(tau) = 1
/// T_n'(tau) = 2*T_{n-1}(tau) + 2*tau*T_{n-1}'(tau) - T_{n-2}'(tau)
/// ```
///
/// Computed jointly with the polynomial values via the "position+velocity"
/// Clenshaw variant.
#[inline]
pub fn evaluate_derivative(coeffs: &[f64], tau: f64) -> f64 {
    let n = coeffs.len();
    if n <= 1 {
        return 0.0;
    }

    // Modified Clenshaw that tracks both value and derivative.
    // We accumulate b_k for the value and db_k for the derivative.
    let two_tau = 2.0 * tau;
    let mut b_kp1 = 0.0;
    let mut b_kp2 = 0.0;
    let mut db_kp1 = 0.0;
    let mut db_kp2 = 0.0;

    for k in (1..n).rev() {
        let b_k = two_tau * b_kp1 - b_kp2 + coeffs[k];
        let db_k = two_tau * db_kp1 - db_kp2 + 2.0 * b_kp1;
        b_kp2 = b_kp1;
        b_kp1 = b_k;
        db_kp2 = db_kp1;
        db_kp1 = db_k;
    }

    // Final step: derivative of (c[0] + tau*b_1 - b_2) w.r.t. tau
    // = b_1 + tau*db_1 - db_2
    b_kp1 + tau * db_kp1 - db_kp2
}

/// Evaluate both the Chebyshev polynomial and its derivative in one pass.
///
/// Returns `(value, d_value_d_tau)`.
#[inline]
pub fn evaluate_both(coeffs: &[f64], tau: f64) -> (f64, f64) {
    let n = coeffs.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    if n == 1 {
        return (coeffs[0], 0.0);
    }

    let two_tau = 2.0 * tau;
    let mut b_kp1 = 0.0;
    let mut b_kp2 = 0.0;
    let mut db_kp1 = 0.0;
    let mut db_kp2 = 0.0;

    for k in (1..n).rev() {
        let b_k = two_tau * b_kp1 - b_kp2 + coeffs[k];
        let db_k = two_tau * db_kp1 - db_kp2 + 2.0 * b_kp1;
        b_kp2 = b_kp1;
        b_kp1 = b_k;
        db_kp2 = db_kp1;
        db_kp1 = db_k;
    }

    let value = coeffs[0] + tau * b_kp1 - b_kp2;
    let deriv = b_kp1 + tau * db_kp1 - db_kp2;
    (value, deriv)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constant() {
        // T_0(x) = 1, so [c] evaluates to c for any tau
        assert!((evaluate(&[3.5], 0.0) - 3.5).abs() < 1e-15);
        assert!((evaluate(&[3.5], 0.7) - 3.5).abs() < 1e-15);
        assert!((evaluate_derivative(&[3.5], 0.7)).abs() < 1e-15);
    }

    #[test]
    fn test_linear() {
        // c[0]*T_0 + c[1]*T_1 = 2 + 3*tau
        let coeffs = [2.0, 3.0];
        assert!((evaluate(&coeffs, 0.0) - 2.0).abs() < 1e-15);
        assert!((evaluate(&coeffs, 1.0) - 5.0).abs() < 1e-15);
        assert!((evaluate(&coeffs, -1.0) - (-1.0)).abs() < 1e-15);
        assert!((evaluate_derivative(&coeffs, 0.0) - 3.0).abs() < 1e-15);
        assert!((evaluate_derivative(&coeffs, 1.0) - 3.0).abs() < 1e-15);
    }

    #[test]
    fn test_quadratic() {
        // c[0]*T_0 + c[1]*T_1 + c[2]*T_2 = 1 + 0*tau + 2*(2*tau^2 - 1)
        // = 1 + 4*tau^2 - 2 = -1 + 4*tau^2
        let coeffs = [1.0, 0.0, 2.0];
        assert!((evaluate(&coeffs, 0.0) - (-1.0)).abs() < 1e-14);
        assert!((evaluate(&coeffs, 1.0) - 3.0).abs() < 1e-14);
        // derivative: 8*tau
        assert!((evaluate_derivative(&coeffs, 0.5) - 4.0).abs() < 1e-14);
    }

    #[test]
    fn test_evaluate_both_matches() {
        let coeffs = [1.0, 2.0, 3.0, 4.0, 5.0];
        let tau = 0.37;
        let (val, deriv) = evaluate_both(&coeffs, tau);
        assert!((val - evaluate(&coeffs, tau)).abs() < 1e-14);
        assert!((deriv - evaluate_derivative(&coeffs, tau)).abs() < 1e-14);
    }
}
