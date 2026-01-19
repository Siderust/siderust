use crate::astro::JulianDate;
use qtty::Days;

/// Finite-difference step for derivative estimation (1 second)
pub const FD_STEP_DAYS: f64 = 1.0 / 86_400.0;

/// Newton-Raphson tolerance (≈ 0.86 µs)
const NEWTON_TOL_DAYS: f64 = 1e-11;

/// Maximum Newton iterations
const NEWTON_MAX_ITERS: usize = 20;

/// Newton-Raphson refinement for a scalar root of `f(x) - threshold = 0`.
///
/// This is a generic routine that takes an initial Julian date guess `jd`, a
/// scalar function `scalar_fn` evaluated at a given `JulianDate`, and the
/// `threshold` to subtract. Returns `Some(jd)` if a root is found, `None` on
/// failure.
pub fn refine_root_newton<F>(mut jd: JulianDate, scalar_fn: F, threshold: f64) -> Option<JulianDate>
where
    F: Fn(JulianDate) -> f64,
{
    for _ in 0..NEWTON_MAX_ITERS {
        let f = scalar_fn(jd) - threshold;

        // Check if already converged
        if f.abs() < 1e-12 {
            return Some(jd);
        }

        // Symmetric finite-difference derivative
        let dt = Days::new(FD_STEP_DAYS);
        let f_plus = scalar_fn(jd + dt);
        let f_minus = scalar_fn(jd - dt);
        let deriv = (f_plus - f_minus) / (2.0 * FD_STEP_DAYS); // units/day

        // Protect against near-zero derivative
        if deriv.abs() < 1e-14 {
            return None;
        }

        // Newton step
        let delta_days = f / deriv;
        jd -= Days::new(delta_days);

        if delta_days.abs() < NEWTON_TOL_DAYS {
            return Some(jd);
        }
    }
    None
}

/// Bisection tolerance (≈ 0.086 ms)
const BISECTION_TOL_DAYS: f64 = 1e-9;

/// Maximum bisection iterations
const BISECTION_MAX_ITERS: usize = 80;

/// Bisection fallback for a scalar root between two Julian dates. Requires
/// `scalar_fn(lo) - threshold` and `scalar_fn(hi) - threshold` to have opposite signs.
pub fn refine_root_bisection<F>(
    mut lo: JulianDate,
    mut hi: JulianDate,
    scalar_fn: F,
    threshold: f64,
) -> Option<JulianDate>
where
    F: Fn(JulianDate) -> f64,
{
    let mut f_lo = scalar_fn(lo) - threshold;
    let f_hi = scalar_fn(hi) - threshold;

    // Check endpoints
    if f_lo.abs() < 1e-12 {
        return Some(lo);
    }
    if f_hi.abs() < 1e-12 {
        return Some(hi);
    }

    // Verify bracket
    if f_lo * f_hi > 0.0 {
        return None;
    }

    for _ in 0..BISECTION_MAX_ITERS {
        let mid = JulianDate::new((lo.value() + hi.value()) * 0.5);
        let f_mid = scalar_fn(mid) - threshold;

        if f_mid.abs() < 1e-12 || (hi.value() - lo.value()).abs() < BISECTION_TOL_DAYS {
            return Some(mid);
        }

        if f_lo * f_mid <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }

    Some(JulianDate::new((lo.value() + hi.value()) * 0.5))
}

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

    refine_root_newton(guess, scalar_fn, threshold)
        .or_else(|| refine_root_bisection(jd_a, jd_b, scalar_fn, threshold))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn newton_returns_guess_when_root_already_exact() {
        let threshold = 25_000.0;
        let guess = JulianDate::new(threshold);

        let scalar = |jd: JulianDate| jd.value();
        let found = refine_root_newton(guess, scalar, threshold);

        assert_eq!(found, Some(guess));
    }

    #[test]
    fn newton_refines_linear_root() {
        let threshold = 2_458_850.0;
        let guess = JulianDate::new(threshold + 0.5);
        let scalar = |jd: JulianDate| jd.value();

        let found = refine_root_newton(guess, scalar, threshold).expect("Newton should converge");
        assert!((found.value() - threshold).abs() < 1e-9);
    }

    #[test]
    fn newton_returns_none_for_constant_function() {
        let threshold = 0.0;
        let guess = JulianDate::J2000;
        let scalar = |_jd: JulianDate| std::f64::consts::PI;

        assert!(refine_root_newton(guess, scalar, threshold).is_none());
    }

    #[test]
    fn newton_returns_none_after_max_iterations() {
        let calls = Cell::new(0usize);
        let scalar = |_jd: JulianDate| {
            let count = calls.get() + 1;
            calls.set(count);
            count as f64
        };

        let guess = JulianDate::J2000;
        assert!(refine_root_newton(guess, scalar, 0.0).is_none());
    }

    #[test]
    fn bisection_returns_lo_when_endpoint_matches_root() {
        let threshold = 12_345.678;
        let lo = JulianDate::new(threshold);
        let hi = JulianDate::new(threshold + 4.0);
        let scalar = |jd: JulianDate| jd.value();

        assert_eq!(refine_root_bisection(lo, hi, scalar, threshold), Some(lo));
    }

    #[test]
    fn bisection_returns_hi_when_endpoint_matches_root() {
        let threshold = 54_321.0;
        let lo = JulianDate::new(threshold - 4.0);
        let hi = JulianDate::new(threshold);
        let scalar = |jd: JulianDate| jd.value();

        assert_eq!(refine_root_bisection(lo, hi, scalar, threshold), Some(hi));
    }

    #[test]
    fn bisection_finds_root_between_brackets() {
        let threshold = 42.0;
        let lo = JulianDate::new(threshold - 1.0);
        let hi = JulianDate::new(threshold + 1.0);
        let scalar = |jd: JulianDate| jd.value();

        let found = refine_root_bisection(lo, hi, scalar, threshold).expect("bisection failed");
        assert!((found.value() - threshold).abs() < 1e-9);
    }

    #[test]
    fn bisection_returns_midpoint_after_max_iterations() {
        let calls = Cell::new(0usize);
        let threshold = 1.0e15 + 1.0;
        let lo = JulianDate::new(0.0);
        let hi = JulianDate::new(2.0e15);
        let scalar = |jd: JulianDate| {
            calls.set(calls.get() + 1);
            jd.value()
        };

        let found = refine_root_bisection(lo, hi, scalar, threshold).expect("bisection failed");
        // Ensure we iterated at least once and terminated early due to tolerance
        // being tighter than the floating precision at this magnitude.
        assert!(calls.get() > 2 && calls.get() <= BISECTION_MAX_ITERS + 2);
        assert!(found.value().is_finite());
    }

    #[test]
    fn bisection_returns_none_for_invalid_bracket() {
        let scalar = |_jd: JulianDate| 99.0;
        let lo = JulianDate::new(0.0);
        let hi = JulianDate::new(1.0);

        assert!(refine_root_bisection(lo, hi, scalar, 0.0).is_none());
    }

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
}
