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
        jd = jd - Days::new(delta_days);

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

    refine_root_newton(guess, |jd| scalar_fn(jd), threshold)
        .or_else(|| refine_root_bisection(jd_a, jd_b, |jd| scalar_fn(jd), threshold))
}
