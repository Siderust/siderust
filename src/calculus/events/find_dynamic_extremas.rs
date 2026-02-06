// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use super::Culmination;
use crate::astro::nutation::corrected_ra_with_nutation;
use crate::astro::precession;
use crate::astro::sidereal::gast_fast;
use crate::astro::JulianDate;
use crate::calculus::math_core::root_finding;
use crate::coordinates::centers::*;
use crate::coordinates::frames;
use crate::coordinates::spherical::*;
use crate::targets::Target;
use qtty::*;

/// Convenience constants.
use core::f64::consts::PI;

/// Scan step: 20 min in days.  Adjust for the usual speed/robustness trade-off.
const STEP_DAYS: Days = Minutes::new(20.0).to::<Day>();

/// Finds every altitude extremum (**upper** and **lower** culmination) of a
/// *time-dependent* target in the half-open interval
/// `[jd_start, jd_end)`.
///
/// `get_equatorial(jd)` must return the target’s **geocentric equatorial**
/// coordinates (RA, Dec) **in degrees** at the supplied JD.
///
/// * `observer` – site coordinates (latitude, longitude) in the geocentric
///   ECEF system, **also in degrees**.
/// * Output is a chronologically ordered `Vec<Culmination>`.
///
/// ### Algorithm
/// 1. Define `H(jd)` – the hour angle of the target at the site (radians).  
/// 2. Scan the interval in `STEP_DAYS` chunks and watch the sign of  
///    `sin H` (upper culm.) and `sin(H−π)` (lower culm.).  
///    A sign change brackets a root.  
/// 3. Refine the midpoint with Newton–Raphson on  
///    `f(jd) = H(jd) − H_target`,    where `H_target` = `0` or `π`.  
///    The derivative `dH/dt` is estimated by a symmetric ±1 s finite-difference.  
/// 4. Collect all roots in `[jd_start, jd_end)` and sort them.
pub fn find_dynamic_extremas<F, U: LengthUnit>(
    get_equatorial: F,
    observer: &position::Geographic,
    jd_start: JulianDate,
    jd_end: JulianDate,
) -> Vec<Culmination>
where
    F: Fn(JulianDate) -> Target<Position<Geocentric, frames::EquatorialMeanJ2000, U>> + Copy,
{
    // ────────────────────────────────────────────────────────────
    // Helper: hour angle H(jd) [rad]
    // ────────────────────────────────────────────────────────────
    let hour_angle = |jd: JulianDate| -> Radians {
        let mean_of_date =
            precession::precess_from_j2000(get_equatorial(jd).get_position().clone(), jd);
        let ra_nut = corrected_ra_with_nutation(&mean_of_date.direction(), jd);
        let ra = ra_nut.to::<Radian>();
        let theta = gast_fast(jd).to::<Radian>() + observer.lon().to::<Radian>(); // local sidereal time
        (theta - ra).wrap_signed() // H in (−π, π]
    };

    // Scalar functions for Brent's method:
    // - For upper culmination (H=0): use H directly (well-behaved near 0)
    // - For lower culmination (H=π): use (H-π) wrapped (crosses 0 at H=π)
    // Using the hour angle directly instead of sin(H) improves convergence rate.
    let h_upper = |jd: JulianDate| -> f64 {
        hour_angle(jd).value() // H in (-π, π], crosses 0 at upper culmination
    };
    let h_lower = |jd: JulianDate| -> f64 {
        (hour_angle(jd) - Radians::new(PI)).wrap_signed().value() // crosses 0 at H=π
    };

    // ────────────────────────────────────────────────────────────
    // Coarse scan + refinement
    // ────────────────────────────────────────────────────────────
    let mut out = Vec::new();
    let mut jd0 = jd_start;

    let h0: Radians = hour_angle(jd0);
    // Pre-compute bracket values for Brent (hour angle based)
    let mut h0_val = h0.value();
    let mut h0_lower = (h0 - Radians::new(PI)).wrap_signed().value();

    while jd0 < jd_end {
        let jd1 = (jd0 + STEP_DAYS).min(jd_end);
        let h1 = hour_angle(jd1);
        let h1_val = h1.value();
        let h1_lower = (h1 - Radians::new(PI)).wrap_signed().value();

        // IMPORTANT: hour angle values returned by `wrap_signed()` have a discontinuity at ±π.
        // A naive sign-change check can therefore detect a false “root” when the value wraps.
        // Brent assumes continuity on the bracket, so we only refine when the function is
        // continuous over the step interval.
        let upper_continuous = (h1_val - h0_val).abs() < PI;
        let lower_continuous = (h1_lower - h0_lower).abs() < PI;

        // Upper culmination: H(jd) crosses 0
        if upper_continuous && h0_val * h1_val < 0.0 {
            if let Some(root_val) =
                root_finding::brent_with_values(jd0.value(), jd1.value(), h0_val, h1_val, |t| {
                    h_upper(JulianDate::new(t))
                })
            {
                let root = JulianDate::new(root_val);
                if root >= jd_start && root < jd_end {
                    out.push(Culmination::Upper { jd: root });
                }
            }
        }

        // Lower culmination: wrap_signed(H(jd) - π) crosses 0
        if lower_continuous && h0_lower * h1_lower < 0.0 {
            if let Some(root_val) = root_finding::brent_with_values(
                jd0.value(),
                jd1.value(),
                h0_lower,
                h1_lower,
                |t| h_lower(JulianDate::new(t)),
            ) {
                let root = JulianDate::new(root_val);
                if root >= jd_start && root < jd_end {
                    out.push(Culmination::Lower { jd: root });
                }
            }
        }

        jd0 = jd1;
        h0_val = h1_val;
        h0_lower = h1_lower;
    }

    // Stable chronological sort
    out.sort_by(|a, b| {
        let jd_a = match a {
            Culmination::Upper { jd } | Culmination::Lower { jd } => jd,
        };
        let jd_b = match b {
            Culmination::Upper { jd } | Culmination::Lower { jd } => jd,
        };
        jd_a.partial_cmp(jd_b).unwrap()
    });

    // Deduplicate: coarse scan can bracket the same root twice when the step
    // boundaries land very close to the crossing.
    const DEDUPE_EPS: f64 = 1e-6; // days (~0.0864 s)
    out.dedup_by(|a, b| match (a, b) {
        (Culmination::Upper { jd: jd_a }, Culmination::Upper { jd: jd_b })
        | (Culmination::Lower { jd: jd_a }, Culmination::Lower { jd: jd_b }) => {
            (jd_a.value() - jd_b.value()).abs() < DEDUPE_EPS
        }
        _ => false,
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::catalog::SIRIUS;
    use crate::observatories::ROQUE_DE_LOS_MUCHACHOS;

    fn approx_eq(a: &Culmination, b: &Culmination, epsilon: f64) -> bool {
        match (a, b) {
            (Culmination::Upper { jd: jd_a }, Culmination::Upper { jd: jd_b })
            | (Culmination::Lower { jd: jd_a }, Culmination::Lower { jd: jd_b }) => {
                (jd_a.value() - jd_b.value()).abs() <= epsilon
            }
            _ => false, // mismatched types: one Upper, one Lower
        }
    }

    #[test]
    fn test_find_sirius_dynamic() {
        let jd_start = JulianDate::new(2_460_677.0);
        let jd_end = jd_start + Days::new(1.0);

        let res =
            find_dynamic_extremas(|_| SIRIUS.target, &ROQUE_DE_LOS_MUCHACHOS, jd_start, jd_end);

        let expected_lower = Culmination::Lower {
            jd: JulianDate::new(2_460_677.050_00),
        };
        let expected_upper = Culmination::Upper {
            jd: JulianDate::new(2_460_677.548_60),
        };

        assert_eq!(res.len(), 2);
        approx_eq(&res[0], &expected_lower, 0.001);
        approx_eq(&res[1], &expected_upper, 0.001);
    }
}
