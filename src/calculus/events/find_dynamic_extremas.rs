// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use super::Culmination;
use crate::astro::nutation::corrected_ra_with_nutation;
use crate::astro::precession;
use crate::astro::JulianDate;
use crate::calculus::root_finding::brent;
use crate::coordinates::centers::*;
use crate::coordinates::frames;
use crate::coordinates::spherical::*;
use crate::targets::Target;
use qtty::*;

/// Convenience constants.
use core::f64::consts::PI;

/// A **quick-and-dirty** JD → GAST approximation
/// (error < 0.1″ for ±50 yr around 2025).  
/// Swap for a rigorous routine if you need sub-arc-second accuracy.
fn gast_fast(jd: JulianDate) -> Degrees {
    // Duffett-Smith & Zwart, *Practical Astronomy*, 4 th ed.
    let t = (jd.value() - 2_451_545.0) / 36_525.0; // centuries since J2000.0
    let gast = 280.460_618_37
        + 360.985_647_366_29 * (jd.value() - 2_451_545.0)   // mean rotation
        + 0.000_387_933 * t * t
        - t * t * t / 38_710.0;
    Degrees::new(gast)
}

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
    // For sign-change detection, we still use sin() which is monotonic in (-π/2, π/2)
    // but for refinement we use the hour angle directly for faster convergence
    let mut s0 = h0.sin(); // for H = 0
    let mut s0_pi = (h0 - Radians::new(PI)).sin(); // for H = π
    // Pre-compute bracket values for Brent (hour angle based)
    let mut h0_val = h0.value();
    let mut h0_lower = (h0 - Radians::new(PI)).wrap_signed().value();

    while jd0 < jd_end {
        let jd1 = (jd0 + STEP_DAYS).min(jd_end);
        let h1 = hour_angle(jd1);
        let s1 = h1.sin();
        let s1_pi = (h1 - Radians::new(PI)).sin();
        let h1_val = h1.value();
        let h1_lower = (h1 - Radians::new(PI)).wrap_signed().value();

        // Sign change in sin(H) ⇒ H crosses 0 (upper culmination)
        // Use Brent with pre-computed hour angle values for faster convergence
        if s0 * s1 < 0.0 {
            if let Some(root) = brent::refine_root_with_values(jd0, jd1, h0_val, h1_val, &h_upper, 0.0) {
                if root >= jd_start && root < jd_end {
                    out.push(Culmination::Upper { jd: root });
                }
            }
        }
        // Sign change in sin(H-π) ⇒ H crosses π (lower culmination)
        if s0_pi * s1_pi < 0.0 {
            if let Some(root) = brent::refine_root_with_values(jd0, jd1, h0_lower, h1_lower, &h_lower, 0.0) {
                if root >= jd_start && root < jd_end {
                    out.push(Culmination::Lower { jd: root });
                }
            }
        }

        jd0 = jd1;
        s0 = s1;
        s0_pi = s1_pi;
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
