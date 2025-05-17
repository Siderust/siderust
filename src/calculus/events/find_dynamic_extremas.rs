use super::Culmination;
use crate::coordinates::{
    centers::Geocentric,
    frames::{ECEF, Equatorial},
    SphericalCoord,
};
use crate::units::{Degrees, Radians, JulianDay, Days};
use crate::astro::nutation::corrected_ra_with_nutation;
use crate::targets::Target;

/// Convenience constants.
use core::f64::consts::PI;

/// A **quick-and-dirty** JD → GAST approximation 
/// (error < 0.1″ for ±50 yr around 2025).  
/// Swap for a rigorous routine if you need sub-arc-second accuracy.
fn gast_fast(jd: JulianDay) -> Degrees {
    // Duffett-Smith & Zwart, *Practical Astronomy*, 4 th ed.
    let t = (jd.value() - 2_451_545.0) / 36_525.0;          // centuries since J2000.0
    let gast = 280.460_618_37
        + 360.985_647_366_29 * (jd.value() - 2_451_545.0)   // mean rotation
        + 0.000_387_933 * t * t
        - t * t * t / 38_710.0;
    Degrees::new(gast)
}

/// Scan step: 20 min in days.  Adjust for the usual speed/robustness trade-off.
const STEP_DAYS: Days        = Days::new(20.0 / 1_440.0); // 20 min
const NEWTON_EPS: f64        = 1e-9;    // ≈ 0.9 ms
const NEWTON_MAX_ITERS: usize = 15;

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
pub fn find_dynamic_extremas<F>(
    get_equatorial: F,
    observer: &SphericalCoord<Geocentric, ECEF>,
    jd_start: JulianDay,
    jd_end: JulianDay,
) -> Vec<Culmination>
where
    F: Fn(JulianDay) -> Target<SphericalCoord<Geocentric, Equatorial>> + Copy,
{
    // ────────────────────────────────────────────────────────────
    // Helper: hour angle H(jd) [rad]
    // ────────────────────────────────────────────────────────────
    let hour_angle = |jd: JulianDay| -> Radians {
        let ra_nut = corrected_ra_with_nutation(get_equatorial(jd).get_position(), jd);
        let ra    = ra_nut.to_radians();
        let theta = gast_fast(jd).to_radians() + observer.lon().to_radians(); // local sidereal time
        (theta - ra).wrap_pi() // H in (−π, π]
    };

    // ────────────────────────────────────────────────────────────
    // Newton-Raphson root-refinement for H(jd) − target = 0
    // ────────────────────────────────────────────────────────────
    let refine = |mut jd: JulianDay, target: Radians| -> Option<JulianDay> {
        for _ in 0..NEWTON_MAX_ITERS {
            let h  = (hour_angle(jd) - target).wrap_pi();
            if h.abs().as_f64() < 1e-12 {
                return Some(jd); // already precise enough
            }

            // Finite-difference dH/dt using ±1 s
            let dt  = Days::new(1.0 / 86_400.0);
            let dh  = (hour_angle(jd + dt) - hour_angle(jd - dt)).wrap_pi();
            let deriv = dh.as_f64() / (2.0 * dt.value()); // rad / day
            if deriv.abs() < 1e-10 {
                return None; // derivative ~ 0, avoid blow-up
            }

            // Newton step
            let delta = Days::new(h.as_f64() / deriv);
            jd -= delta;
            if delta.value().abs() < NEWTON_EPS {
                return Some(jd);
            }
        }
        None
    };

    // ────────────────────────────────────────────────────────────
    // Coarse scan + refinement
    // ────────────────────────────────────────────────────────────
    let mut out = Vec::new();
    let mut jd0 = jd_start;

    let h0      = hour_angle(jd0);
    let mut s0      = h0.sin();                       // for H = 0
    let mut s0_pi   = (h0 - Radians::new(PI)).sin();  // for H = π

    while jd0 < jd_end {
        let jd1 = JulianDay::new((jd0 + STEP_DAYS).value().min(jd_end.value()));
        let h1  = hour_angle(jd1);
        let s1      = h1.sin();
        let s1_pi   = (h1 - Radians::new(PI)).sin();

        // Sign change ⇒ a root lies in (jd0, jd1)
        if s0 * s1 < 0.0 {
            if let Some(root) = refine(
                JulianDay::new((jd0.value() + jd1.value()) * 0.5),
                Radians::new(0.0),
            ) {
                if root >= jd_start && root < jd_end {
                    out.push(Culmination::Upper { jd: root });
                }
            }
        }
        if s0_pi * s1_pi < 0.0 {
            if let Some(root) = refine(
                JulianDay::new((jd0.value() + jd1.value()) * 0.5),
                Radians::new(PI),
            ) {
                if root >= jd_start && root < jd_end {
                    out.push(Culmination::Lower { jd: root });
                }
            }
        }

        jd0   = jd1;
        s0    = s1;
        s0_pi = s1_pi;
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
    out
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::catalog::SIRIUS;
    use crate::observatories::ROQUE_DE_LOS_MUCHACHOS;


    fn approx_eq(a: &Culmination, b: &Culmination, epsilon: f64) -> bool {
        match (a, b) {
            (Culmination::Upper { jd: jd_a }, Culmination::Upper { jd: jd_b }) |
            (Culmination::Lower { jd: jd_a }, Culmination::Lower { jd: jd_b }) => {
                (jd_a.value() - jd_b.value()).abs() <= epsilon
            }
            _ => false, // mismatched types: one Upper, one Lower
        }
    }

    #[test]
    fn test_find_sirius_dynamic() {
        let jd_start = JulianDay::new(2_460_677.0);
        let jd_end = jd_start + Days::new(1.0);

        let res = find_dynamic_extremas(|_| SIRIUS.target, &ROQUE_DE_LOS_MUCHACHOS, jd_start, jd_end);
        println!("\n{:?}", res);

        let expected_lower = Culmination::Lower {jd: JulianDay::new(2_460_677.05000) };
        let expected_upper = Culmination::Upper { jd: JulianDay::new(2_460_677.54860) };

        assert_eq!(res.len(), 2);
        approx_eq(&res[0], &expected_lower, 0.001);
        approx_eq(&res[1], &expected_upper, 0.001);
    }
}
