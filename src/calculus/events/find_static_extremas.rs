use super::Culmination;
use crate::astro::{nutation::corrected_ra_with_nutation, sidereal, JulianDate};
use crate::coordinates::{
    spherical::ext::EcefPositionExt,
    centers::Geocentric, frames::Equatorial, spherical::position::Geographic,
    spherical::Position,
};
use crate::targets::Target;
use qtty::Simplify;
use qtty::*;

/// Returns **all** upper- and lower-culminations (meridian passages) that occur
/// in the interval `[jd_start, jd_end]` for a target whose right-ascension and
/// declination are treated as constant.
///
/// ### Method
/// 1. Let `h₀` be the desired hour angle: `0 °` for an **upper** culmination,
///    `180 °` for a **lower** culmination.  
/// 2. At `jd_start` compute  
///    ```text
///    k₀ ≈ round( (LST − (RA + h₀)) / 360° )
///    ```
///    This gives an initial transit index.  
/// 3. Refine each candidate (`k₀` and `k₀ + 1`) with Newton-Raphson
///    iterations on the equation  
///    `LST(jd) = RA + h₀ + 360 ° · k`  
///    until the update falls below `1 × 10⁻¹¹ days`.  
///    The earliest root ≥ `jd_start` is the first culmination.  
/// 4. From that instant, keep adding one **sidereal day**
///    (≈ 23 h 56 m 4 s) to generate the subsequent culminations until the
///    upper bound `jd_end` is exceeded.  
/// 5. Because RA and Dec are constant, the transit index `k` never needs to be
///    recomputed; a simple time-step suffices.
///
/// ### Parameters
/// * `target`   – Equatorial coordinates (RA/Dec) of the body.  
/// * `observer` – Geodetic longitude/latitude of the observing site.  
/// * `jd_start` – Inclusive lower bound of the search interval (Julian day).  
/// * `jd_end`   – Exclusive upper bound of the search interval (Julian day).
///
/// ### Returns
/// A chronologically ordered `Vec<Culmination>` containing every upper
/// (`Culmination::Upper`) and lower (`Culmination::Lower`) culmination in the
/// specified range.
///
/// ### Accuracy
/// * Newton iterations are capped at **8** and converge to better than
///   `≈ 1 μs` in elapsed time.
/// * The final list is sorted with a stable comparison on the Julian‐day
///   value.
///
/// ---
pub fn find_static_extremas<U: LengthUnit>(
    target: &Target<Position<Geocentric, Equatorial, U>>,
    observer: &Geographic,
    jd_start: JulianDate,
    jd_end: JulianDate,
) -> Vec<Culmination> {
    const MAX_ITER: usize = 8;
    const TOLERANCE: Days = Days::new(1e-11);
    type DegreesPerDay = qtty::Quantity<qtty::Per<Degree, Day>>;
    const D_LST_D_JD: DegreesPerDay = Quantity::new(360.985_647_366_29); // Earth rotation freq.

    let ra: Degrees = corrected_ra_with_nutation(&target.get_position().direction(), jd_start);

    let lon: Degrees = observer.lon();
    let mut out = Vec::new();

    // h₀ = 0° (upper) and 180° (lower)
    for &(h0, is_upper) in &[(Degrees::new(0.0), true), (Degrees::HALF_TURN, false)] {
        // Find first culmination ≥ jd_start
        let lst0: Degrees = sidereal::unmodded_lst(jd_start, lon);
        let raw_k = (lst0 - (ra + h0)) / Degrees::FULL_TURN;
        let k0: i32 = raw_k.simplify().value().round() as i32;

        // Try k₀ and k₀+1, keep the earliest ≥ jd_start
        let mut t_first = JulianDate::new(f64::INFINITY);
        for k in [k0, k0 + 1] {
            let mut t = jd_start;
            // Newton–Raphson refinement
            for _ in 0..MAX_ITER {
                let f: Degrees =
                    sidereal::unmodded_lst(t, lon) - (ra + (h0 + Degrees::FULL_TURN * k as f64));
                let dt: Days = (f / D_LST_D_JD).simplify();
                t -= dt;
                if dt.abs() < TOLERANCE {
                    break;
                }
            }
            if t >= jd_start && t < t_first {
                t_first = t;
            }
        }

        // Step through the remaining culminations
        if t_first.value().is_finite() {
            let mut t = t_first;
            while t <= jd_end + Days::new(1e-12) {
                // numerical guard
                let event = if is_upper {
                    Culmination::Upper { jd: t }
                } else {
                    Culmination::Lower { jd: t }
                };
                out.push(event);
                t += SIDEREAL_DAY.to::<Day>(); // next passage of the same type
            }
        }
    }

    // Chronological, stable sort
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
            (Culmination::Upper { jd: jd_a }, Culmination::Upper { jd: jd_b })
            | (Culmination::Lower { jd: jd_a }, Culmination::Lower { jd: jd_b }) => {
                (jd_a.value() - jd_b.value()).abs() <= epsilon
            }
            _ => false, // mismatched types: one Upper, one Lower
        }
    }

    #[test]
    fn test_find_sirius_static() {
        let jd_start = JulianDate::new(2460677.0);
        let jd_end = jd_start + Days::new(1.0);

        let results =
            find_static_extremas(&SIRIUS.target, &ROQUE_DE_LOS_MUCHACHOS, jd_start, jd_end);
        print!("\n{:?}", results);

        let expected_lower = Culmination::Lower {
            jd: JulianDate::new(2460677.05000),
        };
        let expected_upper = Culmination::Upper {
            jd: JulianDate::new(2460677.54860),
        };

        assert_eq!(results.len(), 2);
        approx_eq(&results[0], &expected_lower, 0.001);
        approx_eq(&results[1], &expected_upper, 0.001);
    }
}
