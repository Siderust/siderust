// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Stellar Altitude Window Periods
//!
//! Routines for finding time intervals where a **fixed star** (static ICRS
//! direction) is above, below, or within a given altitude range.
//!
//! ## Algorithm
//!
//! Unlike the uniform‑scan approach used for Sun and Moon, this module
//! exploits the fact that a star's altitude varies **sinusoidally** with
//! Earth's rotation (period = one sidereal day):
//!
//! 1. Precess J2000 coordinates to the midpoint of the search window
//!    (precession drift < 0.01″/day — one evaluation suffices).
//! 2. Solve `cos(H₀) = (sin(h) − sin(δ)sin(φ)) / (cos(δ)cos(φ))`
//!    to find threshold crossing hour angles analytically.
//! 3. Convert H₀ to Mjd crossing times via the GST rate.
//! 4. Refine each predicted crossing with Brent's method on the
//!    **full‑precision** evaluator (precession + nutation + GAST).
//! 5. Label crossings and assemble periods via [`math_core::intervals`].
//!
//! ## Performance
//!
//! The analytical approach evaluates the full‑precision altitude function
//! ~10–20 times per day (Brent refinement + crossing classification),
//! compared to ~144 for a 10‑minute uniform scan.  For week‑long searches
//! this yields a **7–10× speedup** over the generic scan engine.
//!
//! ## Precision
//!
//! Crossing times are refined to the same tolerance as the generic engine
//! (default ≈ 86 µs).  Precession at the midpoint introduces < 25″ of RA
//! error over a full year, well within the ±15‑minute Brent bracket.

use crate::astro::JulianDate;
use crate::calculus::math_core::{intervals, root_finding};
use crate::coordinates::centers::ObserverSite;
use crate::time::{complement_within, ModifiedJulianDate, Period};
use qtty::*;

use super::star_equations::{StarAltitudeParams, ThresholdResult};

/// Type aliases.
type Mjd = ModifiedJulianDate;

// =============================================================================
// Constants
// =============================================================================

/// Half‑width of the Brent refinement bracket around each analytically
/// predicted crossing.  15 minutes is conservative; the analytical
/// prediction is typically accurate to < 10 seconds.
const BRACKET_HALF: Days = Quantity::new(15.0 / 1440.0);

/// Scan step for the fallback / validation scan path (10 minutes).
const SCAN_STEP_FALLBACK: Days = Quantity::new(10.0 / 1440.0);

// ---------------------------------------------------------------------------
// Fixed Star Altitude
// ---------------------------------------------------------------------------

/// Compute altitude of a fixed RA/Dec object from an observer site.
///
/// Uses precession + nutation-corrected RA, GAST, and the standard
/// equatorial→horizontal formula.
pub(crate) fn fixed_star_altitude_rad(
    mjd: ModifiedJulianDate,
    site: &crate::coordinates::centers::ObserverSite,
    ra_j2000: qtty::Degrees,
    dec_j2000: qtty::Degrees,
) -> qtty::Radians {
    use crate::astro::nutation::corrected_ra_with_nutation;
    use crate::astro::precession;
    use crate::astro::sidereal::{calculate_gst, calculate_lst};
    use crate::coordinates::frames::EquatorialMeanJ2000;
    use crate::coordinates::spherical;
    use qtty::Radian;
    let jd = mjd.to_julian_day();
    // Build a spherical position in EquatorialMeanJ2000
    let pos = spherical::Position::<
        crate::coordinates::centers::Geocentric,
        EquatorialMeanJ2000,
        qtty::LightYear,
    >::new(ra_j2000, dec_j2000, qtty::LightYears::new(1.0));

    // Precess J2000 → mean-of-date
    let mean_of_date = precession::precess_from_j2000(pos, jd);
    // Apply nutation correction to RA
    let ra_corrected = corrected_ra_with_nutation(&mean_of_date.direction(), jd);
    let dec = mean_of_date.polar();

    // Compute hour angle
    let gst = calculate_gst(jd);
    let lst = calculate_lst(gst, site.lon);
    let ha = (lst - ra_corrected).normalize().to::<Radian>();

    // Equatorial → horizontal altitude
    let lat = site.lat.to::<Radian>();
    let dec_rad = dec.to::<Radian>();
    let sin_alt = dec_rad.sin() * lat.sin() + dec_rad.cos() * lat.cos() * ha.cos();
    qtty::Radians::new(sin_alt.asin())
}

// =============================================================================
// Core: analytical bracket + Brent refinement
// =============================================================================

/// Build a full‑precision altitude closure for a fixed star.
#[inline]
fn make_star_fn<'a>(
    ra_j2000: Degrees,
    dec_j2000: Degrees,
    site: &'a ObserverSite,
) -> impl Fn(Mjd) -> Radians + 'a {
    move |t: Mjd| -> Radians { fixed_star_altitude_rad(t, site, ra_j2000, dec_j2000) }
}

/// Find all crossings of a single threshold, refined to full precision.
///
/// Returns chronologically sorted [`LabeledCrossing`]s and whether the
/// function is above threshold at `period.start`.
fn find_crossings_analytical(
    ra_j2000: Degrees,
    dec_j2000: Degrees,
    site: &ObserverSite,
    period: Period<Mjd>,
    threshold_rad: f64,
) -> (Vec<intervals::LabeledCrossing>, bool) {
    let thr = Radians::new(threshold_rad);
    let f = make_star_fn(ra_j2000, dec_j2000, site);

    let start_above = f(period.start).value() > threshold_rad;

    // Build analytical model at the period midpoint
    let mid_jd = JulianDate::new(
        0.5 * (period.start.to_julian_day().value() + period.end.to_julian_day().value()),
    );
    let params = StarAltitudeParams::from_j2000(ra_j2000, dec_j2000, site, mid_jd);

    match params.threshold_ha(threshold_rad) {
        ThresholdResult::AlwaysAbove => {
            // Verify at both endpoints with the full evaluator
            let end_above = f(period.end).value() > threshold_rad;
            if start_above && end_above {
                (Vec::new(), true)
            } else {
                // Precession drift moved the star across the threshold —
                // fall back to uniform scan for this edge case.
                let mut crossings = intervals::find_crossings(period, SCAN_STEP_FALLBACK, &f, thr);
                let labeled = intervals::label_crossings(&mut crossings, &f, thr);
                (labeled, start_above)
            }
        }
        ThresholdResult::NeverAbove => {
            let end_above = f(period.end).value() > threshold_rad;
            if !start_above && !end_above {
                (Vec::new(), false)
            } else {
                let mut crossings = intervals::find_crossings(period, SCAN_STEP_FALLBACK, &f, thr);
                let labeled = intervals::label_crossings(&mut crossings, &f, thr);
                (labeled, start_above)
            }
        }
        ThresholdResult::Crossings { h0_rad } => {
            let predicted = params.predict_crossings(period, h0_rad);

            // Shifted altitude: g(t) = f(t) − threshold
            let g = |t: Mjd| -> Radians { f(t) - thr };
            let g_day = |d: Days| -> Radians { g(Mjd::new(d.value())) };

            let mut refined: Vec<Mjd> = Vec::with_capacity(predicted.len());

            for (t_pred, _dir) in &predicted {
                let lo_v = (t_pred.value() - BRACKET_HALF.value()).max(period.start.value());
                let hi_v = (t_pred.value() + BRACKET_HALF.value()).min(period.end.value());

                if (hi_v - lo_v) < 1e-12 {
                    continue; // degenerate bracket at boundary
                }

                let lo = Mjd::new(lo_v);
                let hi = Mjd::new(hi_v);
                let g_lo = g(lo);
                let g_hi = g(hi);

                if g_lo.value() * g_hi.value() < 0.0 {
                    if let Some(root) = root_finding::brent_with_values(
                        Days::new(lo_v),
                        Days::new(hi_v),
                        g_lo,
                        g_hi,
                        g_day,
                    ) {
                        let rv = root.value();
                        if rv >= period.start.value() && rv <= period.end.value() {
                            refined.push(Mjd::new(rv));
                        }
                    }
                }
            }

            let labeled = intervals::label_crossings(&mut refined, &f, thr);
            (labeled, start_above)
        }
    }
}

// =============================================================================
// Public API
// =============================================================================

/// Finds periods when a fixed star is **above** `threshold` inside `period`.
///
/// Uses the analytical sinusoidal model for O(1) bracket discovery per
/// sidereal cycle, refined by Brent's method on the full‑precision
/// evaluator (precession + nutation + GAST → equatorial → horizontal).
///
/// # Arguments
///
/// * `ra_j2000`  — right ascension in J2000 equatorial coordinates
/// * `dec_j2000` — declination in J2000 equatorial coordinates
/// * `site`      — observer location on Earth
/// * `period`    — time window to search
/// * `threshold` — altitude threshold (e.g. 0° for the geometric horizon)
pub(crate) fn find_star_above_periods(
    ra_j2000: Degrees,
    dec_j2000: Degrees,
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    threshold: Degrees,
) -> Vec<Period<ModifiedJulianDate>> {
    let thr_rad = threshold.to::<Radian>().value();
    let thr = Radians::new(thr_rad);
    let f = make_star_fn(ra_j2000, dec_j2000, &site);

    let (labeled, start_above) =
        find_crossings_analytical(ra_j2000, dec_j2000, &site, period, thr_rad);

    intervals::build_above_periods(&labeled, period, start_above, &f, thr)
}

/// Finds periods when a fixed star is **below** `threshold` inside `period`.
///
/// Complement of [`find_star_above_periods`] within `period`.
pub(crate) fn find_star_below_periods(
    ra_j2000: Degrees,
    dec_j2000: Degrees,
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    threshold: Degrees,
) -> Vec<Period<ModifiedJulianDate>> {
    let above = find_star_above_periods(ra_j2000, dec_j2000, site, period, threshold);
    complement_within(period, &above)
}

/// Finds periods when a fixed star's altitude is within `[min, max]`.
///
/// Computed as `above(min) ∩ complement(above(max))`.
pub(crate) fn find_star_range_periods(
    ra_j2000: Degrees,
    dec_j2000: Degrees,
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    range: (Degrees, Degrees),
) -> Vec<Period<ModifiedJulianDate>> {
    let above_min = find_star_above_periods(ra_j2000, dec_j2000, site, period, range.0);
    let above_max = find_star_above_periods(ra_j2000, dec_j2000, site, period, range.1);
    let below_max = intervals::complement(period, &above_max);
    intervals::intersect(&above_min, &below_max)
}

// =============================================================================
// Scan‑based variants (for comparison / validation)
// =============================================================================

#[cfg(test)]
/// Finds periods where star is above threshold using the **generic
/// 10‑minute scan** + Brent approach.
///
/// Prefer [`find_star_above_periods`] for production use; this function
/// is provided for validation and performance comparison.
fn find_star_above_periods_scan(
    ra_j2000: Degrees,
    dec_j2000: Degrees,
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    threshold: Degrees,
) -> Vec<Period<ModifiedJulianDate>> {
    let thr = threshold.to::<Radian>();
    let f = make_star_fn(ra_j2000, dec_j2000, &site);
    intervals::above_threshold_periods(period, SCAN_STEP_FALLBACK, &f, thr)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn greenwich() -> ObserverSite {
        ObserverSite::new(
            Degrees::new(0.0),
            Degrees::new(51.4769),
            Quantity::<qtty::Meter>::new(0.0),
        )
    }

    fn roque() -> ObserverSite {
        ObserverSite::new(
            Degrees::new(-17.892),
            Degrees::new(28.762),
            Quantity::<qtty::Meter>::new(2396.0),
        )
    }

    fn period_7d() -> Period<ModifiedJulianDate> {
        Period::new(
            ModifiedJulianDate::new(60000.0),
            ModifiedJulianDate::new(60007.0),
        )
    }

    #[test]
    fn polaris_always_above_horizon() {
        let periods = find_star_above_periods(
            Degrees::new(37.95),
            Degrees::new(89.26),
            greenwich(),
            period_7d(),
            Degrees::new(0.0),
        );
        assert_eq!(periods.len(), 1, "Polaris should be continuously above");
        let dur = periods[0].duration_days();
        assert!(
            (dur - 7.0).abs() < 0.01,
            "should span full 7 days, got {}",
            dur
        );
    }

    #[test]
    fn sirius_rises_and_sets() {
        let periods = find_star_above_periods(
            Degrees::new(101.287),
            Degrees::new(-16.716),
            greenwich(),
            period_7d(),
            Degrees::new(0.0),
        );
        assert!(
            periods.len() >= 6 && periods.len() <= 8,
            "expected ~7 above‑horizon periods for Sirius, got {}",
            periods.len()
        );
        for p in &periods {
            let hours = p.duration_days() * 24.0;
            // First/last period may be truncated by the window boundary
            assert!(
                hours > 0.1 && hours < 18.0,
                "unreasonable above‑horizon duration: {} h",
                hours
            );
        }
    }

    #[test]
    fn never_visible_star() {
        let periods = find_star_above_periods(
            Degrees::new(0.0),
            Degrees::new(-80.0),
            greenwich(),
            period_7d(),
            Degrees::new(0.0),
        );
        assert!(periods.is_empty(), "Dec=−80° should never rise at 51°N");
    }

    #[test]
    fn above_plus_below_covers_full_period() {
        let site = greenwich();
        let period = period_7d();
        let ra = Degrees::new(101.287);
        let dec = Degrees::new(-16.716);

        let above = find_star_above_periods(ra, dec, site, period, Degrees::new(0.0));
        let below = find_star_below_periods(ra, dec, site, period, Degrees::new(0.0));

        let total_above: f64 = above.iter().map(|p| p.duration_days()).sum();
        let total_below: f64 = below.iter().map(|p| p.duration_days()).sum();
        assert!(
            (total_above + total_below - 7.0).abs() < 0.01,
            "above + below should cover 7 days, got {}",
            total_above + total_below
        );
    }

    #[test]
    fn range_periods_sirius() {
        let periods = find_star_range_periods(
            Degrees::new(101.287),
            Degrees::new(-16.716),
            roque(),
            period_7d(),
            (Degrees::new(10.0), Degrees::new(30.0)),
        );
        assert!(!periods.is_empty(), "should find range periods for Sirius");
    }

    #[test]
    fn analytical_matches_scan() {
        let site = roque();
        let period = Period::new(
            ModifiedJulianDate::new(60000.0),
            ModifiedJulianDate::new(60003.0),
        );
        let ra = Degrees::new(101.287);
        let dec = Degrees::new(-16.716);
        let thr = Degrees::new(0.0);

        let analytical = find_star_above_periods(ra, dec, site, period, thr);
        let scan = find_star_above_periods_scan(ra, dec, site, period, thr);

        assert_eq!(
            analytical.len(),
            scan.len(),
            "analytical and scan should find same count: {} vs {}",
            analytical.len(),
            scan.len()
        );

        let tolerance = 1.0 / 1440.0; // 1 minute
        for (a, s) in analytical.iter().zip(scan.iter()) {
            assert!(
                (a.start.value() - s.start.value()).abs() < tolerance,
                "start times differ by {} d",
                (a.start.value() - s.start.value()).abs()
            );
            assert!(
                (a.end.value() - s.end.value()).abs() < tolerance,
                "end times differ by {} d",
                (a.end.value() - s.end.value()).abs()
            );
        }
    }
}
