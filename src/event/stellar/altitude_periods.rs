// SPDX-License-Identifier: AGPL-3.0-only
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
//!    (precession drift < 0.01″/day, one evaluation suffices).
//! 2. Solve `cos(H₀) = (sin(h) − sin(δ)sin(φ)) / (cos(δ)cos(φ))`
//!    to find threshold crossing hour angles analytically.
//! 3. Convert H₀ to Mjd crossing times via the GST rate.
//! 4. Refine each predicted crossing with Brent's method on the
//!    **full‑precision** evaluator (precession + nutation + GAST).
//! 5. Label crossings and assemble periods with the shared event-search
//!    primitives.
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

#![allow(dead_code)]

use crate::astro::apparent::CorrectionPolicy;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::event::altitude::search::{
    InternalSearchConfig, CROSSING_DEDUPE_EPS, DEFAULT_SCAN_STEP,
};
use crate::event::altitude::{CrossingDirection, CrossingEvent};
use crate::event::search::{intervals, periods as threshold_periods, scan_fallback};
use crate::qtty::*;
use crate::time::JulianDate;
use crate::time::{Interval, ModifiedJulianDate};
use alloc::vec::Vec;

use super::star_equations::{StarAltitudeParams, ThresholdResult};

/// Type aliases.
type Mjd = ModifiedJulianDate;

#[inline]
fn opposite_sign(a: f64, b: f64) -> bool {
    a.signum() * b.signum() < 0.0
}

// =============================================================================
// Constants
// =============================================================================

/// Half‑width of the Brent refinement bracket around each analytically
/// predicted crossing.  15 minutes is conservative; the analytical
/// prediction is typically accurate to < 10 seconds.
const BRACKET_HALF: Days = Quantity::new(15.0 / 1440.0);

/// Extra span searched on either side of the caller's window. This makes an
/// analytical prediction just outside the window able to refine to a precise
/// crossing on the boundary.
const PREDICTION_GUARD: Days = Quantity::new(30.0 / 1440.0);

// ---------------------------------------------------------------------------
// Fixed Star Altitude
// ---------------------------------------------------------------------------

/// Compute altitude of a fixed RA/Dec object from an observer site.
///
/// Uses IAU 2006 precession + IAU 2000B nutation (via NPB matrix),
/// ERA-based GMST, and the standard equatorial→horizontal formula.
pub(crate) fn fixed_star_altitude_rad(
    mjd: ModifiedJulianDate,
    site: &crate::coordinates::centers::Geodetic<crate::coordinates::frames::ECEF>,
    ra_j2000: crate::qtty::Degrees,
    dec_j2000: crate::qtty::Degrees,
) -> crate::qtty::Radians {
    fixed_star_altitude_rad_with_policy(mjd, site, ra_j2000, dec_j2000, CorrectionPolicy::APPARENT)
}

/// Compute altitude of a fixed RA/Dec object with an explicit correction policy.
pub(crate) fn fixed_star_altitude_rad_with_policy(
    mjd: ModifiedJulianDate,
    site: &crate::coordinates::centers::Geodetic<crate::coordinates::frames::ECEF>,
    ra_j2000: crate::qtty::Degrees,
    dec_j2000: crate::qtty::Degrees,
    policy: CorrectionPolicy,
) -> crate::qtty::Radians {
    use crate::astro::earth_rotation::jd_ut1_from_tt_eop;
    use crate::astro::nutation::nutation_iau2000b;
    use crate::astro::precession::precession_nutation_matrix;
    use crate::astro::sidereal::gast_iau2006;
    use crate::coordinates::transform::AstroContext;
    use crate::qtty::Radian;

    let jd: JulianDate = mjd.to::<crate::JD>();

    // Convert J2000 RA/Dec → unit Cartesian
    let ra_rad = ra_j2000.to::<Radian>().value();
    let dec_rad = dec_j2000.to::<Radian>().value();
    let (sin_dec, cos_dec) = dec_rad.sin_cos();
    let (sin_ra, cos_ra) = ra_rad.sin_cos();
    let x_j2000 = cos_dec * cos_ra;
    let y_j2000 = cos_dec * sin_ra;
    let z_j2000 = sin_dec;

    // Apply full NPB matrix (IAU 2006 precession + IAU 2000B nutation)
    // for the apparent/default pipeline. Geometric-only queries keep the
    // catalogue direction and only rotate it by Earth orientation below.
    let nut = nutation_iau2000b(jd);
    let [x_tod, y_tod, z_tod] = if policy == CorrectionPolicy::GEOMETRIC {
        [x_j2000, y_j2000, z_j2000]
    } else {
        let npb = precession_nutation_matrix(jd, nut.dpsi, nut.deps);
        npb.apply_array([x_j2000, y_j2000, z_j2000])
    };

    // Extract true-of-date RA and Dec
    let ra_tod = y_tod.atan2(x_tod);
    let r_xy = (x_tod * x_tod + y_tod * y_tod).sqrt();
    let dec_tod = z_tod.atan2(r_xy);

    // True-of-date RA/Dec requires apparent sidereal time (GAST).
    let ctx: AstroContext = AstroContext::default();
    let eop = ctx.eop_at_tt(jd);
    let jd_ut1 = jd_ut1_from_tt_eop(jd, &eop);
    let gast = gast_iau2006(jd_ut1, jd, nut.dpsi, nut.mean_obliquity);
    let lst_rad = gast.value() + site.lon.to::<Radian>().value();
    let ha = (lst_rad - ra_tod).rem_euclid(core::f64::consts::TAU);

    // Equatorial → horizontal altitude
    let lat = site.lat.to::<Radian>().value();
    let sin_alt = dec_tod.sin() * lat.sin() + dec_tod.cos() * lat.cos() * ha.cos();
    crate::qtty::Radians::new(sin_alt.asin())
}

// =============================================================================
// Core: analytical bracket + Brent refinement
// =============================================================================

/// Build a full‑precision altitude closure for a fixed star.
#[inline]
fn make_star_fn<'a>(
    ra_j2000: Degrees,
    dec_j2000: Degrees,
    site: &'a Geodetic<ECEF>,
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
    site: &Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Radians,
    opts: InternalSearchConfig,
) -> (Vec<intervals::LabeledCrossing>, bool) {
    let thr = threshold;
    let f = make_star_fn(ra_j2000, dec_j2000, site);
    let signal = |t: Mjd| f(t).sin();
    let threshold_sin = thr.sin();

    if period.end <= period.start {
        return (Vec::new(), false);
    }

    let start_above = signal(period.start) > threshold_sin;

    let generic = || {
        let (crossings, _, _) = crate::event::search::crossings::find_labelled_crossings(
            period,
            DEFAULT_SCAN_STEP,
            &signal,
            threshold_sin,
            opts,
        );
        (crossings, start_above)
    };

    if opts.uses_scan_baseline() {
        return generic();
    }

    // Build analytical model at the period midpoint
    let start_jd: JulianDate = period.start.to::<crate::JD>();
    let end_jd: JulianDate = period.end.to::<crate::JD>();
    let mid_jd = crate::time::JulianDate::new(((start_jd.raw() + end_jd.raw()) / 2.0).value());
    let equatorial_j2000 =
        crate::coordinates::spherical::direction::EquatorialMeanJ2000::new(ra_j2000, dec_j2000);
    let params = StarAltitudeParams::from_j2000(equatorial_j2000, site, mid_jd);

    match params.threshold_ha(thr) {
        // A midpoint model cannot prove that the precise, slowly evolving
        // signal stays on one side of a grazing threshold. Use the generic
        // precise engine for these classifications instead of pruning them.
        ThresholdResult::AlwaysAbove | ThresholdResult::NeverAbove => generic(),
        ThresholdResult::Crossings { h0 } => {
            // Predict across a guarded window so crossings shifted across a
            // query boundary by the approximate model are still considered.
            let guarded = Interval::new(
                ModifiedJulianDate::new((period.start.raw() - PREDICTION_GUARD).value()),
                ModifiedJulianDate::new((period.end.raw() + PREDICTION_GUARD).value()),
            );
            let predicted = params.predict_crossings(guarded, h0);

            // Near a grazing culmination, rising and setting predictions can
            // be closer than their refinement brackets. A same-sign bracket
            // could then hide both roots, so defer to the authoritative
            // generic engine for the complete window.
            if predicted
                .windows(2)
                .any(|pair| (pair[1].0.raw() - pair[0].0.raw()).abs() <= BRACKET_HALF * 2.0)
            {
                return generic();
            }

            let mut refined = Vec::with_capacity(predicted.len());
            let residual_tol = opts.chebyshev.max_residual;
            let time_tol = opts.time_tolerance.value().max(f64::EPSILON);

            for (t_pred, predicted_direction) in &predicted {
                let lo_raw = t_pred.raw() - BRACKET_HALF;
                let lo = crate::time::ModifiedJulianDate::new(
                    (if lo_raw >= period.start.raw() {
                        lo_raw
                    } else {
                        period.start.raw()
                    })
                    .value(),
                );
                let hi_raw = t_pred.raw() + BRACKET_HALF;
                let hi = crate::time::ModifiedJulianDate::new(
                    (if hi_raw <= period.end.raw() {
                        hi_raw
                    } else {
                        period.end.raw()
                    })
                    .value(),
                );

                if hi <= lo {
                    continue;
                }

                let g_lo = signal(lo) - threshold_sin;
                let g_hi = signal(hi) - threshold_sin;
                if !opposite_sign(g_lo, g_hi)
                    && g_lo.abs() > residual_tol
                    && g_hi.abs() > residual_tol
                {
                    return generic();
                }

                let Some(root_days) = scan_fallback::brent_f64(
                    lo.raw().value(),
                    hi.raw().value(),
                    g_lo,
                    g_hi,
                    |days| signal(ModifiedJulianDate::new(days)) - threshold_sin,
                    time_tol,
                    residual_tol,
                ) else {
                    return generic();
                };
                let root = ModifiedJulianDate::new(root_days);
                if root >= period.start && root <= period.end {
                    let direction = if g_lo <= residual_tol && g_hi > residual_tol {
                        1
                    } else if g_lo > residual_tol && g_hi <= residual_tol {
                        -1
                    } else {
                        return generic();
                    };
                    if direction != *predicted_direction {
                        return generic();
                    }
                    refined.push(intervals::LabeledCrossing { t: root, direction });
                }
            }

            refined.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(core::cmp::Ordering::Equal));
            refined.dedup_by(|a, b| (a.t.raw() - b.t.raw()).abs() < CROSSING_DEDUPE_EPS);

            // The precise signal must alternate sides in the direction the
            // analytical model predicted. Any discrepancy means the model was
            // not a safe bracket oracle for this window.
            let mut expected_above = start_above;
            for crossing in &refined {
                let expected_direction = if expected_above { -1 } else { 1 };
                if crossing.direction != expected_direction {
                    return generic();
                }
                expected_above = !expected_above;
            }

            (refined, start_above)
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
/// * `ra_j2000` , right ascension in J2000 equatorial coordinates
/// * `dec_j2000`, declination in J2000 equatorial coordinates
/// * `site`     , observer location on Earth
/// * `period`   , time window to search
/// * `threshold`, altitude threshold (e.g. 0° for the geometric horizon)
pub(crate) fn stellar_above_threshold_impl(
    ra_j2000: Degrees,
    dec_j2000: Degrees,
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: InternalSearchConfig,
) -> Vec<Interval<ModifiedJulianDate>> {
    let thr = threshold.to::<Radian>();
    let (labeled, start_above) =
        find_crossings_analytical(ra_j2000, dec_j2000, &site, period, thr, opts);

    threshold_periods::assemble_above_threshold_periods(&labeled, period, start_above)
}

/// Finds periods when a fixed star is **below** `threshold` inside `period`.
///
/// Complement of [`stellar_above_threshold_impl`] within `period`.
pub(crate) fn stellar_below_threshold_impl(
    ra_j2000: Degrees,
    dec_j2000: Degrees,
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: InternalSearchConfig,
) -> Vec<Interval<ModifiedJulianDate>> {
    let above = stellar_above_threshold_impl(ra_j2000, dec_j2000, site, period, threshold, opts);
    threshold_periods::complement_threshold_periods(period, &above)
}

/// Finds periods when a fixed star's altitude is within `[min, max]`.
///
/// Computed as `above(min) ∩ complement(above(max))`.
pub(crate) fn stellar_altitude_ranges_impl(
    ra_j2000: Degrees,
    dec_j2000: Degrees,
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    range: (Degrees, Degrees),
    opts: InternalSearchConfig,
) -> Vec<Interval<ModifiedJulianDate>> {
    let min = range.0.to::<Radian>();
    let max = range.1.to::<Radian>();
    let (min_crossings, start_above_min) =
        find_crossings_analytical(ra_j2000, dec_j2000, &site, period, min, opts);
    let (max_crossings, start_above_max) =
        find_crossings_analytical(ra_j2000, dec_j2000, &site, period, max, opts);
    threshold_periods::assemble_in_range_periods(
        &min_crossings,
        start_above_min,
        &max_crossings,
        start_above_max,
        period,
    )
}

/// Find precise rising and setting events for a fixed ICRS direction.
pub(crate) fn stellar_crossings_impl(
    ra_j2000: Degrees,
    dec_j2000: Degrees,
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: InternalSearchConfig,
) -> Vec<CrossingEvent> {
    let (crossings, _) = find_crossings_analytical(
        ra_j2000,
        dec_j2000,
        &site,
        period,
        threshold.to::<Radian>(),
        opts,
    );
    crossings
        .into_iter()
        .map(|crossing| CrossingEvent {
            mjd: crossing.t,
            direction: if crossing.direction > 0 {
                CrossingDirection::Rising
            } else {
                CrossingDirection::Setting
            },
        })
        .collect()
}

// =============================================================================
// Scan-based variants (for comparison / validation)
// =============================================================================

#[cfg(any(test, feature = "bench-internals"))]
pub(crate) fn stellar_above_threshold_scan_baseline(
    ra_j2000: Degrees,
    dec_j2000: Degrees,
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: crate::event::altitude::SearchOpts,
) -> Vec<Interval<ModifiedJulianDate>> {
    stellar_above_threshold_impl(
        ra_j2000,
        dec_j2000,
        site,
        period,
        threshold,
        InternalSearchConfig::scan_brent_baseline_config(opts),
    )
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::catalog::SIRIUS;
    use crate::coordinates::spherical::direction;
    use crate::event::altitude::{
        above_threshold, altitude_ranges, below_threshold, crossings, SearchOpts,
    };

    fn greenwich() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(
            Degrees::new(0.0),
            Degrees::new(51.4769),
            Quantity::<crate::qtty::Meter>::new(0.0),
        )
    }

    fn roque() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(
            Degrees::new(-17.892),
            Degrees::new(28.762),
            Quantity::<crate::qtty::Meter>::new(2396.0),
        )
    }

    fn period_7d() -> Interval<ModifiedJulianDate> {
        Interval::new(
            crate::time::ModifiedJulianDate::new(60000.0),
            crate::time::ModifiedJulianDate::new(60007.0),
        )
    }

    fn assert_intervals_close(
        actual: &[Interval<ModifiedJulianDate>],
        expected: &[Interval<ModifiedJulianDate>],
        tolerance: Days,
    ) {
        assert_eq!(actual.len(), expected.len(), "{actual:?} != {expected:?}");
        for (actual, expected) in actual.iter().zip(expected) {
            assert!((actual.start.raw() - expected.start.raw()).abs() <= tolerance);
            assert!((actual.end.raw() - expected.end.raw()).abs() <= tolerance);
        }
    }

    fn assert_crossings_close(
        actual: &[CrossingEvent],
        expected: &[CrossingEvent],
        tolerance: Days,
    ) {
        assert_eq!(actual.len(), expected.len(), "{actual:?} != {expected:?}");
        for (actual, expected) in actual.iter().zip(expected) {
            assert_eq!(actual.direction, expected.direction);
            assert!((actual.mjd.raw() - expected.mjd.raw()).abs() <= tolerance);
        }
    }

    #[test]
    fn polaris_always_above_horizon() {
        let periods = stellar_above_threshold_impl(
            Degrees::new(37.95),
            Degrees::new(89.26),
            greenwich(),
            period_7d(),
            Degrees::new(0.0),
            InternalSearchConfig::default(),
        );
        assert_eq!(periods.len(), 1, "Polaris should be continuously above");
        let dur = periods[0].end.raw() - periods[0].start.raw();
        assert!(
            (dur - Days::new(7.0)).abs() < Days::new(0.01),
            "should span full 7 days, got {}",
            dur
        );
    }

    #[test]
    fn sirius_rises_and_sets() {
        let periods = stellar_above_threshold_impl(
            Degrees::new(101.287),
            Degrees::new(-16.716),
            greenwich(),
            period_7d(),
            Degrees::new(0.0),
            InternalSearchConfig::default(),
        );
        assert!(
            periods.len() >= 6 && periods.len() <= 8,
            "expected ~7 above‑horizon periods for Sirius, got {}",
            periods.len()
        );
        for p in &periods {
            let hours = p.length().to::<Hour>();
            // First/last period may be truncated by the window boundary
            assert!(
                hours > Hours::new(0.1) && hours < Hours::new(18.0),
                "unreasonable above‑horizon duration: {} h",
                hours
            );
        }
    }

    #[test]
    fn never_visible_star() {
        let periods = stellar_above_threshold_impl(
            Degrees::new(0.0),
            Degrees::new(-80.0),
            greenwich(),
            period_7d(),
            Degrees::new(0.0),
            InternalSearchConfig::default(),
        );
        assert!(periods.is_empty(), "Dec=−80° should never rise at 51°N");
    }

    #[test]
    fn above_plus_below_covers_full_period() {
        let site = greenwich();
        let period = period_7d();
        let ra = Degrees::new(101.287);
        let dec = Degrees::new(-16.716);

        let above = stellar_above_threshold_impl(
            ra,
            dec,
            site,
            period,
            Degrees::new(0.0),
            InternalSearchConfig::default(),
        );
        let below = stellar_below_threshold_impl(
            ra,
            dec,
            site,
            period,
            Degrees::new(0.0),
            InternalSearchConfig::default(),
        );

        let total_above: Days = above.iter().map(|p| p.end.raw() - p.start.raw()).sum();
        let total_below: Days = below.iter().map(|p| p.end.raw() - p.start.raw()).sum();
        assert!(
            (total_above + total_below - Days::new(7.0)).abs() < Days::new(0.01),
            "above + below should cover 7 days, got {}",
            total_above + total_below
        );
    }

    #[test]
    fn range_periods_sirius() {
        let periods = stellar_altitude_ranges_impl(
            Degrees::new(101.287),
            Degrees::new(-16.716),
            roque(),
            period_7d(),
            (Degrees::new(10.0), Degrees::new(30.0)),
            InternalSearchConfig::default(),
        );
        assert!(!periods.is_empty(), "should find range periods for Sirius");
    }

    #[test]
    fn analytical_matches_scan() {
        let site = roque();
        let period = Interval::new(
            crate::time::ModifiedJulianDate::new(60000.0),
            crate::time::ModifiedJulianDate::new(60003.0),
        );
        let ra = Degrees::new(101.287);
        let dec = Degrees::new(-16.716);
        let thr = Degrees::new(0.0);

        let opts = crate::event::altitude::SearchOpts::default();
        let analytical = stellar_above_threshold_impl(
            ra,
            dec,
            site,
            period,
            thr,
            InternalSearchConfig::from_public_opts(opts),
        );
        let scan = stellar_above_threshold_scan_baseline(ra, dec, site, period, thr, opts);

        assert_eq!(
            analytical.len(),
            scan.len(),
            "analytical and scan should find same count: {} vs {}",
            analytical.len(),
            scan.len()
        );

        let tolerance = Minutes::new(1.0).to::<Day>();
        for (a, s) in analytical.iter().zip(scan.iter()) {
            assert!(
                (a.start.raw() - s.start.raw()).abs() < tolerance,
                "start times differ by {} d",
                (a.start.raw() - s.start.raw()).abs()
            );
            assert!(
                (a.end.raw() - s.end.raw()).abs() < tolerance,
                "end times differ by {} d",
                (a.end.raw() - s.end.raw()).abs()
            );
        }
    }

    #[test]
    fn public_icrs_and_star_dispatch_match_reference_for_all_event_kinds() {
        let site = roque();
        let target = direction::ICRS::from(&SIRIUS);
        let period = Interval::new(
            ModifiedJulianDate::new(60300.0),
            ModifiedJulianDate::new(60330.0),
        );
        let opts = SearchOpts {
            time_tolerance: Days::new(1e-6),
        };
        let internal = InternalSearchConfig::from_public_opts(opts);
        let reference = InternalSearchConfig::scan_brent_baseline_config(opts);
        let tolerance = opts.time_tolerance * 2.0;

        let icrs_above = above_threshold(&target, &site, period, Degrees::new(0.0), opts);
        let star_above = above_threshold(&SIRIUS, &site, period, Degrees::new(0.0), opts);
        let reference_above = stellar_above_threshold_impl(
            target.ra(),
            target.dec(),
            site,
            period,
            Degrees::new(0.0),
            reference,
        );
        assert_intervals_close(&icrs_above, &reference_above, tolerance);
        assert_intervals_close(&star_above, &icrs_above, tolerance);

        let icrs_below = below_threshold(&target, &site, period, Degrees::new(0.0), opts);
        let star_below = below_threshold(&SIRIUS, &site, period, Degrees::new(0.0), opts);
        let reference_below = stellar_below_threshold_impl(
            target.ra(),
            target.dec(),
            site,
            period,
            Degrees::new(0.0),
            reference,
        );
        assert_intervals_close(&icrs_below, &reference_below, tolerance);
        assert_intervals_close(&star_below, &icrs_below, tolerance);

        let icrs_ranges = altitude_ranges(
            &target,
            &site,
            period,
            Degrees::new(10.0),
            Degrees::new(30.0),
            opts,
        );
        let star_ranges = altitude_ranges(
            &SIRIUS,
            &site,
            period,
            Degrees::new(10.0),
            Degrees::new(30.0),
            opts,
        );
        let reference_ranges = stellar_altitude_ranges_impl(
            target.ra(),
            target.dec(),
            site,
            period,
            (Degrees::new(10.0), Degrees::new(30.0)),
            reference,
        );
        assert_intervals_close(&icrs_ranges, &reference_ranges, tolerance);
        assert_intervals_close(&star_ranges, &icrs_ranges, tolerance);

        let icrs_crossings = crossings(&target, &site, period, Degrees::new(0.0), opts);
        let star_crossings = crossings(&SIRIUS, &site, period, Degrees::new(0.0), opts);
        let reference_crossings = stellar_crossings_impl(
            target.ra(),
            target.dec(),
            site,
            period,
            Degrees::new(0.0),
            reference,
        );
        assert_crossings_close(&icrs_crossings, &reference_crossings, tolerance);
        assert_crossings_close(&star_crossings, &icrs_crossings, tolerance);

        // Ensure this test actually exercises the analytical configuration.
        assert!(!internal.uses_scan_baseline());
    }

    #[test]
    fn long_and_high_latitude_windows_match_reference_boundaries() {
        let target = direction::ICRS::from(&SIRIUS);
        let opts = SearchOpts::default();
        let tolerance = Days::new(5e-8);
        let sites_and_periods = [
            (
                roque(),
                Interval::new(
                    ModifiedJulianDate::new(60000.0),
                    ModifiedJulianDate::new(60365.0),
                ),
            ),
            (
                Geodetic::<ECEF>::new(Degrees::new(15.0), Degrees::new(78.0), Meters::new(20.0)),
                Interval::new(
                    ModifiedJulianDate::new(60300.0),
                    ModifiedJulianDate::new(60330.0),
                ),
            ),
        ];

        for (site, period) in sites_and_periods {
            let actual = above_threshold(&target, &site, period, Degrees::new(0.0), opts);
            let reference = stellar_above_threshold_scan_baseline(
                target.ra(),
                target.dec(),
                site,
                period,
                Degrees::new(0.0),
                opts,
            );
            assert_intervals_close(&actual, &reference, tolerance);
        }
    }

    #[test]
    fn grazing_culminations_are_not_pruned_by_midpoint_model() {
        let site = roque();
        let target = direction::ICRS::from(&SIRIUS);
        let period = Interval::new(
            ModifiedJulianDate::new(60000.0),
            ModifiedJulianDate::new(60001.0),
        );
        let mut sampled_max = Degrees::new(-90.0);
        for sample in 0..=288 {
            let t = ModifiedJulianDate::new(60000.0 + f64::from(sample) / 288.0);
            sampled_max = sampled_max
                .max(fixed_star_altitude_rad(t, &site, target.ra(), target.dec()).to::<Degree>());
        }

        for threshold in [
            sampled_max - Degrees::new(0.02),
            sampled_max + Degrees::new(0.02),
        ] {
            let actual = above_threshold(&target, &site, period, threshold, SearchOpts::default());
            let reference = stellar_above_threshold_scan_baseline(
                target.ra(),
                target.dec(),
                site,
                period,
                threshold,
                SearchOpts::default(),
            );
            assert_intervals_close(&actual, &reference, Days::new(5e-8));
        }
    }
}
