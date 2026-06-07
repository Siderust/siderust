// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Solar daily threshold-crossing predictor with precise Brent validation.
//!
//! Uses a per-day analytic hour-angle model to bracket rise/set candidates,
//! refines each root against `sun_altitude_rad`, and falls back locally to the
//! generic Chebyshev-first engine or scan+Brent when the analytic model is
//! unreliable.

use crate::astro::earth_rotation::gmst_default;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::event::altitude::search::{
    InternalSearchConfig, CROSSING_DEDUPE_EPS, DIURNAL_CHEBY_SCAN_STEP, ONE_DAY, ROOT_INTERVAL_EPS,
};
use crate::event::altitude::{CrossingDirection, CrossingEvent};
use crate::event::search::crossings;
use crate::event::search::intervals::LabeledCrossing;
use crate::event::search::scan_fallback;
use crate::qtty::*;
use crate::time::{Interval, JulianDate, ModifiedJulianDate};

use super::altitude_periods::sun_altitude_rad;

const GRAZE_EPS: f64 = 1e-3;
/// Hour-angle rate used for candidate spacing (rad per mean solar day).
const HA_RATE_RAD_PER_DAY: f64 = std::f64::consts::TAU;

const BRACKET_RADII: [Days; 5] = [
    Minutes::new(15.0).to_const::<Day>(),
    Minutes::new(30.0).to_const::<Day>(),
    Hours::new(1.0).to_const::<Day>(),
    Hours::new(2.0).to_const::<Day>(),
    Hours::new(4.0).to_const::<Day>(),
];

/// Runtime counters for the solar daily predictor (tests/benches only).
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub(crate) struct SolarDailyDiagnostics {
    pub predicted_days: usize,
    pub candidate_crossings: usize,
    pub refined_crossings: usize,
    pub grazing_days: usize,
    pub bracket_failures: usize,
    pub expanded_brackets: usize,
    pub max_bracket_radius_days: f64,
    pub chebyshev_fallback_days: usize,
    pub scan_fallback_days: usize,
    pub precise_evaluations: usize,
}

enum DayThresholdCase {
    AlwaysBelow,
    AlwaysAbove,
    Crossings(f64),
}

/// Analytic daily solar state for one full MJD day (predictor only).
struct SolarDailyState {
    sin_dec: f64,
    cos_dec: f64,
    approx_transit: ModifiedJulianDate,
}

struct DayCrossingInput {
    state: SolarDailyState,
    clipped_day: Interval<ModifiedJulianDate>,
    thr_sin: f64,
    sin_lat: f64,
    cos_lat: f64,
}

/// Find labelled threshold crossings using the solar daily predictor.
pub(crate) fn solar_daily_crossings_impl(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: InternalSearchConfig,
) -> (Vec<LabeledCrossing>, SolarDailyDiagnostics) {
    let mut diagnostics = SolarDailyDiagnostics::default();
    if window.end <= window.start {
        return (Vec::new(), diagnostics);
    }

    if opts.uses_scan_baseline() {
        diagnostics.scan_fallback_days = 1;
        let signal = |t: ModifiedJulianDate| sun_altitude_rad(t, &site).sin();
        let mut search_diag = crossings::SearchDiagnostics::default();
        let labelled = scan_fallback::scan_labelled_crossings(
            window,
            opts.fallback_scan_step(DIURNAL_CHEBY_SCAN_STEP),
            &signal,
            threshold.to::<Radian>().sin(),
            opts.time_tolerance,
            opts.chebyshev.max_residual,
            &mut search_diag,
        );
        diagnostics.precise_evaluations += search_diag.precise_evaluations;
        diagnostics.refined_crossings = labelled.len();
        return (labelled, diagnostics);
    }

    if opts.disable_solar_daily_predictor {
        return chebyshev_labelled_crossings(site, window, threshold, opts, &mut diagnostics);
    }

    let thr_sin = threshold.to::<Radian>().sin();
    let lat = site.lat.to::<Radian>().value();
    let (sin_lat, cos_lat) = lat.sin_cos();

    let mut labelled = Vec::new();
    let mut full_day_start = floor_day(window.start);

    while full_day_start < window.end {
        let full_day_end = add_days(full_day_start, ONE_DAY);
        let full_day = Interval::new(full_day_start, full_day_end);

        let clipped_start = max_mjd(full_day_start, window.start);
        let clipped_end = min_mjd(full_day_end, window.end);
        if clipped_end <= clipped_start {
            full_day_start = full_day_end;
            continue;
        }
        let clipped_day = Interval::new(clipped_start, clipped_end);

        diagnostics.predicted_days += 1;
        let Some(state) = solar_daily_state(site, full_day) else {
            diagnostics.grazing_days += 1;
            labelled.extend(fallback_day(
                site,
                clipped_day,
                thr_sin,
                opts,
                &mut diagnostics,
            ));
            full_day_start = full_day_end;
            continue;
        };
        let day_crossings = predict_day_crossings(
            site,
            DayCrossingInput {
                state,
                clipped_day,
                thr_sin,
                sin_lat,
                cos_lat,
            },
            opts,
            &mut diagnostics,
        );
        labelled.extend(day_crossings);
        full_day_start = full_day_end;
    }

    sort_dedup_labelled(&mut labelled);
    diagnostics.refined_crossings = labelled.len();
    (labelled, diagnostics)
}

/// Public crossing events wrapper around the daily predictor.
pub(crate) fn solar_daily_crossing_events_impl(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: InternalSearchConfig,
) -> Vec<CrossingEvent> {
    let (labelled, _) = solar_daily_crossings_impl(site, window, threshold, opts);
    labelled
        .into_iter()
        .map(|lc| CrossingEvent {
            mjd: lc.t,
            direction: if lc.direction > 0 {
                CrossingDirection::Rising
            } else {
                CrossingDirection::Setting
            },
        })
        .collect()
}

fn predict_day_crossings(
    site: Geodetic<ECEF>,
    input: DayCrossingInput,
    opts: InternalSearchConfig,
    diagnostics: &mut SolarDailyDiagnostics,
) -> Vec<LabeledCrossing> {
    let DayCrossingInput {
        state,
        clipped_day,
        thr_sin,
        sin_lat,
        cos_lat,
    } = input;
    let denom = cos_lat * state.cos_dec;
    if denom.abs() < 1e-12 {
        diagnostics.grazing_days += 1;
        return fallback_day(site, clipped_day, thr_sin, opts, diagnostics);
    }

    let cos_h0 = (thr_sin - sin_lat * state.sin_dec) / denom;
    let case = if cos_h0 > 1.0 + GRAZE_EPS {
        DayThresholdCase::AlwaysBelow
    } else if cos_h0 < -1.0 - GRAZE_EPS {
        DayThresholdCase::AlwaysAbove
    } else if !(-1.0..=1.0).contains(&cos_h0)
        || (cos_h0 + 1.0).abs() < GRAZE_EPS
        || (cos_h0 - 1.0).abs() < GRAZE_EPS
    {
        diagnostics.grazing_days += 1;
        return fallback_day(site, clipped_day, thr_sin, opts, diagnostics);
    } else {
        DayThresholdCase::Crossings(cos_h0.clamp(-1.0, 1.0).acos())
    };

    match case {
        DayThresholdCase::AlwaysBelow => {
            if accept_no_crossing_day(site, clipped_day, thr_sin, false, opts, diagnostics) {
                Vec::new()
            } else {
                fallback_day(site, clipped_day, thr_sin, opts, diagnostics)
            }
        }
        DayThresholdCase::AlwaysAbove => {
            if accept_no_crossing_day(site, clipped_day, thr_sin, true, opts, diagnostics) {
                Vec::new()
            } else {
                fallback_day(site, clipped_day, thr_sin, opts, diagnostics)
            }
        }
        DayThresholdCase::Crossings(h0) => {
            let dt = h0 / HA_RATE_RAD_PER_DAY;
            let candidates = [
                (state.approx_transit.raw() - Days::new(dt), 1_i8),
                (state.approx_transit.raw() + Days::new(dt), -1_i8),
            ];

            let mut expected_candidates = 0usize;
            let mut refined = Vec::new();
            let mut failed = false;

            for (pred_raw, _expected_dir) in candidates {
                let pred = ModifiedJulianDate::new(pred_raw.value());
                if pred < clipped_day.start || pred > clipped_day.end {
                    continue;
                }

                expected_candidates += 1;
                diagnostics.candidate_crossings += 1;

                match refine_candidate(site, clipped_day, pred, thr_sin, opts, diagnostics) {
                    Some((root, direction)) => {
                        refined.push(LabeledCrossing {
                            t: root,
                            direction: direction as i32,
                        });
                    }
                    None => {
                        failed = true;
                        diagnostics.bracket_failures += 1;
                    }
                }
            }

            if failed || refined.len() != expected_candidates {
                return fallback_day(site, clipped_day, thr_sin, opts, diagnostics);
            }

            refined
        }
    }
}

fn solar_daily_state(
    site: Geodetic<ECEF>,
    full_day: Interval<ModifiedJulianDate>,
) -> Option<SolarDailyState> {
    let t_mid = midpoint(full_day);
    let jd: JulianDate = t_mid.to::<crate::JD>();
    let j2000_mjd = crate::J2000.to::<crate::MJD>();
    let n = (t_mid.raw() - j2000_mjd.raw()).value();

    let l_deg = 280.46646 + 0.98564736 * n;
    let g_deg = 357.52911 + 0.98560028 * n;
    let g = deg_to_rad(g_deg);
    let lambda = deg_to_rad(l_deg)
        + deg_to_rad(1.914602) * g.sin()
        + deg_to_rad(0.019993) * (2.0 * g).sin()
        + deg_to_rad(0.000289) * (3.0 * g).sin();
    let epsilon = deg_to_rad(23.439291 - 0.00000036 * n);

    let (sin_lambda, cos_lambda) = lambda.sin_cos();
    let (sin_eps, cos_eps) = epsilon.sin_cos();
    let ra = (cos_eps * sin_lambda).atan2(cos_lambda);
    let dec = (sin_eps * sin_lambda).asin();

    if !ra.is_finite() || !dec.is_finite() {
        return None;
    }

    let gmst = gmst_default(jd).value();
    let lst = normalize_angle(gmst + site.lon.to::<Radian>().value());
    let hour_angle = normalize_angle(lst - ra);
    let transit_offset_days = -hour_angle / HA_RATE_RAD_PER_DAY;
    let mut approx_transit = add_days(t_mid, Days::new(transit_offset_days));

    if approx_transit < full_day.start {
        approx_transit = add_days(approx_transit, ONE_DAY);
    } else if approx_transit > full_day.end {
        approx_transit = add_days(approx_transit, Days::new(-ONE_DAY.value()));
    }

    Some(SolarDailyState {
        sin_dec: dec.sin(),
        cos_dec: dec.cos(),
        approx_transit,
    })
}

fn accept_no_crossing_day(
    site: Geodetic<ECEF>,
    clipped_day: Interval<ModifiedJulianDate>,
    thr_sin: f64,
    expected_above: bool,
    opts: InternalSearchConfig,
    diagnostics: &mut SolarDailyDiagnostics,
) -> bool {
    let tol = residual_tol(opts);
    for t in [clipped_day.start, midpoint(clipped_day), clipped_day.end] {
        let r = solar_residual(site, t, thr_sin, diagnostics);
        if !r.is_finite() || r.abs() <= tol * 10.0 {
            return false;
        }
        if expected_above && r < 0.0 {
            return false;
        }
        if !expected_above && r > 0.0 {
            return false;
        }
    }
    true
}

fn refine_candidate(
    site: Geodetic<ECEF>,
    bounds: Interval<ModifiedJulianDate>,
    predicted: ModifiedJulianDate,
    thr_sin: f64,
    opts: InternalSearchConfig,
    diagnostics: &mut SolarDailyDiagnostics,
) -> Option<(ModifiedJulianDate, i8)> {
    let (lo, hi, f_lo, f_hi) =
        find_candidate_bracket(site, bounds, predicted, thr_sin, opts, diagnostics)?;

    if (hi.raw() - lo.raw()) < ROOT_INTERVAL_EPS {
        let root = lo;
        if root < bounds.start || root > bounds.end {
            return None;
        }
        let direction = classify_direction(site, root, thr_sin, diagnostics);
        return Some((root, direction));
    }

    let root_days = scan_fallback::brent_f64(
        lo.raw().value(),
        hi.raw().value(),
        f_lo,
        f_hi,
        |days| solar_residual(site, ModifiedJulianDate::new(days), thr_sin, diagnostics),
        opts.time_tolerance.value(),
        residual_tol(opts),
    )?;
    let root = ModifiedJulianDate::new(root_days);
    if root < bounds.start || root > bounds.end {
        return None;
    }

    let direction = classify_direction(site, root, thr_sin, diagnostics);
    Some((root, direction))
}

fn find_candidate_bracket(
    site: Geodetic<ECEF>,
    bounds: Interval<ModifiedJulianDate>,
    predicted: ModifiedJulianDate,
    thr_sin: f64,
    opts: InternalSearchConfig,
    diagnostics: &mut SolarDailyDiagnostics,
) -> Option<(ModifiedJulianDate, ModifiedJulianDate, f64, f64)> {
    let tol = residual_tol(opts);
    for (idx, radius) in BRACKET_RADII.iter().enumerate() {
        let lo = max_mjd(
            bounds.start,
            ModifiedJulianDate::new((predicted.raw() - *radius).value()),
        );
        let hi = min_mjd(
            bounds.end,
            ModifiedJulianDate::new((predicted.raw() + *radius).value()),
        );
        if hi <= lo {
            continue;
        }

        let f_lo = solar_residual(site, lo, thr_sin, diagnostics);
        let f_hi = solar_residual(site, hi, thr_sin, diagnostics);
        if !f_lo.is_finite() || !f_hi.is_finite() {
            continue;
        }

        if f_lo.abs() <= tol {
            diagnostics.expanded_brackets += 1;
            diagnostics.max_bracket_radius_days =
                radius.value().max(diagnostics.max_bracket_radius_days);
            return Some((lo, lo, f_lo, f_lo));
        }
        if f_hi.abs() <= tol {
            diagnostics.expanded_brackets += 1;
            diagnostics.max_bracket_radius_days =
                radius.value().max(diagnostics.max_bracket_radius_days);
            return Some((hi, hi, f_hi, f_hi));
        }
        if f_lo.signum() * f_hi.signum() < 0.0 {
            diagnostics.expanded_brackets += 1;
            diagnostics.max_bracket_radius_days =
                radius.value().max(diagnostics.max_bracket_radius_days);
            let _ = idx;
            return Some((lo, hi, f_lo, f_hi));
        }
    }
    None
}

fn classify_direction(
    site: Geodetic<ECEF>,
    root: ModifiedJulianDate,
    thr_sin: f64,
    diagnostics: &mut SolarDailyDiagnostics,
) -> i8 {
    let eps = Minutes::new(1.0).to::<Day>();
    let before = solar_residual(site, add_days(root, -eps), thr_sin, diagnostics);
    let after = solar_residual(site, add_days(root, eps), thr_sin, diagnostics);
    if before <= 0.0 && after > 0.0 {
        1
    } else if before >= 0.0 && after < 0.0 {
        -1
    } else if before < after {
        1
    } else {
        -1
    }
}

fn solar_residual(
    site: Geodetic<ECEF>,
    t: ModifiedJulianDate,
    thr_sin: f64,
    diagnostics: &mut SolarDailyDiagnostics,
) -> f64 {
    diagnostics.precise_evaluations += 1;
    sun_altitude_rad(t, &site).sin() - thr_sin
}

fn residual_tol(opts: InternalSearchConfig) -> f64 {
    opts.chebyshev.max_residual.max(crossings::POLY_ZERO_TOL)
}

fn fallback_day(
    site: Geodetic<ECEF>,
    day: Interval<ModifiedJulianDate>,
    thr_sin: f64,
    opts: InternalSearchConfig,
    diagnostics: &mut SolarDailyDiagnostics,
) -> Vec<LabeledCrossing> {
    diagnostics.chebyshev_fallback_days += 1;
    let (labelled, cheb_diag) = chebyshev_labelled_crossings_raw(site, day, thr_sin, opts);
    diagnostics.precise_evaluations += cheb_diag.precise_evaluations;
    if cheb_diag.fallback_segments > 0 {
        diagnostics.scan_fallback_days += 1;
    }
    labelled
}

fn chebyshev_labelled_crossings(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: InternalSearchConfig,
    diagnostics: &mut SolarDailyDiagnostics,
) -> (Vec<LabeledCrossing>, SolarDailyDiagnostics) {
    let thr_sin = threshold.to::<Radian>().sin();
    let (labelled, search_diag) = chebyshev_labelled_crossings_raw(site, window, thr_sin, opts);
    diagnostics.precise_evaluations += search_diag.precise_evaluations;
    diagnostics.refined_crossings = labelled.len();
    (labelled, *diagnostics)
}

pub(crate) fn chebyshev_labelled_crossings_raw(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    thr_sin: f64,
    opts: InternalSearchConfig,
) -> (Vec<LabeledCrossing>, crossings::SearchDiagnostics) {
    let signal = |t: ModifiedJulianDate| sun_altitude_rad(t, &site).sin();
    let (labelled, _, diag) =
        crossings::find_labelled_crossings(window, DIURNAL_CHEBY_SCAN_STEP, &signal, thr_sin, opts);
    (labelled, diag)
}

#[inline]
fn deg_to_rad(deg: f64) -> f64 {
    deg.to_radians()
}

#[inline]
fn normalize_angle(rad: f64) -> f64 {
    rad.rem_euclid(std::f64::consts::TAU)
}

#[inline]
fn floor_day(t: ModifiedJulianDate) -> ModifiedJulianDate {
    ModifiedJulianDate::new(t.raw().value().floor())
}

#[inline]
fn add_days(t: ModifiedJulianDate, delta: Days) -> ModifiedJulianDate {
    ModifiedJulianDate::new((t.raw() + delta).value())
}

#[inline]
fn max_mjd(a: ModifiedJulianDate, b: ModifiedJulianDate) -> ModifiedJulianDate {
    if a >= b {
        a
    } else {
        b
    }
}

#[inline]
fn min_mjd(a: ModifiedJulianDate, b: ModifiedJulianDate) -> ModifiedJulianDate {
    if a <= b {
        a
    } else {
        b
    }
}

#[inline]
fn midpoint(interval: Interval<ModifiedJulianDate>) -> ModifiedJulianDate {
    ModifiedJulianDate::new((interval.start.raw() + interval.end.raw()) / Days::new(2.0))
}

fn sort_dedup_labelled(crossings: &mut Vec<LabeledCrossing>) {
    crossings.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
    crossings.dedup_by(|a, b| (a.t.raw() - b.t.raw()).abs() < CROSSING_DEDUPE_EPS);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::event::search::intervals;
    use chrono::{TimeZone, Utc};

    fn greenwich() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.4769), Meters::new(0.0))
    }

    fn roque() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(
            Degrees::new(-17.892),
            Degrees::new(28.762),
            Meters::new(2396.0),
        )
    }

    fn utc_window(days: u32) -> Interval<ModifiedJulianDate> {
        let start = Utc.with_ymd_and_hms(2026, 1, 1, 0, 0, 0).single().unwrap();
        let end = start + chrono::Duration::days(days as i64);
        Interval::new(
            ModifiedJulianDate::from(start),
            ModifiedJulianDate::from(end),
        )
    }

    fn scan_labelled(
        site: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
    ) -> Vec<LabeledCrossing> {
        let scan_step = DIURNAL_CHEBY_SCAN_STEP;
        let thr = threshold.to::<Radian>();
        let f = |t: ModifiedJulianDate| sun_altitude_rad(t, &site);
        let mut crossings_t = intervals::find_crossings(window, scan_step, &f, thr);
        intervals::label_crossings(&mut crossings_t, &f, thr)
    }

    fn assert_labelled_close(
        actual: &[LabeledCrossing],
        expected: &[LabeledCrossing],
        threshold: &str,
    ) {
        assert_eq!(
            actual.len(),
            expected.len(),
            "{threshold}: actual={actual:?} expected={expected:?}"
        );
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert!(
                (a.t.raw() - e.t.raw()).abs() < Minutes::new(2.0).to::<Day>(),
                "{threshold}: time mismatch {a:?} vs {e:?}"
            );
        }
    }

    #[test]
    fn daily_predictor_finds_horizon_crossings() {
        let site = greenwich();
        let window = utc_window(1);
        let (labelled, diag) = solar_daily_crossings_impl(
            site,
            window,
            Degrees::new(0.0),
            InternalSearchConfig::default(),
        );
        assert!(labelled.len() >= 2, "expected rise and set: {labelled:?}");
        assert_eq!(diag.scan_fallback_days, 0);
        assert!(diag.precise_evaluations < 200, "too many evals: {diag:?}");
    }

    #[test]
    fn daily_predictor_matches_scan_baseline_week() {
        let site = greenwich();
        let window = utc_window(7);
        for threshold in [0.0, -6.0, -12.0, -18.0] {
            let thr = Degrees::new(threshold);
            let (daily, _) =
                solar_daily_crossings_impl(site, window, thr, InternalSearchConfig::default());
            let scan = scan_labelled(site, window, thr);
            assert_labelled_close(&daily, &scan, &format!("threshold {threshold}"));
        }
    }

    #[test]
    fn analytic_transit_within_hours_of_precise_max() {
        let site = greenwich();
        let full_day = utc_window(1);
        let state = solar_daily_state(site, full_day).expect("daily state");
        let mut best_t = full_day.start;
        let mut best_alt = f64::NEG_INFINITY;
        let step = Minutes::new(10.0).to::<Day>();
        let mut t = full_day.start;
        while t <= full_day.end {
            let alt = sun_altitude_rad(t, &site).value();
            if alt > best_alt {
                best_alt = alt;
                best_t = t;
            }
            t = add_days(t, step);
        }
        let delta = (state.approx_transit.raw() - best_t.raw()).abs();
        assert!(
            delta < Hours::new(3.0).to::<Day>(),
            "transit error {delta:?} approx={:?} best={best_t:?}",
            state.approx_transit
        );
    }

    #[test]
    fn partial_window_before_sunrise_matches_scan() {
        let site = greenwich();
        let day = utc_window(1);
        let scan = scan_labelled(site, day, Degrees::new(0.0));
        let sunrise = scan.iter().find(|c| c.direction > 0).expect("sunrise").t;
        let window = Interval::new(
            add_days(sunrise, -Minutes::new(10.0).to::<Day>()),
            add_days(sunrise, Minutes::new(30.0).to::<Day>()),
        );
        let (daily, _) = solar_daily_crossings_impl(
            site,
            window,
            Degrees::new(0.0),
            InternalSearchConfig::default(),
        );
        let expected = scan_labelled(site, window, Degrees::new(0.0));
        assert_labelled_close(&daily, &expected, "partial before sunrise");
    }

    #[test]
    fn partial_window_after_sunset_matches_scan() {
        let site = greenwich();
        let day = utc_window(1);
        let scan = scan_labelled(site, day, Degrees::new(0.0));
        let sunset = scan.iter().find(|c| c.direction < 0).expect("sunset").t;
        let window = Interval::new(
            add_days(sunset, -Minutes::new(30.0).to::<Day>()),
            add_days(sunset, Minutes::new(10.0).to::<Day>()),
        );
        let (daily, _) = solar_daily_crossings_impl(
            site,
            window,
            Degrees::new(0.0),
            InternalSearchConfig::default(),
        );
        let expected = scan_labelled(site, window, Degrees::new(0.0));
        assert_labelled_close(&daily, &expected, "partial after sunset");
    }

    #[test]
    fn polar_summer_always_above_twilight() {
        let site = Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(80.0), Meters::new(0.0));
        let window = Interval::new(
            ModifiedJulianDate::from(Utc.with_ymd_and_hms(2026, 6, 20, 0, 0, 0).single().unwrap()),
            ModifiedJulianDate::from(Utc.with_ymd_and_hms(2026, 6, 27, 0, 0, 0).single().unwrap()),
        );
        let (daily, diag) = solar_daily_crossings_impl(
            site,
            window,
            Degrees::new(-18.0),
            InternalSearchConfig::default(),
        );
        assert!(daily.is_empty());
        assert_eq!(diag.scan_fallback_days, 0);
    }

    #[test]
    fn mid_latitude_diagnostics_low_fallback() {
        let site = roque();
        let window = utc_window(30);
        let (_, diag) = solar_daily_crossings_impl(
            site,
            window,
            Degrees::new(-18.0),
            InternalSearchConfig::default(),
        );
        assert_eq!(diag.scan_fallback_days, 0);
        assert_eq!(diag.bracket_failures, 0);
    }

    #[test]
    fn expanded_bracket_used_for_normal_crossings() {
        let site = greenwich();
        let window = utc_window(1);
        let (_, diag) = solar_daily_crossings_impl(
            site,
            window,
            Degrees::new(0.0),
            InternalSearchConfig::default(),
        );
        assert!(
            diag.expanded_brackets >= 2,
            "expected bracket expansion: {diag:?}"
        );
    }
}
