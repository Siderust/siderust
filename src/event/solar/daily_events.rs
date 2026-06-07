// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Solar daily threshold-crossing predictor with precise Brent validation.
//!
//! Uses a per-day hour-angle model to bracket rise/set candidates, refines each
//! root against `sun_altitude_rad`, and falls back locally to the generic
//! Chebyshev-first engine or scan+Brent when the analytic model is unreliable.

use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::event::altitude::search::InternalSearchConfig;
use crate::event::altitude::{CrossingDirection, CrossingEvent};
use crate::event::search::crossings;
use crate::event::search::intervals::LabeledCrossing;
use crate::event::search::scan_fallback;
use crate::qtty::*;
use crate::time::{Interval, ModifiedJulianDate};

use super::altitude_periods::sun_altitude_rad;

const GRAZE_EPS: f64 = 1e-5;
const BRACKET_HALF: Days = Minutes::new(15.0).to_const::<Day>();
const SCAN_STEP: Days = Hours::new(2.0).to_const::<Day>();
/// Hour-angle rate: 2π radians per mean solar day.
const HA_RATE_RAD_PER_DAY: f64 = std::f64::consts::TAU;

/// Runtime counters for the solar daily predictor (tests/benches only).
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SolarDailyDiagnostics {
    pub predicted_days: usize,
    pub candidate_crossings: usize,
    pub refined_crossings: usize,
    pub grazing_days: usize,
    pub bracket_failures: usize,
    pub chebyshev_fallback_days: usize,
    pub scan_fallback_days: usize,
    pub precise_evaluations: usize,
}

enum DayThresholdCase {
    AlwaysBelow,
    AlwaysAbove,
    Crossings(f64),
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
            opts.fallback_scan_step(SCAN_STEP),
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
    let mut day_start = floor_day(window.start);

    while day_start < window.end {
        let day_end = day_start.raw() + Days::new(1.0);
        let day_end = ModifiedJulianDate::new(day_end.value().min(window.end.raw().value()));
        let day = Interval::new(
            if day_start < window.start {
                window.start
            } else {
                day_start
            },
            day_end,
        );
        if day.end <= day.start {
            day_start = day_end;
            continue;
        }

        diagnostics.predicted_days += 1;
        let day_crossings =
            predict_day_crossings(site, day, thr_sin, sin_lat, cos_lat, opts, &mut diagnostics);
        labelled.extend(day_crossings);
        day_start = ModifiedJulianDate::new((day_start.raw() + Days::new(1.0)).value());
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
    day: Interval<ModifiedJulianDate>,
    thr_sin: f64,
    sin_lat: f64,
    cos_lat: f64,
    opts: InternalSearchConfig,
    diagnostics: &mut SolarDailyDiagnostics,
) -> Vec<LabeledCrossing> {
    let mjd_mid = ModifiedJulianDate::new((day.start.raw() + day.end.raw()) / Days::new(2.0));
    diagnostics.precise_evaluations += 1;
    let dec = sun_declination_rad(mjd_mid, &site);
    let (sin_dec, cos_dec) = dec.sin_cos();
    let denom = cos_lat * cos_dec;
    if denom.abs() < 1e-12 {
        diagnostics.grazing_days += 1;
        return fallback_day(site, day, thr_sin, opts, diagnostics);
    }

    let cos_h0 = (thr_sin - sin_lat * sin_dec) / denom;
    let case = if cos_h0 > 1.0 + GRAZE_EPS {
        DayThresholdCase::AlwaysBelow
    } else if cos_h0 < -1.0 - GRAZE_EPS {
        DayThresholdCase::AlwaysAbove
    } else if !(-1.0..=1.0).contains(&cos_h0)
        || (cos_h0 + 1.0).abs() < GRAZE_EPS
        || (cos_h0 - 1.0).abs() < GRAZE_EPS
    {
        diagnostics.grazing_days += 1;
        return fallback_day(site, day, thr_sin, opts, diagnostics);
    } else {
        DayThresholdCase::Crossings(cos_h0.clamp(-1.0, 1.0).acos())
    };

    match case {
        DayThresholdCase::AlwaysBelow | DayThresholdCase::AlwaysAbove => Vec::new(),
        DayThresholdCase::Crossings(h0) => {
            let transit = find_transit_mjd(site, day, diagnostics);
            let dt = h0 / HA_RATE_RAD_PER_DAY;
            let candidates = [
                (transit.raw() - Days::new(dt), 1_i8),
                (transit.raw() + Days::new(dt), -1_i8),
            ];
            let mut refined = Vec::new();
            for (pred_raw, expected_dir) in candidates {
                diagnostics.candidate_crossings += 1;
                let pred = ModifiedJulianDate::new(pred_raw.value());
                if pred < day.start || pred > day.end {
                    continue;
                }
                if let Some((root, direction)) =
                    refine_candidate(site, day, pred, thr_sin, opts, diagnostics)
                {
                    let _ = expected_dir;
                    refined.push(LabeledCrossing {
                        t: root,
                        direction: direction as i32,
                    });
                } else {
                    diagnostics.bracket_failures += 1;
                }
            }
            if refined.is_empty() && diagnostics.bracket_failures > 0 {
                return fallback_day(site, day, thr_sin, opts, diagnostics);
            }
            refined
        }
    }
}

fn refine_candidate(
    site: Geodetic<ECEF>,
    day: Interval<ModifiedJulianDate>,
    predicted: ModifiedJulianDate,
    thr_sin: f64,
    opts: InternalSearchConfig,
    diagnostics: &mut SolarDailyDiagnostics,
) -> Option<(ModifiedJulianDate, i8)> {
    let lo = ModifiedJulianDate::new(
        (predicted.raw() - BRACKET_HALF)
            .max(day.start.raw())
            .value(),
    );
    let hi = ModifiedJulianDate::new((predicted.raw() + BRACKET_HALF).min(day.end.raw()).value());
    if (hi.raw() - lo.raw()) < Days::new(1e-12) {
        return None;
    }

    diagnostics.precise_evaluations += 2;
    let f_lo = sun_altitude_rad(lo, &site).sin() - thr_sin;
    let f_hi = sun_altitude_rad(hi, &site).sin() - thr_sin;

    if f_lo.signum() * f_hi.signum() > 0.0 {
        return None;
    }

    let root_days = scan_fallback::brent_f64(
        lo.raw().value(),
        hi.raw().value(),
        f_lo,
        f_hi,
        |days| sun_altitude_rad(ModifiedJulianDate::new(days), &site).sin() - thr_sin,
        opts.time_tolerance.value(),
        opts.chebyshev.max_residual.max(crossings::POLY_ZERO_TOL),
    )?;
    let root = ModifiedJulianDate::new(root_days);
    if root < day.start || root > day.end {
        return None;
    }

    let direction = classify_direction(site, root, thr_sin, diagnostics);
    Some((root, direction))
}

fn classify_direction(
    site: Geodetic<ECEF>,
    root: ModifiedJulianDate,
    thr_sin: f64,
    diagnostics: &mut SolarDailyDiagnostics,
) -> i8 {
    let eps = Minutes::new(1.0).to::<Day>();
    diagnostics.precise_evaluations += 2;
    let before = sun_altitude_rad(ModifiedJulianDate::new((root.raw() - eps).value()), &site).sin()
        - thr_sin;
    let after = sun_altitude_rad(ModifiedJulianDate::new((root.raw() + eps).value()), &site).sin()
        - thr_sin;
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

fn find_transit_mjd(
    site: Geodetic<ECEF>,
    day: Interval<ModifiedJulianDate>,
    diagnostics: &mut SolarDailyDiagnostics,
) -> ModifiedJulianDate {
    let step = Minutes::new(30.0).to::<Day>();
    let mut best_t = day.start;
    let mut best_alt = f64::NEG_INFINITY;
    let mut t = day.start;
    while t <= day.end {
        diagnostics.precise_evaluations += 1;
        let alt = sun_altitude_rad(t, &site).value();
        if alt > best_alt {
            best_alt = alt;
            best_t = t;
        }
        let next = t.raw() + step;
        if next >= day.end.raw() {
            break;
        }
        t = ModifiedJulianDate::new(next.value());
    }
    best_t
}

fn sun_declination_rad(mjd: ModifiedJulianDate, _site: &Geodetic<ECEF>) -> f64 {
    let equ = crate::bodies::solar_system::Sun::get_apparent_geocentric_equ::<AstronomicalUnit>(
        mjd.to::<crate::JD>(),
    );
    equ.dec().to::<Radian>().value()
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

fn chebyshev_labelled_crossings_raw(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    thr_sin: f64,
    opts: InternalSearchConfig,
) -> (Vec<LabeledCrossing>, crossings::SearchDiagnostics) {
    let signal = |t: ModifiedJulianDate| sun_altitude_rad(t, &site).sin();
    let (labelled, _, diag) =
        crossings::find_labelled_crossings(window, SCAN_STEP, &signal, thr_sin, opts);
    (labelled, diag)
}

fn floor_day(t: ModifiedJulianDate) -> ModifiedJulianDate {
    ModifiedJulianDate::new(t.raw().value().floor())
}

fn sort_dedup_labelled(crossings: &mut Vec<LabeledCrossing>) {
    const DEDUPE: Days = Days::new(1e-8);
    crossings.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
    crossings.dedup_by(|a, b| (a.t.raw() - b.t.raw()).abs() < DEDUPE);
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::{TimeZone, Utc};

    fn greenwich() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.4769), Meters::new(0.0))
    }

    fn utc_window(days: u32) -> Interval<ModifiedJulianDate> {
        let start = Utc.with_ymd_and_hms(2026, 1, 1, 0, 0, 0).single().unwrap();
        let end = start + chrono::Duration::days(days as i64);
        Interval::new(
            ModifiedJulianDate::from(start),
            ModifiedJulianDate::from(end),
        )
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
    }

    #[test]
    fn daily_predictor_matches_chebyshev_baseline_week() {
        let site = greenwich();
        let window = utc_window(7);
        for threshold in [0.0, -6.0, -12.0, -18.0] {
            let thr = Degrees::new(threshold);
            let (daily, _) =
                solar_daily_crossings_impl(site, window, thr, InternalSearchConfig::default());
            let (cheb, _) = chebyshev_labelled_crossings_raw(
                site,
                window,
                thr.to::<Radian>().sin(),
                InternalSearchConfig {
                    disable_solar_daily_predictor: true,
                    ..InternalSearchConfig::default()
                },
            );
            assert_eq!(
                daily.len(),
                cheb.len(),
                "threshold {threshold}: daily={daily:?} cheb={cheb:?}"
            );
        }
    }
}
