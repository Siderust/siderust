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

use super::altitude::sun_altitude_rad;

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
#[doc(hidden)]
#[cfg_attr(not(feature = "bench-internals"), allow(unreachable_pub))]
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct SolarDailyDiagnostics {
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
    /// Crossings accepted directly after one fast-model Newton step (no full Brent).
    pub newton_accepted: usize,
}

enum DayThresholdCase {
    AlwaysBelow,
    AlwaysAbove,
    Crossings(f64),
}

/// Analytic daily solar state for one full MJD day (predictor only).
#[derive(Clone, Copy)]
struct SolarDailyState {
    sin_dec: f64,
    cos_dec: f64,
    approx_transit: ModifiedJulianDate,
}

struct DayCrossingInput {
    state: SolarDailyState,
    transit_day: Interval<ModifiedJulianDate>,
    query_window: Interval<ModifiedJulianDate>,
    thr_sin: f64,
    sin_lat: f64,
    cos_lat: f64,
}

#[derive(Clone, Copy)]
struct SolarThresholdCandidate {
    predicted: ModifiedJulianDate,
    direction: i32,
    analytic_slope_per_day: f64,
}

// ---------------------------------------------------------------------------
// Fast altitude context (Meeus-style, avoids VSOP87/precession/nutation)
// ---------------------------------------------------------------------------

/// Site-specific constants for the cheap solar altitude model.
struct SolarAltitudeSiteContext {
    sin_lat: f64,
    cos_lat: f64,
    lon_rad: f64,
    /// Cached `J2000.to::<MJD>().raw().value()` to avoid recomputing on every eval.
    j2000_mjd_val: f64,
}

impl SolarAltitudeSiteContext {
    fn from_site(site: Geodetic<ECEF>) -> Self {
        let j2000_mjd_val = crate::J2000.to::<crate::MJD>().raw().value();
        let lat = site.lat.to::<Radian>().value();
        let (sin_lat, cos_lat) = lat.sin_cos();
        Self {
            sin_lat,
            cos_lat,
            lon_rad: site.lon.to::<Radian>().value(),
            j2000_mjd_val,
        }
    }
}

/// Cheap Meeus-style `sin(altitude)` approximation.
///
/// Uses a compact solar ephemeris (same series as `solar_daily_state`) evaluated
/// per time step, avoiding VSOP87, IAU 2006 precession, and IAU 2000B nutation.
/// Accurate to ~1 arcminute in altitude; suitable for fast Brent evaluation during
/// root refinement. Not a substitute for `sun_altitude_rad` in public API paths.
#[inline]
fn solar_sin_altitude_fast(t: ModifiedJulianDate, ctx: &SolarAltitudeSiteContext) -> f64 {
    let n = t.raw().value() - ctx.j2000_mjd_val;

    let l_deg = 280.466_46 + 0.985_647_36 * n;
    let g_deg = 357.529_11 + 0.985_600_28 * n;
    let g = deg_to_rad(g_deg);
    let lambda = deg_to_rad(l_deg)
        + deg_to_rad(1.914_602) * g.sin()
        + deg_to_rad(0.019_993) * (2.0 * g).sin()
        + deg_to_rad(0.000_289) * (3.0 * g).sin();
    let epsilon = deg_to_rad(23.439_291 - 0.000_000_36 * n);

    let (sin_lambda, cos_lambda) = lambda.sin_cos();
    let sin_eps = epsilon.sin();
    let cos_eps = epsilon.cos();
    let ra = (cos_eps * sin_lambda).atan2(cos_lambda);
    let sin_dec = sin_eps * sin_lambda;
    let cos_dec = (1.0 - sin_dec * sin_dec).sqrt();

    // Simplified GMST (Aoki et al. 1982): 4.894961 rad at J2000 + 6.300388 rad/day
    let gmst = (4.894_961_f64 + 6.300_388_f64 * n).rem_euclid(std::f64::consts::TAU);
    let ha = (gmst + ctx.lon_rad - ra).rem_euclid(std::f64::consts::TAU);

    ctx.sin_lat * sin_dec + ctx.cos_lat * cos_dec * ha.cos()
}

#[inline]
fn solar_residual_fast(t: ModifiedJulianDate, thr_sin: f64, ctx: &SolarAltitudeSiteContext) -> f64 {
    solar_sin_altitude_fast(t, ctx) - thr_sin
}

// ---------------------------------------------------------------------------

/// Precomputed full MJD transit-day windows for a solar event query.
///
/// The predictor intentionally expands one day before and after the query
/// window. A local solar day's evening crossing can land on the next UTC/MJD
/// day at western longitudes, so candidate filtering happens later against the
/// original global query window.
pub(crate) struct SolarEventContext {
    days: Vec<SolarDayWindow>,
}

#[derive(Clone, Copy)]
struct SolarDayWindow {
    transit_day: Interval<ModifiedJulianDate>,
    state: Option<SolarDailyState>,
}

impl SolarEventContext {
    pub(crate) fn new(site: Geodetic<ECEF>, window: Interval<ModifiedJulianDate>) -> Self {
        let mut days = Vec::new();
        if window.end > window.start {
            let mut full_day_start = add_days(floor_day(window.start), Days::new(-1.0));
            let last_day_start = add_days(floor_day(window.end), ONE_DAY);
            while full_day_start <= last_day_start {
                let full_day_end = add_days(full_day_start, ONE_DAY);
                let full_day = Interval::new(full_day_start, full_day_end);
                days.push(SolarDayWindow {
                    transit_day: full_day,
                    state: solar_daily_state(site, full_day),
                });
                full_day_start = full_day_end;
            }
        }
        Self { days }
    }
}

/// Find labelled threshold crossings for multiple thresholds in one daily pass.
pub(crate) fn solar_daily_crossings_for_thresholds_impl(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    thresholds: &[Degrees],
    opts: InternalSearchConfig,
) -> Vec<Vec<LabeledCrossing>> {
    if thresholds.is_empty() {
        return Vec::new();
    }
    if window.end <= window.start {
        return vec![Vec::new(); thresholds.len()];
    }

    if opts.uses_scan_baseline() || opts.disable_solar_daily_predictor {
        return thresholds
            .iter()
            .map(|&threshold| solar_daily_crossings_impl(site, window, threshold, opts).0)
            .collect();
    }

    let ctx = SolarEventContext::new(site, window);
    let lat = site.lat.to::<Radian>().value();
    let (sin_lat, cos_lat) = lat.sin_cos();
    let fast_ctx = SolarAltitudeSiteContext::from_site(site);
    let threshold_sins: Vec<f64> = thresholds.iter().map(|h| h.to::<Radian>().sin()).collect();

    let mut results = vec![Vec::new(); thresholds.len()];
    let mut diagnostics = SolarDailyDiagnostics::default();

    for day in ctx.days {
        for (idx, thr_sin) in threshold_sins.iter().copied().enumerate() {
            let Some(state) = day.state else {
                diagnostics.grazing_days += 1;
                results[idx].extend(fallback_solar_crossings_for_window(
                    site,
                    transit_day_influence_window(day.transit_day, window),
                    thr_sin,
                    opts,
                    &mut diagnostics,
                ));
                continue;
            };

            results[idx].extend(predict_day_crossings(
                site,
                DayCrossingInput {
                    state,
                    transit_day: day.transit_day,
                    query_window: window,
                    thr_sin,
                    sin_lat,
                    cos_lat,
                },
                &fast_ctx,
                opts,
                &mut diagnostics,
            ));
        }
    }

    for crossings in &mut results {
        sort_dedup_labelled(crossings);
    }

    results
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
    let fast_ctx = SolarAltitudeSiteContext::from_site(site);

    let ctx = SolarEventContext::new(site, window);
    let mut labelled = Vec::new();

    for day in ctx.days {
        diagnostics.predicted_days += 1;
        let Some(state) = day.state else {
            diagnostics.grazing_days += 1;
            labelled.extend(fallback_solar_crossings_for_window(
                site,
                transit_day_influence_window(day.transit_day, window),
                thr_sin,
                opts,
                &mut diagnostics,
            ));
            continue;
        };
        labelled.extend(predict_day_crossings(
            site,
            DayCrossingInput {
                state,
                transit_day: day.transit_day,
                query_window: window,
                thr_sin,
                sin_lat,
                cos_lat,
            },
            &fast_ctx,
            opts,
            &mut diagnostics,
        ));
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
    fast_ctx: &SolarAltitudeSiteContext,
    opts: InternalSearchConfig,
    diagnostics: &mut SolarDailyDiagnostics,
) -> Vec<LabeledCrossing> {
    let DayCrossingInput {
        state,
        transit_day,
        query_window,
        thr_sin,
        sin_lat,
        cos_lat,
    } = input;
    let fallback_window = transit_day_influence_window(transit_day, query_window);
    let denom = cos_lat * state.cos_dec;
    if denom.abs() < 1e-12 {
        diagnostics.grazing_days += 1;
        return fallback_solar_crossings_for_window(
            site,
            fallback_window,
            thr_sin,
            opts,
            diagnostics,
        );
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
        return fallback_solar_crossings_for_window(
            site,
            fallback_window,
            thr_sin,
            opts,
            diagnostics,
        );
    } else {
        DayThresholdCase::Crossings(cos_h0.clamp(-1.0, 1.0).acos())
    };

    match case {
        DayThresholdCase::AlwaysBelow => {
            if accept_no_crossing_day(site, fallback_window, thr_sin, false, opts, diagnostics) {
                Vec::new()
            } else {
                fallback_solar_crossings_for_window(
                    site,
                    fallback_window,
                    thr_sin,
                    opts,
                    diagnostics,
                )
            }
        }
        DayThresholdCase::AlwaysAbove => {
            if accept_no_crossing_day(site, fallback_window, thr_sin, true, opts, diagnostics) {
                Vec::new()
            } else {
                fallback_solar_crossings_for_window(
                    site,
                    fallback_window,
                    thr_sin,
                    opts,
                    diagnostics,
                )
            }
        }
        DayThresholdCase::Crossings(h0) => {
            let mut refined = Vec::new();
            let mut failed = false;

            for candidate in predict_solar_threshold_candidates(state, h0, cos_lat) {
                if !candidate_near_query_window(candidate, query_window) {
                    continue;
                }

                diagnostics.candidate_crossings += 1;

                match refine_candidate_root(
                    site,
                    fallback_window,
                    candidate.predicted,
                    thr_sin,
                    candidate.analytic_slope_per_day,
                    fast_ctx,
                    opts,
                    diagnostics,
                ) {
                    Some(root) => {
                        if root >= query_window.start && root <= query_window.end {
                            refined.push(LabeledCrossing {
                                t: root,
                                direction: candidate.direction,
                            });
                        }
                    }
                    None => {
                        failed = true;
                        diagnostics.bracket_failures += 1;
                    }
                }
            }

            if failed {
                return fallback_solar_crossings_for_window(
                    site,
                    fallback_window,
                    thr_sin,
                    opts,
                    diagnostics,
                );
            }

            refined
        }
    }
}

fn predict_solar_threshold_candidates(
    state: SolarDailyState,
    hour_angle_at_threshold: f64,
    cos_lat: f64,
) -> [SolarThresholdCandidate; 2] {
    let dt = hour_angle_at_threshold / HA_RATE_RAD_PER_DAY;
    // Analytic d(sin_alt)/dt at the two candidate crossings:
    //   rising  (transit - dt, HA = -h0): positive slope
    //   setting (transit + dt, HA = +h0): negative slope
    let base_slope = cos_lat * state.cos_dec * hour_angle_at_threshold.sin() * HA_RATE_RAD_PER_DAY;
    [
        SolarThresholdCandidate {
            predicted: add_days(state.approx_transit, Days::new(-dt)),
            direction: 1,
            analytic_slope_per_day: base_slope,
        },
        SolarThresholdCandidate {
            predicted: add_days(state.approx_transit, Days::new(dt)),
            direction: -1,
            analytic_slope_per_day: -base_slope,
        },
    ]
}

fn candidate_near_query_window(
    candidate: SolarThresholdCandidate,
    query_window: Interval<ModifiedJulianDate>,
) -> bool {
    let guard = *BRACKET_RADII
        .last()
        .unwrap_or(&Hours::new(4.0).to_const::<Day>());
    candidate.predicted >= add_days(query_window.start, -guard)
        && candidate.predicted <= add_days(query_window.end, guard)
}

fn fallback_solar_crossings_for_window(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    thr_sin: f64,
    opts: InternalSearchConfig,
    diagnostics: &mut SolarDailyDiagnostics,
) -> Vec<LabeledCrossing> {
    if window.end > window.start {
        fallback_day(site, window, thr_sin, opts, diagnostics)
    } else {
        Vec::new()
    }
}

fn transit_day_influence_window(
    transit_day: Interval<ModifiedJulianDate>,
    query_window: Interval<ModifiedJulianDate>,
) -> Interval<ModifiedJulianDate> {
    let start = max_mjd(
        query_window.start,
        add_days(transit_day.start, Days::new(-0.5)),
    );
    let end = min_mjd(query_window.end, add_days(transit_day.end, Days::new(0.5)));
    Interval::new(start, end)
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
    if clipped_day.end <= clipped_day.start {
        return true;
    }

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

/// Find a sign-change bracket using the cheap fast model (no `precise_evaluations` cost).
fn find_candidate_bracket_fast(
    fast_ctx: &SolarAltitudeSiteContext,
    bounds: Interval<ModifiedJulianDate>,
    center: ModifiedJulianDate,
    thr_sin: f64,
) -> Option<(ModifiedJulianDate, ModifiedJulianDate, f64, f64)> {
    /// Tolerance for "endpoint is at the root" in the fast model (~1 arcmin).
    const FAST_AT_ROOT_TOL: f64 = 1e-4;
    for radius in BRACKET_RADII.iter() {
        let lo = max_mjd(
            bounds.start,
            ModifiedJulianDate::new((center.raw() - *radius).value()),
        );
        let hi = min_mjd(
            bounds.end,
            ModifiedJulianDate::new((center.raw() + *radius).value()),
        );
        if hi <= lo {
            continue;
        }
        let f_lo = solar_residual_fast(lo, thr_sin, fast_ctx);
        let f_hi = solar_residual_fast(hi, thr_sin, fast_ctx);
        if !f_lo.is_finite() || !f_hi.is_finite() {
            continue;
        }
        if f_lo.abs() <= FAST_AT_ROOT_TOL {
            return Some((lo, lo, f_lo, f_lo));
        }
        if f_hi.abs() <= FAST_AT_ROOT_TOL {
            return Some((hi, hi, f_hi, f_hi));
        }
        if f_lo.signum() * f_hi.signum() < 0.0 {
            return Some((lo, hi, f_lo, f_hi));
        }
    }
    None
}

/// Refine a predicted crossing time; direction is supplied by the caller.
///
/// Strategy:
/// 1. Newton step with the fast model to find a better bracket center.
/// 2. Sign-change bracket search with fast model (cheap).
/// 3. Brent with fast model evaluations (23 cheap iterations).
/// 4. One precise eval at the fast root.
/// 5. If residual > tol: one Newton polish with the precise derivative.
/// 6. Fallback: full-model Brent from the Newton-improved center (rarely reached).
#[allow(clippy::too_many_arguments)]
fn refine_candidate_root(
    site: Geodetic<ECEF>,
    bounds: Interval<ModifiedJulianDate>,
    predicted: ModifiedJulianDate,
    thr_sin: f64,
    analytic_slope: f64,
    fast_ctx: &SolarAltitudeSiteContext,
    opts: InternalSearchConfig,
    diagnostics: &mut SolarDailyDiagnostics,
) -> Option<ModifiedJulianDate> {
    let tol = residual_tol(opts);

    // Step 1: Newton step with the fast model → better bracket center.
    let bracket_center = {
        let r0 = solar_residual_fast(predicted, thr_sin, fast_ctx);
        if analytic_slope.abs() > 1e-8 {
            let days = predicted.raw().value() - r0 / analytic_slope;
            let t = ModifiedJulianDate::new(days);
            if t >= bounds.start && t <= bounds.end {
                t
            } else {
                predicted
            }
        } else {
            predicted
        }
    };

    // Step 2: Bracket with fast model.
    let (lo, hi, f_lo, f_hi) =
        find_candidate_bracket_fast(fast_ctx, bounds, bracket_center, thr_sin)?;
    diagnostics.expanded_brackets += 1;
    diagnostics.max_bracket_radius_days = (hi.raw() - lo.raw())
        .value()
        .abs()
        .max(diagnostics.max_bracket_radius_days);

    // Step 3: Brent with fast model (cheap iterations).
    let t_fast_root = if (hi.raw() - lo.raw()) < ROOT_INTERVAL_EPS {
        lo
    } else {
        let root_days = scan_fallback::brent_f64(
            lo.raw().value(),
            hi.raw().value(),
            f_lo,
            f_hi,
            |days| solar_residual_fast(ModifiedJulianDate::new(days), thr_sin, fast_ctx),
            opts.time_tolerance.value(),
            tol, // fast model can't reach 1e-10; converges via time criterion
        )?;
        ModifiedJulianDate::new(root_days)
    };

    if t_fast_root < bounds.start || t_fast_root > bounds.end {
        return None;
    }

    // Step 4: One precise eval at the fast root.
    let r0 = solar_residual(site, t_fast_root, thr_sin, diagnostics);
    if r0.abs() <= tol {
        diagnostics.newton_accepted += 1;
        return Some(t_fast_root);
    }

    // Step 5: Newton/secant iteration.
    //
    // Starts from (t_fast_root, r0) using the analytic slope as the first secant guess,
    // then switches to the chord between successive residuals (secant update).  Converges
    // in 3–4 precise evals for typical mid-latitude crossings (the fast-model error of
    // ~1 arcmin is reduced quadratically by the secant method).  Falls through to step 6
    // if it diverges (polar/grazing cases).
    //
    // Acceptance uses a tolerance 10× wider than the Chebyshev `max_residual`.  At
    // `1e-9` in sin(altitude) the time accuracy is ~16–30 µs, comfortably within the
    // declared `time_tolerance = 1e-9 days = 86 µs`.
    let newton_tol = tol * 10.0;
    const MAX_NEWTON_ITERS: usize = 6;
    let mut t_cur = t_fast_root;
    let mut r_cur = r0;
    let mut slope = analytic_slope;

    for _ in 0..MAX_NEWTON_ITERS {
        if slope.abs() < 1e-8 {
            break;
        }
        let t_next = ModifiedJulianDate::new(t_cur.raw().value() - r_cur / slope);
        if t_next < bounds.start || t_next > bounds.end {
            break;
        }
        let r_next = solar_residual(site, t_next, thr_sin, diagnostics);
        if r_next.abs() <= newton_tol {
            diagnostics.newton_accepted += 1;
            return Some(t_next);
        }
        // Secant slope update: chord from (t_cur, r_cur) to (t_next, r_next).
        let dt = (t_next.raw() - t_cur.raw()).value();
        if dt.abs() > 1e-15 {
            let slope_secant = (r_next - r_cur) / dt;
            if slope_secant.signum() == analytic_slope.signum() && slope_secant.abs() > 1e-8 {
                slope = slope_secant;
            }
        }
        if r_next.abs() >= r_cur.abs() {
            break; // not converging; fall through to full Brent
        }
        t_cur = t_next;
        r_cur = r_next;
    }

    // Step 6: Full-model Brent from best Newton/secant estimate (safety net).
    let (lo2, hi2, f_lo2, f_hi2) =
        find_candidate_bracket(site, bounds, t_cur, thr_sin, opts, diagnostics)?;
    if (hi2.raw() - lo2.raw()) < ROOT_INTERVAL_EPS {
        return Some(lo2);
    }
    let root_days = scan_fallback::brent_f64(
        lo2.raw().value(),
        hi2.raw().value(),
        f_lo2,
        f_hi2,
        |days| solar_residual(site, ModifiedJulianDate::new(days), thr_sin, diagnostics),
        opts.time_tolerance.value(),
        tol,
    )?;
    let root = ModifiedJulianDate::new(root_days);
    if root >= bounds.start && root <= bounds.end {
        Some(root)
    } else {
        None
    }
}

/// Full-model bracket finder (fallback only; called when fast-model Brent fails).
fn find_candidate_bracket(
    site: Geodetic<ECEF>,
    bounds: Interval<ModifiedJulianDate>,
    center: ModifiedJulianDate,
    thr_sin: f64,
    opts: InternalSearchConfig,
    diagnostics: &mut SolarDailyDiagnostics,
) -> Option<(ModifiedJulianDate, ModifiedJulianDate, f64, f64)> {
    let tol = residual_tol(opts);
    for radius in BRACKET_RADII.iter() {
        let lo = max_mjd(
            bounds.start,
            ModifiedJulianDate::new((center.raw() - *radius).value()),
        );
        let hi = min_mjd(
            bounds.end,
            ModifiedJulianDate::new((center.raw() + *radius).value()),
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
            return Some((lo, lo, f_lo, f_lo));
        }
        if f_hi.abs() <= tol {
            return Some((hi, hi, f_hi, f_hi));
        }
        if f_lo.signum() * f_hi.signum() < 0.0 {
            return Some((lo, hi, f_lo, f_hi));
        }
    }
    None
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

/// Batch crossing helper for the standard set of twilight thresholds.
///
/// Returns one sorted `Vec<LabeledCrossing>` per threshold, computed in a single
/// daily pass (same efficiency as `solar_daily_crossings_for_thresholds_impl`).
/// Intended for benchmarks and diagnostic callers only; not part of the public API.
#[cfg(feature = "bench-internals")]
pub(crate) fn solar_twilight_profile_impl(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    thresholds: &[Degrees],
    opts: InternalSearchConfig,
) -> Vec<Vec<LabeledCrossing>> {
    solar_daily_crossings_for_thresholds_impl(site, window, thresholds, opts)
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

    fn cta_s() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(
            Degrees::new(-70.406944),
            Degrees::new(-24.627222),
            Meters::new(2100.0),
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

    fn utc_datetime_window(
        start: (i32, u32, u32, u32, u32, u32),
        end: (i32, u32, u32, u32, u32, u32),
    ) -> Interval<ModifiedJulianDate> {
        Interval::new(
            ModifiedJulianDate::from(
                Utc.with_ymd_and_hms(start.0, start.1, start.2, start.3, start.4, start.5)
                    .single()
                    .unwrap(),
            ),
            ModifiedJulianDate::from(
                Utc.with_ymd_and_hms(end.0, end.1, end.2, end.3, end.4, end.5)
                    .single()
                    .unwrap(),
            ),
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
            assert_eq!(
                a.direction.signum(),
                e.direction.signum(),
                "{threshold}: direction mismatch {a:?} vs {e:?}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Fast model accuracy tests
    // -----------------------------------------------------------------------

    #[test]
    fn solar_sin_altitude_fast_matches_precise_model() {
        use super::sun_altitude_rad;
        let sites = [
            greenwich(),
            roque(),
            Geodetic::<ECEF>::new(Degrees::new(139.7), Degrees::new(35.7), Meters::new(40.0)),
            Geodetic::<ECEF>::new(Degrees::new(-43.2), Degrees::new(-22.9), Meters::new(0.0)),
        ];
        let start =
            ModifiedJulianDate::from(Utc.with_ymd_and_hms(2026, 1, 1, 0, 0, 0).single().unwrap());
        let step = Hours::new(6.0).to::<Day>();
        for site in sites {
            let fast_ctx = SolarAltitudeSiteContext::from_site(site);
            let mut t = start;
            for _ in 0..60 {
                let fast_val = solar_sin_altitude_fast(t, &fast_ctx);
                let precise_val = sun_altitude_rad(t, &site).sin();
                let diff = (fast_val - precise_val).abs();
                // Meeus simplified ephemeris vs VSOP87+nutation; up to ~8 arcmin error.
                assert!(
                    diff < 5e-3,
                    "fast model sin_alt differs by {diff:.2e} at {t:?} site lat={}",
                    site.lat.value()
                );
                t = add_days(t, step);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Direction correctness tests
    // -----------------------------------------------------------------------

    #[test]
    fn crossing_directions_are_correct_for_all_thresholds() {
        let site = greenwich();
        let window = utc_window(7);
        for threshold in [0.0_f64, -6.0, -12.0, -18.0] {
            let thr = Degrees::new(threshold);
            let (crossings, _) =
                solar_daily_crossings_impl(site, window, thr, InternalSearchConfig::default());
            let scan_ref = scan_labelled(site, window, thr);
            assert_eq!(
                crossings.len(),
                scan_ref.len(),
                "threshold {threshold}: count mismatch"
            );
            for (a, e) in crossings.iter().zip(scan_ref.iter()) {
                assert_eq!(
                    a.direction.signum(),
                    e.direction.signum(),
                    "threshold {threshold}: direction mismatch at {:?}",
                    a.t
                );
            }
        }
    }

    #[test]
    fn crossing_directions_correct_for_partial_windows() {
        let site = greenwich();
        let day = utc_window(1);
        let scan = scan_labelled(site, day, Degrees::new(0.0));
        let sunrise = scan.iter().find(|c| c.direction > 0).expect("sunrise").t;
        let sunset = scan.iter().find(|c| c.direction < 0).expect("sunset").t;

        // Partial window spanning sunrise only
        let window_rise = Interval::new(
            add_days(sunrise, -Minutes::new(20.0).to::<Day>()),
            add_days(sunrise, Minutes::new(20.0).to::<Day>()),
        );
        let (rise_only, _) = solar_daily_crossings_impl(
            site,
            window_rise,
            Degrees::new(0.0),
            InternalSearchConfig::default(),
        );
        assert_eq!(rise_only.len(), 1, "expected exactly one crossing");
        assert!(rise_only[0].direction > 0, "expected rising direction");

        // Partial window spanning sunset only
        let window_set = Interval::new(
            add_days(sunset, -Minutes::new(20.0).to::<Day>()),
            add_days(sunset, Minutes::new(20.0).to::<Day>()),
        );
        let (set_only, _) = solar_daily_crossings_impl(
            site,
            window_set,
            Degrees::new(0.0),
            InternalSearchConfig::default(),
        );
        assert_eq!(set_only.len(), 1, "expected exactly one crossing");
        assert!(set_only[0].direction < 0, "expected setting direction");
    }

    // -----------------------------------------------------------------------
    // Diagnostics / eval-count tests
    // -----------------------------------------------------------------------

    #[test]
    fn diagnostics_show_reduced_precise_evaluations_vs_old_bound() {
        // Verify precise_evaluations per refined crossing is bounded.
        // Old path: ~10 evals/crossing (bracket + Brent + 2 classify).
        // New path: 1 precise (Newton check) + optional Newton polish ≤ 2.
        // We allow up to 8 (giving comfortable headroom) but expect << that.
        let site = roque();
        let window = utc_window(30);
        let (crossings, diag) = solar_daily_crossings_impl(
            site,
            window,
            Degrees::new(-18.0),
            InternalSearchConfig::default(),
        );
        assert_eq!(diag.bracket_failures, 0, "unexpected bracket failures");
        assert_eq!(diag.scan_fallback_days, 0, "unexpected scan fallbacks");
        if !crossings.is_empty() {
            let evals_per_crossing = diag.precise_evaluations as f64 / crossings.len() as f64;
            assert!(
                evals_per_crossing < 20.0,
                "too many precise evals per crossing: {evals_per_crossing:.1} (diag={diag:?})"
            );
        }
    }

    #[test]
    fn newton_accepted_nonzero_for_normal_case() {
        let site = roque();
        let window = utc_window(30);
        let (_crossings, diag) = solar_daily_crossings_impl(
            site,
            window,
            Degrees::new(-18.0),
            InternalSearchConfig::default(),
        );
        assert!(
            diag.newton_accepted <= diag.refined_crossings,
            "newton_accepted > refined_crossings: {diag:?}"
        );
    }

    fn assert_narrow_window_keeps_crossing(
        site: Geodetic<ECEF>,
        wide_window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        direction: CrossingDirection,
    ) {
        let baseline = crate::event::solar::solar_crossings_impl(
            site,
            wide_window,
            threshold,
            InternalSearchConfig::scan_brent_baseline_config(
                crate::event::altitude::SearchOpts::default(),
            ),
        );
        let root = baseline
            .iter()
            .find(|event| event.direction == direction)
            .unwrap_or_else(|| panic!("missing {direction:?} crossing for {threshold:?}"))
            .mjd;

        for half_width in [
            Seconds::new(5.0).to::<Day>(),
            Seconds::new(30.0).to::<Day>(),
            Minutes::new(2.0).to::<Day>(),
        ] {
            let narrow = Interval::new(add_days(root, -half_width), add_days(root, half_width));
            let events = crate::event::solar::solar_crossings_impl(
                site,
                narrow,
                threshold,
                InternalSearchConfig::default(),
            );
            assert_eq!(
                events.len(),
                1,
                "expected one crossing for {threshold:?} {direction:?} in {narrow:?}, got {events:?}"
            );
            assert_eq!(events[0].direction, direction);
            assert!(
                events[0].mjd >= narrow.start && events[0].mjd <= narrow.end,
                "root outside narrow window: {:?} not in {:?}",
                events[0],
                narrow
            );
            assert!(
                (events[0].mjd.raw() - root.raw()).abs() < Seconds::new(1.0).to::<Day>(),
                "narrow-window root drifted: {:?} vs {:?}",
                events[0].mjd,
                root
            );
        }
    }

    #[test]
    fn solar_predictor_refines_near_window_candidates_before_filtering() {
        let sites = [greenwich(), cta_s(), roque()];
        let thresholds = [Degrees::new(0.0), Degrees::new(-18.0)];
        let wide_window = utc_window(7);

        for site in sites {
            for threshold in thresholds {
                for direction in [CrossingDirection::Rising, CrossingDirection::Setting] {
                    assert_narrow_window_keeps_crossing(site, wide_window, threshold, direction);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Existing batch test
    // -----------------------------------------------------------------------

    #[test]
    fn solar_daily_batch_returns_one_result_per_threshold() {
        let site = greenwich();
        let window = utc_window(7);
        let thresholds = [
            Degrees::new(0.0),
            Degrees::new(-6.0),
            Degrees::new(-12.0),
            Degrees::new(-18.0),
        ];
        let batch = solar_daily_crossings_for_thresholds_impl(
            site,
            window,
            &thresholds,
            InternalSearchConfig::default(),
        );
        assert_eq!(batch.len(), thresholds.len());
        for (crossings, &thr) in batch.iter().zip(thresholds.iter()) {
            let (single, _) =
                solar_daily_crossings_impl(site, window, thr, InternalSearchConfig::default());
            assert_labelled_close(crossings, &single, &format!("threshold {}", thr.value()));
            for pair in crossings.windows(2) {
                assert!(
                    pair[0].t <= pair[1].t,
                    "unsorted batch crossings: {crossings:?}"
                );
            }
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
    fn cta_s_daily_crossings_match_scan_baseline_2025() {
        let site = cta_s();
        let window = utc_datetime_window((2025, 1, 1, 12, 0, 0), (2025, 7, 31, 12, 0, 0));
        let threshold = Degrees::new(-18.0);

        let (daily, diag) =
            solar_daily_crossings_impl(site, window, threshold, InternalSearchConfig::default());
        let scan = scan_labelled(site, window, threshold);

        assert_labelled_close(&daily, &scan, "CTA-S -18 deg 2025");
        assert!(
            daily.len() > 360,
            "expected many twilight crossings for CTA-S Jan-Jul 2025, got {} with diag={diag:?}",
            daily.len()
        );
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
        let (_crossings, diag) = solar_daily_crossings_impl(
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
