// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Event Finding Functions
//!
//! ## Scientific scope
//!
//! Computes time‑domain altitude events for a topocentric observer:
//! threshold crossings (rises/sets at any specified altitude), upper and
//! lower culminations (local extrema of *h(t)*), and time intervals where
//! the altitude stays inside a user‑defined band. The altitude function
//! itself is delegated to [`AltitudeProvider`], so accuracy and
//! validity inherit from the underlying ephemeris/star model. Atmospheric
//! refraction is not modelled here; observers wanting the standard
//! geometric horizon should pass `−0.833°`.
//!
//! ## Technical scope
//!
//! All `Interval<ModifiedJulianDate>` inputs/outputs are interpreted on the
//! TT axis. Convert UTC timestamps with
//! `tempoch::Time::<tempoch::UTC>::from_chrono(...).to::<tempoch::TT>().into()`
//! into `ModifiedJulianDate` first. Public functions: [`crossings`], [`culminations`],
//! [`altitude_ranges`], [`above_threshold`], [`below_threshold`]. The
//! refinement uses Chebyshev-first crossing discovery with precise validation
//! and local scan+Brent fallback; precision is governed by [`SearchOpts`].
//!
//! ## References
//! None.

use super::provider::AltitudeProvider;
use super::search::{SearchOpts, SearchOptsV2, DEFAULT_SCAN_STEP, EXTREMA_SCAN_STEP};
use super::types::{CrossingDirection, CrossingEvent, CulminationEvent, CulminationKind};
use crate::astro::apparent::CorrectionPolicy;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::event::search::{extrema, intervals};
use crate::qtty::*;
use crate::time::{complement_within, Interval, ModifiedJulianDate};

// ---------------------------------------------------------------------------
// Internal: build altitude function from trait
// ---------------------------------------------------------------------------

/// Build an altitude function from any `AltitudeProvider`.
fn make_altitude_fn<'a, T: AltitudeProvider>(
    target: &'a T,
    site: &'a Geodetic<ECEF>,
    policy: CorrectionPolicy,
) -> impl Fn(ModifiedJulianDate) -> Radians + 'a {
    let site = *site;
    move |t: ModifiedJulianDate| target.altitude_at_with_policy(&site, t, policy)
}

fn scan_step_for_opts<T: AltitudeProvider>(target: &T, opts: &SearchOptsV2) -> Days {
    super::search::resolve_scan_step(target.scan_step_hint(), &opts.legacy(), DEFAULT_SCAN_STEP)
}

fn can_use_provider_search_path(opts: SearchOptsV2, policy: CorrectionPolicy) -> bool {
    policy == CorrectionPolicy::APPARENT && opts.scan_step_days.is_none()
}

fn labelled_crossings_for_altitude<F>(
    window: Interval<ModifiedJulianDate>,
    step: Days,
    altitude: &F,
    threshold: Radians,
    opts: SearchOptsV2,
) -> (Vec<intervals::LabeledCrossing>, bool)
where
    F: Fn(ModifiedJulianDate) -> Radians,
{
    let signal = |t: ModifiedJulianDate| -> f64 { altitude(t).sin() };
    let (labeled, start_above, _) = crate::event::search::crossings::find_labelled_crossings(
        window,
        step,
        &signal,
        threshold.sin(),
        opts,
    );
    (labeled, start_above)
}

fn crossing_events_from_labelled(labeled: &[intervals::LabeledCrossing]) -> Vec<CrossingEvent> {
    labeled
        .iter()
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

// ---------------------------------------------------------------------------
// Crossings
// ---------------------------------------------------------------------------

/// Find all threshold crossings of `target` altitude in the given `window`.
///
/// Returns a chronologically sorted list of [`CrossingEvent`]s.
///
/// # Arguments
/// * `target`, any body implementing [`AltitudeProvider`]
/// * `observer`, site on Earth
/// * `window`, search interval (MJD on TT axis)
/// * `threshold`, altitude threshold
/// * `opts`, search options (tolerances, scan step)
///
/// # Example
/// ```rust
/// use siderust::event::altitude::{crossings, SearchOpts};
/// use siderust::bodies::Sun;
/// use siderust::coordinates::centers::Geodetic;
/// use siderust::coordinates::frames::ECEF;
/// use siderust::time::{ModifiedJulianDate, MJD, Interval};
/// use siderust::qtty::*;
///
/// let site = Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
/// let window = Interval::new(
///     siderust::ModifiedJulianDate::new(60000.0),
///     siderust::ModifiedJulianDate::new(60001.0),
/// );
/// let events = crossings(&Sun, &site, window, Degrees::new(0.0), SearchOpts::default());
/// for e in events {
///     println!("{:?} at MJD {}", e.direction, e.mjd);
/// }
/// ```
///
/// # Returns
///
/// A chronologically sorted `Vec<CrossingEvent>` containing every
/// rising/setting transit found inside `window`.
pub fn crossings<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<CrossingEvent> {
    crossings_with_policy(
        target,
        observer,
        window,
        threshold,
        opts,
        CorrectionPolicy::APPARENT,
    )
}

/// Find threshold crossings using an explicit apparent-position correction
/// policy.
pub fn crossings_with_policy<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
    policy: CorrectionPolicy,
) -> Vec<CrossingEvent> {
    crossings_with_internal_opts_and_policy(
        target,
        observer,
        window,
        threshold,
        SearchOptsV2::from_legacy(opts),
        policy,
    )
}

fn crossings_with_internal_opts_and_policy<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOptsV2,
    policy: CorrectionPolicy,
) -> Vec<CrossingEvent> {
    if can_use_provider_search_path(opts, policy) {
        if let Some(events) = target.crossings_search(*observer, window, threshold, opts.legacy()) {
            return events;
        }
    }

    let f = make_altitude_fn(target, observer, policy);
    let thr_rad = threshold.to::<Radian>();
    let step = scan_step_for_opts(target, &opts);
    let (labeled, _) = labelled_crossings_for_altitude(window, step, &f, thr_rad, opts);
    crossing_events_from_labelled(&labeled)
}

// ---------------------------------------------------------------------------
// Culminations
// ---------------------------------------------------------------------------

/// Find all altitude culminations (local maxima and minima) of `target` in `window`.
///
/// Returns a chronologically sorted list of [`CulminationEvent`]s.
///
/// # Arguments
///
/// * `target`, any body implementing [`AltitudeProvider`]
/// * `observer`, geodetic observer site
/// * `window`, MJD/TT search interval
/// * `opts`, search precision (scan step + refinement tolerance)
///
/// # Returns
///
/// `Vec<CulminationEvent>` sorted by time, mixing upper (`Max`) and lower
/// (`Min`) culminations.
pub fn culminations<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    opts: SearchOpts,
) -> Vec<CulminationEvent> {
    culminations_with_policy(target, observer, window, opts, CorrectionPolicy::APPARENT)
}

/// Find altitude culminations using an explicit apparent-position correction
/// policy.
pub fn culminations_with_policy<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    opts: SearchOpts,
    policy: CorrectionPolicy,
) -> Vec<CulminationEvent> {
    let f = make_altitude_fn(target, observer, policy);
    // For culminations, use a slightly larger step (or the target's hint)
    let step = opts
        .scan_step_days
        .or_else(|| target.scan_step_hint())
        .unwrap_or(EXTREMA_SCAN_STEP);
    let tol = opts.time_tolerance;

    let raw: Vec<extrema::Extremum<Radian>> = extrema::find_extrema_tol(window, step, &f, tol);

    raw.iter()
        .map(|ext| {
            let alt_deg = ext.value.to::<Degree>();
            CulminationEvent {
                mjd: ext.t,
                altitude: alt_deg,
                kind: match ext.kind {
                    extrema::ExtremumKind::Maximum => CulminationKind::Max,
                    extrema::ExtremumKind::Minimum => CulminationKind::Min,
                },
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Altitude Ranges
// ---------------------------------------------------------------------------

/// Find all time intervals where the altitude of `target` is within
/// `[h_min, h_max]`.
///
/// Returns a sorted list of `Interval<ModifiedJulianDate>`.
///
/// # Algorithm
///
/// Uses labelled threshold crossings of `h_min` and `h_max`, then interval
/// algebra: `above(h_min) ∩ complement(above(h_max))`.
///
/// # Example
/// ```ignore
/// // Find astronomical night (Sun between -90° and -18°)
/// let nights = altitude_ranges(
///     &Sun, &site, window,
///     Degrees::new(-90.0), Degrees::new(-18.0),
///     SearchOpts::default(),
/// );
/// ```
///
/// # Arguments
///
/// * `target`, any body implementing [`AltitudeProvider`]
/// * `observer`, geodetic observer site
/// * `window`, MJD/TT search interval
/// * `h_min`, lower altitude bound (inclusive)
/// * `h_max`, upper altitude bound (inclusive)
/// * `opts`, search precision options
///
/// # Returns
///
/// Sorted, non‑overlapping `Vec<Interval<ModifiedJulianDate>>` covering the
/// time intervals where `h_min ≤ altitude(t) ≤ h_max`.
pub fn altitude_ranges<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    h_min: Degrees,
    h_max: Degrees,
    opts: SearchOpts,
) -> Vec<Interval<ModifiedJulianDate>> {
    altitude_ranges_with_policy(
        target,
        observer,
        window,
        h_min,
        h_max,
        opts,
        CorrectionPolicy::APPARENT,
    )
}

/// Find altitude-range periods using an explicit apparent-position correction
/// policy.
pub fn altitude_ranges_with_policy<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    h_min: Degrees,
    h_max: Degrees,
    opts: SearchOpts,
    policy: CorrectionPolicy,
) -> Vec<Interval<ModifiedJulianDate>> {
    altitude_ranges_with_internal_opts_and_policy(
        target,
        observer,
        window,
        h_min,
        h_max,
        SearchOptsV2::from_legacy(opts),
        policy,
    )
}

fn altitude_ranges_with_internal_opts_and_policy<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    h_min: Degrees,
    h_max: Degrees,
    opts: SearchOptsV2,
    policy: CorrectionPolicy,
) -> Vec<Interval<ModifiedJulianDate>> {
    if can_use_provider_search_path(opts, policy) {
        if let Some(periods) =
            target.altitude_range_search(*observer, window, h_min, h_max, opts.legacy())
        {
            return periods;
        }
    }

    let f = make_altitude_fn(target, observer, policy);
    let min_rad = h_min.to::<Radian>();
    let max_rad = h_max.to::<Radian>();
    let step = scan_step_for_opts(target, &opts);

    let (above_min, start_above_min) =
        labelled_crossings_for_altitude(window, step, &f, min_rad, opts);
    let above_min = intervals::build_above_periods_directed(&above_min, window, start_above_min);

    let (above_max, start_above_max) =
        labelled_crossings_for_altitude(window, step, &f, max_rad, opts);
    let above_max = intervals::build_above_periods_directed(&above_max, window, start_above_max);
    let below_max = complement_within(window, &above_max);
    intervals::intersect(&above_min, &below_max)
}

// ---------------------------------------------------------------------------
// Above/Below Threshold
// ---------------------------------------------------------------------------

/// Convenience: find periods where altitude is **above** a threshold.
///
/// Equivalent to `altitude_ranges(target, observer, window, threshold, 90°, opts)`.
///
/// # Arguments
///
/// * `target`, body implementing [`AltitudeProvider`]
/// * `observer`, geodetic site
/// * `window`, MJD/TT search interval
/// * `threshold`, altitude lower bound
/// * `opts`, search precision options
///
/// # Returns
///
/// Sorted, non‑overlapping `Vec<Interval<ModifiedJulianDate>>` covering the
/// times when `altitude(t) ≥ threshold`.
pub fn above_threshold<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<Interval<ModifiedJulianDate>> {
    above_threshold_with_policy(
        target,
        observer,
        window,
        threshold,
        opts,
        CorrectionPolicy::APPARENT,
    )
}

/// Find above-threshold periods using an explicit apparent-position
/// correction policy.
pub fn above_threshold_with_policy<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
    policy: CorrectionPolicy,
) -> Vec<Interval<ModifiedJulianDate>> {
    above_threshold_with_internal_opts_and_policy(
        target,
        observer,
        window,
        threshold,
        SearchOptsV2::from_legacy(opts),
        policy,
    )
}

fn above_threshold_with_internal_opts_and_policy<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOptsV2,
    policy: CorrectionPolicy,
) -> Vec<Interval<ModifiedJulianDate>> {
    if can_use_provider_search_path(opts, policy) {
        if let Some(periods) =
            target.above_threshold_search(*observer, window, threshold, opts.legacy())
        {
            return periods;
        }
    }

    let f = make_altitude_fn(target, observer, policy);
    let thr_rad = threshold.to::<Radian>();
    let step = scan_step_for_opts(target, &opts);
    let (labeled, start_above) = labelled_crossings_for_altitude(window, step, &f, thr_rad, opts);
    intervals::build_above_periods_directed(&labeled, window, start_above)
}

/// Convenience: find periods where altitude is **below** a threshold.
///
/// Equivalent to complement of [`above_threshold`].
///
/// # Arguments
///
/// * `target`, body implementing [`AltitudeProvider`]
/// * `observer`, geodetic site
/// * `window`, MJD/TT search interval
/// * `threshold`, altitude upper bound
/// * `opts`, search precision options
///
/// # Returns
///
/// Sorted, non‑overlapping `Vec<Interval<ModifiedJulianDate>>` covering the
/// times inside `window` when `altitude(t) < threshold`.
pub fn below_threshold<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<Interval<ModifiedJulianDate>> {
    below_threshold_with_policy(
        target,
        observer,
        window,
        threshold,
        opts,
        CorrectionPolicy::APPARENT,
    )
}

/// Find below-threshold periods using an explicit apparent-position
/// correction policy.
pub fn below_threshold_with_policy<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
    policy: CorrectionPolicy,
) -> Vec<Interval<ModifiedJulianDate>> {
    below_threshold_with_internal_opts_and_policy(
        target,
        observer,
        window,
        threshold,
        SearchOptsV2::from_legacy(opts),
        policy,
    )
}

fn below_threshold_with_internal_opts_and_policy<T: AltitudeProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOptsV2,
    policy: CorrectionPolicy,
) -> Vec<Interval<ModifiedJulianDate>> {
    if can_use_provider_search_path(opts, policy) {
        if let Some(periods) =
            target.below_threshold_search(*observer, window, threshold, opts.legacy())
        {
            return periods;
        }
    }

    let above = above_threshold_with_internal_opts_and_policy(
        target, observer, window, threshold, opts, policy,
    );
    complement_within(window, &above)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::solar_system::{Moon, Sun};

    fn greenwich() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(
            Degrees::new(0.0),
            Degrees::new(51.4769),
            Quantity::<Meter>::new(0.0),
        )
    }

    #[test]
    fn crossings_finds_sun_rise_set() {
        let site = greenwich();
        let mjd_start = crate::time::ModifiedJulianDate::new(60000.0);
        let mjd_end = crate::time::ModifiedJulianDate::new(60001.0);
        let window = Interval::new(mjd_start, mjd_end);

        let events = crossings(
            &Sun,
            &site,
            window,
            Degrees::new(0.0),
            SearchOpts::default(),
        );

        // In a normal 24h window at ~51°N, expect 1 rise + 1 set
        assert!(!events.is_empty(), "should find crossings");
        let rises = events
            .iter()
            .filter(|e| e.direction == CrossingDirection::Rising)
            .count();
        let sets = events
            .iter()
            .filter(|e| e.direction == CrossingDirection::Setting)
            .count();
        assert!(
            rises >= 1 || sets >= 1,
            "should find at least one rise or set"
        );
    }

    #[test]
    fn culminations_finds_sun_extrema() {
        let site = greenwich();
        let mjd_start = crate::time::ModifiedJulianDate::new(60000.0);
        let mjd_end = crate::time::ModifiedJulianDate::new(60001.0);
        let window = Interval::new(mjd_start, mjd_end);

        let culms = culminations(&Sun, &site, window, SearchOpts::default());

        // Expect upper and lower culmination in 24h
        assert!(!culms.is_empty(), "should find culminations");
        let maxima = culms
            .iter()
            .filter(|c| c.kind == CulminationKind::Max)
            .count();
        let minima = culms
            .iter()
            .filter(|c| c.kind == CulminationKind::Min)
            .count();
        assert!(maxima >= 1, "should find at least one upper culmination");
        assert!(minima >= 1, "should find at least one lower culmination");
    }

    #[test]
    fn above_threshold_sun_day_periods() {
        let site = greenwich();
        let mjd_start = crate::time::ModifiedJulianDate::new(60000.0);
        let mjd_end = crate::time::ModifiedJulianDate::new(60007.0);
        let window = Interval::new(mjd_start, mjd_end);

        let days = above_threshold(
            &Sun,
            &site,
            window,
            Degrees::new(0.0),
            SearchOpts::default(),
        );

        assert!(!days.is_empty(), "should find daytime periods in 7 days");
        for p in &days {
            assert!(p.length() > Days::new(0.0));
            assert!(p.length() < Days::new(1.0), "each day period < 24h");
        }
    }

    #[test]
    fn below_threshold_sun_night_periods() {
        let site = greenwich();
        let mjd_start = crate::time::ModifiedJulianDate::new(60000.0);
        let mjd_end = crate::time::ModifiedJulianDate::new(60007.0);
        let window = Interval::new(mjd_start, mjd_end);

        let nights = below_threshold(
            &Sun,
            &site,
            window,
            Degrees::new(-18.0), // astronomical twilight
            SearchOpts::default(),
        );

        assert!(!nights.is_empty(), "should find night periods");
    }

    #[test]
    fn public_sun_below_threshold_matches_solar_specialization() {
        let site = greenwich();
        let window = Interval::new(
            crate::time::ModifiedJulianDate::new(60000.0),
            crate::time::ModifiedJulianDate::new(60031.0),
        );

        let below = below_threshold(
            &Sun,
            &site,
            window,
            Degrees::new(-18.0),
            SearchOpts::default(),
        );
        let specialized = crate::event::solar::find_night_periods_with_search_opts(
            site,
            window,
            Degrees::new(-18.0),
            SearchOptsV2::default(),
        );

        assert_eq!(below.len(), specialized.len());
        for (actual, expected) in below.iter().zip(specialized.iter()) {
            assert!((actual.start.raw() - expected.start.raw()).abs() < Days::new(1e-6));
            assert!((actual.end.raw() - expected.end.raw()).abs() < Days::new(1e-6));
        }
    }

    #[test]
    fn altitude_ranges_twilight_band() {
        let site = greenwich();
        let mjd_start = crate::time::ModifiedJulianDate::new(60000.0);
        let mjd_end = crate::time::ModifiedJulianDate::new(60002.0);
        let window = Interval::new(mjd_start, mjd_end);

        let twilight = altitude_ranges(
            &Sun,
            &site,
            window,
            Degrees::new(-18.0),
            Degrees::new(-12.0),
            SearchOpts::default(),
        );

        // Should find nautical-to-astronomical twilight bands
        assert!(!twilight.is_empty(), "should find twilight bands");
    }

    #[test]
    fn moon_above_horizon_7_days() {
        let site = greenwich();
        let mjd_start = crate::time::ModifiedJulianDate::new(60000.0);
        let mjd_end = crate::time::ModifiedJulianDate::new(60007.0);
        let window = Interval::new(mjd_start, mjd_end);

        let periods = above_threshold(
            &Moon,
            &site,
            window,
            Degrees::new(0.0),
            SearchOpts::default(),
        );

        assert!(
            !periods.is_empty(),
            "should find moon-up periods over 7 days"
        );
    }
}
