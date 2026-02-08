// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Event Finding Functions
//!
//! Functions for finding altitude events: crossings, culminations,
//! and periods above/below thresholds or within ranges.

use super::provider::AltitudePeriodsProvider;
use super::search::{SearchOpts, DEFAULT_SCAN_STEP, EXTREMA_SCAN_STEP};
use super::types::{CrossingDirection, CrossingEvent, CulminationEvent, CulminationKind};
use crate::calculus::math_core::{extrema, intervals};
use crate::coordinates::centers::ObserverSite;
use crate::time::{complement_within, ModifiedJulianDate, Period};
use qtty::*;

// ---------------------------------------------------------------------------
// Internal: build altitude function from trait
// ---------------------------------------------------------------------------

/// Build an altitude function from any `AltitudePeriodsProvider`.
fn make_altitude_fn<'a, T: AltitudePeriodsProvider>(
    target: &'a T,
    site: &'a ObserverSite,
) -> impl Fn(ModifiedJulianDate) -> Radians + 'a {
    let site = *site;
    move |t: ModifiedJulianDate| target.altitude_at(&site, t)
}

/// Choose the best scan step for the target.
fn scan_step_for<T: AltitudePeriodsProvider>(target: &T, opts: &SearchOpts) -> Days {
    opts.scan_step_days
        .or_else(|| target.scan_step_hint())
        .unwrap_or(DEFAULT_SCAN_STEP)
}

// ---------------------------------------------------------------------------
// Crossings
// ---------------------------------------------------------------------------

/// Find all threshold crossings of `target` altitude in the given `window`.
///
/// Returns a chronologically sorted list of [`CrossingEvent`]s.
///
/// # Arguments
/// * `target` — any body implementing [`AltitudePeriodsProvider`]
/// * `observer` — site on Earth
/// * `window` — search interval (Mjd)
/// * `threshold` — altitude threshold
/// * `opts` — search options (tolerances, scan step)
///
/// # Example
/// ```rust
/// use siderust::calculus::altitude::{crossings, SearchOpts};
/// use siderust::bodies::Sun;
/// use siderust::coordinates::centers::ObserverSite;
/// use siderust::time::{ModifiedJulianDate, Period};
/// use qtty::*;
///
/// let site = ObserverSite::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
/// let window = Period::new(ModifiedJulianDate::new(60000.0), ModifiedJulianDate::new(60001.0));
/// let events = crossings(&Sun, &site, window, Degrees::new(0.0), SearchOpts::default());
/// for e in events {
///     println!("{:?} at MJD {}", e.direction, e.mjd);
/// }
/// ```
pub fn crossings<T: AltitudePeriodsProvider>(
    target: &T,
    observer: &ObserverSite,
    window: Period<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<CrossingEvent> {
    let f = make_altitude_fn(target, observer);
    let thr_rad = threshold.to::<Radian>();
    let step = scan_step_for(target, &opts);

    // Use the fast scan + label approach from math_core
    let mut raw_crossings = intervals::find_crossings(window, step, &f, thr_rad);
    let labeled = intervals::label_crossings(&mut raw_crossings, &f, thr_rad);

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
// Culminations
// ---------------------------------------------------------------------------

/// Find all altitude culminations (local maxima and minima) of `target` in `window`.
///
/// Returns a chronologically sorted list of [`CulminationEvent`]s.
pub fn culminations<T: AltitudePeriodsProvider>(
    target: &T,
    observer: &ObserverSite,
    window: Period<ModifiedJulianDate>,
    opts: SearchOpts,
) -> Vec<CulminationEvent> {
    let f = make_altitude_fn(target, observer);
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
                jd: ext.t,
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
/// Returns a sorted list of `Period<ModifiedJulianDate>`.
///
/// # Algorithm
///
/// Uses the two‑stage approach:
/// 1. Fast coarse scan to find threshold crossings of `h_min` and `h_max`.
/// 2. Brent refinement for each bracket.
/// 3. Interval algebra: `above(h_min) ∩ complement(above(h_max))`.
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
pub fn altitude_ranges<T: AltitudePeriodsProvider>(
    target: &T,
    observer: &ObserverSite,
    window: Period<ModifiedJulianDate>,
    h_min: Degrees,
    h_max: Degrees,
    opts: SearchOpts,
) -> Vec<Period<ModifiedJulianDate>> {
    let f = make_altitude_fn(target, observer);
    let min_rad = h_min.to::<Radian>();
    let max_rad = h_max.to::<Radian>();
    let step = scan_step_for(target, &opts);

    intervals::in_range_periods(window, step, &f, min_rad, max_rad)
}

// ---------------------------------------------------------------------------
// Above/Below Threshold
// ---------------------------------------------------------------------------

/// Convenience: find periods where altitude is **above** a threshold.
///
/// Equivalent to `altitude_ranges(target, observer, window, threshold, 90°, opts)`.
pub fn above_threshold<T: AltitudePeriodsProvider>(
    target: &T,
    observer: &ObserverSite,
    window: Period<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<Period<ModifiedJulianDate>> {
    let f = make_altitude_fn(target, observer);
    let thr_rad = threshold.to::<Radian>();
    let step = scan_step_for(target, &opts);

    intervals::above_threshold_periods(window, step, &f, thr_rad)
}

/// Convenience: find periods where altitude is **below** a threshold.
///
/// Equivalent to complement of [`above_threshold`].
pub fn below_threshold<T: AltitudePeriodsProvider>(
    target: &T,
    observer: &ObserverSite,
    window: Period<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<Period<ModifiedJulianDate>> {
    let above = above_threshold(target, observer, window, threshold, opts);
    complement_within(window, &above)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::solar_system::{Moon, Sun};

    fn greenwich() -> ObserverSite {
        ObserverSite::new(
            Degrees::new(0.0),
            Degrees::new(51.4769),
            Quantity::<Meter>::new(0.0),
        )
    }

    #[test]
    fn crossings_finds_sun_rise_set() {
        let site = greenwich();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60001.0);
        let window = Period::new(mjd_start, mjd_end);

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
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60001.0);
        let window = Period::new(mjd_start, mjd_end);

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
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);
        let window = Period::new(mjd_start, mjd_end);

        let days = above_threshold(
            &Sun,
            &site,
            window,
            Degrees::new(0.0),
            SearchOpts::default(),
        );

        assert!(!days.is_empty(), "should find daytime periods in 7 days");
        for p in &days {
            assert!(p.duration_days() > 0.0);
            assert!(p.duration_days() < 1.0, "each day period < 24h");
        }
    }

    #[test]
    fn below_threshold_sun_night_periods() {
        let site = greenwich();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);
        let window = Period::new(mjd_start, mjd_end);

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
    fn altitude_ranges_twilight_band() {
        let site = greenwich();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60002.0);
        let window = Period::new(mjd_start, mjd_end);

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
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);
        let window = Period::new(mjd_start, mjd_end);

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
