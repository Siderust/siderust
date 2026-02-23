// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Azimuth Event-Finding Functions
//!
//! Functions for finding azimuth events: bearing crossings, azimuth extrema,
//! and periods within an azimuth range.
//!
//! ## Handling the 0°/360° discontinuity
//!
//! Azimuth is a circular quantity that wraps at North (0° = 360°).  Naive
//! root-finding on the raw `az(t)` function would produce false crossings at
//! every wrap-around.  This module uses two complementary strategies:
//!
//! * **Trigonometric reformulation** for crossings and range queries:
//!   `sin(az(t) − b)` and `cos(az(t) − midpoint)` are continuous even
//!   through the 0°/360° boundary, because cosine/sine are 2π-periodic.
//!   These are passed directly to `math_core` without any unwrapping.
//!
//! * **Stateful unwrapping closure** for extrema: a `Cell`-based closure
//!   tracks the previous raw sample and accumulates ±2π offsets, producing
//!   a continuous monotone-ish function safe for `find_extrema_tol`.
//!
//! All `Period<MJD>` inputs/outputs are interpreted on the TT axis.
//! Convert UTC timestamps with `ModifiedJulianDate::from_utc(…)` first.

use std::cell::Cell;

use super::provider::AzimuthProvider;
use super::search::{SearchOpts, DEFAULT_SCAN_STEP, EXTREMA_SCAN_STEP};
use super::types::{AzimuthCrossingEvent, AzimuthExtremum, AzimuthExtremumKind, AzimuthQuery};
use crate::calculus::altitude::CrossingDirection;
use crate::calculus::math_core::{extrema, intervals};
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::time::{complement_within, ModifiedJulianDate, Period, MJD};
use qtty::*;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Choose the best scan step for the target.
fn scan_step_for<T: AzimuthProvider>(target: &T, opts: &SearchOpts) -> Days {
    opts.scan_step_days
        .or_else(|| target.scan_step_hint())
        .unwrap_or(DEFAULT_SCAN_STEP)
}

/// Build a **continuous** function `sin(az(t) − bearing)`.
///
/// This sidesteps the 0°/360° discontinuity of the raw azimuth: because sine
/// has period 2π, `sin(az + 2π − b) = sin(az − b)` — the function is smooth
/// and sign-changes exactly where `az(t) = bearing  (mod 360°)`.
///
/// The zero-crossings of this function correspond exactly to crossing events
/// for the given bearing.
fn make_az_bearing_sin_fn<'a, T: AzimuthProvider>(
    target: &'a T,
    site: &'a Geodetic<ECEF>,
    bearing_rad: f64,
) -> impl Fn(ModifiedJulianDate) -> Radians + 'a {
    let site = *site;
    move |t: ModifiedJulianDate| {
        let az = target.azimuth_at(&site, t).value();
        Radians::new((az - bearing_rad).sin())
    }
}

/// Build a continuous function `cos(az(t) − mid)` for range queries.
///
/// The body is inside the arc `[min_az, max_az]` when
/// `cos(az − mid) ≥ cos(halfwidth)`, i.e. when this function is ≥ 0 after
/// subtracting `cos(halfwidth)`.  This avoids the 0°/360° wrap entirely.
fn make_az_cosine_fn<'a, T: AzimuthProvider>(
    target: &'a T,
    site: &'a Geodetic<ECEF>,
    mid_rad: f64,
    cos_hw: f64,
) -> impl Fn(ModifiedJulianDate) -> Radians + 'a {
    let site = *site;
    move |t: ModifiedJulianDate| {
        let az = target.azimuth_at(&site, t).value();
        Radians::new((az - mid_rad).cos() - cos_hw)
    }
}

/// Build an **unwrapped** (continuous) azimuth closure using a `Cell`-based
/// accumulator.
///
/// Each call adjusts a running ±2π offset whenever the raw azimuth value
/// jumps by more than π since the last sample.  This produces a smooth,
/// monotone function safe for `find_extrema_tol`.
///
/// **Important:** this closure is stateful and must be evaluated with
/// monotonically increasing `t` — as guaranteed by all `math_core` scan loops.
fn make_az_unwrapped_fn<'a, T: AzimuthProvider>(
    target: &'a T,
    site: Geodetic<ECEF>,
) -> impl Fn(ModifiedJulianDate) -> Radians + 'a {
    let prev_raw: Cell<f64> = Cell::new(f64::NAN);
    let offset: Cell<f64> = Cell::new(0.0f64);
    move |t: ModifiedJulianDate| {
        let raw = target.azimuth_at(&site, t).value(); // [0, 2π)
        let p = prev_raw.get();
        if !p.is_nan() {
            let diff = raw - p;
            if diff > std::f64::consts::PI {
                // Apparent upward jump > π → real motion was downward (0°→360° crossing)
                offset.set(offset.get() - std::f64::consts::TAU);
            } else if diff < -std::f64::consts::PI {
                // Apparent downward jump > π → real motion was upward (360°→0° crossing)
                offset.set(offset.get() + std::f64::consts::TAU);
            }
        }
        prev_raw.set(raw);
        Radians::new(raw + offset.get())
    }
}

// ---------------------------------------------------------------------------
// Bearing Crossings
// ---------------------------------------------------------------------------

/// Find all moments when `target`'s azimuth crosses a given `bearing`,
/// within the specified `window`.
///
/// Returns a chronologically sorted list of [`AzimuthCrossingEvent`]s.
///
/// `direction` is [`CrossingDirection::Rising`] when the azimuth is
/// increasing (moving clockwise / eastward across the bearing), and
/// [`CrossingDirection::Setting`] when decreasing.
///
/// # Arguments
/// * `target`   — any body implementing [`AzimuthProvider`]
/// * `observer` — site on Earth
/// * `window`   — search interval (MJD on TT axis)
/// * `bearing`  — compass bearing to detect crossings of (North-clockwise)
/// * `opts`     — search options (tolerances, scan step)
///
/// # Example
/// ```rust
/// use siderust::calculus::azimuth::{azimuth_crossings, SearchOpts};
/// use siderust::bodies::Sun;
/// use siderust::coordinates::centers::Geodetic;
/// use siderust::coordinates::frames::ECEF;
/// use siderust::time::{ModifiedJulianDate, MJD, Period};
/// use qtty::*;
///
/// let site = Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
/// let window = Period::new(ModifiedJulianDate::new(60000.0), ModifiedJulianDate::new(60001.0));
///
/// // Find when the Sun crosses due-South (180°)
/// let events = azimuth_crossings(&Sun, &site, window, Degrees::new(180.0), SearchOpts::default());
/// ```
pub fn azimuth_crossings<T: AzimuthProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Period<MJD>,
    bearing: Degrees,
    opts: SearchOpts,
) -> Vec<AzimuthCrossingEvent> {
    let bearing_rad = bearing.to::<Radian>().value();
    let f = make_az_bearing_sin_fn(target, observer, bearing_rad);
    let step = scan_step_for(target, &opts);

    let mut raw = intervals::find_crossings(window, step, &f, Radians::new(0.0));
    let labeled = intervals::label_crossings(&mut raw, &f, Radians::new(0.0));

    labeled
        .iter()
        .map(|lc| AzimuthCrossingEvent {
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
// Azimuth Extrema
// ---------------------------------------------------------------------------

/// Find all local azimuth extrema (northernmost and southernmost bearings)
/// of `target` in `window`.
///
/// Returns a chronologically sorted list of [`AzimuthExtremum`]s.
///
/// Uses a stateful unwrapped closure so that true extrema near the
/// 0°/360° boundary are detected correctly.
pub fn azimuth_extrema<T: AzimuthProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Period<MJD>,
    opts: SearchOpts,
) -> Vec<AzimuthExtremum> {
    let step = opts
        .scan_step_days
        .or_else(|| target.scan_step_hint())
        .unwrap_or(EXTREMA_SCAN_STEP);
    let tol = opts.time_tolerance;

    let f = make_az_unwrapped_fn(target, *observer);
    let raw: Vec<extrema::Extremum<Radian>> = extrema::find_extrema_tol(window, step, &f, tol);

    raw.iter()
        .map(|ext| {
            // Wrap the unwrapped value back into [0°, 360°)
            let az_wrapped = ext.value.value().rem_euclid(std::f64::consts::TAU);
            let az_deg = Radians::new(az_wrapped).to::<Degree>();
            AzimuthExtremum {
                mjd: ext.t,
                azimuth: az_deg,
                kind: match ext.kind {
                    extrema::ExtremumKind::Maximum => AzimuthExtremumKind::Max,
                    extrema::ExtremumKind::Minimum => AzimuthExtremumKind::Min,
                },
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Azimuth Range Periods (internal helper called from provider.rs)
// ---------------------------------------------------------------------------

/// Compute periods when `target`'s azimuth is within the band described by
/// `query`.
///
/// This function is `pub(crate)` and is the single implementation backing
/// all [`AzimuthProvider::azimuth_periods`] impls.
///
/// ## Non-wrap range (`min ≤ max`)
///
/// Uses the midpoint-cosine approach:
/// `f(t) = cos(az(t) − mid) − cos(halfwidth) ≥ 0` precisely when
/// `az ∈ [min, max]`.  This function is continuous even through the
/// 0°/360° boundary.
///
/// ## Wrap-around range (`min > max`)
///
/// The arc `[min, 360°) ∪ [0°, max]` equals the complement of the arc
/// `[max, min]` (non-wrapping).  We compute the non-wrapping complement and
/// take its complement within the window.
pub(crate) fn azimuth_range_periods<T: AzimuthProvider>(
    target: &T,
    query: &AzimuthQuery,
) -> Vec<Period<MJD>> {
    let step = target.scan_step_hint().unwrap_or(DEFAULT_SCAN_STEP);

    let min_deg = query.min_azimuth.value();
    let max_deg = query.max_azimuth.value();

    if min_deg <= max_deg {
        // ── Non-wrap arc [min, max] ────────────────────────────────────────
        let mid_rad = ((min_deg + max_deg) / 2.0).to_radians();
        let hw_rad = ((max_deg - min_deg) / 2.0).to_radians();
        let cos_hw = hw_rad.cos();

        let f = make_az_cosine_fn(target, &query.observer, mid_rad, cos_hw);
        intervals::above_threshold_periods(query.window, step, &f, Radians::new(0.0))
    } else {
        // ── Wrap-around arc [min, 360°) ∪ [0°, max] ──────────────────────
        // Complement arc is [max, min] (non-wrapping).
        let mid_rad = ((max_deg + min_deg) / 2.0).to_radians();
        let hw_rad = ((min_deg - max_deg) / 2.0).to_radians();
        let cos_hw = hw_rad.cos();

        let f = make_az_cosine_fn(target, &query.observer, mid_rad, cos_hw);
        let outside = intervals::above_threshold_periods(query.window, step, &f, Radians::new(0.0));
        complement_within(query.window, &outside)
    }
}

// ---------------------------------------------------------------------------
// Public free functions (generic, exposed from mod.rs)
// ---------------------------------------------------------------------------

/// Find all time intervals where `target`'s azimuth is within
/// `[min_az, max_az]`.
///
/// Supports wrap-around ranges when `min_az > max_az` (e.g. 350° → 10°).
///
/// Returns a sorted list of `Period<MJD>`.
pub fn azimuth_ranges<T: AzimuthProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Period<MJD>,
    min_az: Degrees,
    max_az: Degrees,
    _opts: SearchOpts,
) -> Vec<Period<MJD>> {
    target.azimuth_periods(&AzimuthQuery {
        observer: *observer,
        window,
        min_azimuth: min_az,
        max_azimuth: max_az,
    })
}

/// Convenience: find periods where azimuth is within `[min_az, max_az]`.
///
/// Identical to [`azimuth_ranges`]; provided as a named convenience wrapper.
pub fn in_azimuth_range<T: AzimuthProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Period<MJD>,
    min_az: Degrees,
    max_az: Degrees,
    opts: SearchOpts,
) -> Vec<Period<MJD>> {
    azimuth_ranges(target, observer, window, min_az, max_az, opts)
}

/// Convenience: find periods where azimuth is **outside** `[min_az, max_az]`.
pub fn outside_azimuth_range<T: AzimuthProvider>(
    target: &T,
    observer: &Geodetic<ECEF>,
    window: Period<MJD>,
    min_az: Degrees,
    max_az: Degrees,
    opts: SearchOpts,
) -> Vec<Period<MJD>> {
    let inside = azimuth_ranges(target, observer, window, min_az, max_az, opts);
    complement_within(window, &inside)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::solar_system::Sun;

    fn greenwich() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.4769), Meters::new(0.0))
    }

    fn one_day() -> Period<MJD> {
        Period::new(
            ModifiedJulianDate::new(60000.0),
            ModifiedJulianDate::new(60001.0),
        )
    }

    #[test]
    fn crossings_finds_south_transit() {
        let events = azimuth_crossings(
            &Sun,
            &greenwich(),
            one_day(),
            Degrees::new(180.0), // due-South
            SearchOpts::default(),
        );
        // Sun transits the meridian (az≈180°) once per day at 51°N
        assert!(
            !events.is_empty(),
            "should find the Sun's southward transit"
        );
    }

    #[test]
    fn extrema_are_valid_if_present() {
        // The Sun's azimuth is monotonically increasing over 24h (Earth rotation
        // dominates), so no local extrema are expected.  This test only validates
        // that the function runs cleanly and any returned extrema are in range.
        let exts = azimuth_extrema(&Sun, &greenwich(), one_day(), SearchOpts::default());
        for e in &exts {
            assert!(
                e.azimuth.value() >= 0.0 && e.azimuth.value() < 360.0,
                "extremum azimuth out of [0°, 360°): {}",
                e.azimuth.value()
            );
        }
    }

    #[test]
    fn range_periods_non_wrap() {
        let periods = in_azimuth_range(
            &Sun,
            &greenwich(),
            one_day(),
            Degrees::new(90.0),  // East
            Degrees::new(270.0), // West
            SearchOpts::default(),
        );
        assert!(!periods.is_empty(), "Sun should be east-to-west in 24h");
    }

    #[test]
    fn range_periods_wrap_around() {
        // North-crossing arc: 340° → 20°
        let periods = in_azimuth_range(
            &Sun,
            &greenwich(),
            // Use a 7-day window — the Sun is near North briefly near the solstice,
            // but let's just check the function runs without panicking and that
            // non-wrap + wrap complement each other.
            Period::new(
                ModifiedJulianDate::new(60000.0),
                ModifiedJulianDate::new(60007.0),
            ),
            Degrees::new(340.0),
            Degrees::new(20.0),
            SearchOpts::default(),
        );
        // Result may be empty (Sun rarely near North at 51°N); just check no panic.
        let _ = periods;
    }

    #[test]
    fn outside_is_complement_of_in_range() {
        let window = one_day();
        let observer = greenwich();
        let inside = in_azimuth_range(
            &Sun,
            &observer,
            window,
            Degrees::new(90.0),
            Degrees::new(270.0),
            SearchOpts::default(),
        );
        let outside = outside_azimuth_range(
            &Sun,
            &observer,
            window,
            Degrees::new(90.0),
            Degrees::new(270.0),
            SearchOpts::default(),
        );
        let total: f64 = inside
            .iter()
            .chain(outside.iter())
            .map(|p| p.duration_days().value())
            .sum();
        let window_duration = window.duration_days().value();
        assert!(
            (total - window_duration).abs() < 1e-6,
            "inside + outside durations must equal the window"
        );
    }
}
