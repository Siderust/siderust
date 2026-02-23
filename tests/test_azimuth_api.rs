// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Integration tests for the unified azimuth calculus API.
//!
//! These tests verify the behaviour of `AzimuthProvider`, `azimuth_crossings`,
//! `azimuth_extrema`, and interval-query functions across Sun, Moon, and star
//! targets.

use qtty::*;
use siderust::bodies::catalog;
use siderust::bodies::solar_system::{Moon, Sun};
use siderust::calculus::azimuth::{
    azimuth_crossings, azimuth_extrema, azimuth_periods, azimuth_ranges, in_azimuth_range,
    outside_azimuth_range, AzimuthProvider, AzimuthQuery, SearchOpts,
};
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::coordinates::spherical::direction;
use siderust::time::{ModifiedJulianDate, Period, MJD};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn greenwich() -> Geodetic<ECEF> {
    Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.4769), Meters::new(0.0))
}

fn one_day() -> Period<MJD> {
    Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60001.0),
    )
}

fn one_week() -> Period<MJD> {
    Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60007.0),
    )
}

// ===========================================================================
// azimuth_at — valid range
// ===========================================================================

#[test]
fn sun_azimuth_at_in_valid_range() {
    let az = Sun.azimuth_at(&greenwich(), ModifiedJulianDate::new(60000.5));
    assert!(az.value() >= 0.0, "az must be ≥ 0 rad");
    assert!(az.value() < std::f64::consts::TAU, "az must be < 2π rad");
}

#[test]
fn moon_azimuth_at_in_valid_range() {
    let az = Moon.azimuth_at(&greenwich(), ModifiedJulianDate::new(60000.5));
    assert!(az.value() >= 0.0);
    assert!(az.value() < std::f64::consts::TAU);
}

#[test]
fn star_azimuth_at_in_valid_range() {
    let sirius = &catalog::SIRIUS;
    let az = sirius.azimuth_at(&greenwich(), ModifiedJulianDate::new(60000.5));
    assert!(az.value() >= 0.0);
    assert!(az.value() < std::f64::consts::TAU);
}

#[test]
fn star_and_icrs_direction_azimuth_agree() {
    let sirius = &catalog::SIRIUS;
    let dir = direction::ICRS::from(sirius);
    let mjd = ModifiedJulianDate::new(60000.3);
    let az_star = sirius.azimuth_at(&greenwich(), mjd);
    let az_dir = dir.azimuth_at(&greenwich(), mjd);
    assert!(
        (az_star.value() - az_dir.value()).abs() < 1e-10,
        "Star and direction::ICRS must agree on azimuth"
    );
}

// ===========================================================================
// azimuth_crossings
// ===========================================================================

#[test]
fn sun_crosses_south_once_per_day() {
    // At 51°N the Sun transits the meridian (South, az ≈ 180°) once per day.
    let events = azimuth_crossings(
        &Sun,
        &greenwich(),
        one_day(),
        Degrees::new(180.0),
        SearchOpts::default(),
    );
    assert!(
        !events.is_empty(),
        "should detect the Sun's southward meridian transit"
    );
}

#[test]
fn sun_crosses_east_or_west_in_24h() {
    let east = azimuth_crossings(
        &Sun,
        &greenwich(),
        one_day(),
        Degrees::new(90.0),
        SearchOpts::default(),
    );
    let west = azimuth_crossings(
        &Sun,
        &greenwich(),
        one_day(),
        Degrees::new(270.0),
        SearchOpts::default(),
    );
    assert!(
        !east.is_empty() || !west.is_empty(),
        "Sun should cross East or West bearing in 24h"
    );
}

#[test]
fn moon_crossing_south_over_week() {
    let events = azimuth_crossings(
        &Moon,
        &greenwich(),
        one_week(),
        Degrees::new(180.0),
        SearchOpts::default(),
    );
    assert!(
        !events.is_empty(),
        "Moon should cross South bearing in a week"
    );
}

// ===========================================================================
// azimuth_extrema — smoke tests
// ===========================================================================

/// Smoke test: `azimuth_extrema` runs without panicking and any returned
/// extremum has a valid azimuth value.  The Sun's azimuth is monotonically
/// increasing over a 24-hour window (Earth rotation dominates its apparent
/// motion), so the result may well be empty — but the values must be valid.
#[test]
fn sun_azimuth_extrema_smoke() {
    let exts = azimuth_extrema(&Sun, &greenwich(), one_day(), SearchOpts::default());
    for e in &exts {
        assert!(
            e.azimuth.value() >= 0.0 && e.azimuth.value() < 360.0,
            "extremum azimuth must be in [0°, 360°), got {}",
            e.azimuth.value()
        );
    }
}

#[test]
fn star_azimuth_extrema_smoke() {
    let sirius = &catalog::SIRIUS;
    let exts = azimuth_extrema(sirius, &greenwich(), one_day(), SearchOpts::default());
    for e in &exts {
        assert!(
            e.azimuth.value() >= 0.0 && e.azimuth.value() < 360.0,
            "extremum azimuth must be in [0°, 360°), got {}",
            e.azimuth.value()
        );
    }
}

// ===========================================================================
// azimuth_periods / azimuth_ranges
// ===========================================================================

#[test]
fn sun_in_eastern_half_non_empty() {
    let query = AzimuthQuery {
        observer: greenwich(),
        window: one_day(),
        min_azimuth: Degrees::new(90.0),
        max_azimuth: Degrees::new(270.0),
    };
    let periods = Sun.azimuth_periods(&query);
    assert!(
        !periods.is_empty(),
        "Sun should spend time in the eastern half (az 90–270°)"
    );
}

#[test]
fn azimuth_periods_free_fn_matches_trait() {
    let observer = greenwich();
    let window = one_day();
    let query = AzimuthQuery {
        observer,
        window,
        min_azimuth: Degrees::new(90.0),
        max_azimuth: Degrees::new(270.0),
    };
    let via_trait = Sun.azimuth_periods(&query);
    let via_fn = azimuth_periods(&Sun, &query);
    assert_eq!(
        via_trait.len(),
        via_fn.len(),
        "trait and free-function must agree on period count"
    );
}

#[test]
fn in_azimuth_range_equals_azimuth_ranges() {
    let observer = greenwich();
    let window = one_day();
    let via_ranges = azimuth_ranges(
        &Sun,
        &observer,
        window,
        Degrees::new(90.0),
        Degrees::new(270.0),
        SearchOpts::default(),
    );
    let via_in = in_azimuth_range(
        &Sun,
        &observer,
        window,
        Degrees::new(90.0),
        Degrees::new(270.0),
        SearchOpts::default(),
    );
    assert_eq!(
        via_ranges.len(),
        via_in.len(),
        "azimuth_ranges and in_azimuth_range must agree"
    );
}

// ===========================================================================
// outside_azimuth_range — complement property
// ===========================================================================

#[test]
fn outside_plus_inside_equals_window_sun() {
    let observer = greenwich();
    let window = one_day();
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
    let window_len = window.duration_days().value();
    assert!(
        (total - window_len).abs() < 1e-5,
        "inside ({total:.6}) + outside must equal window ({window_len:.6})"
    );
}

#[test]
fn outside_plus_inside_equals_window_moon() {
    let observer = greenwich();
    let window = one_week();
    let inside = in_azimuth_range(
        &Moon,
        &observer,
        window,
        Degrees::new(90.0),
        Degrees::new(270.0),
        SearchOpts::default(),
    );
    let outside = outside_azimuth_range(
        &Moon,
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
    let window_len = window.duration_days().value();
    assert!(
        (total - window_len).abs() < 1e-5,
        "Moon: inside + outside must equal window"
    );
}

// ===========================================================================
// Wrap-around range query
// ===========================================================================

#[test]
fn wrap_range_complement_covers_window() {
    // North-crossing arc 340° → 20° is wrap-around (min > max).
    let observer = greenwich();
    let window = one_week();
    let inside = in_azimuth_range(
        &Sun,
        &observer,
        window,
        Degrees::new(340.0),
        Degrees::new(20.0),
        SearchOpts::default(),
    );
    let outside = outside_azimuth_range(
        &Sun,
        &observer,
        window,
        Degrees::new(340.0),
        Degrees::new(20.0),
        SearchOpts::default(),
    );
    let total: f64 = inside
        .iter()
        .chain(outside.iter())
        .map(|p| p.duration_days().value())
        .sum();
    let window_len = window.duration_days().value();
    assert!(
        (total - window_len).abs() < 1e-5,
        "wrap-around inside + outside must cover the full window"
    );
}
