// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Integration tests for the trait‑based `AltitudePeriodsProvider` API.
//!
//! These tests verify that Sun, Moon, `Star`, and `direction::ICRS` all produce
//! consistent, correct results through the unified [`AltitudePeriodsProvider`]
//! interface while delegating to their respective `calculus` engines.

use siderust::bodies::catalog::{POLARIS, SIRIUS, VEGA};
use siderust::bodies::solar_system::{Moon, Sun};
use siderust::calculus::altitude::{altitude_periods, AltitudePeriodsProvider, AltitudeQuery};
use siderust::coordinates::centers::ObserverSite;
use siderust::coordinates::spherical::direction;
use siderust::time::{ModifiedJulianDate, Period};

use qtty::*;

// ===========================================================================
// Helpers
// ===========================================================================

/// Greenwich, ~51.48°N, 0°E, sea level
fn greenwich() -> ObserverSite {
    ObserverSite::new(Degrees::new(0.0), Degrees::new(51.4769), Meters::new(0.0))
}

/// La Palma (Roque de los Muchachos), ~28.76°N, -17.89°E, 2396 m
fn roque() -> ObserverSite {
    ObserverSite::new(
        Degrees::new(-17.892),
        Degrees::new(28.762),
        Meters::new(2396.0),
    )
}

/// Near the North Pole, ~89°N
fn north_pole() -> ObserverSite {
    ObserverSite::new(Degrees::new(0.0), Degrees::new(89.0), Meters::new(0.0))
}

fn one_day() -> Period<ModifiedJulianDate> {
    Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60001.0),
    )
}

fn one_week() -> Period<ModifiedJulianDate> {
    Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60007.0),
    )
}

/// Generic assertion helper: verifies basic structural invariants of periods.
fn assert_periods_valid(
    periods: &[Period<ModifiedJulianDate>],
    window: Period<ModifiedJulianDate>,
) {
    let win_start = window.start.value();
    let win_end = window.end.value();

    for (i, p) in periods.iter().enumerate() {
        // Positive duration
        assert!(
            p.duration_days() >= 0.0,
            "Period {} has negative duration: {:?}",
            i,
            p
        );
        // Within window
        assert!(
            p.start.value() >= win_start - 1e-9,
            "Period {} starts before window: {} < {}",
            i,
            p.start.value(),
            win_start
        );
        assert!(
            p.end.value() <= win_end + 1e-9,
            "Period {} ends after window: {} > {}",
            i,
            p.end.value(),
            win_end
        );
    }

    // Sorted and non-overlapping
    for w in periods.windows(2) {
        assert!(
            w[0].end.value() <= w[1].start.value() + 1e-9,
            "Periods overlap or out of order: {:?} vs {:?}",
            w[0],
            w[1]
        );
    }
}

// ===========================================================================
// Sun — via trait
// ===========================================================================

#[test]
fn sun_daytime_periods_trait() {
    let periods = Sun.above_threshold(greenwich(), one_day(), Degrees::new(0.0));
    assert_periods_valid(&periods, one_day());
    assert!(!periods.is_empty(), "Sun should rise at 51°N");
    for p in &periods {
        let hours = p.duration_days() * 24.0;
        assert!(
            hours > 5.0 && hours < 20.0,
            "Day length {} hours is unreasonable at 51°N in Feb",
            hours
        );
    }
}

#[test]
fn sun_night_periods_trait() {
    let nights = Sun.below_threshold(greenwich(), one_week(), Degrees::new(-18.0));
    assert_periods_valid(&nights, one_week());
    assert!(!nights.is_empty(), "Should find astronomical night at 51°N");
}

#[test]
fn sun_twilight_band_trait() {
    let query = AltitudeQuery {
        observer: greenwich(),
        window: Period::new(
            ModifiedJulianDate::new(60000.0),
            ModifiedJulianDate::new(60002.0),
        ),
        min_altitude: Degrees::new(-18.0),
        max_altitude: Degrees::new(-12.0),
    };
    let bands = Sun.altitude_periods(&query);
    assert_periods_valid(&bands, query.window);
    assert!(
        bands.len() >= 2,
        "Should find ≥2 twilight bands in 2 days, found {}",
        bands.len()
    );
}

// ===========================================================================
// Moon — via trait
// ===========================================================================

#[test]
fn moon_above_horizon_trait() {
    let periods = Moon.above_threshold(roque(), one_week(), Degrees::new(0.0));
    assert_periods_valid(&periods, one_week());
    assert!(
        !periods.is_empty(),
        "Moon should be above horizon in a week"
    );
}

#[test]
fn moon_below_horizon_trait() {
    let periods = Moon.below_threshold(roque(), one_week(), Degrees::new(0.0));
    assert_periods_valid(&periods, one_week());
    assert!(
        !periods.is_empty(),
        "Moon should be below horizon at some point in a week"
    );
}

// ===========================================================================
// Star — via trait
// ===========================================================================

#[test]
fn star_sirius_above_horizon_trait() {
    let periods = SIRIUS.above_threshold(greenwich(), one_day(), Degrees::new(0.0));
    assert_periods_valid(&periods, one_day());
    assert!(
        !periods.is_empty(),
        "Sirius should be above horizon for part of the day at 51°N"
    );
}

#[test]
fn star_vega_above_horizon_trait() {
    let periods = VEGA.above_threshold(greenwich(), one_day(), Degrees::new(0.0));
    assert_periods_valid(&periods, one_day());
    assert!(
        !periods.is_empty(),
        "Vega should be above horizon for part of the day at 51°N"
    );
}

// ===========================================================================
// direction::ICRS — lightweight path
// ===========================================================================

#[test]
fn icrs_direction_above_horizon() {
    let dir = direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716));
    let periods = dir.above_threshold(greenwich(), one_day(), Degrees::new(0.0));
    assert_periods_valid(&periods, one_day());
    assert!(!periods.is_empty());
}

#[test]
fn icrs_direction_matches_star() {
    let star = &SIRIUS;
    let dir = direction::ICRS::from(star);

    let window = one_day();
    let observer = greenwich();

    let star_periods = star.above_threshold(observer, window, Degrees::new(0.0));
    let dir_periods = dir.above_threshold(observer, window, Degrees::new(0.0));

    assert_eq!(
        star_periods.len(),
        dir_periods.len(),
        "Star and direction::ICRS should yield same count"
    );
    for (sp, dp) in star_periods.iter().zip(dir_periods.iter()) {
        assert!(
            (sp.start.value() - dp.start.value()).abs() < 1e-6,
            "Start mismatch"
        );
        assert!(
            (sp.end.value() - dp.end.value()).abs() < 1e-6,
            "End mismatch"
        );
    }
}

// ===========================================================================
// Generic free function
// ===========================================================================

#[test]
fn free_function_sun() {
    let query = AltitudeQuery {
        observer: greenwich(),
        window: one_day(),
        min_altitude: Degrees::new(0.0),
        max_altitude: Degrees::new(90.0),
    };
    let periods = altitude_periods(&Sun, &query);
    assert!(!periods.is_empty());
}

#[test]
fn free_function_moon() {
    let query = AltitudeQuery {
        observer: roque(),
        window: one_week(),
        min_altitude: Degrees::new(0.0),
        max_altitude: Degrees::new(90.0),
    };
    let periods = altitude_periods(&Moon, &query);
    assert!(!periods.is_empty());
}

#[test]
fn free_function_direction() {
    let dir = direction::ICRS::new(Degrees::new(279.2347), Degrees::new(38.7837)); // Vega
    let query = AltitudeQuery {
        observer: greenwich(),
        window: one_day(),
        min_altitude: Degrees::new(0.0),
        max_altitude: Degrees::new(90.0),
    };
    let periods = altitude_periods(&dir, &query);
    assert!(!periods.is_empty());
}

// ===========================================================================
// altitude_at single-point
// ===========================================================================

#[test]
fn altitude_at_sun_in_range() {
    let alt = Sun.altitude_at(&greenwich(), ModifiedJulianDate::new(51544.5));
    assert!(alt.abs() < std::f64::consts::FRAC_PI_2);
}

#[test]
fn altitude_at_moon_in_range() {
    let alt = Moon.altitude_at(&greenwich(), ModifiedJulianDate::new(51544.5));
    assert!(alt.abs() < std::f64::consts::FRAC_PI_2);
}

#[test]
fn altitude_at_star_in_range() {
    let alt = SIRIUS.altitude_at(&greenwich(), ModifiedJulianDate::new(51544.5));
    assert!(alt.abs() < std::f64::consts::FRAC_PI_2);
}

#[test]
fn altitude_at_icrs_direction_in_range() {
    let dir = direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716));
    let alt = dir.altitude_at(&greenwich(), ModifiedJulianDate::new(51544.5));
    assert!(alt.abs() < std::f64::consts::FRAC_PI_2);
}

// ===========================================================================
// Edge cases
// ===========================================================================

#[test]
fn empty_window_returns_empty() {
    let window = Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60000.0),
    );
    let query = AltitudeQuery {
        observer: greenwich(),
        window,
        min_altitude: Degrees::new(0.0),
        max_altitude: Degrees::new(90.0),
    };
    assert!(Sun.altitude_periods(&query).is_empty());
    assert!(Moon.altitude_periods(&query).is_empty());
    let dir = direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716));
    assert!(dir.altitude_periods(&query).is_empty());
}

#[test]
fn full_sky_range_returns_full_window() {
    let query = AltitudeQuery {
        observer: greenwich(),
        window: one_day(),
        min_altitude: Degrees::new(-90.0),
        max_altitude: Degrees::new(90.0),
    };
    let periods = Sun.altitude_periods(&query);
    assert!(!periods.is_empty());
    let total: f64 = periods.iter().map(|p| p.duration_days().value()).sum();
    assert!(
        (total - 1.0).abs() < 0.01,
        "Full sky range should span ~1 day, got {}",
        total
    );
}

#[test]
fn circumpolar_star_always_above() {
    let periods = POLARIS.above_threshold(greenwich(), one_day(), Degrees::new(0.0));
    assert_eq!(
        periods.len(),
        1,
        "Polaris should be continuously above horizon at 51°N"
    );
    assert!(
        (periods[0].duration_days() - Days::new(1.0)).abs() < 0.01,
        "Polaris period should span full day, got {} days",
        periods[0].duration_days()
    );
}

#[test]
fn circumpolar_star_never_below() {
    let periods = POLARIS.below_threshold(greenwich(), one_day(), Degrees::new(0.0));
    assert!(
        periods.is_empty(),
        "Polaris should never be below 0° at 51°N"
    );
}

#[test]
fn never_visible_star_at_north_pole() {
    // A star at Dec = −60° should never rise above horizon at 89°N
    let dir = direction::ICRS::new(Degrees::new(100.0), Degrees::new(-60.0));
    let periods = dir.above_threshold(north_pole(), one_day(), Degrees::new(0.0));
    assert!(
        periods.is_empty(),
        "Star at Dec=-60° should never be above horizon at 89°N"
    );
}

#[test]
fn periods_at_span_edges_are_clipped() {
    // Use a very short window; periods should not extend beyond it
    let window = Period::new(
        ModifiedJulianDate::new(60000.4),
        ModifiedJulianDate::new(60000.6),
    );
    let periods = Sun.above_threshold(greenwich(), window, Degrees::new(0.0));
    assert_periods_valid(&periods, window);
}

#[test]
fn above_and_below_cover_full_window() {
    let window = one_day();
    let threshold = Degrees::new(0.0);
    let observer = greenwich();

    let above = Sun.above_threshold(observer, window, threshold);
    let below = Sun.below_threshold(observer, window, threshold);

    let total_above: f64 = above.iter().map(|p| p.duration_days().value()).sum();
    let total_below: f64 = below.iter().map(|p| p.duration_days().value()).sum();
    let total = total_above + total_below;

    assert!(
        (total - 1.0).abs() < 0.02,
        "Above + below should cover full day: {} + {} = {}",
        total_above,
        total_below,
        total
    );
}
