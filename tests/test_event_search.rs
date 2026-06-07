// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Chebyshev-first event search integration tests and internal baseline comparisons.

use siderust::bench_internals;
use siderust::bodies::solar_system::{Moon, Sun};
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::event::altitude::{
    above_threshold, altitude_ranges, below_threshold, crossings, SearchOpts,
};
use siderust::qtty::*;
use siderust::time::{Interval, ModifiedJulianDate};

fn roque() -> Geodetic<ECEF> {
    Geodetic::<ECEF>::new(
        Degrees::new(-17.892),
        Degrees::new(28.762),
        Meters::new(2396.0),
    )
}

fn window_days(start: f64, len: f64) -> Interval<ModifiedJulianDate> {
    Interval::new(
        ModifiedJulianDate::try_new(Days::new(start)).unwrap(),
        ModifiedJulianDate::try_new(Days::new(start + len)).unwrap(),
    )
}

fn default_opts() -> SearchOpts {
    SearchOpts::default()
}

fn assert_periods_in_window(
    periods: &[Interval<ModifiedJulianDate>],
    window: Interval<ModifiedJulianDate>,
) {
    for p in periods {
        assert!(p.start >= window.start, "period starts before window");
        assert!(p.end <= window.end, "period ends after window");
        assert!(p.end > p.start, "non-positive period");
    }
    for pair in periods.windows(2) {
        assert!(pair[0].end <= pair[1].start, "overlapping periods");
    }
}

fn assert_period_lists_close(
    actual: &[Interval<ModifiedJulianDate>],
    expected: &[Interval<ModifiedJulianDate>],
) {
    assert_eq!(actual.len(), expected.len(), "period count mismatch");
    for (a, e) in actual.iter().zip(expected.iter()) {
        assert!((a.start.raw() - e.start.raw()).abs() < Days::new(1e-5));
        assert!((a.end.raw() - e.end.raw()).abs() < Days::new(1e-5));
    }
}

#[test]
fn sun_below_threshold_default_matches_scan_baseline() {
    let site = roque();
    let window = window_days(60000.0, 30.0);
    let threshold = Degrees::new(-18.0);
    let opts = default_opts();

    let default = below_threshold(&Sun, &site, window, threshold, opts);
    let scan = bench_internals::solar_below_threshold_scan_baseline(site, window, threshold, opts);
    assert_period_lists_close(&default, &scan);
}

#[test]
fn sun_altitude_ranges_default_matches_scan_baseline() {
    let site = roque();
    let window = window_days(60000.0, 7.0);
    let opts = default_opts();

    let default = altitude_ranges(
        &Sun,
        &site,
        window,
        Degrees::new(-18.0),
        Degrees::new(-12.0),
        opts,
    );
    let scan = bench_internals::solar_altitude_ranges_scan_baseline(
        site,
        window,
        Degrees::new(-18.0),
        Degrees::new(-12.0),
        opts,
    );
    assert_period_lists_close(&default, &scan);
}

#[test]
fn moon_above_threshold_default_matches_scan_baseline() {
    let site = roque();
    let window = window_days(60000.0, 30.0);
    let threshold = Degrees::new(0.0);
    let opts = default_opts();

    let default = above_threshold(&Moon, &site, window, threshold, opts);
    let scan = bench_internals::lunar_above_threshold_scan_baseline(site, window, threshold, opts);
    assert_period_lists_close(&default, &scan);
}

#[test]
fn sun_crossings_finds_events_in_day() {
    let site = roque();
    let window = window_days(60000.0, 1.0);
    let events = crossings(&Sun, &site, window, Degrees::new(0.0), default_opts());
    assert!(!events.is_empty());
}

#[test]
fn periods_respect_window_bounds() {
    let site = roque();
    let window = window_days(60000.0, 30.0);
    let periods = below_threshold(&Sun, &site, window, Degrees::new(-18.0), default_opts());
    assert_periods_in_window(&periods, window);
}
