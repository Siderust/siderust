// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Chebyshev-first event search integration tests and scan baseline comparisons.

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

fn chebyshev_opts() -> SearchOpts {
    SearchOpts::default()
}

fn scan_baseline_opts() -> SearchOpts {
    SearchOpts {
        scan_step_days: Some(Minutes::new(10.0).to::<Day>()),
        ..SearchOpts::default()
    }
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
        assert!(
            pair[0].end <= pair[1].start,
            "overlapping periods: {:?} vs {:?}",
            pair[0],
            pair[1]
        );
    }
}

fn assert_crossings_in_window(
    events: &[siderust::event::altitude::CrossingEvent],
    window: Interval<ModifiedJulianDate>,
) {
    for e in events {
        assert!(e.mjd >= window.start, "crossing before window");
        assert!(e.mjd <= window.end, "crossing after window");
    }
    for pair in events.windows(2) {
        assert!(pair[0].mjd <= pair[1].mjd, "crossings not sorted");
    }
}

#[test]
fn stable_public_api_exports_minimal_search_surface() {
    use siderust::bodies::Sun;
    use siderust::{
        above_threshold, altitude_ranges, below_threshold, crossings, culminations, SearchOpts,
    };
    let site = roque();
    let window = window_days(60000.0, 1.0);
    let opts = SearchOpts::default();
    let _ = crossings(&Sun, &site, window, Degrees::new(0.0), opts);
    let _ = above_threshold(&Sun, &site, window, Degrees::new(0.0), opts);
    let _ = below_threshold(&Sun, &site, window, Degrees::new(0.0), opts);
    let _ = altitude_ranges(
        &Sun,
        &site,
        window,
        Degrees::new(-6.0),
        Degrees::new(0.0),
        opts,
    );
    let _ = culminations(&Sun, &site, window, opts);
}

#[test]
fn solar_thresholds_match_scan_baseline_one_day() {
    let site = roque();
    let window = window_days(60000.0, 1.0);
    for threshold in [0.0, -6.0, -12.0, -18.0] {
        let thr = Degrees::new(threshold);
        let cheb = crossings(&Sun, &site, window, thr, chebyshev_opts());
        let scan = crossings(&Sun, &site, window, thr, scan_baseline_opts());
        assert_crossings_in_window(&cheb, window);
        assert_crossings_in_window(&scan, window);
        assert_eq!(
            cheb.len(),
            scan.len(),
            "threshold {threshold}: cheb={} scan={}",
            cheb.len(),
            scan.len()
        );
    }
}

#[test]
fn solar_above_below_range_match_scan_baseline_week() {
    let site = roque();
    let window = window_days(60000.0, 7.0);
    let cheb = chebyshev_opts();
    let scan = scan_baseline_opts();

    let above_cheb = above_threshold(&Sun, &site, window, Degrees::new(0.0), cheb);
    let above_scan = above_threshold(&Sun, &site, window, Degrees::new(0.0), scan);
    assert_periods_in_window(&above_cheb, window);
    assert_periods_in_window(&above_scan, window);
    assert_eq!(above_cheb.len(), above_scan.len());

    let below_cheb = below_threshold(&Sun, &site, window, Degrees::new(-18.0), cheb);
    let below_scan = below_threshold(&Sun, &site, window, Degrees::new(-18.0), scan);
    assert_periods_in_window(&below_cheb, window);
    assert_periods_in_window(&below_scan, window);
    assert_eq!(below_cheb.len(), below_scan.len());

    let range_cheb = altitude_ranges(
        &Sun,
        &site,
        window,
        Degrees::new(-18.0),
        Degrees::new(-12.0),
        cheb,
    );
    let range_scan = altitude_ranges(
        &Sun,
        &site,
        window,
        Degrees::new(-18.0),
        Degrees::new(-12.0),
        scan,
    );
    assert_periods_in_window(&range_cheb, window);
    assert_periods_in_window(&range_scan, window);
    assert_eq!(range_cheb.len(), range_scan.len());
}

#[test]
fn lunar_above_below_range_match_scan_baseline_week() {
    let site = roque();
    let window = window_days(60000.0, 7.0);
    let cheb = chebyshev_opts();
    let scan = scan_baseline_opts();

    let above_cheb = above_threshold(&Moon, &site, window, Degrees::new(0.0), cheb);
    let above_scan = above_threshold(&Moon, &site, window, Degrees::new(0.0), scan);
    assert_periods_in_window(&above_cheb, window);
    assert_eq!(above_cheb.len(), above_scan.len());

    let below_cheb = below_threshold(&Moon, &site, window, Degrees::new(0.0), cheb);
    let below_scan = below_threshold(&Moon, &site, window, Degrees::new(0.0), scan);
    assert_periods_in_window(&below_cheb, window);
    assert_eq!(below_cheb.len(), below_scan.len());

    let range_cheb = altitude_ranges(
        &Moon,
        &site,
        window,
        Degrees::new(-5.0),
        Degrees::new(45.0),
        cheb,
    );
    let range_scan = altitude_ranges(
        &Moon,
        &site,
        window,
        Degrees::new(-5.0),
        Degrees::new(45.0),
        scan,
    );
    assert_periods_in_window(&range_cheb, window);
    assert_eq!(range_cheb.len(), range_scan.len());
}

#[test]
fn long_windows_solar_below_astronomical_twilight() {
    let site = roque();
    let cheb = chebyshev_opts();
    let scan = scan_baseline_opts();
    for len in [30.0, 184.0] {
        let window = window_days(60000.0, len);
        let cheb_nights = below_threshold(&Sun, &site, window, Degrees::new(-18.0), cheb);
        let scan_nights = below_threshold(&Sun, &site, window, Degrees::new(-18.0), scan);
        assert_periods_in_window(&cheb_nights, window);
        assert_eq!(cheb_nights.len(), scan_nights.len(), "len={len}");
    }
}

#[test]
#[ignore = "slow: one-year window baseline comparison"]
fn long_window_one_year_solar_below_astronomical_twilight() {
    let site = roque();
    let window = window_days(60000.0, 365.0);
    let cheb_nights = below_threshold(&Sun, &site, window, Degrees::new(-18.0), chebyshev_opts());
    let scan_nights = below_threshold(
        &Sun,
        &site,
        window,
        Degrees::new(-18.0),
        scan_baseline_opts(),
    );
    assert_periods_in_window(&cheb_nights, window);
    assert_eq!(cheb_nights.len(), scan_nights.len());
}

#[test]
fn crossing_near_window_start() {
    let site = roque();
    let window = window_days(60000.0, 2.0);
    let events = crossings(&Sun, &site, window, Degrees::new(0.0), chebyshev_opts());
    assert_crossings_in_window(&events, window);
    assert!(!events.is_empty());
}

#[test]
fn explicit_scan_step_forces_scan_path() {
    let site = roque();
    let window = window_days(60000.0, 1.0);
    let scan = scan_baseline_opts();
    let events = crossings(&Sun, &site, window, Degrees::new(0.0), scan);
    assert_crossings_in_window(&events, window);
    assert_eq!(events.len(), 2);
}
