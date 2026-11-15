// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Chebyshev-first event search integration tests and internal baseline comparisons.

use siderust::bench_internals;
use siderust::bodies::solar_system::{Moon, Sun};
use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::event::altitude::{
    above_threshold, altitude_ranges, below_threshold, crossings, CrossingDirection, SearchOpts,
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

#[test]
fn sun_above_threshold_default_matches_scan_baseline() {
    let site = roque();
    let window = window_days(60000.0, 30.0);
    let threshold = Degrees::new(0.0);
    let opts = default_opts();

    let default = above_threshold(&Sun, &site, window, threshold, opts);
    let scan = bench_internals::solar_above_threshold_scan_baseline(site, window, threshold, opts);
    assert_period_lists_close(&default, &scan);
}

#[test]
fn sun_chebyshev_baselines_run() {
    let site = roque();
    let window = window_days(60000.0, 7.0);
    let opts = default_opts();
    let threshold = Degrees::new(-18.0);

    let below =
        bench_internals::solar_below_threshold_chebyshev_baseline(site, window, threshold, opts);
    assert_periods_in_window(&below, window);

    let ranges = bench_internals::solar_altitude_ranges_chebyshev_baseline(
        site,
        window,
        Degrees::new(-18.0),
        Degrees::new(-12.0),
        opts,
    );
    assert_periods_in_window(&ranges, window);
}

#[test]
fn sun_crossings_scan_baseline_finds_events() {
    let site = roque();
    let window = window_days(60000.0, 1.0);
    let events = bench_internals::solar_crossings_scan_baseline(
        site,
        window,
        Degrees::new(0.0),
        default_opts(),
    );
    assert!(!events.is_empty());
}

#[test]
fn solar_daily_diagnostics_and_twilight_profile_run() {
    let site = roque();
    let window = window_days(60000.0, 3.0);
    let opts = default_opts();
    let threshold = Degrees::new(-6.0);

    let diag = bench_internals::solar_daily_diagnostics(site, window, threshold, opts);
    assert!(diag.precise_evaluations > 0);

    let profile = bench_internals::solar_twilight_profile(
        site,
        window,
        &[Degrees::new(-6.0), Degrees::new(-12.0)],
        opts,
    );
    assert_eq!(profile.len(), 2);
}

#[test]
fn lunar_scan_baselines_run() {
    let site = roque();
    let window = window_days(60000.0, 14.0);
    let opts = default_opts();
    let threshold = Degrees::new(0.0);

    let below = bench_internals::lunar_below_threshold_scan_baseline(site, window, threshold, opts);
    assert_periods_in_window(&below, window);

    let ranges = bench_internals::lunar_altitude_ranges_scan_baseline(
        site,
        window,
        Degrees::new(-5.0),
        Degrees::new(45.0),
        opts,
    );
    assert_periods_in_window(&ranges, window);
}

fn direct_moon_altitude(mjd: ModifiedJulianDate, site: Geodetic<ECEF>) -> Radians {
    Moon::get_horizontal::<Kilometer>(mjd.to::<siderust::time::JD>(), site)
        .alt()
        .to::<Radian>()
}

#[test]
fn cached_lunar_crossings_match_direct_altitude() {
    let shallow_site = Geodetic::<ECEF>::new(18.95 * DEG, 69.65 * DEG, 0.0 * M);
    // At the high-latitude site the Moon barely dips below the horizon.
    for (site, start, end) in [
        (ROQUE_DE_LOS_MUCHACHOS.geodetic(), 60_310.0, 60_311.0),
        (shallow_site, 60_335.4, 60_335.6),
    ] {
        let period = Interval::new(ModifiedJulianDate::new(start), ModifiedJulianDate::new(end));
        let events = crossings(&Moon, &site, period, 0.0 * DEG, default_opts());
        assert_eq!(
            events
                .iter()
                .map(|event| event.direction)
                .collect::<Vec<_>>(),
            [CrossingDirection::Setting, CrossingDirection::Rising],
        );

        // At a lunar range of 300,000 km, a one-metre interpolation
        // error subtends 0.00069 arcsec. At the shallow fixture's
        // slope (about 1 arcsec/s), a 2 ms bracket allows for that
        // pointing error and the numerical root tolerance.
        let margin = Seconds::new(0.002).to::<Day>();
        for event in &events {
            let before = ModifiedJulianDate::new((event.mjd.raw() - margin).value());
            let after = ModifiedJulianDate::new((event.mjd.raw() + margin).value());
            let direct_before = direct_moon_altitude(before, site);
            let direct_after = direct_moon_altitude(after, site);
            let sign = if event.direction == CrossingDirection::Rising {
                1.0
            } else {
                -1.0
            };
            assert!(
                direct_before * sign < 0.0 * RAD && direct_after * sign > 0.0 * RAD,
                "{event}: direct altitudes {direct_before} and {direct_after}",
            );
        }
    }
}

#[test]
fn cached_lunar_altitude_matches_direct_across_segment_boundary() {
    let start = ModifiedJulianDate::new(60_310.0);
    let end = ModifiedJulianDate::new(60_340.0);
    let polar = Geodetic::<ECEF>::new(18.95 * DEG, 69.65 * DEG, 0.0 * M);
    // A one-metre vector error at a conservative 300,000 km
    // topocentric range bounds the angular error for this fixture.
    let bound = Radians::new(1.0 / 300_000_000.0);
    for site in [ROQUE_DE_LOS_MUCHACHOS.geodetic(), polar] {
        let context = bench_internals::MoonAltitudeContext::new(start, end, site);
        for t in [
            60_310.0,
            60_312.999_999_99,
            60_313.0,
            60_313.000_000_01,
            60_320.0,
            60_335.0,
        ] {
            let mjd = ModifiedJulianDate::new(t);
            let cached = context.altitude_rad(mjd);
            let direct = direct_moon_altitude(mjd, site);
            assert!(
                (cached - direct).abs() < bound,
                "{mjd}: {cached} vs {direct}"
            );
        }
    }
}
