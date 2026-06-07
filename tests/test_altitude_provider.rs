// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Integration tests for the standardized altitude period API.

use siderust::bodies::catalog::{POLARIS, SIRIUS, VEGA};
use siderust::bodies::solar_system::{Moon, Sun};
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::coordinates::spherical::direction;
use siderust::event::altitude::{
    above_threshold, altitude_ranges, below_threshold, AltitudeProvider, SearchOpts,
};
use siderust::time::{Interval, ModifiedJulianDate};

use siderust::qtty::*;

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

fn north_pole() -> Geodetic<ECEF> {
    Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(89.0), Meters::new(0.0))
}

fn one_day() -> Interval<ModifiedJulianDate> {
    Interval::new(
        ModifiedJulianDate::try_new(Days::new(60000.0)).unwrap(),
        ModifiedJulianDate::try_new(Days::new(60001.0)).unwrap(),
    )
}

fn one_week() -> Interval<ModifiedJulianDate> {
    Interval::new(
        ModifiedJulianDate::try_new(Days::new(60000.0)).unwrap(),
        ModifiedJulianDate::try_new(Days::new(60007.0)).unwrap(),
    )
}

fn assert_periods_valid(
    periods: &[Interval<ModifiedJulianDate>],
    window: Interval<ModifiedJulianDate>,
) {
    let win_start = window.start.raw().value();
    let win_end = window.end.raw().value();

    for (i, p) in periods.iter().enumerate() {
        assert!(
            ((p).end.raw() - (p).start.raw()) >= Days::new(0.0),
            "Period {} has negative duration: {:?}",
            i,
            p
        );
        assert!(
            p.start.raw().value() >= win_start - 1e-9,
            "Period {} starts before window",
            i
        );
        assert!(
            p.end.raw().value() <= win_end + 1e-9,
            "Period {} ends after window",
            i
        );
    }

    for w in periods.windows(2) {
        assert!(
            w[0].end.raw().value() <= w[1].start.raw().value() + 1e-9,
            "Periods overlap or out of order"
        );
    }
}

#[test]
fn sun_daytime_periods() {
    let periods = above_threshold(
        &Sun,
        &greenwich(),
        one_day(),
        Degrees::new(0.0),
        SearchOpts::default(),
    );
    assert_periods_valid(&periods, one_day());
    assert!(!periods.is_empty(), "Sun should rise at 51°N");
}

#[test]
fn sun_night_periods() {
    let nights = below_threshold(
        &Sun,
        &greenwich(),
        one_week(),
        Degrees::new(-18.0),
        SearchOpts::default(),
    );
    assert_periods_valid(&nights, one_week());
    assert!(!nights.is_empty(), "Should find astronomical night at 51°N");
}

#[test]
fn sun_twilight_band() {
    let window = Interval::new(
        ModifiedJulianDate::try_new(Days::new(60000.0)).unwrap(),
        ModifiedJulianDate::try_new(Days::new(60002.0)).unwrap(),
    );
    let bands = altitude_ranges(
        &Sun,
        &greenwich(),
        window,
        Degrees::new(-18.0),
        Degrees::new(-12.0),
        SearchOpts::default(),
    );
    assert_periods_valid(&bands, window);
    assert!(bands.len() >= 2, "Should find ≥2 twilight bands in 2 days");
}

#[test]
fn moon_above_horizon() {
    let periods = above_threshold(
        &Moon,
        &roque(),
        one_week(),
        Degrees::new(0.0),
        SearchOpts::default(),
    );
    assert_periods_valid(&periods, one_week());
    assert!(!periods.is_empty());
}

#[test]
fn moon_below_horizon() {
    let periods = below_threshold(
        &Moon,
        &roque(),
        one_week(),
        Degrees::new(0.0),
        SearchOpts::default(),
    );
    assert_periods_valid(&periods, one_week());
    assert!(!periods.is_empty());
}

#[test]
fn star_sirius_above_horizon() {
    let periods = above_threshold(
        &SIRIUS,
        &greenwich(),
        one_day(),
        Degrees::new(0.0),
        SearchOpts::default(),
    );
    assert_periods_valid(&periods, one_day());
    assert!(!periods.is_empty());
}

#[test]
fn star_vega_above_horizon() {
    let periods = above_threshold(
        &VEGA,
        &greenwich(),
        one_day(),
        Degrees::new(0.0),
        SearchOpts::default(),
    );
    assert_periods_valid(&periods, one_day());
    assert!(!periods.is_empty());
}

#[test]
fn icrs_direction_above_horizon() {
    let dir = direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716));
    let periods = above_threshold(
        &dir,
        &greenwich(),
        one_day(),
        Degrees::new(0.0),
        SearchOpts::default(),
    );
    assert_periods_valid(&periods, one_day());
    assert!(!periods.is_empty());
}

#[test]
fn icrs_direction_matches_star() {
    let star = &SIRIUS;
    let dir = direction::ICRS::from(star);
    let window = one_day();
    let observer = greenwich();
    let opts = SearchOpts::default();

    let star_periods = above_threshold(star, &observer, window, Degrees::new(0.0), opts);
    let dir_periods = above_threshold(&dir, &observer, window, Degrees::new(0.0), opts);
    assert_eq!(star_periods.len(), dir_periods.len());
}

#[test]
fn altitude_at_sun_in_range() {
    let alt = Sun.altitude_at(&greenwich(), ModifiedJulianDate::new(51544.5));
    assert!(alt.abs().value() < std::f64::consts::FRAC_PI_2);
}

#[test]
fn empty_window_returns_empty() {
    let window = Interval::new(
        ModifiedJulianDate::try_new(Days::new(60000.0)).unwrap(),
        ModifiedJulianDate::try_new(Days::new(60000.0)).unwrap(),
    );
    assert!(above_threshold(
        &Sun,
        &greenwich(),
        window,
        Degrees::new(0.0),
        SearchOpts::default()
    )
    .is_empty());
}

#[test]
fn circumpolar_star_always_above() {
    let periods = above_threshold(
        &POLARIS,
        &greenwich(),
        one_day(),
        Degrees::new(0.0),
        SearchOpts::default(),
    );
    assert_eq!(periods.len(), 1);
}

#[test]
fn never_visible_star_at_north_pole() {
    let dir = direction::ICRS::new(Degrees::new(100.0), Degrees::new(-60.0));
    let periods = above_threshold(
        &dir,
        &north_pole(),
        one_day(),
        Degrees::new(0.0),
        SearchOpts::default(),
    );
    assert!(periods.is_empty());
}
