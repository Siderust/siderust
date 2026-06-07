// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Integration tests for stellar altitude periods via the standardized API.

use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::coordinates::spherical::direction;
use siderust::event::altitude::{above_threshold, altitude_ranges, below_threshold, SearchOpts};
use siderust::qtty::*;
use siderust::time::{Interval, ModifiedJulianDate};

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

fn period_7d() -> Interval<ModifiedJulianDate> {
    Interval::new(
        ModifiedJulianDate::try_new(Days::new(60000.0)).unwrap(),
        ModifiedJulianDate::try_new(Days::new(60007.0)).unwrap(),
    )
}

fn period_3d() -> Interval<ModifiedJulianDate> {
    Interval::new(
        ModifiedJulianDate::try_new(Days::new(60000.0)).unwrap(),
        ModifiedJulianDate::try_new(Days::new(60003.0)).unwrap(),
    )
}

fn sirius() -> direction::ICRS {
    direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716))
}

fn polaris() -> direction::ICRS {
    direction::ICRS::new(Degrees::new(37.95), Degrees::new(89.26))
}

#[test]
fn polaris_circumpolar_at_greenwich() {
    let periods = above_threshold(
        &polaris(),
        &greenwich(),
        period_7d(),
        Degrees::new(0.0),
        SearchOpts::default(),
    );
    assert_eq!(periods.len(), 1);
}

#[test]
fn sirius_above_horizon_greenwich_7d() {
    let periods = above_threshold(
        &sirius(),
        &greenwich(),
        period_7d(),
        Degrees::new(0.0),
        SearchOpts::default(),
    );
    assert!(periods.len() >= 6 && periods.len() <= 8);
}

#[test]
fn deep_south_star_never_visible_at_greenwich() {
    let star = direction::ICRS::new(Degrees::new(0.0), Degrees::new(-80.0));
    let periods = above_threshold(
        &star,
        &greenwich(),
        period_7d(),
        Degrees::new(0.0),
        SearchOpts::default(),
    );
    assert!(periods.is_empty());
}

#[test]
fn range_periods_sirius_roque() {
    let periods = altitude_ranges(
        &sirius(),
        &roque(),
        period_7d(),
        Degrees::new(10.0),
        Degrees::new(30.0),
        SearchOpts::default(),
    );
    assert!(!periods.is_empty());
}

#[test]
fn above_plus_below_covers_full_period() {
    let site = roque();
    let period = period_3d();
    let opts = SearchOpts::default();
    let above = above_threshold(&sirius(), &site, period, Degrees::new(5.0), opts);
    let below = below_threshold(&sirius(), &site, period, Degrees::new(5.0), opts);
    let total: f64 = above
        .iter()
        .chain(below.iter())
        .map(|p| (p.end.raw() - p.start.raw()).value())
        .sum();
    assert!((total - 3.0).abs() < 0.01);
}
