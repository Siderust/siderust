// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Integration tests for the stellar altitude engine via the unified
//! [`AltitudePeriodsProvider`] trait API.
//!
//! Validates correctness for circumpolar, rise/set, and never-visible cases.
//! Scan-vs-analytical consistency is tested inside the `calculus::stellar`
//! module's own `#[cfg(test)]` block.

use qtty::*;
use siderust::calculus::altitude::{AltitudePeriodsProvider, AltitudeQuery};
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::coordinates::spherical::direction;
use siderust::time::{ModifiedJulianDate, Period, MJD};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn greenwich() -> Geodetic<ECEF> {
    Geodetic::<ECEF>::new(
        Degrees::new(0.0),
        Degrees::new(51.4769),
        Quantity::<Meter>::new(0.0),
    )
}

fn roque() -> Geodetic<ECEF> {
    Geodetic::<ECEF>::new(
        Degrees::new(-17.892),
        Degrees::new(28.762),
        Quantity::<Meter>::new(2396.0),
    )
}

fn period_7d() -> Period<MJD> {
    Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60007.0),
    )
}

fn period_3d() -> Period<MJD> {
    Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60003.0),
    )
}

// Sirius — J2000 RA/Dec
const SIRIUS_RA: Degrees = Degrees::new(101.287);
const SIRIUS_DEC: Degrees = Degrees::new(-16.716);

// Polaris — J2000 RA/Dec
const POLARIS_RA: Degrees = Degrees::new(37.95);
const POLARIS_DEC: Degrees = Degrees::new(89.26);

fn sirius() -> direction::ICRS {
    direction::ICRS::new(SIRIUS_RA, SIRIUS_DEC)
}

fn polaris() -> direction::ICRS {
    direction::ICRS::new(POLARIS_RA, POLARIS_DEC)
}

// ===========================================================================
// Circumpolar case: Polaris at mid-latitudes
// ===========================================================================

#[test]
fn polaris_circumpolar_at_greenwich() {
    let periods = polaris().above_threshold(greenwich(), period_7d(), Degrees::new(0.0));

    assert_eq!(
        periods.len(),
        1,
        "Polaris should be continuously above horizon at 51°N"
    );
    let dur = periods[0].duration_days();
    assert!(
        (dur - Days::new(7.0)).abs() < 0.01,
        "should span full 7 days, got {}",
        dur
    );
}

#[test]
fn polaris_circumpolar_at_roque() {
    // Roque (28.7°N) — Polaris has Dec ≈ +89.26°, so it's circumpolar
    let periods = polaris().above_threshold(roque(), period_7d(), Degrees::new(0.0));

    assert_eq!(periods.len(), 1, "Polaris circumpolar at 28.7°N");
}

// ===========================================================================
// Rise / set case: Sirius at mid-latitudes
// ===========================================================================

#[test]
fn sirius_above_horizon_greenwich_7d() {
    let periods = sirius().above_threshold(greenwich(), period_7d(), Degrees::new(0.0));

    assert!(
        periods.len() >= 6 && periods.len() <= 8,
        "expected ~7 above-horizon periods for Sirius at 51°N, got {}",
        periods.len()
    );
    for p in &periods {
        let hours = p.duration_days() * 24.0;
        // First/last period may be truncated by the window boundary
        assert!(
            hours > 0.1 && hours < 18.0,
            "unreasonable above-horizon duration: {} h",
            hours
        );
    }
}

#[test]
fn sirius_above_horizon_roque_7d() {
    let periods = sirius().above_threshold(roque(), period_7d(), Degrees::new(0.0));

    assert!(
        periods.len() >= 6 && periods.len() <= 8,
        "expected ~7 above-horizon periods for Sirius at 28°N, got {}",
        periods.len()
    );
}

// ===========================================================================
// Never-visible case
// ===========================================================================

#[test]
fn deep_south_star_never_visible_at_greenwich() {
    // Dec = −80° at 51°N → max altitude ≈ −20° → never above horizon
    let star = direction::ICRS::new(Degrees::new(0.0), Degrees::new(-80.0));
    let periods = star.above_threshold(greenwich(), period_7d(), Degrees::new(0.0));
    assert!(
        periods.is_empty(),
        "Star at Dec=−80° should never be visible at 51°N"
    );
}

// ===========================================================================
// Above + below = full period
// ===========================================================================

#[test]
fn above_plus_below_covers_full_period() {
    let site = greenwich();
    let period = period_7d();
    let thr = Degrees::new(0.0);

    let above = sirius().above_threshold(site, period, thr);
    let below = sirius().below_threshold(site, period, thr);

    let total_above: f64 = above.iter().map(|p| p.duration_days().value()).sum();
    let total_below: f64 = below.iter().map(|p| p.duration_days().value()).sum();
    assert!(
        (total_above + total_below - 7.0).abs() < 0.01,
        "above({}) + below({}) should sum to 7 days",
        total_above,
        total_below,
    );
}

// ===========================================================================
// Range periods
// ===========================================================================

#[test]
fn range_periods_sirius_roque() {
    let periods = sirius().altitude_periods(&AltitudeQuery {
        observer: roque(),
        window: period_7d(),
        min_altitude: Degrees::new(10.0),
        max_altitude: Degrees::new(30.0),
    });
    assert!(
        !periods.is_empty(),
        "should find altitude-range periods for Sirius at Roque"
    );
}

// ===========================================================================
// Consistency: trait API dispatches to the same engine
// ===========================================================================

#[test]
fn trait_api_above_below_consistent() {
    let site = roque();
    let period = period_3d();
    let thr = Degrees::new(5.0);

    let above = sirius().above_threshold(site, period, thr);
    let below = sirius().below_threshold(site, period, thr);

    let total_above: f64 = above.iter().map(|p| p.duration_days().value()).sum();
    let total_below: f64 = below.iter().map(|p| p.duration_days().value()).sum();
    assert!(
        (total_above + total_below - 3.0).abs() < 0.01,
        "above({:.4}) + below({:.4}) should sum to 3.0 days",
        total_above,
        total_below,
    );
}

#[test]
fn trait_api_range_within_above() {
    let site = roque();
    let period = period_3d();

    let above_10 = sirius().above_threshold(site, period, Degrees::new(10.0));
    let range_10_30 = sirius().altitude_periods(&AltitudeQuery {
        observer: site,
        window: period,
        min_altitude: Degrees::new(10.0),
        max_altitude: Degrees::new(30.0),
    });

    // Range [10°, 30°] periods should be subsets of above(10°) periods
    let total_range: f64 = range_10_30.iter().map(|p| p.duration_days().value()).sum();
    let total_above: f64 = above_10.iter().map(|p| p.duration_days().value()).sum();
    assert!(
        total_range <= total_above + 0.01,
        "range time ({:.4}) should not exceed above time ({:.4})",
        total_range,
        total_above,
    );
}
