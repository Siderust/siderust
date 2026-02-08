// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Integration tests for the `calculus::stellar` module.
//!
//! Validates the analytical sinusoidal model against the generic scan engine
//! and checks correctness for circumpolar, rise/set, and never‑visible cases.

use qtty::*;
use siderust::calculus::stellar::{
    find_star_above_periods, find_star_above_periods_scan, find_star_below_periods,
    find_star_range_periods, find_star_range_periods_scan,
};
use siderust::coordinates::centers::ObserverSite;
use siderust::time::{ModifiedJulianDate, Period};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn greenwich() -> ObserverSite {
    ObserverSite::new(
        Degrees::new(0.0),
        Degrees::new(51.4769),
        Quantity::<Meter>::new(0.0),
    )
}

fn roque() -> ObserverSite {
    ObserverSite::new(
        Degrees::new(-17.892),
        Degrees::new(28.762),
        Quantity::<Meter>::new(2396.0),
    )
}

fn period_7d() -> Period<ModifiedJulianDate> {
    Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60007.0),
    )
}

fn period_3d() -> Period<ModifiedJulianDate> {
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

// ===========================================================================
// Circumpolar case: Polaris at mid-latitudes
// ===========================================================================

#[test]
fn polaris_circumpolar_at_greenwich() {
    let periods = find_star_above_periods(
        POLARIS_RA,
        POLARIS_DEC,
        greenwich(),
        period_7d(),
        Degrees::new(0.0),
    );

    assert_eq!(
        periods.len(),
        1,
        "Polaris should be continuously above horizon at 51°N"
    );
    let dur = periods[0].duration_days();
    assert!(
        (dur - 7.0).abs() < 0.01,
        "should span full 7 days, got {}",
        dur
    );
}

#[test]
fn polaris_circumpolar_at_roque() {
    // Roque (28.7°N) — Polaris has Dec ≈ +89.26°, so it's circumpolar
    let periods = find_star_above_periods(
        POLARIS_RA,
        POLARIS_DEC,
        roque(),
        period_7d(),
        Degrees::new(0.0),
    );

    assert_eq!(periods.len(), 1, "Polaris circumpolar at 28.7°N");
}

// ===========================================================================
// Rise / set case: Sirius at mid-latitudes
// ===========================================================================

#[test]
fn sirius_above_horizon_greenwich_7d() {
    let periods = find_star_above_periods(
        SIRIUS_RA,
        SIRIUS_DEC,
        greenwich(),
        period_7d(),
        Degrees::new(0.0),
    );

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
    let periods = find_star_above_periods(
        SIRIUS_RA,
        SIRIUS_DEC,
        roque(),
        period_7d(),
        Degrees::new(0.0),
    );

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
    let periods = find_star_above_periods(
        Degrees::new(0.0),
        Degrees::new(-80.0),
        greenwich(),
        period_7d(),
        Degrees::new(0.0),
    );
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

    let above = find_star_above_periods(SIRIUS_RA, SIRIUS_DEC, site, period, thr);
    let below = find_star_below_periods(SIRIUS_RA, SIRIUS_DEC, site, period, thr);

    let total_above: f64 = above.iter().map(|p| p.duration_days()).sum();
    let total_below: f64 = below.iter().map(|p| p.duration_days()).sum();
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
    let periods = find_star_range_periods(
        SIRIUS_RA,
        SIRIUS_DEC,
        roque(),
        period_7d(),
        (Degrees::new(10.0), Degrees::new(30.0)),
    );
    assert!(
        !periods.is_empty(),
        "should find altitude-range periods for Sirius at Roque"
    );
}

// ===========================================================================
// Analytical ↔ scan consistency
// ===========================================================================

#[test]
fn analytical_matches_scan_above_threshold() {
    let site = roque();
    let period = period_3d();
    let thr = Degrees::new(0.0);

    let analytical = find_star_above_periods(SIRIUS_RA, SIRIUS_DEC, site, period, thr);
    let scan = find_star_above_periods_scan(SIRIUS_RA, SIRIUS_DEC, site, period, thr);

    assert_eq!(
        analytical.len(),
        scan.len(),
        "analytical({}) and scan({}) should agree on period count",
        analytical.len(),
        scan.len(),
    );

    let tol = 1.0 / 1440.0; // 1 minute
    for (a, s) in analytical.iter().zip(scan.iter()) {
        assert!(
            (a.start.value() - s.start.value()).abs() < tol,
            "start times differ by {:.6} d",
            (a.start.value() - s.start.value()).abs()
        );
        assert!(
            (a.end.value() - s.end.value()).abs() < tol,
            "end times differ by {:.6} d",
            (a.end.value() - s.end.value()).abs()
        );
    }
}

#[test]
fn analytical_matches_scan_range() {
    let site = roque();
    let period = period_3d();
    let range = (Degrees::new(10.0), Degrees::new(40.0));

    let analytical = find_star_range_periods(SIRIUS_RA, SIRIUS_DEC, site, period, range);
    let scan = find_star_range_periods_scan(SIRIUS_RA, SIRIUS_DEC, site, period, range);

    assert_eq!(
        analytical.len(),
        scan.len(),
        "analytical({}) and scan({}) range periods should agree",
        analytical.len(),
        scan.len(),
    );

    let tol = 1.0 / 1440.0;
    for (a, s) in analytical.iter().zip(scan.iter()) {
        assert!(
            (a.start.value() - s.start.value()).abs() < tol,
            "range start times differ by {:.6} d",
            (a.start.value() - s.start.value()).abs()
        );
        assert!(
            (a.end.value() - s.end.value()).abs() < tol,
            "range end times differ by {:.6} d",
            (a.end.value() - s.end.value()).abs()
        );
    }
}

// ===========================================================================
// Consistency with the unified altitude API
// ===========================================================================

#[test]
fn stellar_matches_unified_api() {
    use siderust::calculus::altitude::AltitudePeriodsProvider;
    use siderust::coordinates::spherical::direction;

    let site = roque();
    let period = period_3d();
    let thr = Degrees::new(0.0);

    let stellar = find_star_above_periods(SIRIUS_RA, SIRIUS_DEC, site, period, thr);

    // Use direction::ICRS which implements AltitudePeriodsProvider
    // The trait's above_threshold method dispatches to the analytical stellar engine
    let sirius = direction::ICRS::new(SIRIUS_RA, SIRIUS_DEC);
    let unified = sirius.above_threshold(site, period, thr);

    // Since the trait impl dispatches to the stellar module, these should be identical.
    assert_eq!(stellar.len(), unified.len());
    for (s, u) in stellar.iter().zip(unified.iter()) {
        assert!(
            (s.start.value() - u.start.value()).abs() < 1e-12,
            "Start mismatch: stellar={}, unified={}",
            s.start.value(),
            u.start.value()
        );
        assert!(
            (s.end.value() - u.end.value()).abs() < 1e-12,
            "End mismatch: stellar={}, unified={}",
            s.end.value(),
            u.end.value()
        );
    }
}
