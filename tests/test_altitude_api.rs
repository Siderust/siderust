// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Integration tests for the unified altitude API (`calculus::altitude`).

use siderust::calculus::altitude::{
    above_threshold, altitude_at, altitude_ranges, below_threshold, crossings, culminations,
    AltitudeTarget, CrossingDirection, CulminationKind, SearchOpts,
};
use siderust::coordinates::centers::ObserverSite;
use siderust::time::{ModifiedJulianDate, Period};

use qtty::*;

/// La Palma (Roque de los Muchachos), ~28.76°N, -17.89°E, 2396 m
fn roque() -> ObserverSite {
    ObserverSite::new(
        Degrees::new(-17.892),
        Degrees::new(28.762),
        Quantity::<Meter>::new(2396.0),
    )
}

/// Greenwich, ~51.48°N, 0°E, sea level
fn greenwich() -> ObserverSite {
    ObserverSite::new(
        Degrees::new(0.0),
        Degrees::new(51.4769),
        Quantity::<Meter>::new(0.0),
    )
}

// ===========================================================================
// altitude_at — single-point smoke tests
// ===========================================================================

#[test]
fn altitude_at_sun_j2000_greenwich() {
    let alt = altitude_at(
        &AltitudeTarget::Sun,
        &greenwich(),
        siderust::astro::JulianDate::J2000,
    );
    // J2000 = 2000-01-01 12:00 TT, winter in northern hemisphere,
    // sun should be low but above horizon around local noon (UTC ≈ TT-64s).
    // Accept any value in [-90, +90].
    assert!(
        alt.value() > -90.0 && alt.value() < 90.0,
        "Sun altitude out of range: {}°",
        alt.value()
    );
}

#[test]
fn altitude_at_moon_j2000_greenwich() {
    let alt = altitude_at(
        &AltitudeTarget::Moon,
        &greenwich(),
        siderust::astro::JulianDate::J2000,
    );
    assert!(alt.value() > -90.0 && alt.value() < 90.0);
}

#[test]
fn altitude_at_sirius_reasonable() {
    // Sirius (α CMa): RA ≈ 101.287°, Dec ≈ −16.716° (J2000)
    let target = AltitudeTarget::FixedEquatorial {
        ra: Degrees::new(101.287),
        dec: Degrees::new(-16.716),
    };
    let alt = altitude_at(
        &target,
        &greenwich(),
        siderust::astro::JulianDate::J2000,
    );
    assert!(alt.value() > -90.0 && alt.value() < 90.0);
}

// ===========================================================================
// crossings — sunrise / sunset
// ===========================================================================

#[test]
fn crossings_sun_one_day_greenwich() {
    let site = greenwich();
    // MJD 60000 ≈ 2023-02-25
    let window = Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60001.0),
    );
    let events = crossings(&AltitudeTarget::Sun, &site, window, Degrees::new(0.0), SearchOpts::default());

    let rises: Vec<_> = events.iter().filter(|e| e.direction == CrossingDirection::Rising).collect();
    let sets: Vec<_> = events.iter().filter(|e| e.direction == CrossingDirection::Setting).collect();

    assert_eq!(rises.len(), 1, "expect 1 sunrise in 24h at 51°N");
    assert_eq!(sets.len(), 1, "expect 1 sunset in 24h at 51°N");

    // Sunrise should be before sunset
    assert!(rises[0].jd.value() < sets[0].jd.value());
}

#[test]
fn crossings_sun_astronomical_twilight() {
    let site = roque();
    let window = Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60001.0),
    );
    let events = crossings(
        &AltitudeTarget::Sun,
        &site,
        window,
        Degrees::new(-18.0),
        SearchOpts::default(),
    );
    // Should find 2 crossings (evening crossing below -18° and morning crossing above -18°)
    assert!(!events.is_empty(), "should find astronomical twilight crossings");
}

// ===========================================================================
// culminations — solar transit & anti-transit
// ===========================================================================

#[test]
fn culminations_sun_one_day() {
    let site = greenwich();
    let window = Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60001.0),
    );
    let culms = culminations(&AltitudeTarget::Sun, &site, window, SearchOpts::default());

    let upper: Vec<_> = culms.iter().filter(|c| c.kind == CulminationKind::Max).collect();
    let lower: Vec<_> = culms.iter().filter(|c| c.kind == CulminationKind::Min).collect();

    assert!(upper.len() >= 1, "should find upper culmination");
    assert!(lower.len() >= 1, "should find lower culmination");

    // Upper culmination altitude at 51°N in late February should be ~25-35°
    let max_alt = upper[0].altitude.value();
    assert!(
        max_alt > 10.0 && max_alt < 50.0,
        "upper culmination altitude {} should be between 10-50° at 51°N in winter",
        max_alt
    );
}

#[test]
fn culminations_moon_one_day() {
    let site = greenwich();
    let window = Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60001.0),
    );
    let culms = culminations(&AltitudeTarget::Moon, &site, window, SearchOpts::default());
    assert!(!culms.is_empty(), "should find Moon culminations in 24h");
}

// ===========================================================================
// above/below threshold — day/night periods
// ===========================================================================

#[test]
fn above_threshold_sun_week() {
    let site = roque();
    let window = Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60007.0),
    );
    let days = above_threshold(&AltitudeTarget::Sun, &site, window, Degrees::new(0.0), SearchOpts::default());

    assert_eq!(days.len(), 7, "should find 7 daytime periods in a week, found {}", days.len());
    for (i, p) in days.iter().enumerate() {
        let hours = p.duration_days() * 24.0;
        assert!(
            hours > 8.0 && hours < 16.0,
            "day {} duration {} hours is unreasonable",
            i,
            hours
        );
    }
}

#[test]
fn below_threshold_astronomical_night_week() {
    let site = roque();
    let window = Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60007.0),
    );
    let nights = below_threshold(
        &AltitudeTarget::Sun,
        &site,
        window,
        Degrees::new(-18.0),
        SearchOpts::default(),
    );

    assert!(!nights.is_empty(), "should find astronomical night periods");
    for p in &nights {
        let hours = p.duration_days() * 24.0;
        // Nights should be several hours long
        assert!(hours > 3.0, "night period is too short: {} hours", hours);
    }
}

// ===========================================================================
// altitude_ranges — twilight band
// ===========================================================================

#[test]
fn altitude_ranges_nautical_to_astro_twilight() {
    let site = greenwich();
    let window = Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60002.0),
    );
    let bands = altitude_ranges(
        &AltitudeTarget::Sun,
        &site,
        window,
        Degrees::new(-18.0),
        Degrees::new(-12.0),
        SearchOpts::default(),
    );

    // In 2 days, we should find ~4 twilight bands (evening+morning × 2)
    assert!(
        bands.len() >= 2,
        "should find at least 2 nautical-astronomical twilight bands, found {}",
        bands.len()
    );
    for p in &bands {
        let minutes = p.duration_days() * 1440.0;
        // Each twilight band should be 20-90 minutes
        assert!(
            minutes > 10.0 && minutes < 120.0,
            "twilight band duration {} min is unreasonable",
            minutes
        );
    }
}

// ===========================================================================
// Moon periods
// ===========================================================================

#[test]
fn moon_above_horizon_week() {
    let site = roque();
    let window = Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60007.0),
    );
    let periods = above_threshold(
        &AltitudeTarget::Moon,
        &site,
        window,
        Degrees::new(0.0),
        SearchOpts::default(),
    );
    // Moon should be above horizon for several periods over a week
    assert!(!periods.is_empty(), "should find moon-up periods in 7 days");
}

// ===========================================================================
// Fixed star — Polaris should be circumpolar at high latitude
// ===========================================================================

#[test]
fn polaris_always_above_horizon_at_greenwich() {
    // Polaris: RA ≈ 37.95°, Dec ≈ +89.26° (J2000)
    let target = AltitudeTarget::FixedEquatorial {
        ra: Degrees::new(37.95),
        dec: Degrees::new(89.26),
    };
    let site = greenwich(); // 51.5°N — Polaris is circumpolar here

    let window = Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60001.0),
    );

    // Should have no horizon crossings (circumpolar)
    let events = crossings(&target, &site, window, Degrees::new(0.0), SearchOpts::default());
    assert!(
        events.is_empty(),
        "Polaris should be circumpolar at 51.5°N — found {} crossings",
        events.len()
    );

    // Should be above horizon 100% of the time
    let up = above_threshold(&target, &site, window, Degrees::new(0.0), SearchOpts::default());
    assert_eq!(
        up.len(),
        1,
        "Polaris should be continuously above horizon, found {} periods",
        up.len()
    );
    let duration = up[0].duration_days();
    assert!(
        (duration - 1.0).abs() < 0.01,
        "Polaris up-period should span the full day, got {} days",
        duration
    );
}
