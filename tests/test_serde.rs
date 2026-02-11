// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Serde roundtrip tests for coordinate types.
//!
//! These tests verify that all coordinate types can be serialized to JSON
//! and deserialized back without data loss.

#![cfg(feature = "serde")]

use qtty::*;
use serde_json;
use siderust::coordinates::{cartesian, centers, frames, spherical};

// =============================================================================
// Cartesian Direction Tests
// =============================================================================

#[test]
fn test_cartesian_direction_roundtrip() {
    let dir = cartesian::Direction::<frames::ICRS>::new(1.0, 2.0, 2.0);

    let json = serde_json::to_string(&dir).expect("serialize direction");
    let recovered: cartesian::Direction<frames::ICRS> =
        serde_json::from_str(&json).expect("deserialize direction");

    assert!((dir.x() - recovered.x()).abs() < 1e-12);
    assert!((dir.y() - recovered.y()).abs() < 1e-12);
    assert!((dir.z() - recovered.z()).abs() < 1e-12);
}

// =============================================================================
// Cartesian Position Tests
// =============================================================================

#[test]
fn test_cartesian_position_geocentric_roundtrip() {
    let pos = cartesian::Position::<centers::Geocentric, frames::ICRS, AstronomicalUnit>::new(
        1.0, 0.5, -0.25,
    );

    let json = serde_json::to_string(&pos).expect("serialize position");
    let recovered: cartesian::Position<centers::Geocentric, frames::ICRS, AstronomicalUnit> =
        serde_json::from_str(&json).expect("deserialize position");

    assert!((pos.x().value() - recovered.x().value()).abs() < 1e-12);
    assert!((pos.y().value() - recovered.y().value()).abs() < 1e-12);
    assert!((pos.z().value() - recovered.z().value()).abs() < 1e-12);
}

#[test]
fn test_cartesian_position_heliocentric_roundtrip() {
    let pos = cartesian::Position::<centers::Heliocentric, frames::Ecliptic, Kilometer>::new(
        1.496e8, 0.0, 0.0,
    );

    let json = serde_json::to_string(&pos).expect("serialize position");
    let recovered: cartesian::Position<centers::Heliocentric, frames::Ecliptic, Kilometer> =
        serde_json::from_str(&json).expect("deserialize position");

    assert!((pos.x().value() - recovered.x().value()).abs() < 1e-6);
}

// =============================================================================
// Cartesian Vector/Displacement Tests
// =============================================================================

#[test]
fn test_cartesian_displacement_roundtrip() {
    let disp = cartesian::Displacement::<frames::ICRS, Meter>::new(100.0, 200.0, 300.0);

    let json = serde_json::to_string(&disp).expect("serialize displacement");
    let recovered: cartesian::Displacement<frames::ICRS, Meter> =
        serde_json::from_str(&json).expect("deserialize displacement");

    assert!((disp.x().value() - recovered.x().value()).abs() < 1e-12);
    assert!((disp.y().value() - recovered.y().value()).abs() < 1e-12);
    assert!((disp.z().value() - recovered.z().value()).abs() < 1e-12);
}

#[test]
fn test_cartesian_velocity_roundtrip() {
    type KmPerSec = Per<Kilometer, Second>;
    let vel = cartesian::Velocity::<frames::ICRS, KmPerSec>::new(10.0, 20.0, 30.0);

    let json = serde_json::to_string(&vel).expect("serialize velocity");
    let recovered: cartesian::Velocity<frames::ICRS, KmPerSec> =
        serde_json::from_str(&json).expect("deserialize velocity");

    assert!((vel.x().value() - recovered.x().value()).abs() < 1e-12);
    assert!((vel.y().value() - recovered.y().value()).abs() < 1e-12);
    assert!((vel.z().value() - recovered.z().value()).abs() < 1e-12);
}

// =============================================================================
// Spherical Direction Tests
// =============================================================================

#[test]
fn test_spherical_direction_icrs_roundtrip() {
    let dir = spherical::Direction::<frames::ICRS>::new(120.0 * DEG, 45.0 * DEG);

    let json = serde_json::to_string(&dir).expect("serialize spherical direction");
    let recovered: spherical::Direction<frames::ICRS> =
        serde_json::from_str(&json).expect("deserialize spherical direction");

    assert!((dir.ra().value() - recovered.ra().value()).abs() < 1e-12);
    assert!((dir.dec().value() - recovered.dec().value()).abs() < 1e-12);
}

#[test]
fn test_spherical_direction_ecliptic_roundtrip() {
    let dir = spherical::Direction::<frames::Ecliptic>::new(45.0 * DEG, 7.0 * DEG);

    let json = serde_json::to_string(&dir).expect("serialize spherical direction");
    let recovered: spherical::Direction<frames::Ecliptic> =
        serde_json::from_str(&json).expect("deserialize spherical direction");

    assert!((dir.polar.value() - recovered.polar.value()).abs() < 1e-12);
    assert!((dir.azimuth.value() - recovered.azimuth.value()).abs() < 1e-12);
}

// =============================================================================
// Spherical Position Tests
// =============================================================================

#[test]
fn test_spherical_position_heliocentric_roundtrip() {
    let pos =
        spherical::Position::<centers::Heliocentric, frames::Ecliptic, AstronomicalUnit>::new_raw(
            7.0 * DEG,   // latitude
            120.0 * DEG, // longitude
            1.5 * AU,    // distance
        );

    let json = serde_json::to_string(&pos).expect("serialize spherical position");
    let recovered: spherical::Position<centers::Heliocentric, frames::Ecliptic, AstronomicalUnit> =
        serde_json::from_str(&json).expect("deserialize spherical position");

    assert!((pos.polar.value() - recovered.polar.value()).abs() < 1e-12);
    assert!((pos.azimuth.value() - recovered.azimuth.value()).abs() < 1e-12);
    assert!((pos.distance.value() - recovered.distance.value()).abs() < 1e-12);
}

// =============================================================================
// Center Types Tests
// =============================================================================

#[test]
fn test_observer_site_roundtrip() {
    let site = centers::ObserverSite::new(
        -17.8947 * DEG, // longitude (La Palma)
        28.7636 * DEG,  // latitude
        2396.0 * M,     // height
    );

    let json = serde_json::to_string(&site).expect("serialize observer site");
    let recovered: centers::ObserverSite =
        serde_json::from_str(&json).expect("deserialize observer site");

    assert!((site.lon.value() - recovered.lon.value()).abs() < 1e-12);
    assert!((site.lat.value() - recovered.lat.value()).abs() < 1e-12);
    assert!((site.height.value() - recovered.height.value()).abs() < 1e-12);
}

#[test]
fn test_orbit_reference_center_roundtrip() {
    let center = centers::OrbitReferenceCenter::Heliocentric;

    let json = serde_json::to_string(&center).expect("serialize orbit reference center");
    let recovered: centers::OrbitReferenceCenter =
        serde_json::from_str(&json).expect("deserialize orbit reference center");

    assert_eq!(center, recovered);
}

// =============================================================================
// Frame Types Tests
// =============================================================================

#[test]
fn test_frame_types_roundtrip() {
    // Test that frame marker types serialize as unit structs
    let icrs = frames::ICRS;
    let ecliptic = frames::Ecliptic;
    let horizontal = frames::Horizontal;

    let json_icrs = serde_json::to_string(&icrs).expect("serialize ICRS");
    let json_ecliptic = serde_json::to_string(&ecliptic).expect("serialize Ecliptic");
    let json_horizontal = serde_json::to_string(&horizontal).expect("serialize Horizontal");

    assert_eq!(json_icrs, "null"); // unit structs serialize as null
    assert_eq!(json_ecliptic, "null");
    assert_eq!(json_horizontal, "null");

    let _: frames::ICRS = serde_json::from_str(&json_icrs).expect("deserialize ICRS");
    let _: frames::Ecliptic = serde_json::from_str(&json_ecliptic).expect("deserialize Ecliptic");
    let _: frames::Horizontal =
        serde_json::from_str(&json_horizontal).expect("deserialize Horizontal");
}

// =============================================================================
// Time Types Tests
// =============================================================================

#[test]
fn test_julian_date_roundtrip() {
    use siderust::time::JulianDate;

    let jd = JulianDate::new(2451545.0); // J2000

    let json = serde_json::to_string(&jd).expect("serialize julian date");
    let recovered: JulianDate = serde_json::from_str(&json).expect("deserialize julian date");

    assert!((jd.value() - recovered.value()).abs() < 1e-12);
}

#[test]
fn test_modified_julian_date_roundtrip() {
    use siderust::time::ModifiedJulianDate;

    let mjd = ModifiedJulianDate::new(51544.5); // J2000 in MJD

    let json = serde_json::to_string(&mjd).expect("serialize modified julian date");
    let recovered: ModifiedJulianDate =
        serde_json::from_str(&json).expect("deserialize modified julian date");

    assert!((mjd.value() - recovered.value()).abs() < 1e-12);
}
