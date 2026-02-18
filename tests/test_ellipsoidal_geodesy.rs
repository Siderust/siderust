// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Ellipsoidal geodesy tests
//!
//! Tests for the `Ellipsoid` trait, `HasEllipsoid` frame association,
//! the generic `to_cartesian` / `from_cartesian` conversions, and round-trip
//! accuracy across multiple geodetic positions.

use affn::ellipsoid::{Ellipsoid, Grs80, HasEllipsoid, Wgs84};
use siderust::coordinates::centers::Geodetic;
use qtty::*;
use siderust::coordinates::cartesian::Position;
use siderust::coordinates::centers::Geocentric;
use siderust::coordinates::frames::{ECEF, ITRF};

// =============================================================================
// Ellipsoid trait tests
// =============================================================================

#[test]
fn wgs84_ellipsoid_constants() {
    assert!((Wgs84::A - 6_378_137.0).abs() < 1e-6);
    assert!((Wgs84::F - 1.0 / 298.257_223_563).abs() < 1e-15);
    assert!((Wgs84::e2() - 0.006_694_379_990_14).abs() < 1e-12);
    assert!((Wgs84::b() - 6_356_752.314_245_179).abs() < 1e-3);
}

#[test]
fn grs80_ellipsoid_constants() {
    assert!((Grs80::A - 6_378_137.0).abs() < 1e-6);
    assert!((Grs80::F - 1.0 / 298.257_222_101).abs() < 1e-15);
    // GRS80 and WGS84 share the same a, differ slightly in f
    assert!((Grs80::e2() - Wgs84::e2()).abs() < 1e-10);
}

// =============================================================================
// HasEllipsoid derive tests
// =============================================================================

#[test]
fn ecef_has_wgs84_ellipsoid() {
    // ECEF frame should have WGS84 associated via #[frame(ellipsoid = "Wgs84")]
    assert_eq!(<ECEF as HasEllipsoid>::Ellipsoid::A, Wgs84::A);
    assert_eq!(<ECEF as HasEllipsoid>::Ellipsoid::F, Wgs84::F);
}

#[test]
fn itrf_has_grs80_ellipsoid() {
    // ITRF frame should have GRS80 associated via #[frame(ellipsoid = "Grs80")]
    assert_eq!(<ITRF as HasEllipsoid>::Ellipsoid::A, Grs80::A);
    assert_eq!(<ITRF as HasEllipsoid>::Ellipsoid::F, Grs80::F);
}

// =============================================================================
// Generic to_ecef tests
// =============================================================================

/// Absolute tolerance: 1 metre in ECEF XYZ.
const TOL_M: f64 = 1.0;

#[test]
fn to_ecef_greenwich() {
    let coord = Geodetic::<ECEF>::new(0.0 * DEG, 51.4769 * DEG, 65.0 * M);
    let pos: Position<Geocentric, ECEF, Meter> = coord.to_cartesian();
    assert!((pos.x().value() - 3_980_700.035).abs() < TOL_M, "X={}", pos.x().value());
    assert!(pos.y().value().abs() < TOL_M, "Y={}", pos.y().value());
    assert!((pos.z().value() - 4_966_813.011).abs() < TOL_M, "Z={}", pos.z().value());
}

#[test]
fn to_ecef_equator_prime_meridian() {
    let coord = Geodetic::<ECEF>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M);
    let pos: Position<Geocentric, ECEF, Meter> = coord.to_cartesian();
    assert!((pos.x().value() - Wgs84::A).abs() < TOL_M);
    assert!(pos.y().value().abs() < TOL_M);
    assert!(pos.z().value().abs() < TOL_M);
}

#[test]
fn to_ecef_north_pole() {
    let coord = Geodetic::<ECEF>::new(0.0 * DEG, 90.0 * DEG, 0.0 * M);
    let pos: Position<Geocentric, ECEF, Meter> = coord.to_cartesian();
    assert!(pos.x().value().abs() < TOL_M);
    assert!(pos.y().value().abs() < TOL_M);
    assert!((pos.z().value() - Wgs84::b()).abs() < TOL_M);
}

#[test]
fn to_ecef_with_kilometer_output() {
    let coord = Geodetic::<ECEF>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M);
    let pos: Position<Geocentric, ECEF, Kilometer> = coord.to_cartesian();
    assert!((pos.x().value() - Wgs84::A / 1000.0).abs() < 0.001);
}

// =============================================================================
// Generic height unit tests
// =============================================================================

#[test]
fn geodetic_coord_generic_kilometer_height() {
    let coord = Geodetic::<ECEF, Kilometer>::new(0.0 * DEG, 0.0 * DEG, Kilometers::new(10.0));
    assert_eq!(coord.height.value(), 10.0);
}

#[test]
fn to_ecef_from_kilometer_height() {
    // Geodetic::<ECEF> with height in km should produce the same ECEF as one in metres
    let coord_m = Geodetic::<ECEF>::new(0.0 * DEG, 45.0 * DEG, Meters::new(1000.0));
    let coord_km =
        Geodetic::<ECEF, Kilometer>::new(0.0 * DEG, 45.0 * DEG, Kilometers::new(1.0));

    let pos_m: Position<Geocentric, ECEF, Meter> = coord_m.to_cartesian();
    let pos_km: Position<Geocentric, ECEF, Meter> = coord_km.to_cartesian();

    assert!((pos_m.x().value() - pos_km.x().value()).abs() < 0.01);
    assert!((pos_m.y().value() - pos_km.y().value()).abs() < 0.01);
    assert!((pos_m.z().value() - pos_km.z().value()).abs() < 0.01);
}

// =============================================================================
// from_ecef (inverse) tests
// =============================================================================

/// Tolerance for angular round-trip: 0.1 µas ≈ sub-mm on the surface.
const ANG_TOL_DEG: f64 = 1e-10;
/// Tolerance for height round-trip: sub-millimetre.
const H_TOL_M: f64 = 0.001;

fn check_round_trip(lon_deg: f64, lat_deg: f64, h_m: f64) {
    let original = Geodetic::<ECEF>::new(
        Degrees::new(lon_deg),
        Degrees::new(lat_deg),
        Meters::new(h_m),
    );
    let ecef: Position<Geocentric, ECEF, Meter> = original.to_cartesian();
    let back = Geodetic::<ECEF>::from_cartesian::<Meter>(&ecef);

    assert!(
        (back.lon.value() - original.lon.value()).abs() < ANG_TOL_DEG,
        "lon: expected {} got {} (diff {})",
        original.lon.value(),
        back.lon.value(),
        (back.lon.value() - original.lon.value()).abs(),
    );
    assert!(
        (back.lat.value() - original.lat.value()).abs() < ANG_TOL_DEG,
        "lat: expected {} got {} (diff {})",
        original.lat.value(),
        back.lat.value(),
        (back.lat.value() - original.lat.value()).abs(),
    );
    assert!(
        (back.height.value() - original.height.value()).abs() < H_TOL_M,
        "height: expected {} got {} (diff {})",
        original.height.value(),
        back.height.value(),
        (back.height.value() - original.height.value()).abs(),
    );
}

#[test]
fn round_trip_greenwich() {
    check_round_trip(0.0, 51.4769, 65.0);
}

#[test]
fn round_trip_roque_de_los_muchachos() {
    check_round_trip(-17.8925, 28.7543, 2396.0);
}

#[test]
fn round_trip_mauna_kea() {
    check_round_trip(-155.4681, 19.8207, 4205.0);
}

#[test]
fn round_trip_equator() {
    check_round_trip(0.0, 0.0, 0.0);
}

#[test]
fn round_trip_north_pole() {
    check_round_trip(0.0, 90.0, 0.0);
}

#[test]
fn round_trip_south_pole() {
    check_round_trip(0.0, -90.0, 0.0);
}

#[test]
fn round_trip_high_altitude() {
    // ISS orbit altitude
    check_round_trip(45.0, 30.0, 408_000.0);
}

#[test]
fn round_trip_negative_lon() {
    check_round_trip(-120.5, -35.2, 500.0);
}

#[test]
fn round_trip_date_line() {
    check_round_trip(179.9999, 0.0, 0.0);
}

// =============================================================================
// ITRF frame with GRS80 ellipsoid
// =============================================================================

#[test]
fn to_ecef_itrf_uses_grs80() {
    // ITRF has GRS80 ellipsoid — at equator, X should equal GRS80's a
    let coord = Geodetic::<ITRF>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M);
    let pos: Position<Geocentric, ITRF, Meter> = coord.to_cartesian();
    assert!((pos.x().value() - Grs80::A).abs() < TOL_M);
    assert!(pos.y().value().abs() < TOL_M);
    assert!(pos.z().value().abs() < TOL_M);
}

#[test]
fn itrf_vs_ecef_equator_differ_due_to_ellipsoid() {
    // The difference between WGS84 and GRS80 is tiny but nonzero
    let coord_ecef = Geodetic::<ECEF>::new(0.0 * DEG, 45.0 * DEG, 1000.0 * M);
    let coord_itrf = Geodetic::<ITRF>::new(0.0 * DEG, 45.0 * DEG, 1000.0 * M);

    let pos_ecef: Position<Geocentric, ECEF, Meter> = coord_ecef.to_cartesian();
    let pos_itrf: Position<Geocentric, ITRF, Meter> = coord_itrf.to_cartesian();

    // At 45° latitude, the difference should be in the sub-metre range
    let diff_x = (pos_ecef.x().value() - pos_itrf.x().value()).abs();
    let diff_z = (pos_ecef.z().value() - pos_itrf.z().value()).abs();

    // WGS84 and GRS80 differ in the 9th decimal of f, so total diff is tiny
    assert!(diff_x < 0.01, "x diff = {diff_x}");
    assert!(diff_z < 0.01, "z diff = {diff_z}");
    // But they should not be exactly equal (they use different ellipsoids)
    let total_diff = diff_x + diff_z;
    assert!(total_diff > 1e-10, "total diff should be >0, got {total_diff}");
}
