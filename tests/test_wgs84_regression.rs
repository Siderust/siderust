// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # WGS84 geodetic → ECEF regression tests
//!
//! These tests verify that [`ObserverSite::geocentric_itrf`] produces correct
//! Earth-Centred Earth-Fixed (ECEF) Cartesian coordinates for well-known
//! observatory locations.
//!
//! ## Reference values
//!
//! Ground-truth ECEF coordinates were computed from the WGS84 definition:
//!
//! - Semi-major axis:  a = 6 378 137.0 m
//! - Flattening:       f = 1 / 298.257 223 563
//! - e² = 2f − f²     ≈ 0.006 694 379 990 14
//!
//! using the standard geodetic-to-ECEF formulae (Bowring 1985 / IERS Conventions 2010 §4):
//!
//! ```text
//! N(φ) = a / √(1 − e² sin²φ)
//! X = (N + h) cos φ cos λ
//! Y = (N + h) cos φ sin λ
//! Z = (N(1 − e²) + h) sin φ
//! ```
//!
//! ## Tolerance
//!
//! All assertions use a 1 m tolerance, which is well within the accuracy of
//! the WGS84 formula and far tighter than any pointing/observing need.

use affn::geodesy::GeodeticCoord;
use qtty::*;
use siderust::coordinates::centers::ObserverSite;
use siderust::coordinates::cartesian::Position;
use siderust::coordinates::centers::Geocentric;
use siderust::coordinates::frames::ECEF;

/// Absolute tolerance: 1 metre in ECEF XYZ.
const TOL_M: f64 = 1.0;

fn check_ecef(site: ObserverSite, expected_x: f64, expected_y: f64, expected_z: f64, label: &str) {
    let pos: Position<Geocentric, ECEF, Meter> = site.geocentric_itrf();
    let (x, y, z) = (pos.x().value(), pos.y().value(), pos.z().value());
    assert!(
        (x - expected_x).abs() < TOL_M,
        "{label}: X expected {expected_x:.3} m, got {x:.3} m (diff {:.3} m)",
        (x - expected_x).abs()
    );
    assert!(
        (y - expected_y).abs() < TOL_M,
        "{label}: Y expected {expected_y:.3} m, got {y:.3} m (diff {:.3} m)",
        (y - expected_y).abs()
    );
    assert!(
        (z - expected_z).abs() < TOL_M,
        "{label}: Z expected {expected_z:.3} m, got {z:.3} m (diff {:.3} m)",
        (z - expected_z).abs()
    );
}

/// Royal Observatory Greenwich (prime meridian, ~51.5°N, 65 m elevation).
///
/// Expected ECEF (WGS84):
///   X =  3 980 700.035 m
///   Y =          0.000 m  (on prime meridian → Y ≈ 0)
///   Z =  4 966 813.011 m
#[test]
fn greenwich_ecef() {
    let coord = GeodeticCoord::new(0.0 * DEG, 51.4769 * DEG, 65.0 * M);
    let site = ObserverSite::from_geodetic(&coord);
    check_ecef(site, 3_980_700.035, 0.0, 4_966_813.011, "Greenwich");
}

/// Roque de los Muchachos Observatory, La Palma, Spain.
///
/// Expected ECEF (WGS84):
///   X =  5 327 336.209 m
///   Y = -1 719 912.648 m
///   Z =  3 051 208.060 m
#[test]
fn roque_de_los_muchachos_ecef() {
    let coord = GeodeticCoord::new(-17.8925 * DEG, 28.7543 * DEG, 2396.0 * M);
    let site = ObserverSite::from_geodetic(&coord);
    check_ecef(site, 5_327_336.209, -1_719_912.648, 3_051_208.060, "Roque");
}

/// W. M. Keck Observatory (Mauna Kea, Hawaiʻi).
///
/// Expected ECEF (WGS84):
///   X = -5 464 341.898 m
///   Y = -2 493 919.182 m
///   Z =  2 150 459.985 m
#[test]
fn mauna_kea_ecef() {
    let coord = GeodeticCoord::new(-155.4681 * DEG, 19.8207 * DEG, 4205.0 * M);
    let site = ObserverSite::from_geodetic(&coord);
    check_ecef(site, -5_464_341.898, -2_493_919.182, 2_150_459.985, "MaunaKea");
}

/// Sanity check: a point at the equator on the prime meridian with zero height
/// should have X = a (semi-major axis), Y = 0, Z = 0.
#[test]
fn equator_prime_meridian_ecef() {
    let coord = GeodeticCoord::new(0.0 * DEG, 0.0 * DEG, 0.0 * M);
    let site = ObserverSite::from_geodetic(&coord);
    let pos: Position<Geocentric, ECEF, Meter> = site.geocentric_itrf();
    let a = 6_378_137.0_f64; // WGS84 semi-major axis
    assert!(
        (pos.x().value() - a).abs() < TOL_M,
        "X should equal semi-major axis a={a}m, got {}m",
        pos.x().value()
    );
    assert!(
        pos.y().value().abs() < TOL_M,
        "Y should be ≈0 on prime meridian, got {}m",
        pos.y().value()
    );
    assert!(
        pos.z().value().abs() < TOL_M,
        "Z should be ≈0 on equator, got {}m",
        pos.z().value()
    );
}

/// Sanity check: north pole (lat=90°, h=0) should have X=Y=0, Z≈b (semi-minor axis).
#[test]
fn north_pole_ecef() {
    let coord = GeodeticCoord::new(0.0 * DEG, 90.0 * DEG, 0.0 * M);
    let site = ObserverSite::from_geodetic(&coord);
    let pos: Position<Geocentric, ECEF, Meter> = site.geocentric_itrf();
    // b = a*(1-f) = 6 356 752.314 m
    let b = 6_356_752.314_f64;
    assert!(
        pos.x().value().abs() < TOL_M,
        "X should be ≈0 at north pole, got {}m",
        pos.x().value()
    );
    assert!(
        pos.y().value().abs() < TOL_M,
        "Y should be ≈0 at north pole, got {}m",
        pos.y().value()
    );
    assert!(
        (pos.z().value() - b).abs() < TOL_M,
        "Z should equal semi-minor axis b={b}m, got {}m",
        pos.z().value()
    );
}
