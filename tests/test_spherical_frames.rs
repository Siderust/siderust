// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use affn::geodesy::GeodeticCoord;
use qtty::*;
use siderust::coordinates::centers::ObserverSite;
use siderust::coordinates::spherical::{direction, position};

const EPS: f64 = 1e-6;

#[test]
fn ecef_normalization_and_altitude() {
    // Input angles outside their canonical ranges
    // new(lon, lat) normalizes lat via wrap_quarter_fold to [-90, 90]
    // and lon via normalize to [0, 360)
    // Note: IAU convention - lon first, lat second
    let dir = direction::EcefDir::new(190.0 * DEG, 95.0 * DEG);
    // 95° lat after wrap_quarter_fold: clamped to 85° (90 - 5)
    // 190° lon after normalize stays as 190°
    assert!(
        (dir.lat().value() - 85.0).abs() < EPS,
        "lat mismatch: {}",
        dir.lat().value()
    );
    assert!(
        (dir.lon().value() - 190.0).abs() < EPS,
        "lon mismatch: {}",
        dir.lon().value()
    );

    // GeodeticCoord stores geodetic (lon, lat, height) without normalisation —
    // validation happens in ObserverSite::from_geodetic.
    let coord = GeodeticCoord::new(190.0 * DEG, 85.0 * DEG, 10_000.0 * M);
    assert!(
        (coord.lat.value() - 85.0).abs() < EPS,
        "geodetic lat mismatch: {}",
        coord.lat.value()
    );
    assert!(
        (coord.lon.value() - 190.0).abs() < EPS,
        "geodetic lon mismatch: {}",
        coord.lon.value()
    );
}

#[test]
fn ecliptic_normalization() {
    // Direction::new(lon, lat) normalizes both
    let dir = direction::EclipticMeanJ2000::new(120.0 * DEG, -45.0 * DEG);
    // After normalization: lat = -45° (in [-90, 90]), lon = 120° (in [0, 360))
    assert!(
        (dir.lon().value() - 120.0).abs() < EPS,
        "lon mismatch: {}",
        dir.lon().value()
    );
    assert!(
        (dir.lat().value() - (-45.0)).abs() < EPS,
        "lat mismatch: {}",
        dir.lat().value()
    );

    // Position::new(lon, lat, distance) also normalizes
    let pos =
        position::EclipticMeanJ2000::<AstronomicalUnit>::new(120.0 * DEG, -45.0 * DEG, 2.0 * AU);
    assert!(
        (pos.lon().value() - 120.0).abs() < EPS,
        "pos lon mismatch: {}",
        pos.lon().value()
    );
    assert!(
        (pos.lat().value() - (-45.0)).abs() < EPS,
        "pos lat mismatch: {}",
        pos.lat().value()
    );
    assert!((pos.distance - 2.0 * AU).abs() < EPS * AU);
}

#[test]
fn horizontal_normalization() {
    // Direction is now frame-only (no site parameter)
    // Note: new(alt, az) - IAU Alt-Az convention (altitude first)
    let dir = direction::Horizontal::new(120.0 * DEG, -30.0 * DEG);
    // Altitude 120° wraps via wrap_quarter_fold: 90 - |120-90| = 60°, azimuth -30° normalizes to 330°
    assert!((dir.alt().value() - 60.0).abs() < EPS, "alt={}", dir.alt());
    assert!((dir.az().value() - 330.0).abs() < EPS, "az={}", dir.az());

    // Positions use new_with_site for Topocentric center
    let site = ObserverSite::default();
    // Note: new_with_site takes (site, alt, az, dist) - IAU Alt-Az convention (altitude first)
    let pos = siderust::coordinates::centers::Topocentric::horizontal(
        site,
        120.0 * DEG, // alt - wraps to 60
        -30.0 * DEG, // az - normalizes to 330
        2.0 * AU,
    );
    assert!(
        (pos.alt().value() - 60.0).abs() < EPS,
        "pos alt={}",
        pos.alt()
    );
    assert!(
        (pos.az().value() - 330.0).abs() < EPS,
        "pos az={}",
        pos.az()
    );
    assert!((pos.distance - 2.0 * AU).abs() < EPS * AU);
}
