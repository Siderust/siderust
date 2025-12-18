use qtty::*;
use siderust::bodies::EARTH;
use siderust::coordinates::centers::ObserverSite;
use siderust::coordinates::spherical::{direction, position};

const EPS: f64 = 1e-6;

#[test]
fn ecef_normalization_and_altitude() {
    // Input angles outside their canonical ranges
    // new_geographic(lon, lat) normalizes lat via wrap_quarter_fold to [-90, 90]
    // and lon via wrap_signed_lo to [-180, 180]
    let dir = direction::Geographic::new_geographic(95.0 * DEG, 190.0 * DEG);
    // 190° lat after wrap_quarter_fold: 90 - |280 - 180| = 90 - 100 = -10°
    // 95° lon after wrap_signed_lo stays as 95° (within [-180, 180])
    assert!((dir.lat().value() - (-10.0)).abs() < EPS, "lat mismatch: {}", dir.lat().value());
    assert!((dir.lon().value() - 95.0).abs() < EPS, "lon mismatch: {}", dir.lon().value());

    // Position::new(lon, lat, alt) uses the same normalization
    let pos = position::Geographic::new(95.0 * DEG, 190.0 * DEG, 10.0 * KM);
    assert!((pos.lat().value() - (-10.0)).abs() < EPS, "pos lat mismatch: {}", pos.lat().value());
    assert!((pos.lon().value() - 95.0).abs() < EPS, "pos lon mismatch: {}", pos.lon().value());
    assert!((pos.distance - (EARTH.radius + 10.0 * KM)).abs() < EPS * KM);
}

#[test]
fn ecliptic_normalization() {
    // Direction::new_ecliptic(lon, lat) normalizes both
    let dir = direction::Ecliptic::new_ecliptic(120.0 * DEG, -45.0 * DEG);
    // After normalization: lat = -45° (in [-90, 90]), lon = 120° (in [0, 360))
    assert!((dir.lon().value() - 120.0).abs() < EPS, "lon mismatch: {}", dir.lon().value());
    assert!((dir.lat().value() - (-45.0)).abs() < EPS, "lat mismatch: {}", dir.lat().value());

    // Position::new(lon, lat, distance) also normalizes
    let pos = position::Ecliptic::<AstronomicalUnit>::new(120.0 * DEG, -45.0 * DEG, 2.0 * AU);
    assert!((pos.lon().value() - 120.0).abs() < EPS, "pos lon mismatch: {}", pos.lon().value());
    assert!((pos.lat().value() - (-45.0)).abs() < EPS, "pos lat mismatch: {}", pos.lat().value());
    assert!((pos.distance - 2.0 * AU).abs() < EPS * AU);
}

#[test]
fn horizontal_normalization() {
    // Direction is now frame-only (no site parameter)
    // Note: new_horizontal takes (alt, az) order
    let dir = direction::Horizontal::new_horizontal(120.0 * DEG, -30.0 * DEG);
    // Altitude 120° wraps to 60° (quarter fold), azimuth -30° normalizes to 330°
    assert!((dir.alt().value() - 60.0).abs() < EPS);
    assert!((dir.az().value() - 330.0).abs() < EPS);

    // Positions still support observer site via with_site
    let site = ObserverSite::default();
    let pos = position::Horizontal::<AstronomicalUnit>::with_site(
        site,
        120.0 * DEG,
        -30.0 * DEG,
        2.0 * AU,
    );
    assert!((pos.alt().value() - 60.0).abs() < EPS);
    assert!((pos.az().value() - 330.0).abs() < EPS);
    assert!((pos.distance - 2.0 * AU).abs() < EPS * AU);
}
