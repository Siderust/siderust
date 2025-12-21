use qtty::*;
use siderust::bodies::EARTH;
use siderust::coordinates::spherical::{
    EcefDirectionExt, EcefPositionExt, EclipticDirectionExt, EclipticPositionExt,
    HorizontalDirectionExt, HorizontalPositionExt, HorizontalPositionReadExt,
};
use siderust::coordinates::centers::ObserverSite;
use siderust::coordinates::spherical::{direction, position};

const EPS: f64 = 1e-6;

#[test]
fn ecef_normalization_and_altitude() {
    // Input angles outside their canonical ranges
    // new_geographic(lat, lon) normalizes lat via wrap_quarter_fold to [-90, 90]
    // and lon via normalize to [0, 360)
    let dir = direction::Geographic::new_geographic(95.0 * DEG, 190.0 * DEG);
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

    // Position::new_geographic(lat, lon, alt) uses the same normalization
    let pos = position::Geographic::new_geographic(95.0 * DEG, 190.0 * DEG, 10.0 * KM);
    assert!(
        (pos.lat().value() - 85.0).abs() < EPS,
        "pos lat mismatch: {}",
        pos.lat().value()
    );
    assert!(
        (pos.lon().value() - 190.0).abs() < EPS,
        "pos lon mismatch: {}",
        pos.lon().value()
    );
    // Note: distance is now just the third coordinate, not radius + altitude
    assert!((pos.distance - 10.0 * KM).abs() < EPS * KM);
}

#[test]
fn ecliptic_normalization() {
    // Direction::new_ecliptic(lon, lat) normalizes both
    let dir = direction::Ecliptic::new_ecliptic(120.0 * DEG, -45.0 * DEG);
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

    // Position::new_ecliptic(lon, lat, distance) also normalizes
    let pos = position::Ecliptic::<AstronomicalUnit>::new_ecliptic(120.0 * DEG, -45.0 * DEG, 2.0 * AU);
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
    // Note: new_horizontal takes (alt, az) order
    let dir = direction::Horizontal::new_horizontal(120.0 * DEG, -30.0 * DEG);
    // Altitude 120° wraps via wrap_quarter_fold: 90 - |120-90| = 60°, azimuth -30° normalizes to 330°
    assert!((dir.alt().value() - 60.0).abs() < EPS, "alt={}", dir.alt());
    assert!((dir.az().value() - 330.0).abs() < EPS, "az={}", dir.az());

    // Positions use new_raw_with_params for Topocentric center
    let site = ObserverSite::default();
    // Note: new_raw takes (polar, azimuth, distance) = (alt, az, dist)
    // wrap_quarter_fold(120°) = 60°, not 90°
    let pos = position::Horizontal::<AstronomicalUnit>::new_raw_with_params(
        site,
        (120.0 * DEG).wrap_quarter_fold(),  // alt (polar) - wraps to 60
        (-30.0 * DEG).normalize(),          // az (azimuth) - normalizes to 330
        2.0 * AU,
    );
    assert!((pos.alt().value() - 60.0).abs() < EPS, "pos alt={}", pos.alt());
    assert!((pos.az().value() - 330.0).abs() < EPS, "pos az={}", pos.az());
    assert!((pos.distance - 2.0 * AU).abs() < EPS * AU);
}
