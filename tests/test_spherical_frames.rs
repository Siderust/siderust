use qtty::*;
use siderust::bodies::EARTH;
use siderust::coordinates::spherical::{direction, position};

const EPS: f64 = 1e-6;

#[test]
fn ecef_normalization_and_altitude() {
    // Input angles outside their canonical ranges
    let dir = direction::Geographic::new(190.0 * DEG, 95.0 * DEG);
    assert!((dir.lat().value() - (-170.0)).abs() < EPS);
    assert!((dir.lon().value() - 85.0).abs() < EPS);

    let pos = position::Geographic::new(190.0 * DEG, 95.0 * DEG, 10.0 * KM);
    assert!((pos.lat().value() - (-170.0)).abs() < EPS);
    assert!((pos.lon().value() - 85.0).abs() < EPS);
    assert!((pos.distance - (EARTH.radius + 10.0 * KM)).abs() < EPS * KM);
}

#[test]
fn ecliptic_normalization() {
    let dir = direction::Ecliptic::new(-45.0 * DEG, 120.0 * DEG);
    assert!((dir.lon().value() - 315.0).abs() < EPS);
    assert!((dir.lat().value() - 60.0).abs() < EPS);

    let pos = position::Ecliptic::<AstronomicalUnit>::new(-45.0 * DEG, 120.0 * DEG, 2.0 * AU);
    assert!((pos.lon().value() - 315.0).abs() < EPS);
    assert!((pos.lat().value() - 60.0).abs() < EPS);
    assert!((pos.distance - 2.0 * AU).abs() < EPS * AU);
}

#[test]
fn horizontal_normalization() {
    let dir = direction::Horizontal::new(120.0 * DEG, -30.0 * DEG);
    assert!((dir.alt().value() - 60.0).abs() < EPS);
    assert!((dir.az().value() - 330.0).abs() < EPS);

    let pos = position::Horizontal::<AstronomicalUnit>::new(120.0 * DEG, -30.0 * DEG, 2.0 * AU);
    assert!((pos.alt().value() - 60.0).abs() < EPS);
    assert!((pos.az().value() - 330.0).abs() < EPS);
    assert!((pos.distance - 2.0 * AU).abs() < EPS * AU);
}
