use siderust::coordinates::cartesian::Direction;
use siderust::coordinates::centers::Heliocentric;
use siderust::coordinates::frames::Ecliptic;
use siderust::units::*;

#[test]
fn direction_position_scales_with_magnitude() {
    let dir = Direction::<Heliocentric, Ecliptic>::normalize(1.0, 2.0, 2.0);
    let magnitude = 3.5 * AU;
    let pos = dir.position(magnitude);
    let expected_x = magnitude.value() / 3.0;
    let expected_y = 2.0 * magnitude.value() / 3.0;
    let expected_z = 2.0 * magnitude.value() / 3.0;
    assert!((pos.x().value() - expected_x).abs() < 1e-12);
    assert!((pos.y().value() - expected_y).abs() < 1e-12);
    assert!((pos.z().value() - expected_z).abs() < 1e-12);
}

#[test]
fn direction_normalize_returns_unit_vector() {
    let dir = Direction::<Heliocentric, Ecliptic>::normalize(3.0, 4.0, 0.0);
    let norm = (dir.x().value().powi(2) + dir.y().value().powi(2) + dir.z().value().powi(2)).sqrt();
    assert!((norm - 1.0).abs() < 1e-12);
    assert!((dir.x().value() - 0.6).abs() < 1e-12);
    assert!((dir.y().value() - 0.8).abs() < 1e-12);
    assert!(dir.z().value().abs() < 1e-12);
}

#[test]
fn direction_display_includes_center_frame_and_components() {
    let dir = Direction::<Heliocentric, Ecliptic>::normalize(1.0, 2.0, 2.0);
    let formatted = format!("{}", dir);
    assert!(formatted.starts_with("Center: Heliocentric, Frame: Ecliptic"));
    assert!(formatted.contains("X: 0.333333"));
    assert!(formatted.contains("Y: 0.666666"));
    assert!(formatted.contains("Z: 0.666666"));
}
