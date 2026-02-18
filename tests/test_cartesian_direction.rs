// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use qtty::*;
use siderust::coordinates::cartesian::Direction;
use siderust::coordinates::centers::Heliocentric;
use siderust::coordinates::frames::EclipticMeanJ2000;

#[test]
fn direction_position_scales_with_magnitude() {
    // Direction is now frame-only (no center parameter)
    let dir = Direction::<EclipticMeanJ2000>::normalize(1.0, 2.0, 2.0);
    let magnitude = 3.5 * AU;
    // position() now requires explicit center type
    let pos = dir.position::<Heliocentric, _>(magnitude);
    let expected_x = magnitude.value() / 3.0;
    let expected_y = 2.0 * magnitude.value() / 3.0;
    let expected_z = 2.0 * magnitude.value() / 3.0;
    assert!((pos.x().value() - expected_x).abs() < 1e-12);
    assert!((pos.y().value() - expected_y).abs() < 1e-12);
    assert!((pos.z().value() - expected_z).abs() < 1e-12);
}

#[test]
fn direction_normalize_returns_unit_vector() {
    let dir = Direction::<EclipticMeanJ2000>::normalize(3.0, 4.0, 0.0);
    let norm = (dir.x().powi(2) + dir.y().powi(2) + dir.z().powi(2)).sqrt();
    assert!((norm - 1.0).abs() < 1e-12);
    assert!((dir.x() - 0.6).abs() < 1e-12);
    assert!((dir.y() - 0.8).abs() < 1e-12);
    assert!(dir.z().abs() < 1e-12);
}

#[test]
fn direction_display_includes_frame_and_components() {
    let dir = Direction::<EclipticMeanJ2000>::normalize(1.0, 2.0, 2.0);
    let formatted = dir.display();
    // Direction no longer has center, only frame
    assert!(formatted.starts_with("Frame: EclipticMeanJ2000"));
    assert!(formatted.contains("X: 0.333333"));
    assert!(formatted.contains("Y: 0.666667"));
    assert!(formatted.contains("Z: 0.666667"));
}
