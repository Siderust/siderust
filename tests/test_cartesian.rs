// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use qtty::*;
use siderust::coordinates::cartesian::{Displacement, Position};
use siderust::coordinates::centers::Heliocentric;
use siderust::coordinates::frames::Ecliptic;

type DispAu = Displacement<Ecliptic, Au>;
type PosAu = Position<Heliocentric, Ecliptic, Au>;

#[test]
fn displacement_basic_operations() {
    let v1 = DispAu::new(1.0, 0.0, 0.0);
    let v2 = DispAu::new(0.0, 1.0, 0.0);

    let sum = v1 + v2;
    assert_eq!(sum.x(), 1.0);
    assert_eq!(sum.y(), 1.0);

    let diff = v1 - v2;
    assert_eq!(diff.x(), 1.0);
    assert_eq!(diff.y(), -1.0);

    let dist = v1.magnitude();
    assert!((dist.value() - 1.0).abs() < 1e-12);
}

#[test]
fn position_affine_operations() {
    let p1 = PosAu::new(1.0, 0.0, 0.0);
    let p2 = PosAu::new(0.0, 1.0, 0.0);

    // Position - Position = Displacement
    let disp = p1 - p2;
    assert_eq!(disp.x(), 1.0);
    assert_eq!(disp.y(), -1.0);

    let dist = p1.distance();
    assert!((dist.value() - 1.0).abs() < 1e-12);

    let between = p1.distance_to(&p2);
    assert!((between.value() - 2.0_f64.sqrt()).abs() < 1e-12);
}
