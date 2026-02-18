// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use qtty::*;
use siderust::coordinates::cartesian::{Direction, Displacement, Position};
use siderust::coordinates::centers::Heliocentric;
use siderust::coordinates::frames::EclipticMeanJ2000;

fn approx_eq_pos<C, F, U>(a: Position<C, F, U>, b: Position<C, F, U>, tol: f64)
where
    C: siderust::coordinates::centers::ReferenceCenter,
    F: siderust::coordinates::frames::ReferenceFrame,
    U: LengthUnit,
    Quantity<U>: std::cmp::PartialOrd + std::fmt::Display,
{
    let eps: Quantity<U> = tol.into();
    assert!(
        (a.x() - b.x()).abs() < eps,
        "x mismatch: {} vs {}",
        a.x(),
        b.x()
    );
    assert!(
        (a.y() - b.y()).abs() < eps,
        "y mismatch: {} vs {}",
        a.y(),
        b.y()
    );
    assert!(
        (a.z() - b.z()).abs() < eps,
        "z mismatch: {} vs {}",
        a.z(),
        b.z()
    );
}

fn approx_eq_disp<F, U>(a: Displacement<F, U>, b: Displacement<F, U>, tol: f64)
where
    F: siderust::coordinates::frames::ReferenceFrame,
    U: LengthUnit,
    Quantity<U>: std::cmp::PartialOrd + std::fmt::Display,
{
    let eps: Quantity<U> = tol.into();
    assert!(
        (a.x() - b.x()).abs() < eps,
        "x mismatch: {} vs {}",
        a.x(),
        b.x()
    );
    assert!(
        (a.y() - b.y()).abs() < eps,
        "y mismatch: {} vs {}",
        a.y(),
        b.y()
    );
    assert!(
        (a.z() - b.z()).abs() < eps,
        "z mismatch: {} vs {}",
        a.z(),
        b.z()
    );
}

#[test]
fn position_sub_distance() {
    let a = Position::<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>::new(
        1.0 * AU,
        2.0 * AU,
        2.0 * AU,
    );
    let b = Position::<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>::new(
        0.5 * AU,
        -1.0 * AU,
        1.0 * AU,
    );

    // Position - Position = Displacement
    let diff = a - b;
    approx_eq_disp(diff, Displacement::new(0.5 * AU, 3.0 * AU, 1.0 * AU), 1e-12);

    let dist = a.distance();
    assert!((dist.value() - (1.0_f64 * 1.0 + 2.0 * 2.0 + 2.0 * 2.0).sqrt()).abs() < 1e-12);
    let dist_to = a.distance_to(&b);
    assert!(
        (dist_to.value() - ((0.5_f64).powi(2) + (3.0_f64).powi(2) + (1.0_f64).powi(2)).sqrt())
            .abs()
            < 1e-12
    );
}

#[test]
fn displacement_add_sub() {
    let v1 = Displacement::<EclipticMeanJ2000, AstronomicalUnit>::new(1.0 * AU, 2.0 * AU, 2.0 * AU);
    let v2 =
        Displacement::<EclipticMeanJ2000, AstronomicalUnit>::new(0.5 * AU, -1.0 * AU, 1.0 * AU);

    let sum = v1 + v2;
    approx_eq_disp(sum, Displacement::new(1.5 * AU, 1.0 * AU, 3.0 * AU), 1e-12);

    let diff = v1 - v2;
    approx_eq_disp(diff, Displacement::new(0.5 * AU, 3.0 * AU, 1.0 * AU), 1e-12);
}

#[test]
fn direction_normalize_position() {
    let p = Position::<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>::new(
        3.0 * AU,
        0.0 * AU,
        4.0 * AU,
    );
    let dir = p.direction().expect("position should have a direction");
    // direction() now returns Direction<F> (frame-only), position() requires explicit center
    let scaled = dir.position::<Heliocentric, _>(2.5 * AU);
    approx_eq_pos(scaled, Position::new(1.5 * AU, 0.0 * AU, 2.0 * AU), 1e-12);
    // Direction is now frame-only (no center parameter)
    let norm = Direction::<EclipticMeanJ2000>::normalize(2.0, 0.0, 0.0);
    let pos = norm.position::<Heliocentric, _>(1.0 * AU);
    approx_eq_pos(pos, Position::new(1.0 * AU, 0.0 * AU, 0.0 * AU), 1e-12);
}
