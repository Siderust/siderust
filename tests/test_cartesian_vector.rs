use siderust::coordinates::cartesian::{Direction, Position};
use siderust::coordinates::centers::Heliocentric;
use siderust::coordinates::frames::Ecliptic;
use siderust::units::*;
use siderust::coordinates::*;
use siderust::coordinates::centers::*;
use siderust::coordinates::frames::*;

fn approx_eq<C, F, U>(
    a: cartesian::Vector<C, F, U>,
    b: cartesian::Vector<C, F, U>,
    tol: f64,
)
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
    Quantity<U>: std::cmp::PartialOrd + std::fmt::Display,
{
    let eps: Quantity<U> = tol.into();
    assert!((a.x() - b.x()).abs() < eps, "x mismatch: {} vs {}", a.x(), b.x());
    assert!((a.y() - b.y()).abs() < eps, "y mismatch: {} vs {}", a.y(), b.y());
    assert!((a.z() - b.z()).abs() < eps, "z mismatch: {} vs {}", a.z(), b.z());
}


#[test]
fn vector_add_sub_distance() {
    let a = Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(1.0*AU, 2.0*AU, 2.0*AU);
    let b = Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(0.5*AU, -1.0*AU, 1.0*AU);
    let sum = a + b;
    approx_eq(sum, Position::new(1.5*AU, 1.0*AU, 3.0*AU), 1e-12);
    let diff = sum.sub(&a);
    approx_eq(diff, b, 1e-12);
    let dist = a.distance();
    assert!((dist.value() - (1.0_f64*1.0 + 2.0*2.0 + 2.0*2.0).sqrt()).abs() < 1e-12);
    let dist_to = a.distance_to(&b);
    assert!((dist_to.value() - ((0.5_f64).powi(2) + (3.0_f64).powi(2) + (1.0_f64).powi(2)).sqrt()).abs() < 1e-12);
}

#[test]
fn direction_normalize_position() {
    let p = Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(3.0*AU, 0.0*AU, 4.0*AU);
    let dir = p.direction();
    let scaled = dir.position(2.5*AU);
    approx_eq(scaled, Position::new(1.5*AU, 0.0*AU, 2.0*AU), 1e-12);
    let norm = Direction::<Heliocentric, Ecliptic>::normalize(2.0, 0.0, 0.0);
    let pos = norm.position(1.0*AU);
    approx_eq(pos, Position::new(1.0*AU, 0.0*AU, 0.0*AU), 1e-12);
}
