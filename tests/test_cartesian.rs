use siderust::coordinates::cartesian::Vector;
use siderust::coordinates::centers::Heliocentric;
use siderust::coordinates::frames::Ecliptic;
use qtty::*;

type VecAu = Vector<Heliocentric, Ecliptic, Au>;

#[test]
fn vector_basic_operations() {
    let v1 = VecAu::new(1.0, 0.0, 0.0);
    let v2 = VecAu::new(0.0, 1.0, 0.0);

    let sum = v1 + v2;
    assert_eq!(sum.x().value(), 1.0);
    assert_eq!(sum.y().value(), 1.0);

    let diff = v1 - v2;
    assert_eq!(diff.x().value(), 1.0);
    assert_eq!(diff.y().value(), -1.0);

    let dist = v1.distance();
    assert!((dist.value() - 1.0).abs() < 1e-12);

    let between = v1.distance_to(&v2);
    assert!((between.value() - 2.0_f64.sqrt()).abs() < 1e-12);
}
