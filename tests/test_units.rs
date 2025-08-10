use siderust::units::angular::{HourAngles, Degrees};

#[test]
fn hour_angles_from_hms() {
    let ha = HourAngles::from_hms(5, 30, 0.0);
    assert!((ha.value() - 5.5).abs() < 1e-12);

    let neg = HourAngles::from_hms(-5, 15, 30.0);
    assert!((neg.value() + 5.258333333333333).abs() < 1e-12);
}

#[test]
fn degrees_from_dms() {
    let deg = Degrees::from_dms(-33, 52, 0.0);
    assert!((deg.value() + 33.86666666666667).abs() < 1e-12);

    let with_sign = Degrees::from_dms_sign(-1, 33, 52, 0.0);
    assert!((with_sign.value() + 33.86666666666667).abs() < 1e-12);
}
