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


use siderust::units::{Meter, Meters, Km, AstronomicalUnits, LightYears, Seconds, Unitless, Quantity, Simplify};

#[test]
fn meters_kilometers_conversion() {
    let m = Meters::new(5000.0);
    let km = m.to::<Km>();
    assert!((km.value() - 5.0).abs() < 1e-12);

    let m_back = km.to::<Meter>();
    assert!((m_back.value() - 5000.0).abs() < 1e-12);
}

#[test]
fn astronomical_unit_lightyear_roundtrip() {
    let au = AstronomicalUnits::new(2.0);
    let ly: LightYears = au.into();
    let au_back: AstronomicalUnits = ly.into();
    assert!((au_back.value() - au.value()).abs() < 1e-10);
}

#[test]
fn quantity_arithmetic_and_simplify() {
    let distance = Meters::new(120.0);
    let extra = Meters::new(5.0);
    let sum = distance + extra;
    assert_eq!(sum, 125.0);

    let diff = sum - extra;
    assert!((diff.value() - distance.value()).abs() < 1e-12);

    let time = Seconds::new(10.0);
    let velocity = distance / time; // Quantity<Per<Meter, Second>>
    let dist_back = velocity * time;
    assert!((dist_back.value() - distance.value()).abs() < 1e-12);

    let time_back = (distance / velocity).simplify();
    assert!((time_back.value() - time.value()).abs() < 1e-12);

    let unitless_ratio = distance / distance;
    assert!((unitless_ratio.asin() - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
}

#[test]
fn unitless_from_length_and_display() {
    let meters = Meters::new(42.0);
    let unitless: Quantity<Unitless> = meters.into();
    assert_eq!(unitless.value(), 42.0);
    assert_eq!(format!("{}", unitless), "42");
}
