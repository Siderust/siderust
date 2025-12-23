use qtty::*;
use siderust::astro::orbit::Orbit;
use siderust::astro::JulianDate;
use siderust::bodies::Satellite;
use std::borrow::Cow;

#[test]
fn satellite_new_const() {
    let orbit = Orbit::new(
        AstronomicalUnits::new(1.0),
        0.01,
        Degrees::new(5.0),
        Degrees::new(10.0),
        Degrees::new(20.0),
        Degrees::new(30.0),
        JulianDate::J2000,
    );
    let sat = Satellite::new_const(
        "TestSat",
        Kilograms::new(1000.0),
        Kilometers::new(10.0),
        orbit,
    );
    assert_eq!(sat.name, "TestSat");
    assert!(matches!(sat.name, Cow::Borrowed("TestSat")));
    assert_eq!(sat.mass, Kilograms::new(1000.0));
    assert_eq!(sat.radius, Kilometers::new(10.0));
    assert_eq!(sat.orbit.semi_major_axis, AstronomicalUnits::new(1.0));
    assert!((sat.orbit.eccentricity - 0.01).abs() < 1e-10);
}

#[test]
fn satellite_new_owned() {
    let orbit = Orbit::new(
        AstronomicalUnits::new(2.0),
        0.1,
        Degrees::new(1.0),
        Degrees::new(2.0),
        Degrees::new(3.0),
        Degrees::new(4.0),
        JulianDate::J2000,
    );
    let name = String::from("OwnedSat");
    let sat = Satellite::new(
        name.clone(),
        Kilograms::new(2000.0),
        Kilometers::new(20.0),
        orbit,
    );
    assert_eq!(sat.name, "OwnedSat");
    assert!(matches!(sat.name, Cow::Owned(ref s) if s == &name));
    assert_eq!(sat.mass, Kilograms::new(2000.0));
    assert_eq!(sat.radius, Kilometers::new(20.0));
    assert_eq!(sat.orbit.semi_major_axis, AstronomicalUnits::new(2.0));
    assert!((sat.orbit.eccentricity - 0.1).abs() < 1e-10);
}
