use siderust::bodies::{EARTH, MARS, MOON};
use siderust::bodies::planets::OrbitExt;
use siderust::units::*;

#[test]
fn earth_constants() {
    assert!((EARTH.mass - Kilograms::new(5.97237e24)).abs() < Kilograms::new(1e16));
    assert!((EARTH.radius - Kilometers::new(6371.0)).abs() < Kilometers::new(1e-6));
    let orbit = &EARTH.orbit;
    assert!((orbit.semi_major_axis - AstronomicalUnits::new(1.00000011)).abs() < AstronomicalUnits::new(1e-8));
    assert!((orbit.eccentricity - 0.01671022).abs() < 1e-8);
    assert!((orbit.inclination - Degrees::new(0.00005)).abs().value() < 1e-8);
}

#[test]
fn moon_constants() {
    assert!((MOON.mass - Kilograms::new(7.346e22)).abs() < Kilograms::new(1e16));
    assert!((MOON.radius - Kilometers::new(1737.4)).abs() < Kilometers::new(1e-4));
    let orbit = &MOON.orbit;
    assert!((orbit.semi_major_axis - AstronomicalUnits::new(2.566881e-6)).abs() < AstronomicalUnits::new(1e-12));
    assert!((orbit.eccentricity - 0.0549).abs() < 1e-6);
}

#[test]
fn orbital_period() {
    let earth_period = EARTH.orbit.period();
    // sidereal year ~ 365.25 days
    let days = earth_period.to::<Day>().value();
    assert!((days - 365.25).abs() < 0.1);

    let mars_period = MARS.orbit.period();
    let days = mars_period.to::<Day>().value();
    assert!((days - 686.98).abs() < 0.5);
}
