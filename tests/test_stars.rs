// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use qtty::length::nominal::SolarRadiuses;
use qtty::*;
use siderust::bodies::Star;
use siderust::coordinates::spherical::position::EquatorialMeanJ2000;
use siderust::targets::CoordinateWithPM;
use siderust::time::JulianDate;
use std::borrow::Cow;

#[test]
fn star_new_const() {
    let position =
        EquatorialMeanJ2000::<LightYear>::new(Degrees::new(10.0), Degrees::new(20.0), 4.0);
    let expected_ra = position.ra();
    let expected_dec = position.dec();
    let coord = CoordinateWithPM::new_static(position, JulianDate::J2000);
    let star = Star::new_const(
        "Alpha",
        LightYears::new(4.0),
        SolarMasses::new(1.1),
        SolarRadiuses::new(1.2),
        SolarLuminosities::new(1.3),
        coord.clone(),
    );
    assert_eq!(star.name, "Alpha");
    assert!(matches!(star.name, Cow::Borrowed("Alpha")));
    assert_eq!(star.distance, LightYears::new(4.0));
    assert_eq!(star.mass, SolarMasses::new(1.1));
    assert_eq!(star.radius, SolarRadiuses::new(1.2));
    assert_eq!(star.luminosity, SolarLuminosities::new(1.3));
    assert_eq!(star.coordinate.position.ra(), expected_ra);
    assert_eq!(star.coordinate.position.dec(), expected_dec);
    assert_eq!(star.coordinate.time, JulianDate::J2000);
}

#[test]
fn star_new_owned() {
    let position =
        EquatorialMeanJ2000::<LightYear>::new(Degrees::new(30.0), Degrees::new(-10.0), 10.0);
    let expected_ra = position.ra();
    let expected_dec = position.dec();
    let coord = CoordinateWithPM::new_static(position, JulianDate::J2000);
    let name = String::from("Beta");
    let star = Star::new(
        name.clone(),
        LightYears::new(10.0),
        SolarMasses::new(2.0),
        SolarRadiuses::new(3.0),
        SolarLuminosities::new(4.0),
        coord.clone(),
    );
    assert_eq!(star.name, "Beta");
    assert!(matches!(star.name, Cow::Owned(ref s) if s == &name));
    assert_eq!(star.distance, LightYears::new(10.0));
    assert_eq!(star.mass, SolarMasses::new(2.0));
    assert_eq!(star.radius, SolarRadiuses::new(3.0));
    assert_eq!(star.luminosity, SolarLuminosities::new(4.0));
    assert_eq!(star.coordinate.position.ra(), expected_ra);
    assert_eq!(star.coordinate.position.dec(), expected_dec);
    assert_eq!(star.coordinate.time, JulianDate::J2000);
}
