// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use qtty::*;
use siderust::astro::orbit::Orbit;
use siderust::astro::JulianDate;
use siderust::bodies::asteroid::{Asteroid, AsteroidClass};
use siderust::bodies::comet::{Comet, OrbitFrame, HALLEY};
use siderust::bodies::planets::{OrbitExt, Planet, PlanetBuilderError};
use siderust::bodies::{EARTH, MARS, MOON};

#[test]
fn earth_constants() {
    assert!((EARTH.mass - Kilograms::new(5.97237e24)).abs() < Kilograms::new(1e16));
    assert!((EARTH.radius - Kilometers::new(6371.0)).abs() < Kilometers::new(1e-6));
    let orbit = &EARTH.orbit;
    assert!(
        (orbit.semi_major_axis - AstronomicalUnits::new(1.00000011)).abs()
            < AstronomicalUnits::new(1e-8)
    );
    assert!((orbit.eccentricity - 0.01671022).abs() < 1e-8);
    assert!((orbit.inclination - Degrees::new(0.00005)).abs().value() < 1e-8);
}

#[test]
fn moon_constants() {
    assert!((MOON.mass - Kilograms::new(7.346e22)).abs() < Kilograms::new(1e16));
    assert!((MOON.radius - Kilometers::new(1737.4)).abs() < Kilometers::new(1e-4));
    let orbit = &MOON.orbit;
    assert!(
        (orbit.semi_major_axis - AstronomicalUnits::new(2.566881e-6)).abs()
            < AstronomicalUnits::new(1e-12)
    );
    assert!((orbit.eccentricity - 0.0549).abs() < 1e-6);
}

#[test]
fn orbital_period() {
    let earth_period = EARTH.orbit.period();
    // sidereal year ~ 365.25 days
    let days = earth_period.to::<Day>().value();
    assert!((days - 365.25).abs() < 0.01);

    let mars_period = MARS.orbit.period();
    let days = mars_period.to::<Day>().value();
    assert!((days - 686.98).abs() < 0.05);
}

#[test]
fn asteroid_builder_defaults() {
    let orbit = Orbit::new(
        AstronomicalUnits::new(1.0),
        0.1,
        Degrees::new(1.0),
        Degrees::new(2.0),
        Degrees::new(3.0),
        Degrees::new(4.0),
        JulianDate::J2000,
    );
    let asteroid = Asteroid::builder()
        .name("Test")
        .designation("T1")
        .orbit(orbit)
        .build();
    assert_eq!(asteroid.composition, "Unknown");
    assert_eq!(asteroid.class, AsteroidClass::MainBelt);
}

#[test]
fn comet_builder_defaults_and_period() {
    let orbit = HALLEY.orbit;
    let comet = Comet::builder()
        .name("X/1 Test")
        .tail_length(Kilometers::new(1.0e6))
        .orbit(orbit)
        .build();
    assert_eq!(comet.reference, OrbitFrame::Heliocentric);
    let period = HALLEY.period_years();
    assert!((period - 75.0).abs() < 5.0);
}

#[test]
fn planet_builder_errors_and_period() {
    use std::f64::consts::PI;

    let orbit = Orbit::new(
        AstronomicalUnits::new(1.0),
        0.0,
        Degrees::new(0.0),
        Degrees::new(0.0),
        Degrees::new(0.0),
        Degrees::new(0.0),
        JulianDate::J2000,
    );
    let err = Planet::builder()
        .mass(Kilograms::new(1.0))
        .orbit(orbit)
        .try_build()
        .unwrap_err();
    assert!(matches!(err, PlanetBuilderError::MissingRadius));

    let planet = Planet::builder()
        .mass(Kilograms::new(1.0))
        .radius(Kilometers::new(1.0))
        .orbit(orbit)
        .build();
    let k = 0.017_202_098_95_f64; // rad/day
    let period_days = planet.orbit.period().value() / 86_400.0;
    let expected_days = (2.0 * PI / k) * (1.0f64).powf(1.5);
    assert!((period_days - expected_days).abs() < 1e-6);
}
