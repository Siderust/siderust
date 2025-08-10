use siderust::astro::orbit::Orbit;
use siderust::astro::JulianDate;
use siderust::units::{AstronomicalUnits, Degrees, Kilometers, Kilograms};
use siderust::bodies::asteroid::{Asteroid, AsteroidClass};
use siderust::bodies::comet::{Comet, OrbitFrame, HALLEY};
use siderust::bodies::planets::{Planet, PlanetBuilderError, OrbitExt};

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
    let period_days = planet.orbit.period().value() / 86_400.0;
    let expected_days = 2.0 * std::f64::consts::PI / 0.986 * (1.0f64).powf(1.5);
    assert!((period_days - expected_days).abs() < 1e-6);
}
