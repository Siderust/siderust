use siderust::bodies::planets::{Planet, PlanetBuilderError, OrbitExt};
use siderust::units::{Kilograms, Kilometers, AstronomicalUnits, Degrees};
use siderust::astro::{orbit::Orbit, JulianDate};

#[test]
fn planet_builder_errors() {
    let builder = Planet::builder().radius(Kilometers::new(1.0)).orbit(
        Orbit::new(
            AstronomicalUnits::new(1.0),
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        )
    );
    let err = builder.clone().try_build().unwrap_err();
    assert!(matches!(err, PlanetBuilderError::MissingMass));

    let err = Planet::builder().mass(Kilograms::new(1.0)).try_build().unwrap_err();
    assert!(matches!(err, PlanetBuilderError::MissingRadius));
}

#[test]
fn orbit_period_computation() {
    let orbit = Orbit::new(
        AstronomicalUnits::new(1.0),
        0.0,
        Degrees::new(0.0),
        Degrees::new(0.0),
        Degrees::new(0.0),
        Degrees::new(0.0),
        JulianDate::J2000,
    );
    let planet = Planet::builder()
        .mass(Kilograms::new(1.0))
        .radius(Kilometers::new(1.0))
        .orbit(orbit.clone())
        .build();

    let p = planet.orbit.period().value();
    let expected = 2.0 * std::f64::consts::PI / 0.986 * 1.0_f64.powf(1.5) * 86400.0;
    assert!((p - expected).abs() < 1e-6);
}
