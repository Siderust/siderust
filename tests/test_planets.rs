// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use siderust::astro::orbit::KeplerianOrbit;
use siderust::bodies::planets::{OrbitExt, Planet, PlanetBuilderError};
use siderust::qtty::*;

#[test]
fn planet_builder_errors() {
    let builder = Planet::builder()
        .radius(Kilometers::new(1.0))
        .orbit(KeplerianOrbit::new(
            AstronomicalUnits::new(1.0),
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            siderust::time::J2000,
        ));
    let err = builder.clone().try_build().unwrap_err();
    assert!(matches!(err, PlanetBuilderError::MissingMass));

    let err = Planet::builder()
        .mass(Kilograms::new(1.0))
        .try_build()
        .unwrap_err();
    assert!(matches!(err, PlanetBuilderError::MissingRadius));
}

#[test]
fn orbit_period_computation() {
    let k = 0.017_202_098_95_f64; // rad/day

    let orbit = KeplerianOrbit::new(
        AstronomicalUnits::new(1.0),
        0.0,
        Degrees::new(0.0),
        Degrees::new(0.0),
        Degrees::new(0.0),
        Degrees::new(0.0),
        siderust::time::J2000,
    );
    let planet = Planet::builder()
        .mass(Kilograms::new(1.0))
        .radius(Kilometers::new(1.0))
        .orbit(orbit)
        .build();

    let p = planet.orbit.period().value();
    let expected = 2.0 * std::f64::consts::PI / k * 1.0_f64.powf(1.5) * 86400.0;
    // The inherent rounding between the 11-digit k and the stored Gaussian-year
    // constant (365.256898326 days) is ~2.8×10⁻⁵ s; use 1e-3 s tolerance.
    assert!(
        (p - expected).abs() < 1e-3,
        "period {p} vs expected {expected}, diff = {}",
        (p - expected).abs()
    );
}

#[test]
fn const_constructor_and_unchecked_builder() {
    let orbit = KeplerianOrbit::new(
        AstronomicalUnits::new(1.0),
        0.01,
        Degrees::new(1.0),
        Degrees::new(2.0),
        Degrees::new(3.0),
        Degrees::new(4.0),
        siderust::time::J2000,
    );

    let planet = Planet::new_const(Kilograms::new(5.0e24), Kilometers::new(6_000.0), orbit);
    assert!((planet.mass.value() - 5.0e24).abs() < 1e-6);
    assert!((planet.radius.value() - 6_000.0).abs() < 1e-6);

    let unchecked = Planet::builder()
        .mass(Kilograms::new(1.0e24))
        .radius(Kilometers::new(3_000.0))
        .orbit(orbit)
        .build_unchecked();
    assert!((unchecked.mass.value() - 1.0e24).abs() < 1e-6);
    assert!((unchecked.radius.value() - 3_000.0).abs() < 1e-6);
}

#[test]
#[should_panic(expected = "PlanetBuilder::build_unchecked called with missing fields")]
fn build_unchecked_missing_fields_panics() {
    let _ = Planet::builder().build_unchecked();
}
