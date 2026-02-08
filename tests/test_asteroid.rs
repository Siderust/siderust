// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use qtty::*;
use siderust::astro::orbit::Orbit;
use siderust::bodies::asteroid::{
    Asteroid, AsteroidClass, APOPHIS, ASTEROID_PRESETS, BENNU, CERES_AST,
};
use siderust::time::JulianDate;

const TEST_ORBIT: Orbit = Orbit::new(
    AstronomicalUnits::new(1.0),
    0.0,
    Degrees::new(0.0),
    Degrees::new(0.0),
    Degrees::new(0.0),
    Degrees::new(0.0),
    JulianDate::J2000,
);

#[test]
fn builder_sets_all_fields() {
    let orbit = TEST_ORBIT;
    let asteroid = Asteroid::builder()
        .name("MyAst")
        .designation("M-01")
        .composition("Carbonaceous")
        .class(AsteroidClass::NearEarth)
        .orbit(orbit)
        .build();

    assert_eq!(asteroid.name, "MyAst");
    assert_eq!(asteroid.designation, "M-01");
    assert_eq!(asteroid.composition, "Carbonaceous");
    assert_eq!(asteroid.class, AsteroidClass::NearEarth);
    assert_eq!(asteroid.orbit.semi_major_axis, AstronomicalUnits::new(1.0));
}

#[test]
#[should_panic(expected = "missing name")]
fn builder_missing_name_panics() {
    let orbit = TEST_ORBIT;
    let _ = Asteroid::builder().designation("D-1").orbit(orbit).build();
}

#[test]
#[should_panic(expected = "missing designation")]
fn builder_missing_designation_panics() {
    let orbit = TEST_ORBIT;
    let _ = Asteroid::builder().name("Test").orbit(orbit).build();
}

#[test]
#[should_panic(expected = "missing orbit")]
fn builder_missing_orbit_panics() {
    let _ = Asteroid::builder().name("Test").designation("D-1").build();
}

const CONST_ASTEROID: Asteroid = Asteroid::new_const(
    "ConstAst",
    "C-1",
    "Silicate",
    AsteroidClass::MainBelt,
    TEST_ORBIT,
);

#[test]
fn const_constructor_and_presets() {
    assert_eq!(CONST_ASTEROID.name, "ConstAst");
    assert_eq!(CONST_ASTEROID.designation, "C-1");
    assert_eq!(CONST_ASTEROID.composition, "Silicate");
    assert_eq!(CONST_ASTEROID.class, AsteroidClass::MainBelt);
    assert_eq!(CONST_ASTEROID.orbit.eccentricity, 0.0);

    assert_eq!(ASTEROID_PRESETS.len(), 3);
    assert!(std::ptr::eq(ASTEROID_PRESETS[0], &CERES_AST));
    assert!(std::ptr::eq(ASTEROID_PRESETS[1], &BENNU));
    assert!(std::ptr::eq(ASTEROID_PRESETS[2], &APOPHIS));
}
