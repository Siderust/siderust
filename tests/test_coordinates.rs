// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use qtty::*;
use siderust::astro::JulianDate;
use siderust::bodies::solar_system::Mars;
use siderust::coordinates::centers::*;
use siderust::coordinates::frames::*;
use siderust::coordinates::transform::TransformFrame;
use siderust::coordinates::*;

fn approx_eq_pos<C, F, U>(a: &cartesian::Position<C, F, U>, b: &cartesian::Position<C, F, U>)
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
    Quantity<U>: std::cmp::PartialOrd,
{
    assert!(
        (a.x() - b.x()).abs() < (1e-6).into(),
        "x mismatch: {} vs {}",
        a.x().value(),
        b.x().value()
    );
    assert!(
        (a.y() - b.y()).abs() < (1e-6).into(),
        "y mismatch: {} vs {}",
        a.y().value(),
        b.y().value()
    );
    assert!(
        (a.z() - b.z()).abs() < (1e-6).into(),
        "z mismatch: {} vs {}",
        a.z().value(),
        b.z().value()
    );
}

fn approx_eq_dir<F>(a: &cartesian::Direction<F>, b: &cartesian::Direction<F>)
where
    F: ReferenceFrame,
{
    assert!(
        (a.x() - b.x()).abs() < 1e-6,
        "x mismatch: {} vs {}",
        a.x(),
        b.x()
    );
    assert!(
        (a.y() - b.y()).abs() < 1e-6,
        "y mismatch: {} vs {}",
        a.y(),
        b.y()
    );
    assert!(
        (a.z() - b.z()).abs() < 1e-6,
        "z mismatch: {} vs {}",
        a.z(),
        b.z()
    );
}

fn sph_approx_eq<C, F, U>(a: &spherical::Position<C, F, U>, b: &spherical::Position<C, F, U>)
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
    Quantity<U>: std::cmp::PartialOrd + std::fmt::Display,
{
    assert!(
        (a.polar - b.polar).abs().value() < 1e-6,
        "polar mismatch: {} vs {}",
        a.polar,
        b.polar
    );
    assert!(
        (a.azimuth - b.azimuth).abs().value() < 1e-6,
        "polar mismatch: {} vs {}",
        a.azimuth,
        b.azimuth
    );
    assert!(
        (a.distance - b.distance).abs().value() < 1e-6,
        "polar mismatch: {} vs {}",
        a.distance,
        b.distance
    );
}

/// Test position coordinate transformations (positions support both center and frame transforms)
#[test]
fn test_position_transformations() {
    use siderust::coordinates::transform::Transform;

    let original = *Mars::vsop87a(JulianDate::J2000).get_position();

    // Heliocentric Ecliptic -> Heliocentric EquatorialMeanJ2000 -> back
    let helio_eq: cartesian::Position<Heliocentric, EquatorialMeanJ2000, _> =
        original.transform(JulianDate::J2000);
    let helio_from_eq: cartesian::Position<Heliocentric, Ecliptic, _> =
        helio_eq.transform(JulianDate::J2000);
    approx_eq_pos(&original, &helio_from_eq);
}

/// Test direction frame transformations (directions only support frame transforms, not center transforms)
#[test]
fn test_direction_frame_transformations() {
    use siderust::coordinates::cartesian::Direction;

    // Create a direction from Mars position
    let original: Direction<Ecliptic> = Mars::vsop87a(JulianDate::J2000)
        .get_position()
        .clone()
        .direction()
        .expect("Mars position should have a direction");

    // Ecliptic -> EquatorialMeanJ2000 -> back (frame rotation only)
    let equatorial: Direction<EquatorialMeanJ2000> = TransformFrame::to_frame(&original);
    let ecliptic_back: Direction<Ecliptic> = TransformFrame::to_frame(&equatorial);
    approx_eq_dir(&original, &ecliptic_back);

    // Verify directions are still unit vectors after transformation
    let norm = (equatorial.x().powi(2) + equatorial.y().powi(2) + equatorial.z().powi(2)).sqrt();
    assert!(
        (norm - 1.0).abs() < 1e-12,
        "direction should be unit vector"
    );
}

/// Test spherical direction transformations
#[test]
fn test_spherical_direction_transformations() {
    // Create a direction from Mars position
    let cart_original: cartesian::Direction<Ecliptic> = Mars::vsop87a(JulianDate::J2000)
        .get_position()
        .clone()
        .direction()
        .expect("Mars position should have a direction");

    // Convert to spherical and back
    let sph = cart_original.to_spherical();
    let back = sph.to_cartesian();
    approx_eq_dir(&cart_original, &back);
}

/// Test spherical position round-trip
#[test]
fn serialize_cartesian_spherical() {
    let sph_orig = spherical::Position::<Barycentric, ICRS, AstronomicalUnit>::new(
        Degrees::new(101.28715533),
        Degrees::new(-16.71611586),
        1.0,
    );
    let cart = sph_orig.to_cartesian();
    let sph_rec = spherical::Position::from_cartesian(&cart);
    sph_approx_eq(&sph_orig, &sph_rec);
}

/// Test line_of_sight function for computing direction from observer to target
#[test]
fn test_line_of_sight() {
    use siderust::coordinates::cartesian::line_of_sight;

    // Create two positions (Earth and Mars at some point)
    let observer = cartesian::Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(
        AstronomicalUnits::new(1.0),
        AstronomicalUnits::new(0.0),
        AstronomicalUnits::new(0.0),
    );

    let target = cartesian::Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(
        AstronomicalUnits::new(2.0),
        AstronomicalUnits::new(0.0),
        AstronomicalUnits::new(0.0),
    );

    let los = line_of_sight(&observer, &target);

    // Direction should point in +X direction
    assert!((los.x() - 1.0).abs() < 1e-12);
    assert!(los.y().abs() < 1e-12);
    assert!(los.z().abs() < 1e-12);
}
