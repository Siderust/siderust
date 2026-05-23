// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

#![allow(missing_docs)]

use siderust::bodies::solar_system::Mars;
use siderust::coordinates::centers::*;
use siderust::coordinates::frames::*;
use siderust::coordinates::transform::TransformFrame;
use siderust::coordinates::*;
use siderust::qtty::*;

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
        (a.distance - b.distance).abs() < (1e-6).into(),
        "polar mismatch: {} vs {}",
        a.distance,
        b.distance
    );
}

/// Test position coordinate transformations (positions support both center and frame transforms)
#[test]
fn test_position_transformations() {
    use siderust::coordinates::transform::Transform;

    let original = Mars::vsop87a(siderust::time::J2000);

    // Heliocentric EclipticMeanJ2000 -> Heliocentric EquatorialMeanJ2000 -> back
    let helio_eq: cartesian::Position<Heliocentric, EquatorialMeanJ2000, _> =
        original.transform(siderust::time::J2000);
    let helio_from_eq: cartesian::Position<Heliocentric, EclipticMeanJ2000, _> =
        helio_eq.transform(siderust::time::J2000);
    approx_eq_pos(&original, &helio_from_eq);
}

/// Test direction frame transformations (directions only support frame transforms, not center transforms)
#[test]
fn test_direction_frame_transformations() {
    use siderust::coordinates::cartesian::Direction;

    // Create a direction from Mars position
    let original: Direction<EclipticMeanJ2000> = Mars::vsop87a(siderust::time::J2000)
        .direction()
        .expect("Mars position should have a direction");

    // EclipticMeanJ2000 -> EquatorialMeanJ2000 -> back (frame rotation only)
    let equatorial: Direction<EquatorialMeanJ2000> = TransformFrame::to_frame(&original);
    let ecliptic_back: Direction<EclipticMeanJ2000> = TransformFrame::to_frame(&equatorial);
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
    let cart_original: cartesian::Direction<EclipticMeanJ2000> =
        Mars::vsop87a(siderust::time::J2000)
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
    let observer = cartesian::Position::<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>::new(
        AstronomicalUnits::new(1.0),
        AstronomicalUnits::new(0.0),
        AstronomicalUnits::new(0.0),
    );

    let target = cartesian::Position::<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>::new(
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

#[test]
fn transform_frame_trait_preserves_expected_metadata() {
    use siderust::coordinates::centers::Heliocentric;

    let vec_ecl = cartesian::position::EclipticMeanJ2000::<AstronomicalUnit>::new(
        0.1 * AU,
        0.2 * AU,
        0.3 * AU,
    );
    let vec_same: cartesian::position::EclipticMeanJ2000<AstronomicalUnit> =
        TransformFrame::to_frame(&vec_ecl);
    assert_eq!(vec_same.x(), vec_ecl.x());

    let sph_ecl =
        spherical::direction::EclipticMeanJ2000::new(Degrees::new(10.0), Degrees::new(5.0));
    let cart_ecl = sph_ecl.to_cartesian();
    let cart_equatorial: cartesian::direction::EquatorialMeanJ2000 =
        TransformFrame::to_frame(&cart_ecl);
    assert!(cart_equatorial.x().is_finite());

    let vec_equatorial: cartesian::Position<Heliocentric, EquatorialMeanJ2000, AstronomicalUnit> =
        TransformFrame::to_frame(&vec_ecl);
    assert!(vec_equatorial.x().value().is_finite());
}

#[test]
fn topocentric_horizontal_position_transform_preserves_distance() {
    use siderust::coordinates::transform::{
        providers::frame_rotation, AstroContext, Transform, TransformCenter,
    };

    let observer = Geodetic::<ECEF>::new(
        Degrees::new(-17.89),
        Degrees::new(28.76),
        Meters::new(2400.0),
    );
    let jd = siderust::time::J2000;
    let site = observer;

    let eq_pos = spherical::position::EquatorialMeanJ2000::<AstronomicalUnit>::new(
        Degrees::new(83.0),
        Degrees::new(-5.0),
        1.0,
    );
    let cart_pos = eq_pos.to_cartesian();
    let topo_cart_j2000: cartesian::Position<Topocentric, EquatorialMeanJ2000, AstronomicalUnit> =
        cart_pos.to_center((site, jd));

    let ctx = AstroContext::default();
    let rot = frame_rotation::<EquatorialMeanJ2000, EquatorialMeanOfDate>(jd, &ctx);
    let [x, y, z] = rot
        * [
            topo_cart_j2000.x(),
            topo_cart_j2000.y(),
            topo_cart_j2000.z(),
        ];
    let topo_cart_mod =
        cartesian::Position::<Topocentric, EquatorialMeanOfDate, AstronomicalUnit>::new_with_params(
            *topo_cart_j2000.center_params(),
            x,
            y,
            z,
        );

    let horiz_cart_pos: cartesian::position::Horizontal<AstronomicalUnit> =
        topo_cart_mod.transform(jd);
    let horiz_pos = horiz_cart_pos.to_spherical();

    assert!((horiz_pos.distance - eq_pos.distance).abs().value() < 1e-4);
    assert!(horiz_cart_pos.z().value().is_finite());
}
