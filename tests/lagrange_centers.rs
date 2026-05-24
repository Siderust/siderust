//! Integration tests for feature-gated Sun-Earth Lagrange centers.

#![cfg(feature = "lagrange-centers")]

use siderust::coordinates::cartesian::Position;
use siderust::coordinates::centers::{
    Barycentric, ReferenceCenter, SunEarthL1, SunEarthL2, SunEarthL3, SunEarthL4, SunEarthL5,
};
use siderust::coordinates::frames::EclipticMeanJ2000;
use siderust::coordinates::transform::Transform;
use siderust::ephemeris::lagrange::{try_sun_earth_lagrange_barycentric, SunEarthLagrangePoint};
use siderust::qtty::AstronomicalUnit;

#[test]
fn sun_earth_lagrange_centers_are_unit_parameterized() {
    fn assert_unit_params<C: ReferenceCenter<Params = ()>>() {}
    assert_unit_params::<SunEarthL1>();
    assert_unit_params::<SunEarthL2>();
    assert_unit_params::<SunEarthL3>();
    assert_unit_params::<SunEarthL4>();
    assert_unit_params::<SunEarthL5>();
    let names = [
        SunEarthL1::center_name(),
        SunEarthL2::center_name(),
        SunEarthL3::center_name(),
        SunEarthL4::center_name(),
        SunEarthL5::center_name(),
    ];
    for (idx, name) in names.iter().enumerate() {
        assert!(names[(idx + 1)..].iter().all(|other| other != name));
    }
}

#[test]
fn placeholder_archive_returns_out_of_range() {
    let err = try_sun_earth_lagrange_barycentric(SunEarthLagrangePoint::L1, siderust::J2000)
        .expect_err("placeholder records are intentionally empty");
    assert!(matches!(
        err,
        siderust::ephemeris::EphemerisError::OutOfRange { .. }
    ));
}

#[test]
fn barycentric_lagrange_roundtrips() {
    let bary = Position::<Barycentric, EclipticMeanJ2000, AstronomicalUnit>::new(0.1, -0.2, 0.3);
    roundtrip::<SunEarthL1>(&bary);
    roundtrip::<SunEarthL2>(&bary);
    roundtrip::<SunEarthL3>(&bary);
    roundtrip::<SunEarthL4>(&bary);
    roundtrip::<SunEarthL5>(&bary);
}

fn roundtrip<C>(bary: &Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>)
where
    C: ReferenceCenter<Params = ()>,
    Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>:
        Transform<Position<C, EclipticMeanJ2000, AstronomicalUnit>>,
    Position<C, EclipticMeanJ2000, AstronomicalUnit>:
        Transform<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>>,
{
    let shifted: Position<C, EclipticMeanJ2000, AstronomicalUnit> = bary.transform(siderust::J2000);
    let back: Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> =
        shifted.transform(siderust::J2000);
    assert!((back.x().value() - bary.x().value()).abs() < 1.0e-8);
    assert!((back.y().value() - bary.y().value()).abs() < 1.0e-8);
    assert!((back.z().value() - bary.z().value()).abs() < 1.0e-8);
}
