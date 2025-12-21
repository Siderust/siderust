use qtty::{AstronomicalUnit, AstronomicalUnits, Days, Degrees, Kilograms, Kilometers, Years, AU};
use siderust::astro::orbit::Orbit;
use siderust::astro::JulianDate;
use siderust::bodies::asteroid::{Asteroid, AsteroidClass};
use siderust::bodies::comet::{Comet, CometBuilder, OrbitFrame};
use siderust::bodies::planets::{Planet, PlanetBuilder};
use siderust::coordinates::transform::centers::position::to_topocentric::ToTopocentricExt;
use siderust::coordinates::{
    cartesian,
    centers::ObserverSite,
    frames, spherical,
    transform::{Transform, TransformFrame},
};

#[test]
fn julian_date_arithmetic_and_display_branches() {
    let mut jd = JulianDate::new(2_450_000.0);
    let printed = format!("{jd}");
    assert!(printed.contains("Julian Day"));

    jd += Days::new(2.0);
    jd -= Days::new(0.5);

    let with_years = jd + Years::new(1.0);
    let day_span: Days = with_years - jd;
    assert!((day_span.value() - JulianDate::JULIAN_YEAR.value()).abs() < 1e-9);

    let ratio = with_years / Days::new(2.0);
    assert!((ratio - with_years.value() / 2.0).abs() < 1e-12);
    let ratio_plain = with_years / 2.0;
    assert!((ratio_plain - with_years.value() / 2.0).abs() < 1e-12);

    let min = with_years.min(jd);
    assert_eq!(min, jd);

    let utc = jd.to_utc().expect("UTC conversion should succeed");
    let roundtrip = JulianDate::from_utc(utc);
    assert!((roundtrip.value() - jd.value()).abs() < 1e-6);
}

#[test]
fn cartesian_vector_display_includes_metadata() {
    let v = cartesian::position::ICRS::<AstronomicalUnit>::new(1.0 * AU, 2.0 * AU, 3.0 * AU);
    let rendered = format!("{v}");
    assert!(rendered.contains("Center: Barycentric"));
    assert!(rendered.contains("Frame: ICRS"));
    assert!(rendered.contains("X: 1"));
}

#[test]
fn horizontal_conversion_variants_cover_all_impls() {
    let observer = spherical::position::Geographic::new(
        Degrees::new(-17.89),
        Degrees::new(28.76),
        Kilometers::new(2.4),
    );
    let jd = JulianDate::J2000;

    // Convert Geographic to ObserverSite for the new API
    let site = ObserverSite::from_geographic(&observer);

    // Test position (with distance) conversion - positions still support center transforms
    let eq_pos = spherical::position::Equatorial::<AstronomicalUnit>::new(
        Degrees::new(83.0),
        Degrees::new(-5.0),
        1.0,
    );
    let cart_pos = eq_pos.to_cartesian();
    let topo_cart_pos = cart_pos.to_topocentric(site, jd);
    let horiz_cart_pos: cartesian::position::Horizontal<AstronomicalUnit> =
        topo_cart_pos.transform(jd);
    let horiz_pos = horiz_cart_pos.to_spherical();
    // Distance changes slightly due to real topocentric parallax (observer is ~6000 km from Earth center)
    // For an object at 1 AU, this is a very small fractional change (Earth radius / 1 AU ≈ 4e-5)
    assert!((horiz_pos.distance - eq_pos.distance).abs().value() < 1e-4);
    assert!(horiz_cart_pos.z().value().is_finite());

    // Note: Directions no longer support center transforms (to_topocentric).
    // Directions are free vectors - they can only undergo frame transformations.
}

#[test]
fn frame_transform_traits_exercised() {
    use siderust::coordinates::centers::Heliocentric;

    let vec_ecl =
        cartesian::position::Ecliptic::<AstronomicalUnit>::new(0.1 * AU, 0.2 * AU, 0.3 * AU);
    let vec_same: cartesian::position::Ecliptic<AstronomicalUnit> =
        TransformFrame::to_frame(&vec_ecl);
    assert_eq!(vec_same.x(), vec_ecl.x());

    // Spherical direction is now frame-only (no center parameter)
    let sph_ecl = spherical::direction::Ecliptic::new(Degrees::new(10.0), Degrees::new(5.0));
    // Convert to cartesian, transform frame
    let cart_ecl = sph_ecl.to_cartesian();
    let cart_equatorial: cartesian::direction::Equatorial = TransformFrame::to_frame(&cart_ecl);
    assert!(cart_equatorial.x().is_finite());

    // Test frame transform on position (must preserve center type)
    let vec_equatorial: cartesian::Position<Heliocentric, frames::Equatorial, AstronomicalUnit> =
        TransformFrame::to_frame(&vec_ecl);
    assert!(vec_equatorial.x().value().is_finite());
}

#[test]
fn body_const_constructors_and_builders() {
    let orbit = Orbit::new(
        AstronomicalUnits::new(1.0),
        0.01,
        Degrees::new(1.0),
        Degrees::new(2.0),
        Degrees::new(3.0),
        Degrees::new(4.0),
        JulianDate::J2000,
    );

    let asteroid = Asteroid::new_const("Test", "T-1", "Rock", AsteroidClass::NearEarth, orbit);
    assert_eq!(asteroid.designation, "T-1");

    let comet = Comet::new_const(
        "TestComet",
        Kilometers::new(1_234.0),
        orbit,
        OrbitFrame::Barycentric,
    );
    assert!(comet.period_years() > 0.0);

    let comet_from_builder = CometBuilder::default()
        .name("Builder")
        .tail_length(Kilometers::new(1.0))
        .reference(OrbitFrame::Heliocentric)
        .orbit(orbit)
        .build();
    assert_eq!(comet_from_builder.reference, OrbitFrame::Heliocentric);

    let planet = Planet::new_const(Kilograms::new(5.0e24), Kilometers::new(6_000.0), orbit);
    assert!((planet.radius.value() - 6_000.0).abs() < 1e-6);

    let unchecked = PlanetBuilder::default()
        .mass(Kilograms::new(1.0e24))
        .radius(Kilometers::new(3_000.0))
        .orbit(orbit)
        .build_unchecked();
    assert!((unchecked.mass.value() - 1.0e24).abs() < 1e-6);
}

#[test]
#[should_panic(expected = "PlanetBuilder::build_unchecked called with missing fields")]
fn planet_build_unchecked_missing_fields_panics() {
    let _ = PlanetBuilder::default().build_unchecked();
}
