use siderust::bodies::solar_system::Mars;
use siderust::coordinates::*;
use siderust::coordinates::centers::*;
use siderust::coordinates::frames::*;
use siderust::units::*;
use siderust::astro::JulianDate;

fn approx_eq<C, F, U>(a: &cartesian::Vector<C, F, U>, b: &cartesian::Vector<C, F, U>)
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
    Quantity<U>: std::cmp::PartialOrd + std::fmt::Display,
{
    assert!((a.x() - b.x()).abs() < (1e-6).into(), "x mismatch: {} vs {}", a.x(), b.x());
    assert!((a.y() - b.y()).abs() < (1e-6).into(), "y mismatch: {} vs {}", a.y(), b.y());
    assert!((a.z() - b.z()).abs() < (1e-6).into(), "z mismatch: {} vs {}", a.z(), b.z());
}

fn sph_approx_eq<C, F, U>(a: &spherical::Position<C, F, U>, b: &spherical::Position<C, F, U>)
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
    Quantity<U>: std::cmp::PartialOrd + std::fmt::Display,
{
    assert!((a.polar    - b.polar   ).abs().value() < 1e-6, "polar mismatch: {} vs {}", a.polar, b.polar);
    assert!((a.azimuth  - b.azimuth ).abs().value() < 1e-6, "polar mismatch: {} vs {}", a.azimuth, b.azimuth);
    assert!((a.distance - b.distance).abs().value() < 1e-6, "polar mismatch: {} vs {}", a.distance, b.distance);
}


#[test]
fn test_coord_transformations() {
    use siderust::coordinates::cartesian::Direction;

    let original: Direction<Heliocentric, Ecliptic> = Mars::vsop87a(JulianDate::J2000).get_position().clone().direction(); // Heliocentric, Ecliptic

    // Heliocentric Ecliptic -> Geocentric Ecliptic -> back
    let geo_ecl = Direction::<Geocentric, Ecliptic>::from(&original);
    let helio_ecl = Direction::<Heliocentric, Ecliptic>::from(&geo_ecl);
    approx_eq(&original, &helio_ecl);

    // Heliocentric Ecliptic -> Barycentric Ecliptic -> back
    let bary_ecl = Direction::<Barycentric, Ecliptic>::from(&original);
    let helio_from_bary = Direction::<Heliocentric, Ecliptic>::from(&bary_ecl);
    approx_eq(&original, &helio_from_bary);

    // Heliocentric Ecliptic -> Heliocentric Equatorial -> back
    let helio_eq = Direction::<Heliocentric, Equatorial>::from(&original);
    let helio_from_eq = Direction::<Heliocentric, Ecliptic>::from(&helio_eq);
    approx_eq(&original, &helio_from_eq);

    // Heliocentric Equatorial -> ICRS -> back
    let icrs = Direction::<Heliocentric, frames::ICRS>::from(&helio_eq);
    let helio_eq_back = Direction::<Heliocentric, Equatorial>::from(&icrs);
    approx_eq(&helio_eq, &helio_eq_back);
}


#[test]
fn test_spherical_transformations() {
    let original = Mars::vsop87a(JulianDate::J2000).get_position().clone().direction(); // Heliocentric, Ecliptic

    let geo_eq = cartesian::Direction::<Geocentric, Equatorial>::from(&original);
    let helio_icrs = cartesian::Direction::<Heliocentric, frames::ICRS>::from(&original);
    let bary_eq = cartesian::Direction::<Barycentric, Equatorial>::from(&original);

    let sph_helio_ecl = original.to_spherical();
    let back_helio_ecl = sph_helio_ecl.to_cartesian();
    approx_eq(&original, &back_helio_ecl);

    let sph_geo_eq = spherical::Direction::from_cartesian(&geo_eq);
    let back_geo_eq = sph_geo_eq.to_cartesian();
    approx_eq(&geo_eq, &back_geo_eq);

    let sph_helio_icrs = spherical::Direction::from_cartesian(&helio_icrs);
    let back_helio_icrs = sph_helio_icrs.to_cartesian();
    approx_eq(&helio_icrs, &back_helio_icrs);

    let sph_bary_eq = spherical::Direction::from_cartesian(&bary_eq);
    let back_bary_eq = sph_bary_eq.to_cartesian();
    approx_eq(&bary_eq, &back_bary_eq);
}


#[test]
fn serialize_cartesian_spherical() {
    let sph_orig = spherical::Position::<Barycentric, frames::ICRS, AstronomicalUnit>::new(Degrees::new(101.28715533), Degrees::new(-16.71611586), 1.0);
    let cart = sph_orig.to_cartesian();
    let sph_rec = cart.to_spherical();
    sph_approx_eq(&sph_orig, &sph_rec);
}