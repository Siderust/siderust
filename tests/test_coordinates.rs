use siderust::bodies::solar_system::Mars;
use siderust::coordinates::*;
use siderust::coordinates::centers::*;
use siderust::coordinates::frames::*;
use siderust::units::{JulianDay, Degrees};

fn approx_eq<C, F>(a: &CartesianCoord<C, F>, b: &CartesianCoord<C, F>)
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    assert!((a.x() - b.x()).abs() < 1e-6, "x mismatch: {} vs {}", a.x(), b.x());
    assert!((a.y() - b.y()).abs() < 1e-6, "y mismatch: {} vs {}", a.y(), b.y());
    assert!((a.z() - b.z()).abs() < 1e-6, "z mismatch: {} vs {}", a.z(), b.z());
}

fn sph_approx_eq<C, F>(a: &SphericalCoord<C, F>, b: &SphericalCoord<C, F>)
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    assert!((a.polar.as_f64()   - b.polar.as_f64()).abs()   < 1e-6, "polar mismatch: {} vs {}", a.polar, b.polar);
    assert!((a.azimuth.as_f64() - b.azimuth.as_f64()).abs() < 1e-6, "polar mismatch: {} vs {}", a.azimuth, b.azimuth);
    assert!((a.radial_distance  - b.radial_distance).abs()  < 1e-6, "polar mismatch: {} vs {}", a.radial_distance, b.radial_distance);
}


#[test]
fn test_coord_transformations() {
    let original = Mars::vsop87a(JulianDay::J2000).get_position().clone(); // Heliocentric, Ecliptic (implícito)

    // Heliocentric Ecliptic -> Geocentric Ecliptic -> back
    let geo_ecl = CartesianCoord::<Geocentric, Ecliptic>::from(&original);
    let helio_ecl = CartesianCoord::<Heliocentric, Ecliptic>::from(&geo_ecl);
    approx_eq(&original, &helio_ecl);

    // Heliocentric Ecliptic -> Barycentric Ecliptic -> back
    let bary_ecl = CartesianCoord::<Barycentric, Ecliptic>::from(&original);
    let helio_from_bary = CartesianCoord::<Heliocentric, Ecliptic>::from(&bary_ecl);
    approx_eq(&original, &helio_from_bary);

    // Heliocentric Ecliptic -> Heliocentric Equatorial -> back
    let helio_eq = CartesianCoord::<Heliocentric, Equatorial>::from(&original);
    let helio_from_eq = CartesianCoord::<Heliocentric, Ecliptic>::from(&helio_eq);
    approx_eq(&original, &helio_from_eq);

    // Heliocentric Equatorial -> ICRS -> back
    let icrs = CartesianCoord::<Heliocentric, frames::ICRS>::from(&helio_eq);
    let helio_eq_back = CartesianCoord::<Heliocentric, Equatorial>::from(&icrs);
    approx_eq(&helio_eq, &helio_eq_back);
}


#[test]
fn test_spherical_transformations() {
    let original = Mars::vsop87a(JulianDay::J2000).get_position().clone(); // Heliocentric, Ecliptic

    // Convertir a otros marcos cartesianos
    let geo_eq = CartesianCoord::<Geocentric, Equatorial>::from(&original);
    let helio_icrs = CartesianCoord::<Heliocentric, frames::ICRS>::from(&original);
    let bary_eq = CartesianCoord::<Barycentric, Equatorial>::from(&original);

    // Probar conversiones Cartesiano -> Esférico -> Cartesiano
    let sph_helio_ecl = SphericalCoord::from_cartesian(&original);
    let back_helio_ecl = sph_helio_ecl.to_cartesian();
    approx_eq(&original, &back_helio_ecl);

    let sph_geo_eq = SphericalCoord::from_cartesian(&geo_eq);
    let back_geo_eq = sph_geo_eq.to_cartesian();
    approx_eq(&geo_eq, &back_geo_eq);

    let sph_helio_icrs = SphericalCoord::from_cartesian(&helio_icrs);
    let back_helio_icrs = sph_helio_icrs.to_cartesian();
    approx_eq(&helio_icrs, &back_helio_icrs);

    let sph_bary_eq = SphericalCoord::from_cartesian(&bary_eq);
    let back_bary_eq = sph_bary_eq.to_cartesian();
    approx_eq(&bary_eq, &back_bary_eq);
}


#[test]
fn serialize_cartesian_spherical() {
    let sph_orig = SphericalCoord::<Barycentric, frames::ICRS>::new(Degrees::new(101.28715533), Degrees::new(-16.71611586), 1.0);
    let cart = sph_orig.to_cartesian();
    let sph_rec = cart.to_spherical();
    sph_approx_eq(&sph_orig, &sph_rec);
}