use qtty::*;
use siderust::astro::JulianDate;
use siderust::bodies::solar_system::Mars;
use siderust::coordinates::{cartesian, centers::*, frames::*, spherical, transform::Transform};
use siderust::targets::Target;

const EPS: f64 = 1e-9;

#[test]
fn target_cartesian_position_transform() {
    let jd = JulianDate::J2000;
    let orig: Target<cartesian::Position<Heliocentric, Ecliptic, AstronomicalUnit>> =
        Mars::vsop87a(jd);

    let converted: Target<cartesian::Position<Geocentric, Equatorial, AstronomicalUnit>> =
        Target::from(&orig);

    let step: cartesian::Position<Heliocentric, Equatorial, AstronomicalUnit> =
        orig.position.transform(jd);
    let expected: cartesian::Position<Geocentric, Equatorial, AstronomicalUnit> =
        step.transform(jd);

    assert!(converted.position.distance_to(&expected).value() < EPS);
    assert_eq!(converted.time, orig.time);
    assert_eq!(
        converted.proper_motion.is_none(),
        orig.proper_motion.is_none()
    );
}

#[test]
fn target_spherical_position_transform() {
    let jd = JulianDate::J2000;
    let cart_orig: Target<cartesian::Position<Heliocentric, Ecliptic, AstronomicalUnit>> =
        Mars::vsop87a(jd);
    let sph_pos = cart_orig.position.to_spherical();
    let orig: Target<spherical::Position<Heliocentric, Ecliptic, AstronomicalUnit>> =
        Target::new_static(sph_pos, jd);

    let converted: Target<spherical::Position<Geocentric, Equatorial, AstronomicalUnit>> =
        Target::from(&orig);

    let step_cart: cartesian::Position<Heliocentric, Equatorial, AstronomicalUnit> =
        orig.position.to_cartesian().transform(jd);
    let expected_cart: cartesian::Position<Geocentric, Equatorial, AstronomicalUnit> =
        step_cart.transform(jd);
    let converted_cart = converted.position.to_cartesian();

    assert!(converted_cart.distance_to(&expected_cart).value() < EPS);
    assert_eq!(converted.time, orig.time);
    assert_eq!(
        converted.proper_motion.is_none(),
        orig.proper_motion.is_none()
    );
}

#[test]
fn target_cartesian_direction_transform() {
    let jd = JulianDate::J2000;
    let pos = cartesian::Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(
        AstronomicalUnits::new(1.0),
        AstronomicalUnits::new(0.5),
        AstronomicalUnits::new(0.2),
    );
    let dir = pos.direction();
    let orig: Target<cartesian::Direction<Heliocentric, Ecliptic>> = Target::new_static(dir, jd);

    let converted: Target<cartesian::Direction<Geocentric, Equatorial>> = Target::from(&orig);

    let step: cartesian::Direction<Heliocentric, Equatorial> = orig.position.transform(jd);
    let expected: cartesian::Direction<Geocentric, Equatorial> = step.transform(jd);

    assert!(converted.position.distance_to(&expected).value() < EPS);
    assert_eq!(converted.time, orig.time);
    assert_eq!(
        converted.proper_motion.is_none(),
        orig.proper_motion.is_none()
    );
}

#[test]
fn target_spherical_direction_transform() {
    let jd = JulianDate::J2000;
    let sph_dir =
        spherical::Direction::<Heliocentric, Ecliptic>::new(Degrees::new(10.0), Degrees::new(20.0));
    let orig: Target<spherical::Direction<Heliocentric, Ecliptic>> =
        Target::new_static(sph_dir, jd);

    let converted: Target<spherical::Direction<Geocentric, Equatorial>> = Target::from(&orig);

    let step_cart: cartesian::Direction<Heliocentric, Equatorial> =
        orig.position.to_cartesian().transform(jd);
    let expected_cart: cartesian::Direction<Geocentric, Equatorial> = step_cart.transform(jd);
    let converted_cart = converted.position.to_cartesian();

    assert!(converted_cart.distance_to(&expected_cart).value() < EPS);
    assert_eq!(converted.time, orig.time);
    assert_eq!(
        converted.proper_motion.is_none(),
        orig.proper_motion.is_none()
    );
}
