use qtty::*;
use siderust::astro::JulianDate;
use siderust::bodies::solar_system::Mars;
use siderust::coordinates::{
    cartesian,
    centers::*,
    frames::*,
    spherical,
    transform::{Transform, TransformFrame},
};
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
    let sph_pos: spherical::Position<Heliocentric, Ecliptic, AstronomicalUnit> = 
        spherical::Position::from_cartesian(&cart_orig.position);
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
fn cartesian_direction_frame_transform() {
    // Directions are now frame-only (no center parameter).
    // They can only undergo frame transformations, not center transformations.
    let dir = cartesian::Direction::<Ecliptic>::normalize(1.0, 0.5, 0.2);

    // Frame transform from Ecliptic to Equatorial (rotation only)
    let dir_equatorial: cartesian::Direction<Equatorial> = TransformFrame::to_frame(&dir);

    // Verify it's still a unit vector
    let norm = (dir_equatorial.x().powi(2)
        + dir_equatorial.y().powi(2)
        + dir_equatorial.z().powi(2))
    .sqrt();
    assert!((norm - 1.0).abs() < 1e-12);

    // X component is unchanged in ecliptic/equatorial rotation
    assert!((dir_equatorial.x() - dir.x()).abs() < 1e-12);
}

#[test]
fn spherical_direction_frame_transform() {
    // Directions are now frame-only (no center parameter).
    let sph_dir = spherical::Direction::<Ecliptic>::new(Degrees::new(10.0), Degrees::new(20.0));

    // Convert to cartesian, transform frame, then back
    let cart_dir = sph_dir.to_cartesian();
    let cart_equatorial: cartesian::Direction<Equatorial> = TransformFrame::to_frame(&cart_dir);

    // Verify the Cartesian direction is still unit vector
    let norm = (cart_equatorial.x().powi(2)
        + cart_equatorial.y().powi(2)
        + cart_equatorial.z().powi(2))
    .sqrt();
    assert!((norm - 1.0).abs() < 1e-12);
}
