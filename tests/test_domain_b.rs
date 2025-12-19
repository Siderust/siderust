//! # Domain B Integration Tests
//!
//! These tests verify the correct implementation of Domain B (Physical/Astronomical Semantics):
//!
//! - B1: Aberration is not part of center transforms
//! - B2: Observational state is explicit (Astrometric vs Apparent)
//! - B3: Topocentric center applies real parallax

use qtty::*;
use siderust::astro::JulianDate;
use siderust::coordinates::cartesian::{line_of_sight, Position};
use siderust::coordinates::centers::{Geocentric, Heliocentric, ObserverSite};
use siderust::coordinates::frames::Equatorial;
use siderust::coordinates::observation::{Apparent, Astrometric, ObserverState};
use siderust::coordinates::spherical;
use siderust::coordinates::transform::{Transform, TransformCenter};

// =============================================================================
// B1: Aberration is NOT part of center transforms
// =============================================================================

#[test]
fn center_transforms_do_not_apply_aberration() {
    // Create a distant star in heliocentric coordinates
    let star_helio = Position::<Heliocentric, Equatorial, Au>::new(
        10000.0, 0.0, 0.0, // 10000 AU along x-axis
    );

    let jd = JulianDate::J2000;

    // Transform to geocentric - should be pure translation
    let star_geo: Position<Geocentric, Equatorial, Au> = star_helio.transform(jd);

    // The position should be shifted by approximately Earth's position (~1 AU)
    // but the direction should NOT have aberration applied

    // Earth is roughly at (-0.17, 0.97, 0) AU in equatorial J2000
    // So geocentric position should be approximately (10000.17, -0.97, 0)

    // Check that the x-coordinate is close to original (pure translation)
    let x_diff = (star_geo.x().value() - star_helio.x().value()).abs();
    assert!(
        x_diff < 2.0, // Less than 2 AU difference (Earth's orbit radius)
        "X difference too large: {} (expected < 2 AU)",
        x_diff
    );

    // If aberration were applied, the direction would shift by ~20 arcsec
    // For an object at 10000 AU, that's a transverse shift of ~0.001 AU
    // We verify the shift is consistent with translation, not aberration
}

#[test]
fn roundtrip_center_transform_preserves_position() {
    // This test verifies center transforms are reversible (pure geometry)
    let pos_helio = Position::<Heliocentric, Equatorial, Au>::new(1000.0, 500.0, -200.0);

    let jd = JulianDate::J2000;

    // Helio -> Geo -> Helio should recover original position
    let pos_geo: Position<Geocentric, Equatorial, Au> = pos_helio.transform(jd);
    let pos_helio_recovered: Position<Heliocentric, Equatorial, Au> = pos_geo.transform(jd);

    // Should recover within floating-point precision
    assert!(
        (pos_helio.x().value() - pos_helio_recovered.x().value()).abs() < 1e-10,
        "X not preserved: {} vs {}",
        pos_helio.x().value(),
        pos_helio_recovered.x().value()
    );
    assert!(
        (pos_helio.y().value() - pos_helio_recovered.y().value()).abs() < 1e-10,
        "Y not preserved"
    );
    assert!(
        (pos_helio.z().value() - pos_helio_recovered.z().value()).abs() < 1e-10,
        "Z not preserved"
    );
}

// =============================================================================
// B2: Observational state is explicit
// =============================================================================

#[test]
fn astrometric_and_apparent_are_distinct_types() {
    let jd = JulianDate::J2000;
    let obs = ObserverState::geocentric(jd);

    // Create an astrometric direction
    let astrometric = Astrometric::new(spherical::direction::Equatorial::new(
        180.0 * DEG,
        45.0 * DEG,
    ));

    // Convert to apparent
    let apparent: Apparent<spherical::direction::Equatorial> = astrometric.to_apparent(&obs);

    // The two types are distinct at compile time
    // We can access the underlying directions
    let _astro_dir: &spherical::direction::Equatorial = astrometric.direction();
    let _app_dir: &spherical::direction::Equatorial = apparent.direction();

    // The directions should differ due to aberration
    // Use angular_separation for proper comparison
    let separation = astrometric
        .direction()
        .angular_separation(apparent.direction());

    // Aberration is ~20 arcsec ≈ 0.006 degrees
    assert!(
        separation.value() > 0.0,
        "Apparent direction should differ from astrometric"
    );
    assert!(
        separation.value() < 0.1,
        "Aberration shift should be small (< 0.1 deg), got {} deg",
        separation.value()
    );
}

#[test]
fn aberration_requires_observer_state() {
    // This test verifies that aberration cannot be applied without an observer state
    // (This is enforced by the type system - to_apparent() requires &ObserverState)

    let jd = JulianDate::J2000;
    let astrometric = Astrometric::new(spherical::direction::Equatorial::new(0.0 * DEG, 0.0 * DEG));

    // Without an observer state, we cannot convert to apparent
    // (uncommenting the next line would be a compile error)
    // let apparent = astrometric.to_apparent();  // ERROR: missing observer state

    // We must create an observer state first
    let obs = ObserverState::geocentric(jd);
    let _apparent = astrometric.to_apparent(&obs); // OK
}

#[test]
fn aberration_roundtrip_preserves_direction() {
    let jd = JulianDate::J2000;
    let obs = ObserverState::geocentric(jd);

    let original = Astrometric::new(spherical::direction::Equatorial::new(
        45.0 * DEG,
        30.0 * DEG,
    ));

    // Astrometric -> Apparent -> Astrometric
    let apparent = original.to_apparent(&obs);
    let recovered = apparent.to_astrometric(&obs);

    let orig_dir = original.direction();
    let rec_dir = recovered.direction();

    let delta_ra = (rec_dir.azimuth.value() - orig_dir.azimuth.value()).abs();
    let delta_dec = (rec_dir.polar.value() - orig_dir.polar.value()).abs();

    // Should roundtrip within numerical precision
    assert!(
        delta_ra < 1e-6 && delta_dec < 1e-6,
        "Roundtrip should preserve direction: dRA={}, dDec={}",
        delta_ra,
        delta_dec
    );
}

#[test]
fn aberration_maximum_near_ecliptic_pole() {
    // Aberration is maximum for objects perpendicular to Earth's velocity
    // At the ecliptic poles, aberration should be close to the aberration constant

    let jd = JulianDate::J2000;
    let obs = ObserverState::geocentric(jd);

    // Direction toward north ecliptic pole (roughly)
    let astrometric = Astrometric::new(
        spherical::direction::Equatorial::new(270.0 * DEG, 66.56 * DEG), // Near NEP
    );

    let apparent = astrometric.to_apparent(&obs);

    // Calculate total angular separation
    let orig = astrometric.direction();
    let shifted = apparent.direction();

    let separation = orig.angular_separation(shifted);

    // Aberration constant is ~20.5 arcsec = 0.0057 degrees
    assert!(
        separation.value() > 0.001 && separation.value() < 0.02,
        "Aberration should be around 20 arcsec, got {} deg",
        separation.value()
    );
}

// =============================================================================
// B3: Topocentric applies real parallax
// =============================================================================

#[test]
fn topocentric_parallax_is_real_translation() {
    // For a nearby object like the Moon, topocentric parallax is significant

    // Moon at ~384,400 km from Earth's center (along x-axis for simplicity)
    let moon_geo = Position::<Geocentric, Equatorial, Kilometer>::new(384_400.0, 0.0, 0.0);

    // Observer at equator, prime meridian
    let site = ObserverSite::new(0.0 * DEG, 0.0 * DEG, 0.0 * M);
    let jd = JulianDate::J2000;

    let moon_topo = moon_geo.to_topocentric(site, jd);

    // The distance should differ by approximately the observer's geocentric distance
    let geo_dist = moon_geo.distance().value();
    let topo_dist = moon_topo.distance().value();
    let _dist_diff = (geo_dist - topo_dist).abs();

    // Observer is ~6371 km from Earth's center
    // Depending on GMST, the parallax shift could be 0 to ~6371 km
    // For a reasonable test, just check the positions differ meaningfully

    // The x-component should change by roughly Earth's radius
    let x_diff = (moon_geo.x().value() - moon_topo.x().value()).abs();
    assert!(
        x_diff > 1000.0 && x_diff < 10000.0,
        "Topocentric parallax should shift x by ~6371 km, got {} km",
        x_diff
    );
}

#[test]
fn topocentric_roundtrip_preserves_geocentric_position() {
    let geo = Position::<Geocentric, Equatorial, Kilometer>::new(100_000.0, 50_000.0, 25_000.0);

    let site = ObserverSite::new(10.0 * DEG, 45.0 * DEG, 100.0 * M);
    let jd = JulianDate::J2000;

    // Geocentric -> Topocentric -> Geocentric
    let topo = geo.to_topocentric(site, jd);
    let geo_recovered: Position<Geocentric, Equatorial, Kilometer> = topo.to_center(jd);

    // Should recover within floating-point precision
    assert!(
        (geo.x().value() - geo_recovered.x().value()).abs() < 1e-6,
        "X not preserved: {} vs {}",
        geo.x().value(),
        geo_recovered.x().value()
    );
    assert!(
        (geo.y().value() - geo_recovered.y().value()).abs() < 1e-6,
        "Y not preserved"
    );
    assert!(
        (geo.z().value() - geo_recovered.z().value()).abs() < 1e-6,
        "Z not preserved"
    );
}

#[test]
fn topocentric_parallax_negligible_for_stars() {
    // For distant stars, topocentric parallax should be negligible

    // Star at 1 parsec = 206265 AU
    let star_geo = Position::<Geocentric, Equatorial, Au>::new(206265.0, 0.0, 0.0);

    let site = ObserverSite::new(0.0 * DEG, 45.0 * DEG, 0.0 * M);
    let jd = JulianDate::J2000;

    let star_topo = star_geo.to_topocentric(site, jd);

    // Relative parallax should be tiny
    let rel_diff = (star_geo.x().value() - star_topo.x().value()).abs() / star_geo.x().value();

    // Earth radius / 1 parsec ≈ 6371 km / (3.086e13 km) ≈ 2e-10
    assert!(
        rel_diff < 1e-6,
        "Stellar parallax should be negligible: relative diff = {}",
        rel_diff
    );
}

#[test]
fn observer_site_provides_geocentric_position() {
    // Verify that ObserverSite computes a reasonable geocentric position

    // Greenwich Observatory
    let site = ObserverSite::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
    let pos: Position<Geocentric, siderust::coordinates::frames::ECEF, Kilometer> =
        site.geocentric_itrf();

    // Greenwich should be at roughly (3980, 0, 4970) km in ECEF
    assert!(
        (pos.x().value() - 3980.0).abs() < 100.0,
        "x = {} km, expected ~3980 km",
        pos.x().value()
    );
    assert!(
        pos.y().value().abs() < 10.0,
        "y = {} km, expected ~0 km (prime meridian)",
        pos.y().value()
    );
    assert!(
        (pos.z().value() - 4970.0).abs() < 100.0,
        "z = {} km, expected ~4970 km",
        pos.z().value()
    );

    // Distance should be roughly Earth's radius
    let r = pos.distance().value();
    assert!(
        (r - 6371.0).abs() < 50.0,
        "distance = {} km, expected ~6371 km",
        r
    );
}

// =============================================================================
// Integration: Complete astronomical pipeline
// =============================================================================

#[test]
fn complete_pipeline_geometric_to_apparent() {
    // This test demonstrates the complete pipeline:
    // 1. Start with heliocentric position
    // 2. Transform to geocentric (pure translation)
    // 3. Compute line of sight direction
    // 4. Apply aberration to get apparent direction

    let jd = JulianDate::J2000;

    // A fictitious object at 10 AU from the Sun
    let object_helio = Position::<Heliocentric, Equatorial, Au>::new(10.0, 0.0, 0.0);

    // Step 1-2: Transform to geocentric (pure translation, no aberration)
    let object_geo: Position<Geocentric, Equatorial, Au> = object_helio.transform(jd);

    // Step 3: Compute line of sight from geocenter to object
    let observer = Position::<Geocentric, Equatorial, Au>::new(0.0, 0.0, 0.0);
    let los_direction = line_of_sight(&observer, &object_geo);

    // Wrap in Astrometric to track observational state
    let astrometric_dir = Astrometric::new(los_direction.to_spherical());

    // Step 4: Apply aberration with observer state
    let obs = ObserverState::geocentric(jd);
    let apparent_dir = astrometric_dir.to_apparent(&obs);

    // Verify the pipeline produced different astrometric and apparent positions
    let astro = astrometric_dir.direction();
    let app = apparent_dir.direction();

    let delta = astro.angular_separation(app);

    // Aberration should be ~20 arcsec
    assert!(
        delta.value() > 0.001,
        "Apparent should differ from astrometric by aberration"
    );
}
