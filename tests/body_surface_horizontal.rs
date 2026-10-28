// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Scientific validation of generic body-surface geometric Alt/Az.
//!
//! ## Provenance
//!
//! Reference cases below are computed from first-principles ENU geometry
//! using siderust's published IAU body-fixed orientation and reference
//! surfaces. They do **not** require network access.
//!
//! Where a case claims agreement with an external tool (e.g. SPICE), the
//! body-fixed frame realisation must match. Siderust's lunar marker
//! `MoonIauMeanEarth` (legacy alias `MoonPrincipalAxes`) uses the NAIF
//! IAU_MOON analytical Mean-Earth model, not the JPL DE `MOON_PA` or
//! binary-PCK `MOON_ME` realizations.
//!
//! Geometric scope only: no refraction, aberration, light-time, or terrain.

use siderust::bodies::solar_system::{Earth, Mars, Moon};
use siderust::coordinates::cartesian::{Direction, Position};
use siderust::coordinates::centers::{Geocentric, Marscentric, Selenocentric};
use siderust::coordinates::frames::{MarsFixed, MoonIauMeanEarth, ECEF, ICRS};
use siderust::coordinates::surface::{
    BodySurfaceModel, EnuBasis, LatitudeConvention, SurfaceSite, EARTH_WGS84_SURFACE,
};
use siderust::coordinates::transform::context::{DefaultEop, DefaultEphemeris};
use siderust::coordinates::transform::{
    direction_to_horizontal, horizontal_to_direction, position_to_horizontal, AstroContext,
    DirectionAstroExt, SurfaceHorizontalExt,
};
use siderust::qtty::*;
use siderust::time::JulianDate;

type Ctx = AstroContext<DefaultEphemeris, DefaultEop>;

fn ctx() -> Ctx {
    AstroContext::default()
}

fn chord<F: siderust::coordinates::frames::ReferenceFrame>(
    a: &Direction<F>,
    b: &Direction<F>,
) -> f64 {
    ((a.x() - b.x()).powi(2) + (a.y() - b.y()).powi(2) + (a.z() - b.z()).powi(2)).sqrt()
}

fn assert_near(actual: f64, expected: f64, tol: f64, label: &str) {
    assert!(
        (actual - expected).abs() < tol,
        "{label}: got {actual}, expected {expected} (tol {tol})"
    );
}

// =============================================================================
// ENU orthonormality & poles
// =============================================================================

#[test]
fn enu_basis_orthonormal_earth_mars_moon() {
    for (lon, lat) in [
        (0.0, 0.0),
        (45.0, 30.0),
        (-120.0, -60.0),
        (0.0, 89.9),
        (10.0, -89.9),
    ] {
        let moon = SurfaceSite::<Moon>::new(lon * DEG, lat * DEG, 0.0 * M).unwrap();
        moon.enu_basis().unwrap().assert_orthonormal(1e-12).unwrap();

        let mars = SurfaceSite::<Mars>::new(lon * DEG, lat * DEG, 0.0 * M).unwrap();
        mars.enu_basis().unwrap().assert_orthonormal(1e-12).unwrap();

        let earth = SurfaceSite::<Earth>::new(lon * DEG, lat * DEG, 0.0 * M).unwrap();
        earth
            .enu_basis()
            .unwrap()
            .assert_orthonormal(1e-12)
            .unwrap();
    }
}

#[test]
fn north_pole_azimuth_follows_longitude_meridian() {
    let site = SurfaceSite::<Moon>::new(30.0 * DEG, 90.0 * DEG, 0.0 * M).unwrap();
    let enu = site.enu_basis().unwrap();
    // up ≈ +Z
    assert_near(enu.up[2], 1.0, 1e-12, "up.z");
    // east = (−sin λ, cos λ, 0)
    let lon = 30.0_f64.to_radians();
    assert_near(enu.east[0], -lon.sin(), 1e-12, "east.x");
    assert_near(enu.east[1], lon.cos(), 1e-12, "east.y");
}

// =============================================================================
// Cardinal directions (body-fixed, no orientation) — Moon sphere
// =============================================================================

#[test]
fn moon_equator_cardinals_body_fixed() {
    let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let jd = siderust::J2000;
    let ctx = ctx();

    // Body-fixed unit vectors at (0,0): +X = up, +Y = east, +Z = north
    let cases = [
        ([1.0, 0.0, 0.0], 90.0, 0.0, "up/+X"),
        ([0.0, 1.0, 0.0], 0.0, 90.0, "east/+Y"),
        ([0.0, 0.0, 1.0], 0.0, 0.0, "north/+Z"),
        ([0.0, -1.0, 0.0], 0.0, 270.0, "west/−Y"),
        ([0.0, 0.0, -1.0], 0.0, 180.0, "south/−Z"),
        ([-1.0, 0.0, 0.0], -90.0, 0.0, "nadir/−X"),
    ];
    for (v, alt_deg, az_deg, label) in cases {
        let fixed = Direction::<MoonIauMeanEarth>::new_unchecked(v[0], v[1], v[2]);
        let icrs: Direction<ICRS> = fixed.to_frame_with(&jd, &ctx);
        let horiz = direction_to_horizontal(&icrs, &site, &jd, &ctx).unwrap();
        let sph = horiz.to_spherical();
        assert_near(sph.alt().value(), alt_deg, 1e-6, &format!("{label} alt"));
        // Azimuth undefined at zenith/nadir (we adopt 0°); skip az check there.
        if alt_deg.abs() < 89.0 {
            let daz = (sph.az().value() - az_deg).abs();
            let daz = daz.min(360.0 - daz);
            assert!(
                daz < 1e-4,
                "{label} az: got {}, expected {az_deg}",
                sph.az()
            );
        }
    }
}

// =============================================================================
// Absolute IAU pole-orientation reference at J2000 (not a round-trip test)
// =============================================================================

#[test]
fn iau_body_poles_match_published_j2000_right_ascension_declination() {
    let jd = siderust::J2000;
    let ctx = ctx();

    // NAIF pck00011.tpc IAU_MOON (ME) including all 13 periodic terms;
    // The +Z body-fixed axis is the spin pole, independently of W and of
    // the observer's longitude. These reference vectors are computed from
    // the tabulated angles, not by inverting the same provider rotation.
    let pole_reference = |ra_deg: f64, dec_deg: f64| -> [f64; 3] {
        let ra = ra_deg.to_radians();
        let dec = dec_deg.to_radians();
        [dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin()]
    };

    let moon_pole: Direction<ICRS> =
        Direction::<MoonIauMeanEarth>::new_unchecked(0.0, 0.0, 1.0).to_frame_with(&jd, &ctx);
    let moon_expected = pole_reference(266.85773344495135, 65.64110274784532);
    assert_near(moon_pole.x(), moon_expected[0], 1e-10, "Moon pole ICRS x");
    assert_near(moon_pole.y(), moon_expected[1], 1e-10, "Moon pole ICRS y");
    assert_near(moon_pole.z(), moon_expected[2], 1e-10, "Moon pole ICRS z");

    let mars_pole: Direction<ICRS> =
        Direction::<MarsFixed>::new_unchecked(0.0, 0.0, 1.0).to_frame_with(&jd, &ctx);
    let mars_expected = pole_reference(317.269, 54.432);
    assert_near(mars_pole.x(), mars_expected[0], 1e-10, "Mars pole ICRS x");
    assert_near(mars_pole.y(), mars_expected[1], 1e-10, "Mars pole ICRS y");
    assert_near(mars_pole.z(), mars_expected[2], 1e-10, "Mars pole ICRS z");
}

// =============================================================================
// Absolute IAU prime-meridian orientation (independent of the provider)
// =============================================================================

#[test]
fn iau_body_prime_meridians_match_j2000_cartographic_axes() {
    // IAU body-fixed +X is the reference-equator node rotated eastward by W.
    // This is a direct published-angle reference, not an inverse round-trip.
    fn expected_x(ra_deg: f64, dec_deg: f64, w_deg: f64) -> [f64; 3] {
        let (sa, ca) = ra_deg.to_radians().sin_cos();
        let (sd, cd) = dec_deg.to_radians().sin_cos();
        let (sw, cw) = w_deg.to_radians().sin_cos();
        [-sa * cw - sd * ca * sw, ca * cw - sd * sa * sw, cd * sw]
    }

    let jd = siderust::J2000;
    let ctx = ctx();
    let moon_x: Direction<ICRS> =
        Direction::<MoonIauMeanEarth>::new_unchecked(1.0, 0.0, 0.0).to_frame_with(&jd, &ctx);
    let mars_x: Direction<ICRS> =
        Direction::<MarsFixed>::new_unchecked(1.0, 0.0, 0.0).to_frame_with(&jd, &ctx);

    for (actual, expected, label) in [
        (
            moon_x,
            expected_x(266.85773344495135, 65.64110274784532, 41.19526396965),
            "Moon",
        ),
        (mars_x, expected_x(317.269, 54.432, 176.049), "Mars"),
    ] {
        assert_near(actual.x(), expected[0], 1e-10, &format!("{label} X-axis x"));
        assert_near(actual.y(), expected[1], 1e-10, &format!("{label} X-axis y"));
        assert_near(actual.z(), expected[2], 1e-10, &format!("{label} X-axis z"));
    }
}

#[test]
fn mars_altaz_invariant_under_equivalent_latitude_conventions() {
    let s = Mars::default_surface();
    let graphic = 45.0_f64;
    let centric = ((1.0 - s.e2()) * graphic.to_radians().tan())
        .atan()
        .to_degrees();
    let site_g = SurfaceSite::<Mars>::with_surface(
        137.0 * DEG,
        graphic * DEG,
        100.0 * M,
        s,
        LatitudeConvention::Planetographic,
    )
    .unwrap();
    let site_c = SurfaceSite::<Mars>::with_surface(
        137.0 * DEG,
        centric * DEG,
        100.0 * M,
        s,
        LatitudeConvention::Planetocentric,
    )
    .unwrap();

    let target = Direction::<ICRS>::new_unchecked(0.0, 0.0, 1.0);
    let jd = siderust::J2000;
    let ctx = ctx();
    let g = direction_to_horizontal(&target, &site_g, &jd, &ctx)
        .unwrap()
        .to_spherical();
    let c = direction_to_horizontal(&target, &site_c, &jd, &ctx)
        .unwrap()
        .to_spherical();

    assert_near(
        g.alt().value(),
        c.alt().value(),
        1e-10,
        "equivalent Mars alt",
    );
    let daz = (g.az().value() - c.az().value()).abs();
    assert!(daz.min(360.0 - daz) < 1e-10);
}

// =============================================================================
// Finite-distance parallax
// =============================================================================

#[test]
fn moon_finite_parallax_shifts_altitude() {
    let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let jd = siderust::J2000;
    let ctx = ctx();

    // Same sky direction (+Z body-fixed): direction-only vs finite near target.
    let dir_fixed = Direction::<MoonIauMeanEarth>::new_unchecked(0.0, 0.0, 1.0);
    let dir_icrs: Direction<ICRS> = dir_fixed.to_frame_with(&jd, &ctx);
    let dir_alt = direction_to_horizontal(&dir_icrs, &site, &jd, &ctx)
        .unwrap()
        .to_spherical()
        .alt()
        .value();

    let near = Position::<Selenocentric, MoonIauMeanEarth, Kilometer>::new(0.0, 0.0, 2_000.0);
    let near_alt = position_to_horizontal(&near, &site, &jd, &ctx)
        .unwrap()
        .alt()
        .value();

    assert_near(dir_alt, 0.0, 1e-6, "direction alt");
    assert!(
        near_alt < dir_alt - 10.0,
        "finite near target must show large negative parallax vs direction-only: dir={dir_alt}, near={near_alt}"
    );
}

#[test]
fn direction_does_not_invent_parallax() {
    let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let jd = siderust::J2000;
    let ctx = ctx();
    let fixed = Direction::<MoonIauMeanEarth>::new_unchecked(0.0, 0.0, 1.0);
    let icrs: Direction<ICRS> = fixed.to_frame_with(&jd, &ctx);
    let horiz = site.observe_direction(&icrs, &jd, &ctx).unwrap();
    assert_near(horiz.to_spherical().alt().value(), 0.0, 1e-8, "dir alt");
    assert_eq!(horiz.site_params(), &site.to_params());
    let other = SurfaceSite::<Moon>::new(15.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let bad: Result<Direction<ICRS>, _> =
        siderust::coordinates::transform::horizontal_observation_to_direction(
            &horiz, &other, &jd, &ctx,
        );
    assert!(matches!(
        bad,
        Err(siderust::coordinates::surface::SurfaceError::MismatchedObserver)
    ));
    let recovered: Direction<ICRS> =
        siderust::coordinates::transform::horizontal_observation_to_direction(
            &horiz, &site, &jd, &ctx,
        )
        .unwrap();
    assert!(chord(&icrs, &recovered) < 1e-8);
}

// =============================================================================
// Forward / inverse round-trips
// =============================================================================

#[test]
fn moon_direction_roundtrip_random_sky() {
    let site = SurfaceSite::<Moon>::new((-15.0) * DEG, 42.0 * DEG, 100.0 * M).unwrap();
    let jd = JulianDate::new(2_460_000.5);
    let ctx = ctx();

    let samples = [
        (0.1, 0.2, 0.97),
        (-0.5, 0.3, 0.81),
        (0.7, -0.2, 0.68),
        (0.0, 0.0, 1.0),
    ];
    for (x, y, z) in samples {
        let n = (x * x + y * y + z * z).sqrt();
        let icrs = Direction::<ICRS>::new_unchecked(x / n, y / n, z / n);
        let horiz = direction_to_horizontal(&icrs, &site, &jd, &ctx).unwrap();
        let back: Direction<ICRS> = horizontal_to_direction(&horiz, &site, &jd, &ctx).unwrap();
        let residual = chord(&icrs, &back);
        assert!(residual < 1e-10, "roundtrip chord residual {residual}");
    }
}

#[test]
fn mars_direction_roundtrip() {
    let site = SurfaceSite::<Mars>::new(75.0 * DEG, (-15.0) * DEG, 0.0 * M).unwrap();
    let jd = siderust::J2000;
    let ctx = ctx();
    let n = (0.2_f64.powi(2) + 0.5_f64.powi(2) + 0.84_f64.powi(2)).sqrt();
    let icrs = Direction::<ICRS>::new_unchecked(0.2 / n, -0.5 / n, 0.84 / n);
    let horiz = direction_to_horizontal(&icrs, &site, &jd, &ctx).unwrap();
    let back: Direction<ICRS> = horizontal_to_direction(&horiz, &site, &jd, &ctx).unwrap();
    assert!(chord(&icrs, &back) < 1e-10);
}

// =============================================================================
// Earth regression vs ENU geometry (WGS84 site)
// =============================================================================

#[test]
fn earth_wgs84_zenith_matches_geodetic_cartesian() {
    use siderust::coordinates::centers::Geodetic;

    let geo = Geodetic::<ECEF>::new(0.0 * DEG, 45.0 * DEG, 0.0 * M);
    let site = SurfaceSite::<Earth>::from_geodetic(geo).unwrap();
    let cart: Position<Geocentric, ECEF, Meter> = geo.to_cartesian();
    let enu = site.enu_basis().unwrap();
    // Target along the ellipsoid normal (local up), not the radius vector.
    let target = Position::<Geocentric, ECEF, Meter>::new(
        Quantity::<Meter>::new(cart.x().value() + enu.up[0] * 1_000.0),
        Quantity::<Meter>::new(cart.y().value() + enu.up[1] * 1_000.0),
        Quantity::<Meter>::new(cart.z().value() + enu.up[2] * 1_000.0),
    );
    let jd = siderust::J2000;
    let ctx = ctx();
    let obs = position_to_horizontal(&target, &site, &jd, &ctx).unwrap();
    assert_near(obs.alt().value(), 90.0, 1e-6, "earth zenith alt");
}

#[test]
fn earth_surface_matches_wgs84_constants() {
    assert_eq!(Earth::default_surface(), EARTH_WGS84_SURFACE);
    assert_eq!(
        Earth::default_latitude_convention(),
        LatitudeConvention::Planetographic
    );
}

// =============================================================================
// Mars planetographic vs planetocentric difference
// =============================================================================

#[test]
fn mars_lat_conventions_differ_on_ellipsoid() {
    let lon = 0.0 * DEG;
    let lat = 45.0 * DEG;
    let h = 0.0 * M;
    let surf = Mars::default_surface();

    let site_g =
        SurfaceSite::<Mars>::with_surface(lon, lat, h, surf, LatitudeConvention::Planetographic)
            .unwrap();
    let site_c =
        SurfaceSite::<Mars>::with_surface(lon, lat, h, surf, LatitudeConvention::Planetocentric)
            .unwrap();
    let pg = site_g.body_fixed_cartesian().unwrap();
    let pc = site_c.body_fixed_cartesian().unwrap();
    let dx = (pg[0] - pc[0]).value();
    let dy = (pg[1] - pc[1]).value();
    let dz = (pg[2] - pc[2]).value();
    let dist = (dx * dx + dy * dy + dz * dz).sqrt();
    // At 45° on Mars ellipsoid the two conventions differ by kilometres.
    assert!(
        dist > 1_000.0,
        "expected multi-km difference between conventions, got {dist} m"
    );
}

// =============================================================================
// Body rotation changes Alt/Az across epochs
// =============================================================================

#[test]
fn mars_azimuth_changes_with_epoch() {
    let site = SurfaceSite::<Mars>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let ctx = ctx();
    // Fixed ICRS direction
    let icrs = Direction::<ICRS>::new_unchecked(1.0, 0.0, 0.0);
    let jd1 = siderust::J2000;
    let jd2 = JulianDate::new(siderust::J2000.value() + 0.5); // half day later
    let h1 = direction_to_horizontal(&icrs, &site, &jd1, &ctx).unwrap();
    let h2 = direction_to_horizontal(&icrs, &site, &jd2, &ctx).unwrap();
    let daz = (h1.to_spherical().az().value() - h2.to_spherical().az().value()).abs();
    let daz = daz.min(360.0 - daz);
    // Mars rotates ~350.9°/day → half day ≈ 175°; expect a large az change.
    assert!(daz > 90.0, "expected large az change over 0.5 d, got {daz}");
}

// =============================================================================
// Invalid inputs
// =============================================================================

#[test]
fn rejects_non_finite_site() {
    let err = SurfaceSite::<Moon>::new(Degrees::new(f64::NAN), 0.0 * DEG, 0.0 * M);
    assert!(err.is_err());
}

#[test]
fn zero_range_is_error() {
    let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let [x, y, z] = site.body_fixed_cartesian().unwrap();
    let target = Position::<Selenocentric, MoonIauMeanEarth, Meter>::new(x, y, z);
    let err = position_to_horizontal(&target, &site, &siderust::J2000, &ctx());
    // A coincident observer and target has no meaningful sky direction.
    // Numerical frame round-trips must not turn it into a success.
    assert!(
        matches!(
            err,
            Err(siderust::coordinates::surface::SurfaceError::DegenerateGeometry { .. })
        ),
        "coincident position must fail; got {err:?}"
    );
}

// =============================================================================
// Unit-independent coincident-target tolerance
// =============================================================================

#[test]
fn finite_observation_range_is_invariant_under_length_units() {
    let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let [x, y, z] = site.body_fixed_cartesian().unwrap();
    // 1 km above the surface: 1 km = ~6.7e-9 AU. Comparing the
    // raw numeric range to 1e-6 would incorrectly reject the AU target.
    let target_m = Position::<Selenocentric, MoonIauMeanEarth, Meter>::new(
        Quantity::<Meter>::new(x.value() + 1_000.0),
        y,
        z,
    );
    let target_km = target_m.to_unit::<Kilometer>();
    let target_au = target_m.to_unit::<AstronomicalUnit>();
    let jd = siderust::J2000;
    let ctx = ctx();

    let obs_m = position_to_horizontal(&target_m, &site, &jd, &ctx).unwrap();
    let obs_km = position_to_horizontal(&target_km, &site, &jd, &ctx).unwrap();
    let obs_au = position_to_horizontal(&target_au, &site, &jd, &ctx).unwrap();

    assert_near(obs_m.range().value(), 1_000.0, 1e-3, "metre slant range");
    assert_near(
        obs_km.range().to::<Meter>().value(),
        1_000.0,
        1e-3,
        "km slant range",
    );
    assert_near(
        obs_au.range().to::<Meter>().value(),
        1_000.0,
        1e-2,
        "AU slant range",
    );
    assert_near(
        obs_m.alt().value(),
        obs_au.alt().value(),
        1e-6,
        "unit invariant altitude",
    );
}

// =============================================================================
// Independent analytical reference (Moon sphere ENU, no frame rotation)
// =============================================================================

#[test]
fn analytical_enu_projection_reference() {
    // Pure geometry check of los_to_horizontal via site ENU, independent of
    // orientation providers: observer at lunar equator, LOS = east+up.
    let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let EnuBasis { east, north: _, up } = site.enu_basis().unwrap();
    let los = [
        (east[0] + up[0]) / 2.0_f64.sqrt(),
        (east[1] + up[1]) / 2.0_f64.sqrt(),
        (east[2] + up[2]) / 2.0_f64.sqrt(),
    ];
    // Reconstruct via body-fixed direction → ICRS → horizontal at J2000.
    // Because we use the same ENU, alt should be 45°, az 90°.
    let n = (los[0] * los[0] + los[1] * los[1] + los[2] * los[2]).sqrt();
    let fixed = Direction::<MoonIauMeanEarth>::new_unchecked(los[0] / n, los[1] / n, los[2] / n);
    let jd = siderust::J2000;
    let ctx = ctx();
    let icrs: Direction<ICRS> = fixed.to_frame_with(&jd, &ctx);
    let horiz = direction_to_horizontal(&icrs, &site, &jd, &ctx).unwrap();
    let sph = horiz.to_spherical();
    assert_near(sph.alt().value(), 45.0, 1e-6, "alt 45");
    assert_near(sph.az().value(), 90.0, 1e-6, "az 90");
}

#[test]
fn mars_zenith_along_ellipsoid_normal() {
    let site = SurfaceSite::<Mars>::new(10.0 * DEG, 25.0 * DEG, 0.0 * M).unwrap();
    let [x, y, z] = site.body_fixed_cartesian().unwrap();
    // Radial extension is *not* the zenith at nonzero planetographic
    // latitude on an oblate ellipsoid.
    let radial_target = Position::<Marscentric, MarsFixed, Meter>::new(
        Quantity::<Meter>::new(x.value() * 2.0),
        Quantity::<Meter>::new(y.value() * 2.0),
        Quantity::<Meter>::new(z.value() * 2.0),
    );
    let radial = position_to_horizontal(&radial_target, &site, &siderust::J2000, &ctx()).unwrap();
    assert!(
        radial.alt().value() < 89.95,
        "radial direction must differ from planetographic zenith"
    );

    // The ellipsoid normal is the physical local up for this reference surface.
    let enu = site.enu_basis().unwrap();
    let normal_target = Position::<Marscentric, MarsFixed, Meter>::new(
        Quantity::<Meter>::new(x.value() + enu.up[0] * 1_000.0),
        Quantity::<Meter>::new(y.value() + enu.up[1] * 1_000.0),
        Quantity::<Meter>::new(z.value() + enu.up[2] * 1_000.0),
    );
    let normal = position_to_horizontal(&normal_target, &site, &siderust::J2000, &ctx()).unwrap();
    assert_near(
        normal.alt().value(),
        90.0,
        1e-6,
        "mars ellipsoid-normal zenith",
    );
    let nadir_target = Position::<Marscentric, MarsFixed, Meter>::new(
        Quantity::<Meter>::new(x.value() - enu.up[0] * 1_000.0),
        Quantity::<Meter>::new(y.value() - enu.up[1] * 1_000.0),
        Quantity::<Meter>::new(z.value() - enu.up[2] * 1_000.0),
    );
    let nadir = position_to_horizontal(&nadir_target, &site, &siderust::J2000, &ctx()).unwrap();
    assert_near(nadir.alt().value(), -90.0, 1e-6, "mars normal nadir");
}

// Compile-time: SurfaceSite is parameterized by body (Moon ≠ Mars).
#[test]
fn site_params_body_name() {
    let moon = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let mars = SurfaceSite::<Mars>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    assert_eq!(moon.body_name(), "Moon");
    assert_eq!(mars.body_name(), "Mars");
    // Params equality is within the same body centre type only.
    assert_eq!(moon.to_params(), moon.to_params());
    // The default params are body-specific, never implicitly lunar for Mars.
    let moon_default = siderust::coordinates::surface::SurfaceObserverParams::<Moon>::default();
    let mars_default = siderust::coordinates::surface::SurfaceObserverParams::<Mars>::default();
    assert_eq!(*moon_default.surface(), Moon::default_surface());
    assert_eq!(*mars_default.surface(), Mars::default_surface());
    assert_eq!(
        mars_default.latitude_convention(),
        LatitudeConvention::Planetographic
    );
}

#[test]
fn finite_east_matches_direction_only_and_affn_horizontal_axes() {
    // On the lunar equator at longitude zero, +Y is local east. A target
    // 1000 metres along +Y must report azimuth +90°, not 270°.
    let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let [x, y, z] = site.body_fixed_cartesian().unwrap();
    let target = Position::<Selenocentric, MoonIauMeanEarth, Meter>::new(x, y + 1_000.0 * M, z);
    let jd = siderust::J2000;
    let ctx = ctx();
    let finite = position_to_horizontal(&target, &site, &jd, &ctx).unwrap();
    let fixed_east = Direction::<MoonIauMeanEarth>::new_unchecked(0.0, 1.0, 0.0);
    let icrs_east: Direction<ICRS> = fixed_east.to_frame_with(&jd, &ctx);
    let direction = direction_to_horizontal(&icrs_east, &site, &jd, &ctx).unwrap();
    let sph = direction.to_spherical();

    assert_near(finite.alt().value(), 0.0, 1e-6, "finite east altitude");
    assert_near(finite.az().value(), 90.0, 1e-6, "finite east azimuth");
    assert_near(
        sph.az().value(),
        finite.az().value(),
        1e-6,
        "east azimuth parity",
    );
    assert_near(
        sph.alt().value(),
        finite.alt().value(),
        1e-6,
        "east altitude parity",
    );
    assert!(
        finite.position().y().value() > 0.0,
        "affn Horizontal +Y must mean east"
    );
}

#[test]
fn lunar_iau_mean_earth_orientation_matches_independent_pck_vectors() {
    // CSPICE N0067 via SpiceyPy 8.3.0, pck00011.tpc + naif0012.tls:
    // et = unitim((jd_tt - 2451545) * 86400, "TDT", "TDB")
    // axes = pxform("IAU_MOON", "J2000", et).T
    // All three axes, 1900--2050. The library uses tempoch's TT->TDB
    // series, rather than SPICE's DELTET approximation. 1e-10 in direction
    // cosine allows <40 microseconds of timescale-model disagreement.
    let samples = [
        (
            2451545.0,
            [
                [0.7842270522119748, 0.5578471123206701, 0.2716514860154575],
                [-0.6200619150990117, 0.720556665574856, 0.31035675139972596],
                [
                    -0.022608671404191448,
                    -0.4118309009425287,
                    0.9109797785934671,
                ],
            ],
        ),
        (
            2455197.5,
            [
                [
                    0.25823924831040457,
                    -0.8876830435656354,
                    -0.3812234315963898,
                ],
                [0.9657683594503107, 0.24724448879202301, 0.07849610593288507],
                [
                    0.024575730238036195,
                    -0.3884443035082652,
                    0.9211444276307847,
                ],
            ],
        ),
        (
            2447892.5,
            [
                [
                    -0.8802151956698334,
                    0.43320727597351794,
                    0.19378510096372514,
                ],
                [-0.4742629440570647, -0.817761083811191, -0.3261003981878464],
                [
                    0.017200848997720463,
                    -0.3789436182963777,
                    0.9252598904881844,
                ],
            ],
        ),
        (
            2415020.5,
            [
                [-0.1462181824757145, 0.904463719684885, 0.4007064048493364],
                [
                    -0.9889031878073268,
                    -0.14439998435015727,
                    -0.03491603734935458,
                ],
                [0.02628170957159834, -0.4013652006507987, 0.9155409589136568],
            ],
        ),
        (
            2469807.5,
            [
                [
                    -0.9795554098249627,
                    -0.19121626094134914,
                    -0.06251032422135915,
                ],
                [
                    0.1999606555147074,
                    -0.8913653381783955,
                    -0.40679671844822346,
                ],
                [0.02206661117566709, -0.41097953166273926, 0.911377468026008],
            ],
        ),
    ];
    let ctx = ctx();
    for (jd_tt, expected_axes) in samples {
        let jd = JulianDate::new(jd_tt);
        for (axis, expected) in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            .into_iter()
            .zip(expected_axes)
        {
            let actual: Direction<ICRS> =
                Direction::<MoonIauMeanEarth>::new_unchecked(axis[0], axis[1], axis[2])
                    .to_frame_with(&jd, &ctx);
            for (measured, reference) in [actual.x(), actual.y(), actual.z()]
                .into_iter()
                .zip(expected)
            {
                assert_near(measured, reference, 1e-10, "IAU_MOON SPICE axis");
            }
        }
    }
}

#[test]
fn negative_height_domain_for_spheres_and_ellipsoids() {
    use siderust::coordinates::surface::SurfaceError;

    let moon = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, -100.0 * M).unwrap();
    assert_near(
        moon.body_fixed_cartesian().unwrap()[0].value(),
        1_737_300.0,
        1e-6,
        "valid Moon subsurface",
    );
    assert!(matches!(
        SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, -1_737_400.0 * M),
        Err(SurfaceError::DegenerateGeometry { .. })
    ));
    assert!(SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, -1_737_401.0 * M).is_err());

    let ellipsoid = Mars::default_surface();
    let a = ellipsoid.equatorial_radius().value();
    let f = ellipsoid.flattening();
    let meridional_radius = a * (1.0 - (2.0 * f - f * f));
    let safe = SurfaceSite::<Mars>::new(0.0 * DEG, 0.0 * DEG, (-meridional_radius + 100.0) * M);
    assert!(safe.is_ok(), "valid shallow side of normal singularity");
    assert!(matches!(
        SurfaceSite::<Mars>::new(0.0 * DEG, 0.0 * DEG, (-meridional_radius - 100.0) * M),
        Err(SurfaceError::DegenerateGeometry { .. })
    ));
}

#[cfg(feature = "serde")]
#[test]
fn invalid_subsurface_site_deserialization_is_rejected() {
    let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let mut payload = serde_json::to_value(site).unwrap();
    payload["height"] = serde_json::to_value(-1_737_400.0 * M).unwrap();
    assert!(serde_json::from_value::<SurfaceSite<Moon>>(payload).is_err());
}

#[test]
fn microscopic_finite_range_preserves_altaz_in_au() {
    // affn spherical Cartesian->spherical historically treats radii smaller
    // than machine epsilon in the *current unit* as zero. A 10-µm LOS
    // expressed in AU is below that numeric threshold but not physically zero.
    let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let [x, y, z] = site.body_fixed_cartesian().unwrap();
    let target = Position::<Selenocentric, MoonIauMeanEarth, Meter>::new(x, y + 0.00001 * M, z);
    let target_au = target.to_unit::<AstronomicalUnit>();
    let jd = siderust::J2000;
    let ctx = ctx();
    let obs_m = position_to_horizontal(&target, &site, &jd, &ctx).unwrap();
    let obs_au = position_to_horizontal(&target_au, &site, &jd, &ctx).unwrap();
    assert_near(obs_m.az().value(), 90.0, 0.25, "microscopic metre azimuth");
    assert_near(obs_au.az().value(), 90.0, 0.25, "microscopic AU azimuth");
    assert_near(
        obs_m.alt().value(),
        obs_au.alt().value(),
        0.25,
        "micro alt parity",
    );
    assert_near(
        obs_m.az().value(),
        obs_au.az().value(),
        0.25,
        "micro az parity",
    );
    assert_near(
        obs_au.direction().to_spherical().az().value(),
        obs_au.az().value(),
        1e-6,
        "micro direction consistency",
    );
}

#[test]
fn finite_cardinals_and_extreme_ranges_match_analytical_geometry() {
    let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let jd = siderust::J2000;
    let ctx = ctx();
    for (los, altitude, azimuth) in [
        ([0.0, 0.0, 1.0], 0.0, 0.0),
        ([0.0, 1.0, 0.0], 0.0, 90.0),
        ([0.0, 0.0, -1.0], 0.0, 180.0),
        ([0.0, -1.0, 0.0], 0.0, 270.0),
        ([1.0, 0.0, 0.0], 90.0, 0.0),
        ([-1.0, 0.0, 0.0], -90.0, 0.0),
    ] {
        for range in [1_000.0, 1e160] {
            let target = Position::<Selenocentric, MoonIauMeanEarth, Meter>::new(
                1_737_400.0 + range * los[0],
                range * los[1],
                range * los[2],
            );
            let finite = position_to_horizontal(&target, &site, &jd, &ctx).unwrap();
            assert_near(finite.range().value() / range, 1.0, 1e-12, "scaled range");
            assert_near(finite.alt().value(), altitude, 1e-8, "finite altitude");
            if altitude.abs() < 90.0 {
                let delta = (finite.az().value() - azimuth + 180.0).rem_euclid(360.0) - 180.0;
                assert_near(delta, 0.0, 1e-8, "finite cardinal azimuth");
            }
        }
    }
}

#[test]
fn invalid_horizontal_directions_and_targets_return_errors() {
    use siderust::coordinates::frames::Horizontal;
    let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    for xyz in [[0.0; 3], [f64::NAN, 0.0, 0.0], [0.0, f64::INFINITY, 0.0]] {
        let horizontal = Direction::<Horizontal>::new_unchecked(xyz[0], xyz[1], xyz[2]);
        let inverse: Result<Direction<ICRS>, _> =
            horizontal_to_direction(&horizontal, &site, &siderust::J2000, &ctx());
        assert!(inverse.is_err());
        let target = Position::<Selenocentric, MoonIauMeanEarth, Meter>::new(
            xyz[0] + 1_737_400.0,
            xyz[1],
            xyz[2],
        );
        assert!(position_to_horizontal(&target, &site, &siderust::J2000, &ctx()).is_err());
    }
}

#[test]
fn ellipsoid_geometry_is_scale_independent() {
    use siderust::coordinates::surface::ReferenceSurface;
    // Squared dimensional axes used to overflow/underflow here, even though
    // each position and its unit normal are representable.
    for radius in [1e-160, 1.0, 1e160] {
        let surface = ReferenceSurface::ellipsoid(radius * M, 0.25).unwrap();
        let graphic_lat = 35.0_f64.to_radians();
        let centric_lat = ((1.0 - surface.e2()) * graphic_lat.tan()).atan();
        let graphic = SurfaceSite::<Mars>::with_surface(
            40.0 * DEG,
            graphic_lat.to_degrees() * DEG,
            radius * 0.1 * M,
            surface,
            LatitudeConvention::Planetographic,
        )
        .unwrap();
        let centric = SurfaceSite::<Mars>::with_surface(
            40.0 * DEG,
            centric_lat.to_degrees() * DEG,
            radius * 0.1 * M,
            surface,
            LatitudeConvention::Planetocentric,
        )
        .unwrap();
        for (g, c) in graphic
            .body_fixed_cartesian()
            .unwrap()
            .into_iter()
            .zip(centric.body_fixed_cartesian().unwrap())
        {
            assert_near(
                g.value() / radius,
                c.value() / radius,
                1e-14,
                "scaled ellipsoid",
            );
        }
        for (g, c) in graphic
            .enu_basis()
            .unwrap()
            .up
            .into_iter()
            .zip(centric.enu_basis().unwrap().up)
        {
            assert_near(g, c, 1e-14, "scale-free normal");
        }
    }
}

#[test]
fn legacy_earth_cartesian_horizontal_matches_spherical_east() {
    use siderust::coordinates::centers::{Geodetic, Topocentric};
    use siderust::coordinates::frames::{EquatorialMeanOfDate, Horizontal};
    use siderust::coordinates::transform::Transform;
    let jd = siderust::J2000;
    let site = Geodetic::<ECEF>::new(0.0 * DEG, 45.0 * DEG, 0.0 * M);
    let ra =
        siderust::astro::earth_rotation::gmst_from_tt(jd).value() + core::f64::consts::FRAC_PI_2;
    let target = Position::<Topocentric, EquatorialMeanOfDate, Meter>::from_array(
        site,
        [ra.cos() * M, ra.sin() * M, 0.0 * M],
    );
    let horizontal: Position<Topocentric, Horizontal, Meter> = target.transform(jd);
    assert_near(
        horizontal.to_spherical().az().value(),
        90.0,
        1e-10,
        "legacy east azimuth",
    );
    assert_near(horizontal.y().value(), 1.0, 1e-14, "legacy east component");
    let inverse: Position<Topocentric, EquatorialMeanOfDate, Meter> = horizontal.transform(jd);
    assert_near(
        inverse.x().value(),
        target.x().value(),
        1e-14,
        "legacy inverse x",
    );
    assert_near(
        inverse.y().value(),
        target.y().value(),
        1e-14,
        "legacy inverse y",
    );
}

#[test]
fn shared_planetary_rotation_matches_independent_axes_at_multiple_epochs() {
    use siderust::coordinates::frames::*;
    // Explicit scalar IAU basis: node=(-sin RA, cos RA, 0),
    // equatorial tangent=(-sin DEC cos RA,-sin DEC sin RA,cos DEC).
    // This checks every user of the corrected shared ZXZ provider. Existing
    // planetary models retain their linear elements (including omissions of
    // periodic terms); this test validates the convention, not full PCK accuracy.
    macro_rules! check {
        ($frame:ty, $ra:expr, $rar:expr, $dec:expr, $decr:expr, $w:expr, $wr:expr) => {
            for days in [-36525.0, 0.0, 18262.5] {
                let ra = ($ra + $rar * days / 36525.0_f64).to_radians();
                let dec = ($dec + $decr * days / 36525.0_f64).to_radians();
                let w = (($w + $wr * days).rem_euclid(360.0_f64)).to_radians();
                let expected = [
                    [
                        -ra.sin() * w.cos() - dec.sin() * ra.cos() * w.sin(),
                        ra.cos() * w.cos() - dec.sin() * ra.sin() * w.sin(),
                        dec.cos() * w.sin(),
                    ],
                    [
                        ra.sin() * w.sin() - dec.sin() * ra.cos() * w.cos(),
                        -ra.cos() * w.sin() - dec.sin() * ra.sin() * w.cos(),
                        dec.cos() * w.cos(),
                    ],
                    [dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin()],
                ];
                let jd = JulianDate::new(2451545.0 + days);
                for (axis, expected) in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
                    .into_iter()
                    .zip(expected)
                {
                    let actual: Direction<ICRS> =
                        Direction::<$frame>::new_unchecked(axis[0], axis[1], axis[2])
                            .to_frame_with(&jd, &ctx());
                    for (value, reference) in [actual.x(), actual.y(), actual.z()]
                        .into_iter()
                        .zip(expected)
                    {
                        assert_near(value, reference, 1e-10, stringify!($frame));
                    }
                }
            }
        };
    }
    check!(
        MercuryFixed,
        281.0103,
        -0.0328,
        61.4155,
        -0.0049,
        329.5988,
        6.1385108
    );
    check!(VenusFixed, 272.76, 0.0, 67.16, 0.0, 160.20, -1.4813688);
    check!(
        MarsFixed,
        317.269,
        -0.10927,
        54.432,
        -0.05827,
        176.049,
        350.891982443
    );
    check!(
        JupiterSystemIII,
        268.057,
        -0.006,
        64.495,
        0.002,
        284.95,
        870.536
    );
    check!(
        SaturnFixed,
        40.589,
        -0.036,
        83.537,
        -0.004,
        38.90,
        810.7939024
    );
    check!(
        UranusFixed,
        257.311,
        0.0,
        -15.175,
        0.0,
        203.81,
        -501.1600928
    );
    check!(NeptuneFixed, 299.36, 0.0, 43.46, 0.0, 249.978, 541.1397757);
    check!(PlutoFixed, 132.993, 0.0, -6.163, 0.0, 302.695, 56.3625225);
}

#[test]
fn requested_length_unit_must_represent_finite_geometry() {
    use siderust::coordinates::surface::SurfaceError;
    // The site and LOS are finite in metres, but not in millimetres.
    let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 1e307 * M).unwrap();
    assert!(site.observer_position::<Meter>().is_ok());
    assert!(matches!(
        site.observer_position::<Millimeter>(),
        Err(SurfaceError::NonFiniteInput { .. })
    ));
    let target = Position::<Selenocentric, MoonIauMeanEarth, Millimeter>::new(0.0, 0.0, 0.0);
    assert!(matches!(
        position_to_horizontal(&target, &site, &siderust::J2000, &ctx()),
        Err(SurfaceError::NonFiniteInput { .. })
    ));
    let target_m = target.to_unit::<Meter>();
    let observation = position_to_horizontal(&target_m, &site, &siderust::J2000, &ctx()).unwrap();
    assert_near(
        observation.range().value() / 1e307,
        1.0,
        1e-14,
        "finite metre geometry",
    );
    let tiny = SurfaceSite::<Moon>::with_surface(
        0.0 * DEG,
        0.0 * DEG,
        0.0 * M,
        siderust::coordinates::surface::ReferenceSurface::sphere(1e-320 * M).unwrap(),
        LatitudeConvention::Planetocentric,
    )
    .unwrap();
    assert!(tiny.observer_position::<Meter>().is_ok());
    assert!(matches!(
        tiny.observer_position::<AstronomicalUnit>(),
        Err(SurfaceError::DegenerateGeometry { .. })
    ));
}

// =============================================================================
// End-to-end geometric Alt/Az against offline CSPICE references
// =============================================================================
//
// Provenance (regenerate with tools/cspice_surface_references.py):
//   CSPICE N0067 via SpiceyPy 8.3.0
//   Kernels: pck00011.tpc, naif0012.tls (NAIF generic kernels)
//   Epoch: JD TT 2451545.0 (J2000.0); ET = unitim((jd_tt-2451545)*86400,"TDT","TDB")
//   Lunar frame: IAU_MOON (analytical text-PCK Mean-Earth; same model as siderust)
//   Surface: IAU mean sphere R = 1737.4 km; lon east-positive, planetocentric lat
//   Target: body-fixed point offset +50 km east, +10 km north, +5 km up from site
//   Procedure: pxform("IAU_MOON","J2000",et) rotates body-fixed → J2000/ICRS;
//              Alt/Az/range from independent ENU projection of the body-fixed LOS
//   Tolerances: 1e-9 deg / 1e-9 km — composition of CSPICE orientation (already
//              checked to 1e-10 direction cosine) with exact sphere ENU geometry.
//
// Mars case below is body-fixed only: siderust's linear IAU_MARS model differs
// from pck00011 periodic terms by ~0.02 in direction cosine at J2000, so an
// ICRS-mediated Mars comparison would mix inequivalent orientation models.

#[test]
fn lunar_e2e_altaz_matches_offline_cspice_iau_moon_reference() {
    // Site: lon=20°E, lat=−3.5°, h=0 on IAU sphere.
    let site = SurfaceSite::<Moon>::new(20.0 * DEG, (-3.5) * DEG, 0.0 * M).unwrap();
    // CSPICE-derived body-fixed target [km], then converted to metres:
    // obs + 50*east + 10*north + 5*up on R=1737.4 km sphere.
    let target = Position::<Selenocentric, MoonIauMeanEarth, Kilometer>::new(
        1617.7391512763131,
        642.0177864954426,
        -96.38962730131595,
    );
    let jd = JulianDate::new(2_451_545.0);
    let ctx = ctx();

    // Path A: observe in the body-fixed frame (orientation cancelled).
    let obs_fixed = position_to_horizontal(&target, &site, &jd, &ctx).unwrap();
    assert_near(
        obs_fixed.alt().value(),
        5.60040918481664,
        1e-9,
        "CSPICE alt",
    );
    assert_near(obs_fixed.az().value(), 78.69006752597979, 1e-9, "CSPICE az");
    assert_near(
        obs_fixed.range().to::<Kilometer>().value(),
        51.23475382979799,
        1e-9,
        "CSPICE slant range",
    );

    // Path B: same target expressed in ICRS using the CSPICE pxform matrix
    // (columns = body axes in J2000). Exercises TT→TDB, IAU_MOON orientation,
    // observer position, LOS, ENU, Alt/Az, and slant range together.
    // M = pxform("IAU_MOON","J2000", et) at JD TT 2451545.0:
    let m = [
        [
            0.7842270522119748,
            -0.6200619150990117,
            -0.022608671404191448,
        ],
        [0.5578471123206701, 0.720556665574856, -0.4118309009425287],
        [0.2716514860154575, 0.31035675139972596, 0.9109797785934671],
    ];
    let tf = [1617.7391512763131, 642.0177864954426, -96.38962730131595];
    let target_icrs = Position::<Selenocentric, ICRS, Kilometer>::new(
        m[0][0] * tf[0] + m[0][1] * tf[1] + m[0][2] * tf[2],
        m[1][0] * tf[0] + m[1][1] * tf[1] + m[1][2] * tf[2],
        m[2][0] * tf[0] + m[2][1] * tf[1] + m[2][2] * tf[2],
    );
    let obs_icrs = position_to_horizontal(&target_icrs, &site, &jd, &ctx).unwrap();
    // ICRS mediation recomposes the CSPICE matrix through siderust's
    // TT→TDB + IAU_MOON provider. Orientation alone is held to 1e-10 in
    // direction cosine (~2e-8 deg); Alt/Az inherit that floor.
    assert_near(
        obs_icrs.alt().value(),
        5.60040918481664,
        5e-8,
        "ICRS path alt",
    );
    assert_near(
        obs_icrs.az().value(),
        78.69006752597979,
        5e-8,
        "ICRS path az",
    );
    assert_near(
        obs_icrs.range().to::<Kilometer>().value(),
        51.23475382979799,
        1e-6, // ~1 mm; matrix round-trip through ICRS
        "ICRS path range",
    );
}

#[test]
fn mars_e2e_body_fixed_cardinal_and_spice_longitude_interop() {
    // SPICE pgrrec (west-positive 137.4°, lat −4.6°, a=3396.19, b=3376.20 km)
    // agrees with Siderust when longitude is converted to east-positive −137.4°.
    let site =
        SurfaceSite::<Mars>::from_planetographic_longitude(137.4 * DEG, (-4.6) * DEG, 0.0 * M)
            .unwrap();
    let [x, y, z] = site.body_fixed_cartesian().unwrap();
    assert_near(x.value() / 1000.0, -2491.967067275659, 1e-9, "Mars x");
    assert_near(y.value() / 1000.0, -2291.481184471932, 1e-9, "Mars y");
    assert_near(z.value() / 1000.0, -269.1840291055111, 1e-9, "Mars z");

    // 1 km along local east in the body-fixed frame → az=90°, alt≈0, range=1 km.
    let enu = site.enu_basis().unwrap();
    let target = Position::<Marscentric, MarsFixed, Meter>::new(
        Quantity::<Meter>::new(x.value() + enu.east[0] * 1_000.0),
        Quantity::<Meter>::new(y.value() + enu.east[1] * 1_000.0),
        Quantity::<Meter>::new(z.value() + enu.east[2] * 1_000.0),
    );
    let obs = position_to_horizontal(&target, &site, &siderust::J2000, &ctx()).unwrap();
    assert_near(obs.alt().value(), 0.0, 1e-8, "Mars east alt");
    assert_near(obs.az().value(), 90.0, 1e-8, "Mars east az");
    assert_near(obs.range().value(), 1_000.0, 1e-6, "Mars east range");
}
