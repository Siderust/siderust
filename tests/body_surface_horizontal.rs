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
//! body-fixed frame realisation must match. Siderust's `MoonPrincipalAxes`
//! is the IAU 2015 analytical `(α₀, δ₀, W)` model, which is **not** identical
//! to JPL DE `MOON_PA` or the Mean-Earth (`MOON_ME`) cartographic frame.
//!
//! Geometric scope only: no refraction, aberration, light-time, or terrain.

use siderust::coordinates::cartesian::{Direction, Position};
use siderust::coordinates::centers::{Geocentric, Marscentric, Selenocentric};
use siderust::coordinates::frames::{MarsFixed, MoonPrincipalAxes, ECEF, ICRS};
use siderust::coordinates::surface::{
    BodySurfaceModel, EarthSurface, EnuBasis, LatitudeConvention, MarsSurface, MoonSurface,
    SurfaceSite, EARTH_WGS84_SURFACE,
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
        let moon = SurfaceSite::<MoonSurface>::new(lon * DEG, lat * DEG, 0.0 * M).unwrap();
        moon.enu_basis().unwrap().assert_orthonormal(1e-12).unwrap();

        let mars = SurfaceSite::<MarsSurface>::new(lon * DEG, lat * DEG, 0.0 * M).unwrap();
        mars.enu_basis().unwrap().assert_orthonormal(1e-12).unwrap();

        let earth = SurfaceSite::<EarthSurface>::new(lon * DEG, lat * DEG, 0.0 * M).unwrap();
        earth
            .enu_basis()
            .unwrap()
            .assert_orthonormal(1e-12)
            .unwrap();
    }
}

#[test]
fn north_pole_azimuth_follows_longitude_meridian() {
    let site = SurfaceSite::<MoonSurface>::new(30.0 * DEG, 90.0 * DEG, 0.0 * M).unwrap();
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
    let site = SurfaceSite::<MoonSurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
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
        let fixed = Direction::<MoonPrincipalAxes>::new_unchecked(v[0], v[1], v[2]);
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
// Finite-distance parallax
// =============================================================================

#[test]
fn moon_finite_parallax_shifts_altitude() {
    let site = SurfaceSite::<MoonSurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let jd = siderust::J2000;
    let ctx = ctx();

    // Same sky direction (+Z body-fixed): direction-only vs finite near target.
    let dir_fixed = Direction::<MoonPrincipalAxes>::new_unchecked(0.0, 0.0, 1.0);
    let dir_icrs: Direction<ICRS> = dir_fixed.to_frame_with(&jd, &ctx);
    let dir_alt = direction_to_horizontal(&dir_icrs, &site, &jd, &ctx)
        .unwrap()
        .to_spherical()
        .alt()
        .value();

    let near = Position::<Selenocentric, MoonPrincipalAxes, Kilometer>::new(0.0, 0.0, 2_000.0);
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
    let site = SurfaceSite::<MoonSurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let jd = siderust::J2000;
    let ctx = ctx();
    let fixed = Direction::<MoonPrincipalAxes>::new_unchecked(0.0, 0.0, 1.0);
    let icrs: Direction<ICRS> = fixed.to_frame_with(&jd, &ctx);
    let horiz = site.observe_direction(&icrs, &jd, &ctx).unwrap();
    assert_near(horiz.to_spherical().alt().value(), 0.0, 1e-8, "dir alt");
    assert_eq!(horiz.site_params(), &site.to_params());
    let other = SurfaceSite::<MoonSurface>::new(15.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
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
            &horiz, &site, &jd, &ctx
        ).unwrap();
    assert!(icrs.angle_to(&recovered).value() < 1e-8);
}

// =============================================================================
// Forward / inverse round-trips
// =============================================================================

#[test]
fn moon_direction_roundtrip_random_sky() {
    let site = SurfaceSite::<MoonSurface>::new((-15.0) * DEG, 42.0 * DEG, 100.0 * M).unwrap();
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
        let ang = icrs.angle_to(&back).value();
        // Orientation round-trips through ICRS at ~1e-8 rad numerical noise.
        assert!(ang < 1e-7, "roundtrip angle {ang} rad");
    }
}

#[test]
fn mars_direction_roundtrip() {
    let site = SurfaceSite::<MarsSurface>::new(75.0 * DEG, (-15.0) * DEG, 0.0 * M).unwrap();
    let jd = siderust::J2000;
    let ctx = ctx();
    let n = (0.2_f64.powi(2) + 0.5_f64.powi(2) + 0.84_f64.powi(2)).sqrt();
    let icrs = Direction::<ICRS>::new_unchecked(0.2 / n, -0.5 / n, 0.84 / n);
    let horiz = direction_to_horizontal(&icrs, &site, &jd, &ctx).unwrap();
    let back: Direction<ICRS> = horizontal_to_direction(&horiz, &site, &jd, &ctx).unwrap();
    assert!(icrs.angle_to(&back).value() < 1e-10);
}

// =============================================================================
// Earth regression vs ENU geometry (WGS84 site)
// =============================================================================

#[test]
fn earth_wgs84_zenith_matches_geodetic_cartesian() {
    use siderust::coordinates::centers::Geodetic;

    let geo = Geodetic::<ECEF>::new(0.0 * DEG, 45.0 * DEG, 0.0 * M);
    let site = SurfaceSite::<EarthSurface>::from_geodetic(geo).unwrap();
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
    assert_eq!(EarthSurface::default_surface(), EARTH_WGS84_SURFACE);
    assert_eq!(
        EarthSurface::default_latitude_convention(),
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
    let surf = MarsSurface::default_surface();

    let site_g = SurfaceSite::<MarsSurface>::with_surface(
        lon,
        lat,
        h,
        surf,
        LatitudeConvention::Planetographic,
    )
    .unwrap();
    let site_c = SurfaceSite::<MarsSurface>::with_surface(
        lon,
        lat,
        h,
        surf,
        LatitudeConvention::Planetocentric,
    )
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
    let site = SurfaceSite::<MarsSurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
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
    let err = SurfaceSite::<MoonSurface>::new(Degrees::new(f64::NAN), 0.0 * DEG, 0.0 * M);
    assert!(err.is_err());
}

#[test]
fn zero_range_is_error() {
    let site = SurfaceSite::<MoonSurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let [x, y, z] = site.body_fixed_cartesian().unwrap();
    let target = Position::<Selenocentric, MoonPrincipalAxes, Meter>::new(x, y, z);
    let err = position_to_horizontal(&target, &site, &siderust::J2000, &ctx());
    // Frame round-trip MoonPA→ICRS→MoonPA can leave a sub-micron residual;
    // treat either an explicit zero-range error or a sub-millimetre range as OK.
    match err {
        Err(siderust::coordinates::surface::SurfaceError::DegenerateGeometry { .. }) => {}
        Ok(obs) => assert!(
            obs.range().value() < 1e-3,
            "expected near-zero range, got {}",
            obs.range()
        ),
        Err(other) => panic!("unexpected error: {other}"),
    }
}


// =============================================================================
// Unit-independent coincident-target tolerance
// =============================================================================

#[test]
fn finite_observation_range_is_invariant_under_length_units() {
    let site = SurfaceSite::<MoonSurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let [x, y, z] = site.body_fixed_cartesian().unwrap();
    // 1 km above the surface: 1 km = ~6.7e-9 AU. Comparing the
    // raw numeric range to 1e-6 would incorrectly reject the AU target.
    let target_m = Position::<Selenocentric, MoonPrincipalAxes, Meter>::new(
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
    assert_near(obs_km.range().to::<Meter>().value(), 1_000.0, 1e-3, "km slant range");
    assert_near(obs_au.range().to::<Meter>().value(), 1_000.0, 1e-2, "AU slant range");
    assert_near(obs_m.alt().value(), obs_au.alt().value(), 1e-6, "unit invariant altitude");
}

// =============================================================================
// Independent analytical reference (Moon sphere ENU, no frame rotation)
// =============================================================================

#[test]
fn analytical_enu_projection_reference() {
    // Pure geometry check of los_to_horizontal via site ENU, independent of
    // orientation providers: observer at lunar equator, LOS = east+up.
    let site = SurfaceSite::<MoonSurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let EnuBasis { east, north: _, up } = site.enu_basis().unwrap();
    let los = [
        (east[0] + up[0]) / 2.0_f64.sqrt(),
        (east[1] + up[1]) / 2.0_f64.sqrt(),
        (east[2] + up[2]) / 2.0_f64.sqrt(),
    ];
    // Reconstruct via body-fixed direction → ICRS → horizontal at J2000.
    // Because we use the same ENU, alt should be 45°, az 90°.
    let n = (los[0] * los[0] + los[1] * los[1] + los[2] * los[2]).sqrt();
    let fixed = Direction::<MoonPrincipalAxes>::new_unchecked(los[0] / n, los[1] / n, los[2] / n);
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
    let site = SurfaceSite::<MarsSurface>::new(10.0 * DEG, 25.0 * DEG, 0.0 * M).unwrap();
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
    assert_near(normal.alt().value(), 90.0, 1e-6, "mars ellipsoid-normal zenith");
}

// Compile-time: SurfaceSite is parameterized by body (Moon ≠ Mars).
#[test]
fn site_params_body_name() {
    let moon = SurfaceSite::<MoonSurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    let mars = SurfaceSite::<MarsSurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
    assert_eq!(moon.body_name(), "Moon");
    assert_eq!(mars.body_name(), "Mars");
    // Params equality is within the same body centre type only.
    assert_eq!(moon.to_params(), moon.to_params());
    // The default params are body-specific, never implicitly lunar for Mars.
    let moon_default = siderust::coordinates::surface::SurfaceObserverParams::<MoonSurface>::default();
    let mars_default = siderust::coordinates::surface::SurfaceObserverParams::<MarsSurface>::default();
    assert_eq!(moon_default.surface, MoonSurface::default_surface());
    assert_eq!(mars_default.surface, MarsSurface::default_surface());
    assert_eq!(mars_default.latitude_convention, LatitudeConvention::Planetographic);
}
