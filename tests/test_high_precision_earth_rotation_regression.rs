// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Regression tests for high-precision Earth-rotation usage in topocentric paths.
//!
//! Reference numbers in this file were precomputed with SOFA routines:
//! - sidereal angle: GMST06 + Nut00B equation-of-equinoxes form
//! - terrestrial to celestial chain: transpose of ERFA c2t06a, followed by
//!   the J2000 frame bias from bp06

use siderust::astro::eop::{EopProvider, EopValues, NullEop};
#[cfg(feature = "full-nutation")]
use siderust::astro::nutation::Iau2006A;
use siderust::astro::nutation::{Iau2000B, NutationModel};
use siderust::coordinates::cartesian::position;
use siderust::coordinates::centers::{Geocentric, Geodetic, Topocentric};
use siderust::coordinates::frames::{self, ECEF};
use siderust::coordinates::spherical;
use siderust::coordinates::transform::context::{AstroContext, DefaultEphemeris};
use siderust::coordinates::transform::{
    horizontal::{
        project_direction_with, project_position_with, unproject_direction_with,
        unproject_position_with,
    },
    PositionAstroExt, TransformCenter,
};
use siderust::qtty::*;
use siderust::time::JulianDate;

#[inline]
fn wrapped_arcsec_error(actual_deg: f64, expected_deg: f64) -> f64 {
    let d = (actual_deg - expected_deg).abs();
    d.min(360.0 - d) * 3600.0
}

fn tod_unit_position(
    site: Geodetic<ECEF>,
    ra_deg: f64,
    dec_deg: f64,
) -> spherical::Position<Topocentric, frames::EquatorialTrueOfDate, AstronomicalUnit> {
    affn::spherical::Position::<Topocentric, frames::EquatorialTrueOfDate, AstronomicalUnit>::new_unchecked_with_params(
        site,
        Degrees::new(dec_deg),
        Degrees::new(ra_deg),
        AstronomicalUnits::new(1.0),
    )
}

#[test]
fn horizontal_true_of_date_matches_sofa_roque_sirius_2020() {
    let jd_tt = JulianDate::try_new(Days::new(2_459_015.5)).unwrap();
    let site = Geodetic::<ECEF>::new(-17.8925 * DEG, 28.7543 * DEG, 2396.0 * M);
    let eq = tod_unit_position(site, 101.287, -16.716);
    let ctx: AstroContext = AstroContext::default();

    let horiz = project_position_with(&eq, jd_tt, &ctx);

    // SOFA reference (degrees)
    let expected_alt = -55.077_584_745_179_26;
    let expected_az = 282.286_824_123_029_6;

    let alt_err = wrapped_arcsec_error(horiz.alt().value(), expected_alt);
    let az_err = wrapped_arcsec_error(horiz.az().value(), expected_az);

    assert!(
        alt_err < 0.5,
        "altitude error {}\" exceeds tolerance",
        alt_err
    );
    assert!(az_err < 0.5, "azimuth error {}\" exceeds tolerance", az_err);
}

#[test]
fn horizontal_true_of_date_matches_sofa_greenwich_2024() {
    let jd_tt = JulianDate::try_new(Days::new(2_460_310.25)).unwrap();
    let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
    let eq = tod_unit_position(site, 210.1234, 35.6789);
    let ctx: AstroContext = AstroContext::default();

    let horiz = project_position_with(&eq, jd_tt, &ctx);

    // SOFA reference (degrees)
    let expected_alt = -1.006_037_970_004_032_3;
    let expected_az = 343.464_018_807_157_86;

    let alt_err = wrapped_arcsec_error(horiz.alt().value(), expected_alt);
    let az_err = wrapped_arcsec_error(horiz.az().value(), expected_az);

    assert!(
        alt_err < 0.5,
        "altitude error {}\" exceeds tolerance",
        alt_err
    );
    assert!(az_err < 0.5, "azimuth error {}\" exceeds tolerance", az_err);
}

/// Compare a mean-equator J2000 site vector against ERFA 2.0.1 (km).
/// Reproduce with PyERFA using TT split as (jd, 0), UT1 split as
/// (jd, (DUT1 - 69.184) / 86400) and xp/yp in radians:
/// `bp06(2451545, 0)[0] @ c2t06a(jd, 0, *ut1, xp, yp).T @ gd2gc(1, lon, lat, height) / 1000`.
/// Pass WGS84 east-positive longitude and geodetic latitude in radians, and
/// ellipsoidal height in meters.
/// TT–UTC is 69.184 s at both epochs, with zero celestial pole offsets.
fn assert_site_vector<Eop: EopProvider + Default>(
    jd: f64,
    site: Geodetic<ECEF>,
    expected_km: [f64; 3],
) {
    let jd_tt = JulianDate::try_new(Days::new(jd)).unwrap();
    let origin = position::EquatorialMeanJ2000::<Kilometer, Geocentric>::new(0.0, 0.0, 0.0);
    let ctx = AstroContext::<DefaultEphemeris, Eop>::default();
    let topo: position::EquatorialMeanJ2000<Kilometer, Topocentric> =
        origin.to_center_with(site, jd_tt, &ctx);

    // Target at geocenter yields -r_site. Compare the physical 3D error.
    let actual_km = topo.as_array().map(|axis| -axis.value());
    let error_m = actual_km
        .into_iter()
        .zip(expected_km)
        .map(|(actual, expected)| ((actual - expected) * 1000.0).powi(2))
        .sum::<f64>()
        .sqrt();
    // A 0.5 m displacement at Earth radius allows ~16 mas of rotation.
    // This covers the mas-scale 2000B pole error, omitted periodic CIO terms
    // (the leading s06 term is 2.64 mas), and <2 cm from one ~40 µs UT1 JD ulp.
    // Both pipelines receive identical EOP, including for the NullEop case.
    let tolerance_m = 0.5;
    assert!(
        error_m < tolerance_m,
        "site {site:?} at JD(TT) {jd}: mean J2000 vector error {error_m:.6} m exceeds {tolerance_m} m (actual {actual_km:?} km, ERFA {expected_km:?} km)"
    );
}

#[test]
fn topocentric_site_vector_matches_sofa_chain_roque_2020() {
    assert_site_vector::<NullEop>(
        2_459_015.5,
        Geodetic::new(-17.8925 * DEG, 28.7543 * DEG, 2396.0 * M),
        [-2_333.676_983_87, -5_085.751_188_61, 3_055.743_414_41],
    );
}

#[derive(Default)]
struct ParanalEop;

impl EopProvider for ParanalEop {
    fn eop_at(&self, _: tempoch::JulianDate<tempoch::UTC>) -> EopValues {
        EopValues {
            dut1: Seconds::new(0.4),
            xp: Arcseconds::new(0.8),
            yp: Arcseconds::new(-0.6),
            ..Default::default()
        }
    }
}

#[test]
fn topocentric_site_vector_matches_sofa_chain_paranal_with_eop() {
    // Synthetic EOP activate the UT1 and both terrestrial pole corrections.
    // Southern latitude and a different ERA probe the rotation composition.
    assert_site_vector::<ParanalEop>(
        2_460_310.25,
        Geodetic::new(-70.403 * DEG, -24.625 * DEG, 2635.0 * M),
        [2_799.501_427_33, -5_080.725_942_18, -2_648.816_162_52],
    );
}

#[test]
fn topocentric_composition_preserves_site_range_and_context() {
    let jd = JulianDate::new(2_460_310.25);
    let site = Geodetic::new(-70.403 * DEG, -24.625 * DEG, 2635.0 * M);
    let other_site = Geodetic::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
    let target =
        position::EquatorialMeanJ2000::<Kilometer, Geocentric>::new(384_400.0, 52_000.0, -23_000.0);
    let context = AstroContext::<DefaultEphemeris, OffsetEop>::default();
    let model = context.with_model::<Iau2000B>();
    let relative: position::EquatorialMeanJ2000<Kilometer, Topocentric> =
        target.to_center_with(site, jd, &model);
    let equatorial = relative
        .to_frame_with::<frames::EquatorialTrueOfDate, _>(&jd, &model)
        .to_spherical();
    let icrs = relative.to_frame_with::<frames::ICRS, _>(&jd, &model);
    let via_icrs = icrs.to_frame_with::<frames::EquatorialTrueOfDate, _>(&jd, &model);
    assert!(
        equatorial
            .direction()
            .angular_separation(&via_icrs.to_spherical().direction())
            < 1.0e-10 * DEG
    );

    let horizon = project_position_with(&equatorial, jd, &model);
    assert_eq!(horizon.center_params(), &site);
    assert_eq!(horizon.distance, equatorial.distance);
    let direction = project_direction_with(&equatorial.direction(), &site, jd, &model);
    assert_eq!(
        (horizon.alt(), horizon.az()),
        (direction.alt(), direction.az())
    );

    let recovered = unproject_position_with(&horizon, jd, &model);
    assert_eq!(recovered.center_params(), &site);
    assert_eq!(recovered.distance, equatorial.distance);
    assert!(
        recovered
            .direction()
            .angular_separation(&equatorial.direction())
            < 1.0e-10 * DEG
    );

    let recentered: position::EquatorialMeanJ2000<Kilometer, Topocentric> =
        target.to_center_with(other_site, jd, &model);
    let other = recentered.to_frame_with::<frames::EquatorialTrueOfDate, _>(&jd, &model);
    let other_horizon = project_position_with(&other.to_spherical(), jd, &model);
    assert_eq!(other_horizon.center_params(), &other_site);
    let reprojected = project_direction_with(&equatorial.direction(), &other_site, jd, &model);
    assert!(other_horizon.direction().angular_separation(&reprojected) > 0.01 * DEG);
}

#[derive(Default)]
struct OffsetEop;

impl EopProvider for OffsetEop {
    fn eop_at(&self, utc: tempoch::JulianDate<tempoch::UTC>) -> EopValues {
        EopValues {
            dx: MilliArcseconds::new(300.0),
            dy: MilliArcseconds::new(-200.0),
            ..ParanalEop.eop_at(utc)
        }
    }
}

// ERFA 2.0.1 reference recipe:
// B = Pnm06a, or Fw2m(Pfw06 angles with Nut00b added to psib and epsa).
// Take X,Y from B's last row. Q = C2ixys(X+dX, Y+dY, S06(TT,X+dX,Y+dY)).
// W = Pom00(xp,yp,Sp00(TT)). T = C2tcio(Q,Era00(UT1),W).
// Project T * transpose(B) * direction onto the WGS84 local north/east/up basis.
// TT split is (jd,0) and UT1 is (jd,(0.4-69.184)/86400).
// With zero dX,dY, C2teqx(B, ERA-Eors(B,S06), W) agrees with T to 1e-15.
fn assert_erfa_horizon<Nut: NutationModel>(
    jd: JulianDate,
    site: Geodetic<ECEF>,
    equatorial: spherical::Direction<frames::EquatorialTrueOfDate>,
    nominal: spherical::Direction<frames::Horizontal>,
    offset: spherical::Direction<frames::Horizontal>,
) {
    let context = AstroContext::<DefaultEphemeris, ParanalEop>::default();
    let model = context.with_model::<Nut>();
    let actual_nominal = project_direction_with(&equatorial, &site, jd, &model);

    let context = AstroContext::<DefaultEphemeris, OffsetEop>::default();
    let model = context.with_model::<Nut>();
    let actual_offset = project_direction_with(&equatorial, &site, jd, &model);
    let recovered = unproject_direction_with(&offset, &site, jd, &model);

    // Allow 12 mas for the truncated s06 series and rounding UT1 to a
    // single JD. Synthetic terrestrial and celestial pole offsets exceed it.
    let tolerance = Arcseconds::new(0.012).to::<Degree>();
    assert!(
        actual_nominal.angular_separation(&nominal) < tolerance,
        "{actual_nominal:?}"
    );
    assert!(
        actual_offset.angular_separation(&offset) < tolerance,
        "{actual_offset:?}"
    );
    assert!(actual_nominal.angular_separation(&actual_offset) > tolerance);
    assert!(
        recovered.angular_separation(&equatorial) < tolerance,
        "{recovered:?}"
    );
}

#[test]
fn horizontal_with_eop_matches_erfa_2000b() {
    assert_erfa_horizon::<Iau2000B>(
        JulianDate::new(2_459_015.5),
        Geodetic::new(-17.8925 * DEG, 28.7543 * DEG, 2396.0 * M),
        spherical::Direction::<frames::EquatorialTrueOfDate>::new(101.287 * DEG, -16.716 * DEG),
        spherical::Direction::<frames::Horizontal>::new(
            -55.079_775_062_185_17 * DEG,
            282.289_250_395_418 * DEG,
        ),
        spherical::Direction::<frames::Horizontal>::new(
            -55.079_824_889_797_85 * DEG,
            282.289_156_280_150_9 * DEG,
        ),
    );
}

#[cfg(feature = "full-nutation")]
#[test]
fn horizontal_with_eop_matches_erfa_2006a() {
    assert_erfa_horizon::<Iau2006A>(
        JulianDate::new(2_460_310.25),
        Geodetic::new(-70.403 * DEG, -24.625 * DEG, 2635.0 * M),
        spherical::Direction::<frames::EquatorialTrueOfDate>::new(210.1234 * DEG, 35.6789 * DEG),
        spherical::Direction::<frames::Horizontal>::new(
            -13.374_342_598_225_061 * DEG,
            303.401_387_018_247_36 * DEG,
        ),
        spherical::Direction::<frames::Horizontal>::new(
            -13.374_375_514_929_3 * DEG,
            303.401_324_952_877_9 * DEG,
        ),
    );
}

fn assert_round_trips<Nut: NutationModel>() {
    let context = AstroContext::<DefaultEphemeris, OffsetEop>::default();
    let model = context.with_model::<Nut>();
    let tolerance = 1.0e-10 * DEG;
    for jd in [JulianDate::new(2_459_015.5), JulianDate::new(2_460_310.25)] {
        for (lon, lat) in [(-70.403, -24.625), (45.0, 90.0), (-120.0, -90.0)] {
            let site = Geodetic::new(lon * DEG, lat * DEG, 0.0 * M);
            for (ra, dec) in [
                (0.0, 90.0),
                (123.0, -90.0),
                (359.999, 89.999_999),
                (75.0, 0.0),
            ] {
                let eq =
                    spherical::Direction::<frames::EquatorialTrueOfDate>::new(ra * DEG, dec * DEG);
                let hz = project_direction_with(&eq, &site, jd, &model);
                let back = unproject_direction_with(&hz, &site, jd, &model);
                assert!(back.angular_separation(&eq) < tolerance);
            }
            for (az, alt) in [
                (127.0, 90.0),
                (223.0, -90.0),
                (359.999, 0.000_001),
                (0.001, -0.000_001),
            ] {
                let hz = spherical::Direction::<frames::Horizontal>::new(alt * DEG, az * DEG);
                let eq = unproject_direction_with(&hz, &site, jd, &model);
                let back = project_direction_with(&eq, &site, jd, &model);
                assert!(back.angular_separation(&hz) < tolerance);
                if alt == 90.0 || alt == -90.0 {
                    assert_eq!(back.az(), 0.0 * DEG);
                }
            }
        }
    }
}

#[test]
fn horizontal_round_trips_at_singularities() {
    assert_round_trips::<Iau2000B>();
    #[cfg(feature = "full-nutation")]
    assert_round_trips::<Iau2006A>();
}
