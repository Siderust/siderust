// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Regression tests for high-precision Earth-rotation usage in topocentric paths.
//!
//! Reference numbers in this file were precomputed with SOFA routines:
//! - sidereal angle: GMST06 + Nut00B equation-of-equinoxes form
//! - terrestrial to celestial chain: transpose of ERFA c2t06a, followed by
//!   the J2000 frame bias from bp06

use siderust::astro::eop::{EopProvider, EopValues, NullEop};
use siderust::coordinates::cartesian::position;
use siderust::coordinates::centers::{Geocentric, Geodetic, Topocentric};
use siderust::coordinates::frames::{self, ECEF};
use siderust::coordinates::spherical;
use siderust::coordinates::transform::context::{AstroContext, DefaultEphemeris};
use siderust::coordinates::transform::TransformCenter;
use siderust::event::horizontal::equatorial_to_horizontal_true_of_date_with_ctx;
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

    let horiz = equatorial_to_horizontal_true_of_date_with_ctx(&eq, site, jd_tt, &ctx);

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

    let horiz = equatorial_to_horizontal_true_of_date_with_ctx(&eq, site, jd_tt, &ctx);

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
