// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Regression tests for high-precision Earth-rotation usage in topocentric paths.
//!
//! Reference numbers in this file were precomputed with ERFA/SOFA routines:
//! - sidereal angle: GMST06 + Nut00B equation-of-equinoxes form
//! - terrestrial->celestial chain: W(xp,yp,s') · R3(-ERA) · Q(X,Y,s)
//!   with frame-bias to EquatorialMeanJ2000

use qtty::*;
use siderust::astro::eop::NullEop;
use siderust::calculus::horizontal::equatorial_to_horizontal_true_of_date_with_ctx;
use siderust::coordinates::cartesian::position;
use siderust::coordinates::centers::{Geocentric, Geodetic, Topocentric};
use siderust::coordinates::frames::{self, ECEF};
use siderust::coordinates::spherical;
use siderust::coordinates::transform::centers::position::to_topocentric::ToTopocentricExt;
use siderust::coordinates::transform::context::{
    AstroContext, DefaultEphemeris, DefaultNutationModel,
};
use siderust::time::JulianDate;

#[inline]
fn wrapped_arcsec_error(actual_deg: f64, expected_deg: f64) -> f64 {
    let d = (actual_deg - expected_deg).abs();
    d.min(360.0 - d) * 3600.0
}

fn tod_unit_position(
    site: Geodetic::<ECEF>,
    ra_deg: f64,
    dec_deg: f64,
) -> spherical::Position<Topocentric, frames::EquatorialTrueOfDate, AstronomicalUnit> {
    affn::spherical::Position::<Topocentric, frames::EquatorialTrueOfDate, AstronomicalUnit>::new_raw_with_params(
        site,
        Degrees::new(dec_deg),
        Degrees::new(ra_deg),
        AstronomicalUnits::new(1.0),
    )
}

#[test]
fn horizontal_true_of_date_matches_erfa_roque_sirius_2020() {
    let jd_tt = JulianDate::new(2_459_015.5);
    let site = Geodetic::<ECEF>::new(-17.8925 * DEG, 28.7543 * DEG, 2396.0 * M);
    let eq = tod_unit_position(site, 101.287, -16.716);
    let ctx: AstroContext = AstroContext::default();

    let horiz = equatorial_to_horizontal_true_of_date_with_ctx(&eq, site, jd_tt, &ctx);

    // ERFA reference (degrees)
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
fn horizontal_true_of_date_matches_erfa_greenwich_2024() {
    let jd_tt = JulianDate::new(2_460_310.25);
    let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
    let eq = tod_unit_position(site, 210.1234, 35.6789);
    let ctx: AstroContext = AstroContext::default();

    let horiz = equatorial_to_horizontal_true_of_date_with_ctx(&eq, site, jd_tt, &ctx);

    // ERFA reference (degrees)
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

#[test]
fn topocentric_site_vector_matches_erfa_chain_roque_2020() {
    let jd_tt = JulianDate::new(2_459_015.5);
    let site = Geodetic::<ECEF>::new(-17.8925 * DEG, 28.7543 * DEG, 2396.0 * M);
    let origin = position::EquatorialMeanJ2000::<Kilometer, Geocentric>::new(0.0, 0.0, 0.0);

    let topo_default = origin.to_topocentric(site, jd_tt);

    // If target is at geocenter, topocentric vector is -r_site.
    let site_eq_x = -topo_default.x();
    let site_eq_y = -topo_default.y();
    let site_eq_z = -topo_default.z();

    // ERFA reference (km) for this epoch/site using full IAU chain.
    let ex = 1_081.752_851_23;
    let ey = 5_493.726_969_09;
    let ez = 3_049.140_205_07;

    assert!(
        (site_eq_x.value() - ex).abs() < 0.05,
        "x mismatch: {} vs {}",
        site_eq_x,
        ex
    );
    assert!(
        (site_eq_y.value() - ey).abs() < 0.05,
        "y mismatch: {} vs {}",
        site_eq_y,
        ey
    );
    assert!(
        (site_eq_z.value() - ez).abs() < 0.05,
        "z mismatch: {} vs {}",
        site_eq_z,
        ez
    );

    // Null-EOP path remains selectable and should differ measurably.
    let null_ctx: AstroContext<DefaultEphemeris, NullEop, DefaultNutationModel> =
        AstroContext::default();
    let topo_null = origin.to_topocentric_with_ctx(site, jd_tt, &null_ctx);
    let dx = topo_default.x() - topo_null.x();
    let dy = topo_default.y() - topo_null.y();
    let dz = topo_default.z() - topo_null.z();
    let delta = (dx * dx + dy * dy + dz * dz).sqrt();
    assert!(
        delta > 0.001,
        "default vs NullEop delta too small: {}",
        delta
    );
}
