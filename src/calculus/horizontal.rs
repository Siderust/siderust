// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Shared Horizontal Coordinate Helpers
//!
//! ## Scientific scope
//!
//! This module implements the apparent-topocentric pipeline that converts a
//! geocentric mean-J2000 position to an observed altitude / azimuth for an
//! Earth-surface observer.  The pipeline applies, in order:
//!
//! 1. **Topocentric parallax** — shifts the apparent direction from the
//!    geocentre to the observer site (important for the Moon, ~57′ maximum).
//! 2. **Precession** (IAU 2006 / P03) — rotates the mean J2000 frame to
//!    the mean equator / equinox of date.
//! 3. **Nutation** (IAU 2000B) — adds short-period oscillations of the
//!    Earth's rotation axis to reach the true-of-date frame.
//! 4. **Apparent sidereal time** (GAST, IAU 2006) — projects the equatorial
//!    position onto the local meridian to obtain the hour angle.
//! 5. **Equatorial → horizontal** — converts hour angle + declination to
//!    altitude + azimuth for the observer's latitude.
//!
//! The routines are shared between the Sun and Moon modules to avoid
//! duplicating identical transform chains.
//!
//! ## Technical scope
//!
//! - [`geocentric_j2000_to_apparent_topocentric`] — intermediate stage;
//!   takes a geocentric J2000 Cartesian position and returns a topocentric
//!   true-of-date spherical position (`RA`, `Dec`, distance).
//! - [`apparent_topocentric_to_horizontal`] — final stage; converts the
//!   true-of-date equatorial spherical position to altitude / azimuth using
//!   GAST and the observer latitude.
//! - [`geocentric_j2000_to_horizontal`] — convenience wrapper composing the
//!   full pipeline in one call.
//!
//! All public functions accept and return typed `qtty` quantities (degrees,
//! radians, AU, metres).
//!
//! ## References
//!
//! - Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed. ch. 13.
//!   Willmann-Bell.
//! - IAU SOFA (2023). *Standards of Fundamental Astronomy*, release 18.
//!   <https://www.iausofa.org/>
//! - Capitaine, N., et al. (2003). "Expressions for IAU 2000 precession
//!   quantities". *Astronomy and Astrophysics* 412, 567–586.
//!   <https://doi.org/10.1051/0004-6361:20031539>

use crate::astro::earth_rotation::jd_ut1_from_tt_eop;
use crate::astro::nutation::{nutation_iau2000b, NutationModel};
use crate::astro::precession;
use crate::astro::sidereal::gast_iau2006;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::coordinates::transform::centers::position::to_topocentric::to_topocentric_with;
use crate::coordinates::transform::{AstroContext, TransformContext};
use crate::coordinates::{cartesian, centers::*, frames, spherical};
use crate::qtty::{AstronomicalUnits, Degree, LengthUnit, Meter, Quantity, Radian, Radians};
use crate::time::JulianDate;

// =============================================================================
// Apparent Topocentric Pipeline
// =============================================================================

/// Transforms a **geocentric, mean-J2000 equatorial** Cartesian position into
/// **topocentric, true-of-date equatorial** spherical coordinates.
///
/// This pipeline is shared between Sun and Moon; the only prerequisite is
/// that the caller provides the body's geocentric J2000 Cartesian position
/// (obtained from VSOP87 for the Sun or ELP2000 for the Moon).
///
/// ## Pipeline Steps
///
/// 1. **Topocentric parallax**: Geocentric → Topocentric in J2000 frame
/// 2. **Precession**: J2000 → Mean-of-Date rotation
/// 3. **Nutation**: Mean-of-Date → True-of-Date rotation
/// 4. **Spherical conversion**: Cartesian → Spherical (RA, Dec, distance)
pub fn geocentric_j2000_to_apparent_topocentric<U: LengthUnit>(
    geo_cart_j2000: &cartesian::Position<Geocentric, frames::EquatorialMeanJ2000, U>,
    site: Geodetic<frames::ECEF>,
    jd: JulianDate,
) -> spherical::Position<Topocentric, frames::EquatorialTrueOfDate, U>
where
    Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
{
    let ctx: AstroContext = AstroContext::default();
    geocentric_j2000_to_apparent_topocentric_with_ctx(geo_cart_j2000, site, jd, &ctx)
}

/// Context-aware variant of [`geocentric_j2000_to_apparent_topocentric`].
pub fn geocentric_j2000_to_apparent_topocentric_with_ctx<U: LengthUnit, Ctx>(
    geo_cart_j2000: &cartesian::Position<Geocentric, frames::EquatorialMeanJ2000, U>,
    site: Geodetic<frames::ECEF>,
    jd: JulianDate,
    ctx: &Ctx,
) -> spherical::Position<Topocentric, frames::EquatorialTrueOfDate, U>
where
    Ctx: TransformContext,
    Ctx::Eph: crate::calculus::ephemeris::Ephemeris,
    Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
{
    // 1) Translate geocentric → topocentric in J2000 frame (applies parallax)
    let topo_cart_j2000 = to_topocentric_with(geo_cart_j2000, site, jd, ctx);

    // 2) Rotate J2000 → mean-of-date using IAU 2006 precession
    let rot_prec = precession::precession_matrix_iau2006(jd);
    let [x_m, y_m, z_m] = rot_prec
        * [
            topo_cart_j2000.x(),
            topo_cart_j2000.y(),
            topo_cart_j2000.z(),
        ];

    let topo_cart_mod =
        cartesian::Position::<Topocentric, frames::EquatorialMeanOfDate, U>::new_with_params(
            *topo_cart_j2000.center_params(),
            x_m,
            y_m,
            z_m,
        );

    // 3) Apply the selected nutation rotation (mean-of-date → true-of-date)
    let nut = Ctx::Nut::nutation(jd);
    let rot_nut = affn::Rotation3::fused_rx_rz_rx(
        nut.mean_obliquity + nut.deps,
        nut.dpsi,
        -nut.mean_obliquity,
    );
    let [x_t, y_t, z_t] = rot_nut * [topo_cart_mod.x(), topo_cart_mod.y(), topo_cart_mod.z()];

    let topo_cart_true =
        cartesian::Position::<Topocentric, frames::EquatorialTrueOfDate, U>::new_with_params(
            *topo_cart_mod.center_params(),
            x_t,
            y_t,
            z_t,
        );

    // 4) Convert to spherical topocentric equatorial (true of date)
    spherical::Position::from_cartesian(&topo_cart_true)
}

// =============================================================================
// Equatorial → Horizontal Transform
// =============================================================================

/// Converts an **apparent topocentric equatorial** position to **horizontal**
/// coordinates (altitude, azimuth, distance).
///
/// This is the standard equatorial-to-horizontal transform using the observer's
/// latitude and the local apparent sidereal time (GAST → LAST → HA).
///
/// Both Sun and Moon `get_horizontal` methods delegate here after computing
/// their body-specific apparent topocentric equatorial coordinates.
pub fn equatorial_to_horizontal<U: LengthUnit>(
    eq_position: &spherical::Position<Topocentric, frames::EquatorialTrueOfDate, U>,
    site: Geodetic<frames::ECEF>,
    jd: JulianDate,
) -> spherical::Position<Topocentric, frames::Horizontal, U> {
    equatorial_to_horizontal_true_of_date(eq_position, site, jd)
}

/// Converts apparent topocentric **true-of-date** equatorial coordinates to
/// horizontal coordinates using IAU 2006 GAST.
pub fn equatorial_to_horizontal_true_of_date<U: LengthUnit>(
    eq_position: &spherical::Position<Topocentric, frames::EquatorialTrueOfDate, U>,
    site: Geodetic<frames::ECEF>,
    jd: JulianDate,
) -> spherical::Position<Topocentric, frames::Horizontal, U> {
    let ctx: AstroContext = AstroContext::default();
    equatorial_to_horizontal_true_of_date_with_ctx(eq_position, site, jd, &ctx)
}

/// Context-aware variant of [`equatorial_to_horizontal_true_of_date`].
///
/// Uses the context EOP provider to derive UT1 for GAST.
pub fn equatorial_to_horizontal_true_of_date_with_ctx<U: LengthUnit, Ctx>(
    eq_position: &spherical::Position<Topocentric, frames::EquatorialTrueOfDate, U>,
    site: Geodetic<ECEF>,
    jd: JulianDate,
    ctx: &Ctx,
) -> spherical::Position<Topocentric, frames::Horizontal, U>
where
    Ctx: TransformContext,
{
    // Extract RA, Dec, distance from the equatorial position
    let ra = eq_position.azimuth;
    let dec = eq_position.polar;
    let distance = eq_position.distance;

    // True-of-date RA/Dec must use apparent sidereal time (GAST).
    let eop = ctx.astro_context().eop_at(jd);
    let jd_ut1 = jd_ut1_from_tt_eop(jd, &eop);
    let nut = Ctx::Nut::nutation(jd);
    let gast = gast_iau2006(jd_ut1, jd, nut.dpsi, nut.true_obliquity());
    let lst = gast + site.lon.to::<Radian>();
    let ha = (lst - ra.to::<Radian>())
        .value()
        .rem_euclid(std::f64::consts::TAU);
    let ha = Radians::new(ha);

    // Convert equatorial to horizontal using standard spherical trig
    let lat = site.lat.to::<Radian>();
    let dec_rad = dec.to::<Radian>();

    // Altitude: sin(alt) = sin(dec)·sin(lat) + cos(dec)·cos(lat)·cos(HA)
    let sin_alt = dec_rad.sin() * lat.sin() + dec_rad.cos() * lat.cos() * ha.cos();
    let alt = Radians::new(sin_alt.asin()).to::<Degree>();

    // Azimuth: tan(az) = -sin(HA) / (cos(lat)·tan(dec) - sin(lat)·cos(HA))
    let az_rad = (-dec_rad.cos() * ha.sin())
        .atan2(dec_rad.sin() * lat.cos() - dec_rad.cos() * ha.cos() * lat.sin());
    let az = Radians::new(az_rad).normalize().to::<Degree>();

    Topocentric::horizontal(site, alt, az, distance)
}

/// Computes the **horizontal direction** (altitude, azimuth) for a fixed-star
/// RA/Dec (J2000) from an observer site at a given Julian Date.
///
/// This is the stellar equivalent of `Sun::get_horizontal` / `Moon::get_horizontal`.
/// Since stars lack a distance, only a `Direction<Horizontal>` is returned.
///
/// ## Pipeline
///
/// 1. Precess J2000 → mean-of-date
/// 2. Apply nutation correction to RA
/// 3. Compute GAST → LAST → HA
/// 4. Standard equatorial→horizontal altitude/azimuth formula
pub fn star_horizontal(
    ra_j2000: crate::qtty::Degrees,
    dec_j2000: crate::qtty::Degrees,
    site: &Geodetic<frames::ECEF>,
    jd: JulianDate,
) -> spherical::Direction<frames::Horizontal> {
    // Full IAU 2006/2000B NPB matrix-based approach:
    // 1. Convert J2000 direction to Cartesian unit vector
    // 2. Apply full precession-nutation-bias (NPB) matrix → true-of-date
    // 3. Convert back to spherical → RA/Dec of date
    // 4. Compute GAST (IAU 2006) → LAST → HA → altitude/azimuth

    let ra_rad = ra_j2000.to::<Radian>();
    let dec_rad = dec_j2000.to::<Radian>();

    // J2000 direction as unit vector
    let (sin_ra, cos_ra) = ra_rad.sin_cos();
    let (sin_dec, cos_dec) = dec_rad.sin_cos();
    let x0 = cos_dec * cos_ra;
    let y0 = cos_dec * sin_ra;
    let z0 = sin_dec;

    // Full NPB matrix: GCRS → true equator/equinox of date
    let nut = nutation_iau2000b(jd);
    let npb = precession::precession_nutation_matrix(jd, nut.dpsi, nut.deps);
    let [x_t, y_t, z_t] = npb.apply_array([x0, y0, z0]);

    // True-of-date RA/Dec
    let ra_tod = Radians::new(y_t.atan2(x_t));
    let dec_tod = Radians::new(z_t.asin());

    // True-of-date RA/Dec requires apparent sidereal time (GAST).
    let ctx: AstroContext = AstroContext::default();
    let eop = ctx.eop_at(jd);
    let jd_ut1 = jd_ut1_from_tt_eop(jd, &eop);
    let gast = gast_iau2006(jd_ut1, jd, nut.dpsi, nut.true_obliquity());
    let lst = gast + site.lon.to::<Radian>();
    let ha = (lst - ra_tod).wrap_pos();

    // Equatorial → horizontal
    let lat = site.lat.to::<Radian>();
    let sin_alt = dec_tod.sin() * lat.sin() + dec_tod.cos() * lat.cos() * ha.cos();
    let alt = Radians::new(sin_alt.asin()).to::<Degree>();

    let az_rad = (-dec_tod.cos() * ha.sin())
        .atan2(dec_tod.sin() * lat.cos() - dec_tod.cos() * ha.cos() * lat.sin());
    let az = Radians::new(az_rad).normalize().to::<Degree>();

    spherical::Direction::<frames::Horizontal>::new(alt, az)
}
