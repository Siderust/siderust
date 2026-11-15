// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Compatibility helpers for Earth-based horizon and topocentric coordinates.
//!
//! Compose [`TransformCenter`], [`PositionAstroExt`] and
//! [`project_position_with`](crate::coordinates::transform::horizontal::project_position_with)
//! for new geometric workflows. The legacy names use “apparent” for true-of-date axes and apparent sidereal time. These
//! geometric operations apply neither light propagation nor refraction.

use crate::astro::apparent::CorrectionPolicy;
use crate::astro::earth_rotation::jd_ut1_from_tt_eop;
use crate::astro::nutation::nutation_iau2000b;
use crate::astro::precession;
use crate::astro::sidereal::gast_iau2006;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::coordinates::transform::horizontal::{project, project_direction_with};
use crate::coordinates::transform::{
    AstroContext, PositionAstroExt, TransformCenter, TransformContext,
};
use crate::coordinates::{cartesian, centers::*, frames, spherical};
use crate::qtty::{AstronomicalUnits, Degree, LengthUnit, Meter, Quantity, Radian, Radians};
#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Transcendental};
use crate::time::JulianDate;

// =============================================================================
// Topocentric geometry
// =============================================================================

/// Transforms a **geocentric, mean-J2000 equatorial** Cartesian position into
/// **topocentric, true-of-date equatorial** spherical coordinates.
///
/// Applies observer parallax, followed by precession and nutation from the
/// default context, then expresses the position as right ascension,
/// declination, and distance. To compose the same geometry directly, use
/// [`TransformCenter::to_center_with`] and [`PositionAstroExt::to_frame_with`].
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
    Ctx::Eph: crate::ephemeris::Ephemeris,
    Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
{
    let topocentric: cartesian::Position<Topocentric, _, _> =
        geo_cart_j2000.to_center_with(site, jd, ctx);
    topocentric
        .to_frame_with::<frames::EquatorialTrueOfDate, _>(&jd, ctx)
        .to_spherical()
}

// =============================================================================
// Equatorial → Horizontal Transform
// =============================================================================

/// Converts a **true-of-date topocentric equatorial** position to **horizontal**
/// coordinates (altitude, azimuth, distance).
///
/// This is the standard equatorial-to-horizontal transform using the observer's
/// latitude and the local apparent sidereal time (GAST → LAST → HA).
///
/// The Sun and Moon horizon helpers delegate here after computing their
/// topocentric equatorial coordinates.
pub fn equatorial_to_horizontal<U: LengthUnit>(
    eq_position: &spherical::Position<Topocentric, frames::EquatorialTrueOfDate, U>,
    site: Geodetic<frames::ECEF>,
    jd: JulianDate,
) -> spherical::Position<Topocentric, frames::Horizontal, U> {
    equatorial_to_horizontal_true_of_date(eq_position, site, jd)
}

/// Converts topocentric **true-of-date** equatorial coordinates to
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
/// Uses the context EOP provider to derive UT1 for GAST. The explicit
/// `site` controls the projection and output metadata, even if it differs
/// from the input center.
/// [`project_position_with`](crate::coordinates::transform::horizontal::project_position_with)
/// uses the observer encoded in the input position.
pub fn equatorial_to_horizontal_true_of_date_with_ctx<U: LengthUnit, Ctx>(
    eq_position: &spherical::Position<Topocentric, frames::EquatorialTrueOfDate, U>,
    site: Geodetic<ECEF>,
    jd: JulianDate,
    ctx: &Ctx,
) -> spherical::Position<Topocentric, frames::Horizontal, U>
where
    Ctx: TransformContext,
{
    // Preserve the legacy explicit-site override. The composable position API
    // instead obtains its horizon from the observer encoded in the position.
    let direction = project_direction_with(&eq_position.direction(), &site, jd, ctx);
    Topocentric::horizontal(site, direction.alt(), direction.az(), eq_position.distance)
}

/// Computes the **horizontal direction** (altitude, azimuth) for a fixed-star
/// RA/Dec (J2000) from an observer site at a given Julian Date.
///
/// Applies the legacy IAU 2006/2000B precession-nutation-bias model, then
/// projects the resulting true-of-date direction using GAST. A direction-only
/// input carries no parallax information.
pub fn star_horizontal(
    ra_j2000: crate::qtty::Degrees,
    dec_j2000: crate::qtty::Degrees,
    site: &Geodetic<frames::ECEF>,
    jd: JulianDate,
) -> spherical::Direction<frames::Horizontal> {
    star_horizontal_with_policy(ra_j2000, dec_j2000, site, jd, CorrectionPolicy::APPARENT)
}

/// Computes the horizontal direction for a fixed-star RA/Dec with an explicit
/// legacy mode selector. `GEOMETRIC` interprets the supplied angles as
/// true-of-date, while every nonempty policy applies IAU 2006/2000B
/// precession, nutation and frame bias. Individual policy bits do not
/// select physical corrections here.
pub fn star_horizontal_with_policy(
    ra_j2000: crate::qtty::Degrees,
    dec_j2000: crate::qtty::Degrees,
    site: &Geodetic<frames::ECEF>,
    jd: JulianDate,
    policy: CorrectionPolicy,
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

    let nut = nutation_iau2000b(jd);
    let [x_t, y_t, z_t] = if policy == CorrectionPolicy::GEOMETRIC {
        [x0, y0, z0]
    } else {
        let npb = precession::precession_nutation_matrix(jd, nut.dpsi, nut.deps);
        npb.apply_array([x0, y0, z0])
    };

    // True-of-date RA/Dec
    let ra_tod = Radians::new(y_t.atan2(x_t));
    let dec_tod = Radians::new(z_t.asin());

    let ctx: AstroContext = AstroContext::default();
    let eop = ctx.eop_at_tt(jd);
    let ut1 = jd_ut1_from_tt_eop(jd, &eop);
    let gast = gast_iau2006(ut1, jd, nut.dpsi, nut.mean_obliquity);
    // The legacy policy can interpret the supplied angles directly as of-date.
    project(
        spherical::Direction::new_unchecked(dec_tod.to::<Degree>(), ra_tod.to::<Degree>()),
        site,
        gast,
    )
}
