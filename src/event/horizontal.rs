// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Legacy fixed-star horizon helpers.
//!
//! The policy on these functions selects a historical input-frame
//! interpretation. For typed geometric projection of an existing sightline,
//! use [`project_direction_with`](crate::coordinates::transform::horizontal::project_direction_with).

use crate::astro::apparent::CorrectionPolicy;
use crate::astro::earth_rotation::jd_ut1_from_tt_eop;
use crate::astro::nutation::nutation_iau2000b;
use crate::astro::precession;
use crate::astro::sidereal::gast_iau2006;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::transform::horizontal::project;
use crate::coordinates::transform::AstroContext;
use crate::coordinates::{frames, spherical};
use crate::qtty::{Degree, Radian, Radians};
#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Transcendental};
use crate::time::JulianDate;

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
    let dec_tod = Radians::new(z_t.atan2((x_t * x_t + y_t * y_t).sqrt()));

    let ctx: AstroContext = AstroContext::default();
    let eop = ctx.eop_at_tt(jd);
    let ut1 = jd_ut1_from_tt_eop(jd, &eop);
    let gast = gast_iau2006(ut1, jd, nut.dpsi, nut.mean_obliquity);
    // The legacy policy can interpret the supplied angles directly as of-date.
    project(
        &spherical::Direction::new_unchecked(dec_tod.to::<Degree>(), ra_tod.to::<Degree>()),
        site,
        affn::Rotation3::rz(-gast),
    )
}
