// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Shared Horizontal Coordinate Helpers
//!
//! Common routines for converting equatorial coordinates to horizontal
//! (altitude/azimuth), shared across Sun, Moon, and stellar modules.
//!
//! These helpers factor out the duplicated equatorial→horizontal transform
//! that previously existed independently in `solar::sun_equations` and
//! `lunar::moon_equations`.

use crate::astro::earth_rotation::gmst_from_tt;
use crate::astro::nutation::nutation_iau2000b;
use crate::astro::precession;
use crate::coordinates::centers::ObserverSite;
use crate::coordinates::{cartesian, centers::*, frames, spherical};
use crate::time::JulianDate;
use qtty::{AstronomicalUnits, Degree, LengthUnit, Meter, Quantity, Radian, Radians};

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
    site: ObserverSite,
    jd: JulianDate,
) -> spherical::Position<Topocentric, frames::EquatorialTrueOfDate, U>
where
    Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
{
    use crate::coordinates::transform::centers::ToTopocentricExt;

    // 1) Translate geocentric → topocentric in J2000 frame (applies parallax)
    let topo_cart_j2000 = geo_cart_j2000.to_topocentric(site, jd);

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

    // 3) Apply nutation rotation (mean-of-date → true-of-date) — IAU 2000B
    let nut = nutation_iau2000b(jd);
    let rot_nut = crate::astro::nutation::nutation_rotation_iau2000b(jd);
    let _ = nut; // nut used by equatorial_to_horizontal if we need GAST later
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
/// latitude and the local sidereal time (GST → LST → HA).
///
/// Both Sun and Moon `get_horizontal` methods delegate here after computing
/// their body-specific apparent topocentric equatorial coordinates.
pub fn equatorial_to_horizontal<U: LengthUnit>(
    eq_position: &spherical::Position<Topocentric, frames::EquatorialTrueOfDate, U>,
    site: ObserverSite,
    jd: JulianDate,
) -> spherical::Position<Topocentric, frames::Horizontal, U> {
    // Extract RA, Dec, distance from the equatorial position
    let ra = eq_position.azimuth;
    let dec = eq_position.polar;
    let distance = eq_position.distance;

    // Compute hour angle: HA = LST - RA (using IAU 2006 GMST)
    // Uses tempoch ΔT for proper TT→UT1 conversion.
    let gmst = gmst_from_tt(jd);
    let lst = gmst + site.lon.to::<Radian>();
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
/// 3. Compute GST → LST → HA
/// 4. Standard equatorial→horizontal altitude/azimuth formula
pub fn star_horizontal(
    ra_j2000: qtty::Degrees,
    dec_j2000: qtty::Degrees,
    site: &ObserverSite,
    jd: JulianDate,
) -> spherical::Direction<frames::Horizontal> {
    // Full IAU 2006/2000B NPB matrix-based approach:
    // 1. Convert J2000 direction to Cartesian unit vector
    // 2. Apply full precession-nutation-bias (NPB) matrix → true-of-date
    // 3. Convert back to spherical → RA/Dec of date
    // 4. Compute GMST (IAU 2006) → LST → HA → altitude/azimuth

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

    // GMST (IAU 2006 ERA-based) → LST → HA
    // Uses tempoch ΔT for proper TT→UT1 conversion.
    let gmst = gmst_from_tt(jd);
    let lst = gmst + site.lon.to::<Radian>();
    let ha = Radians::new((lst - ra_tod).value().rem_euclid(std::f64::consts::TAU));

    // Equatorial → horizontal
    let lat = site.lat.to::<Radian>();
    let sin_alt = dec_tod.sin() * lat.sin() + dec_tod.cos() * lat.cos() * ha.cos();
    let alt = Radians::new(sin_alt.asin()).to::<Degree>();

    let az_rad = (-dec_tod.cos() * ha.sin())
        .atan2(dec_tod.sin() * lat.cos() - dec_tod.cos() * ha.cos() * lat.sin());
    let az = Radians::new(az_rad).normalize().to::<Degree>();

    spherical::Direction::<frames::Horizontal>::new(alt, az)
}
