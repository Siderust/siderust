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

use crate::astro::nutation::nutation_rotation;
use crate::astro::precession;
use crate::astro::sidereal::{calculate_gst, calculate_lst};
use crate::time::JulianDate;
use crate::coordinates::centers::ObserverSite;
use crate::coordinates::{cartesian, centers::*, frames, spherical};
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

    // 2) Rotate J2000 → mean-of-date using precession rotation
    let rot_prec = precession::precession_rotation_from_j2000(jd);
    let [x_m, y_m, z_m] = rot_prec.apply_array([
        topo_cart_j2000.x().value(),
        topo_cart_j2000.y().value(),
        topo_cart_j2000.z().value(),
    ]);

    let topo_cart_mod =
        cartesian::Position::<Topocentric, frames::EquatorialMeanOfDate, U>::new_with_params(
            *topo_cart_j2000.center_params(),
            Quantity::<U>::new(x_m),
            Quantity::<U>::new(y_m),
            Quantity::<U>::new(z_m),
        );

    // 3) Apply nutation rotation (mean-of-date → true-of-date)
    let rot_nut = nutation_rotation(jd);
    let [x_t, y_t, z_t] = rot_nut.apply_array([
        topo_cart_mod.x().value(),
        topo_cart_mod.y().value(),
        topo_cart_mod.z().value(),
    ]);

    let topo_cart_true =
        cartesian::Position::<Topocentric, frames::EquatorialTrueOfDate, U>::new_with_params(
            *topo_cart_mod.center_params(),
            Quantity::<U>::new(x_t),
            Quantity::<U>::new(y_t),
            Quantity::<U>::new(z_t),
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
    let ra = eq_position.azimuth();
    let dec = eq_position.polar();
    let distance = eq_position.distance();

    // Compute hour angle: HA = LST - RA
    let gst = calculate_gst(jd);
    let lst = calculate_lst(gst, site.lon);
    let ha = (lst - ra).normalize().to::<Radian>();

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

    spherical::Position::<Topocentric, frames::Horizontal, U>::new_with_site(
        site, alt, az, distance,
    )
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
    use crate::astro::nutation::corrected_ra_with_nutation;
    use crate::coordinates::frames::EquatorialMeanJ2000;

    // Build a spherical position in EquatorialMeanJ2000 (unit distance, irrelevant)
    let pos = spherical::Position::<Geocentric, EquatorialMeanJ2000, qtty::LightYear>::new(
        ra_j2000,
        dec_j2000,
        qtty::LightYears::new(1.0),
    );

    // Precess J2000 → mean-of-date
    let mean_of_date = precession::precess_from_j2000(pos, jd);
    // Apply nutation correction to RA
    let ra_corrected = corrected_ra_with_nutation(&mean_of_date.direction(), jd);
    let dec = mean_of_date.polar();

    // Compute hour angle
    let gst = calculate_gst(jd);
    let lst = calculate_lst(gst, site.lon);
    let ha = (lst - ra_corrected).normalize().to::<Radian>();

    // Equatorial → horizontal
    let lat = site.lat.to::<Radian>();
    let dec_rad = dec.to::<Radian>();
    let sin_alt = dec_rad.sin() * lat.sin() + dec_rad.cos() * lat.cos() * ha.cos();
    let alt = Radians::new(sin_alt.asin()).to::<Degree>();

    let az_rad = (-dec_rad.cos() * ha.sin())
        .atan2(dec_rad.sin() * lat.cos() - dec_rad.cos() * ha.cos() * lat.sin());
    let az = Radians::new(az_rad).normalize().to::<Degree>();

    spherical::Direction::<frames::Horizontal>::new(alt, az)
}
