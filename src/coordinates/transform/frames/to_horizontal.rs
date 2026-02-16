// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame transformation: EquatorialMeanOfDate → Horizontal
//!
//! This module implements the transformation from the equatorial frame to the
//! horizontal (alt-az) frame for topocentric coordinates. The transformation
//! uses the observer's site information embedded in the coordinate's center params.

use crate::astro::sidereal::gmst_iau2006;
use crate::coordinates::centers::{ObserverSite, Topocentric};
use crate::coordinates::frames::{EquatorialMeanOfDate, Horizontal};
use crate::coordinates::{cartesian, spherical};
use crate::time::JulianDate;
use qtty::{Deg, Degrees, LengthUnit, Quantity, Radian, Radians};

/// Precomputed trigonometric values for an observer's latitude.
///
/// Computing `sin(lat)` and `cos(lat)` is expensive relative to the
/// rest of the horizontal transform. This struct caches these values
/// for repeated use with the same site.
#[derive(Debug, Clone, Copy)]
struct SiteTrig {
    sin_lat: f64,
    cos_lat: f64,
}

impl SiteTrig {
    /// Precompute sin/cos of latitude from an observer site.
    #[inline]
    fn from_site(site: &ObserverSite) -> Self {
        let lat_rad: Radians = site.lat.to::<Radian>();
        let (sin_lat, cos_lat) = lat_rad.sin_cos();
        Self { sin_lat, cos_lat }
    }
}

/// Performs the equatorial to horizontal coordinate transformation.
///
/// # Arguments
/// - `ra`: Right ascension
/// - `dec`: Declination
/// - `site`: Observer's geographic location
/// - `site_trig`: Precomputed sin/cos of observer latitude
/// - `jd`: Julian Date of observation
///
/// # Returns
/// Tuple of (altitude, azimuth) as `Radians`
fn equatorial_to_horizontal_angles(
    ra: Degrees,
    dec: Degrees,
    site: &ObserverSite,
    site_trig: &SiteTrig,
    jd: JulianDate,
) -> (Radians, Radians) {
    let gmst = gmst_iau2006(jd, jd);
    let lst_rad = gmst.value() + site.lon.to::<Radian>().value();

    let ra_rad: Radians = ra.to::<Radian>();
    let ha = Quantity::<Radian>::new((lst_rad - ra_rad.value()).rem_euclid(std::f64::consts::TAU));
    let dec_rad: Radians = dec.to::<Radian>();

    let (sin_dec, cos_dec) = dec_rad.sin_cos();
    let cos_ha = ha.cos();

    // Use precomputed sin/cos of latitude
    let alt_val = (sin_dec * site_trig.sin_lat + cos_dec * site_trig.cos_lat * cos_ha).asin();
    let az_val = (-cos_dec * ha.sin())
        .atan2(sin_dec * site_trig.cos_lat - cos_dec * cos_ha * site_trig.sin_lat);

    (
        Quantity::<Radian>::new(alt_val),
        Quantity::<Radian>::new(az_val),
    )
}

/// Performs the horizontal to equatorial coordinate transformation.
///
/// # Arguments
/// - `alt`: Altitude
/// - `az`: Azimuth (from North, clockwise)
/// - `site`: Observer's geographic location
/// - `site_trig`: Precomputed sin/cos of observer latitude
/// - `jd`: Julian Date of observation
///
/// # Returns
/// Tuple of (right ascension, declination) as `Radians`
fn horizontal_to_equatorial_angles(
    alt: Radians,
    az: Radians,
    site: &ObserverSite,
    site_trig: &SiteTrig,
    jd: JulianDate,
) -> (Radians, Radians) {
    // Use precomputed sin/cos of latitude
    let (sin_alt, cos_alt) = alt.sin_cos();

    // Calculate declination
    let sin_dec = sin_alt * site_trig.sin_lat + cos_alt * site_trig.cos_lat * az.cos();
    let dec_val = sin_dec.asin();

    // Calculate hour angle
    let cos_ha =
        (sin_alt - site_trig.sin_lat * dec_val.sin()) / (site_trig.cos_lat * dec_val.cos());
    let sin_ha = -cos_alt * az.sin() / dec_val.cos();
    let ha_val = sin_ha.atan2(cos_ha);

    // Convert hour angle to right ascension
    let gmst = gmst_iau2006(jd, jd);
    let lst = gmst + site.lon.to::<Radian>();
    let ha: Radians = Quantity::<Radian>::new(ha_val);
    let ra = (lst - ha).normalize();

    (ra, Quantity::<Radian>::new(dec_val))
}

// =============================================================================
// EquatorialMeanOfDate → Horizontal (for Topocentric center)
// =============================================================================

use crate::coordinates::transform::Transform;

/// Transform from EquatorialMeanOfDate to Horizontal frame for Topocentric coordinates.
///
/// This transformation requires the Julian Date to compute the local sidereal time.
/// The observer's site information is taken from the coordinate's center params.
impl<U: LengthUnit> Transform<cartesian::Position<Topocentric, Horizontal, U>>
    for cartesian::Position<Topocentric, EquatorialMeanOfDate, U>
{
    fn transform(&self, jd: JulianDate) -> cartesian::Position<Topocentric, Horizontal, U> {
        let site = self.center_params();
        let r = self.distance();

        let dec: Radians = if r > 0.0 {
            Quantity::<Radian>::new((self.z() / r).asin())
        } else {
            Quantity::<Radian>::new(0.0)
        };

        let ra: Radians = Quantity::<Radian>::new(self.y().value().atan2(self.x().value()));

        let site_trig = SiteTrig::from_site(site);
        let (alt, az) =
            equatorial_to_horizontal_angles(ra.to::<Deg>(), dec.to::<Deg>(), site, &site_trig, jd);

        // Convert back to Cartesian in horizontal frame
        // In horizontal: x = North, y = West, z = Zenith
        // Trig functions return f64, multiply by Quantity<U> to get Quantity<U>
        let new_x = r * alt.cos() * az.cos();
        let new_y = -r * alt.cos() * az.sin(); // negative for East-positive azimuth
        let new_z = r * alt.sin();

        cartesian::Position::from_vec3(*site, nalgebra::Vector3::new(new_x, new_y, new_z))
    }
}

// =============================================================================
// Horizontal → EquatorialMeanOfDate (for Topocentric center)
// =============================================================================

/// Transform from Horizontal to EquatorialMeanOfDate frame for Topocentric coordinates.
impl<U: LengthUnit> Transform<cartesian::Position<Topocentric, EquatorialMeanOfDate, U>>
    for cartesian::Position<Topocentric, Horizontal, U>
{
    fn transform(
        &self,
        jd: JulianDate,
    ) -> cartesian::Position<Topocentric, EquatorialMeanOfDate, U> {
        let site = self.center_params();

        // Get distance and angles from Cartesian vector
        let r = self.distance();
        // Division of same units gives Per<U,U> which has .asin() returning f64
        let alt: Radians = if r > 0.0 {
            Quantity::<Radian>::new((self.z() / r).asin())
        } else {
            Quantity::<Radian>::new(0.0)
        };
        // atan2 on Quantity<U> values - extract raw values for atan2
        let az: Radians = Quantity::<Radian>::new((-self.y()).value().atan2(self.x().value()));

        let site_trig = SiteTrig::from_site(site);
        let (ra, dec) = horizontal_to_equatorial_angles(alt, az, site, &site_trig, jd);

        // Convert back to Cartesian in equatorial frame
        // Trig functions return f64, multiply by Quantity<U> to get Quantity<U>
        let new_x = r * dec.cos() * ra.cos();
        let new_y = r * dec.cos() * ra.sin();
        let new_z = r * dec.sin();

        cartesian::Position::from_vec3(*site, nalgebra::Vector3::new(new_x, new_y, new_z))
    }
}

// =============================================================================
// Spherical position implementations
// =============================================================================

impl<U: LengthUnit> Transform<spherical::Position<Topocentric, Horizontal, U>>
    for spherical::Position<Topocentric, EquatorialMeanOfDate, U>
{
    fn transform(&self, jd: JulianDate) -> spherical::Position<Topocentric, Horizontal, U> {
        let cart: cartesian::Position<Topocentric, Horizontal, U> =
            self.to_cartesian().transform(jd);
        spherical::Position::from_cartesian(&cart)
    }
}

impl<U: LengthUnit> Transform<spherical::Position<Topocentric, EquatorialMeanOfDate, U>>
    for spherical::Position<Topocentric, Horizontal, U>
{
    fn transform(
        &self,
        jd: JulianDate,
    ) -> spherical::Position<Topocentric, EquatorialMeanOfDate, U> {
        let cart: cartesian::Position<Topocentric, EquatorialMeanOfDate, U> =
            self.to_cartesian().transform(jd);
        spherical::Position::from_cartesian(&cart)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::{DEG, M, RAD};

    #[test]
    fn test_equatorial_to_horizontal_zenith() {
        // An object at the zenith should have altitude 90°
        let site = ObserverSite::new(0.0 * DEG, 45.0 * DEG, 0.0 * M);
        let jd = JulianDate::J2000;

        // At the zenith, the object's declination equals the latitude
        // and hour angle is 0 (on the meridian)
        let gmst = gmst_iau2006(jd, jd);
        let lst_rad = gmst.value() + site.lon.to::<Radian>().value();

        // RA = LST means HA = 0 (on meridian)
        let ra: Degrees = Quantity::<Deg>::new(lst_rad.to_degrees());
        let dec: Degrees = site.lat;

        let (alt, _az) =
            equatorial_to_horizontal_angles(ra, dec, &site, &SiteTrig::from_site(&site), jd);

        let expected_alt = std::f64::consts::FRAC_PI_2 * RAD;
        assert!(
            (alt - expected_alt).abs() < 1e-6 * RAD,
            "Expected altitude ~90°, got {}",
            alt.to::<Deg>()
        );
    }

    #[test]
    fn test_roundtrip_angles() {
        // Test that the angle conversion functions are inverses of each other
        let site = ObserverSite::new(-17.8925 * DEG, 28.7543 * DEG, 2396.0 * M);
        let jd = JulianDate::new(2460677.0);

        let ra = 101.29 * DEG;
        let dec = -16.72 * DEG;

        let site_trig = SiteTrig::from_site(&site);
        let (alt, az) = equatorial_to_horizontal_angles(ra, dec, &site, &site_trig, jd);
        let (ra_back, dec_back) = horizontal_to_equatorial_angles(alt, az, &site, &site_trig, jd);

        let ra_rad: Radians = ra.to::<Radian>();
        let dec_rad: Radians = dec.to::<Radian>();

        // RA wraps around 2π, so compare with modulo
        let two_pi = 2.0 * std::f64::consts::PI * RAD;
        let ra_diff = (ra_rad - ra_back).abs().wrap_pos();
        let ra_diff = ra_diff.min(two_pi - ra_diff);

        assert!(
            ra_diff < 1e-6 * RAD,
            "RA mismatch: {} vs {} (diff: {})",
            ra,
            ra_back.to::<Deg>(),
            ra_diff.to::<Deg>()
        );
        assert!(
            (dec_rad - dec_back).abs() < 1e-6 * RAD,
            "Dec mismatch: {} vs {}",
            dec,
            dec_back.to::<Deg>()
        );
    }
}
