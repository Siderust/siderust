//! Frame transformation: Equatorial → Horizontal
//!
//! This module implements the transformation from the equatorial frame to the
//! horizontal (alt-az) frame for topocentric coordinates. The transformation
//! uses the observer's site information embedded in the coordinate's center params.

use crate::astro::sidereal::{calculate_gst, calculate_lst};
use crate::astro::JulianDate;
use crate::coordinates::cartesian::Vector;
use crate::coordinates::centers::{ObserverSite, Topocentric};
use crate::coordinates::frames::{Equatorial, Horizontal};
use qtty::{Deg, Degrees, Quantity, Radian, Radians, Unit};

/// Performs the equatorial to horizontal coordinate transformation.
///
/// # Arguments
/// - `ra`: Right ascension
/// - `dec`: Declination
/// - `site`: Observer's geographic location
/// - `jd`: Julian Date of observation
///
/// # Returns
/// Tuple of (altitude, azimuth) as `Radians`
fn equatorial_to_horizontal_angles(
    ra: Degrees,
    dec: Degrees,
    site: &ObserverSite,
    jd: JulianDate,
) -> (Radians, Radians) {
    let gst = calculate_gst(jd);
    let lst = calculate_lst(gst, site.lon);

    let ha: Radians = (lst - ra).normalize().to::<Radian>();
    let lat: Radians = site.lat.to::<Radian>();
    let dec_rad: Radians = dec.to::<Radian>();

    // Trig functions return f64, wrap results back into Radians
    let alt_val = (dec_rad.sin() * lat.sin() + dec_rad.cos() * lat.cos() * ha.cos()).asin();
    let az_val = (-dec_rad.cos() * ha.sin())
        .atan2(dec_rad.sin() * lat.cos() - dec_rad.cos() * ha.cos() * lat.sin());

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
/// - `jd`: Julian Date of observation
///
/// # Returns
/// Tuple of (right ascension, declination) as `Radians`
fn horizontal_to_equatorial_angles(
    alt: Radians,
    az: Radians,
    site: &ObserverSite,
    jd: JulianDate,
) -> (Radians, Radians) {
    let lat: Radians = site.lat.to::<Radian>();

    // Calculate declination - trig functions return f64
    let sin_dec = alt.sin() * lat.sin() + alt.cos() * lat.cos() * az.cos();
    let dec_val = sin_dec.asin();

    // Calculate hour angle
    let cos_ha = (alt.sin() - lat.sin() * dec_val.sin()) / (lat.cos() * dec_val.cos());
    let sin_ha = -alt.cos() * az.sin() / dec_val.cos();
    let ha_val = sin_ha.atan2(cos_ha);

    // Convert hour angle to right ascension
    let gst = calculate_gst(jd);
    let lst: Radians = calculate_lst(gst, site.lon).to::<Radian>();
    let ra_val = (lst.value() - ha_val).rem_euclid(2.0 * std::f64::consts::PI);

    (
        Quantity::<Radian>::new(ra_val),
        Quantity::<Radian>::new(dec_val),
    )
}

// =============================================================================
// Equatorial → Horizontal (for Topocentric center)
// =============================================================================

use crate::coordinates::transform::Transform;

/// Transform from Equatorial to Horizontal frame for Topocentric coordinates.
///
/// This transformation requires the Julian Date to compute the local sidereal time.
/// The observer's site information is taken from the coordinate's center params.
impl<U: Unit> Transform<Vector<Topocentric, Horizontal, U>> for Vector<Topocentric, Equatorial, U> {
    fn transform(&self, jd: JulianDate) -> Vector<Topocentric, Horizontal, U> {
        let site = self.center_params();
        let r = self.distance();

        let dec: Radians = if r.value() > 0.0 {
            Quantity::<Radian>::new((self.z() / r).asin())
        } else {
            Quantity::<Radian>::new(0.0)
        };

        let ra: Radians = Quantity::<Radian>::new(self.y().value().atan2(self.x().value()));

        let (alt, az) = equatorial_to_horizontal_angles(ra.to::<Deg>(), dec.to::<Deg>(), site, jd);

        // Convert back to Cartesian in horizontal frame
        // In horizontal: x = North, y = West, z = Zenith
        // Trig functions return f64, multiply by Quantity<U> to get Quantity<U>
        let new_x = r * alt.cos() * az.cos();
        let new_y = -r * alt.cos() * az.sin(); // negative for East-positive azimuth
        let new_z = r * alt.sin();

        Vector::from_vec3(*site, nalgebra::Vector3::new(new_x, new_y, new_z))
    }
}

// =============================================================================
// Horizontal → Equatorial (for Topocentric center)
// =============================================================================

/// Transform from Horizontal to Equatorial frame for Topocentric coordinates.
impl<U: Unit> Transform<Vector<Topocentric, Equatorial, U>> for Vector<Topocentric, Horizontal, U> {
    fn transform(&self, jd: JulianDate) -> Vector<Topocentric, Equatorial, U> {
        let site = self.center_params();

        // Get distance and angles from Cartesian vector
        let r = self.distance();
        // Division of same units gives Per<U,U> which has .asin() returning f64
        let alt: Radians = if r.value() > 0.0 {
            Quantity::<Radian>::new((self.z() / r).asin())
        } else {
            Quantity::<Radian>::new(0.0)
        };
        // atan2 on Quantity<U> values - extract raw values for atan2
        let az: Radians = Quantity::<Radian>::new((-self.y()).value().atan2(self.x().value()));

        let (ra, dec) = horizontal_to_equatorial_angles(alt, az, site, jd);

        // Convert back to Cartesian in equatorial frame
        // Trig functions return f64, multiply by Quantity<U> to get Quantity<U>
        let new_x = r * dec.cos() * ra.cos();
        let new_y = r * dec.cos() * ra.sin();
        let new_z = r * dec.sin();

        Vector::from_vec3(*site, nalgebra::Vector3::new(new_x, new_y, new_z))
    }
}

// =============================================================================
// SphericalCoord implementations
// =============================================================================

use crate::coordinates::spherical::SphericalCoord;

impl<U: Unit> Transform<SphericalCoord<Topocentric, Horizontal, U>>
    for SphericalCoord<Topocentric, Equatorial, U>
{
    fn transform(&self, jd: JulianDate) -> SphericalCoord<Topocentric, Horizontal, U> {
        let cart: Vector<Topocentric, Horizontal, U> = self.to_cartesian().transform(jd);
        cart.to_spherical()
    }
}

impl<U: Unit> Transform<SphericalCoord<Topocentric, Equatorial, U>>
    for SphericalCoord<Topocentric, Horizontal, U>
{
    fn transform(&self, jd: JulianDate) -> SphericalCoord<Topocentric, Equatorial, U> {
        let cart: Vector<Topocentric, Equatorial, U> = self.to_cartesian().transform(jd);
        cart.to_spherical()
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
        let gst = calculate_gst(jd);
        let lst = calculate_lst(gst, site.lon);

        // RA = LST means HA = 0 (on meridian)
        let ra: Degrees = lst;
        let dec: Degrees = site.lat;

        let (alt, _az) = equatorial_to_horizontal_angles(ra, dec, &site, jd);

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

        let (alt, az) = equatorial_to_horizontal_angles(ra, dec, &site, jd);
        let (ra_back, dec_back) = horizontal_to_equatorial_angles(alt, az, &site, jd);

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
