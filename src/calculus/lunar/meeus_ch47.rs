// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Simplified Meeus Ch.47 Lunar Position
//!
//! Implements the reduced-term geocentric Moon model from
//! *Astronomical Algorithms*, 2nd ed., Ch. 47 (Jean Meeus).
//!
//! This is **intentionally simplified** (6 longitude, 4 latitude, 4 distance
//! terms) for use as a lightweight benchmark reference. For production
//! accuracy see [`crate::calculus::lunar`] (ELP2000-82B).
//!
//! ## References
//! * Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed., Ch. 47.

use crate::astro::precession_iau2006;
use crate::time::JulianDate;
use qtty::{Kilometers, Radians};

/// Geocentric apparent equatorial + ecliptic coordinates of the Moon.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MoonMeeusCh47 {
    /// Right ascension (equatorial).
    pub ra: Radians,
    /// Declination (equatorial).
    pub dec: Radians,
    /// Geocentric distance.
    pub dist: Kilometers,
    /// Ecliptic longitude (of date).
    pub ecl_lon: Radians,
    /// Ecliptic latitude (of date).
    pub ecl_lat: Radians,
}

/// Compute simplified Meeus Ch.47 geocentric lunar position.
///
/// Returns both ecliptic-of-date and equatorial coordinates. The ecliptic→
/// equatorial conversion uses the IAU 2006 mean obliquity of date.
pub fn moon_position_meeus_ch47(jd_tt: JulianDate) -> MoonMeeusCh47 {
    let date2 = jd_tt.value() - 2_451_545.0;
    let t = date2 / 36_525.0;

    // Mean elements (degrees), Meeus Ch.47
    let lp = (218.316_447_7 + 481_267.881_234_21 * t - 0.001_578_6 * t * t + t * t * t / 538_841.0
        - t * t * t * t / 65_194_000.0)
        % 360.0;
    let d = (297.850_192_1 + 445_267.111_403_4 * t - 0.001_881_9 * t * t + t * t * t / 545_868.0
        - t * t * t * t / 113_065_000.0)
        % 360.0;
    let m_sun = (357.529_109_2 + 35_999.050_290_9 * t - 0.000_153_6 * t * t
        + t * t * t / 24_490_000.0)
        % 360.0;
    let mp = (134.963_396_4 + 477_198.867_505_5 * t + 0.008_741_4 * t * t + t * t * t / 69_699.0
        - t * t * t * t / 14_712_000.0)
        % 360.0;
    let f = (93.272_095_0 + 483_202.017_523_3 * t - 0.003_653_9 * t * t - t * t * t / 3_526_000.0
        + t * t * t * t / 863_310_000.0)
        % 360.0;

    let d_r = d.to_radians();
    let m_r = m_sun.to_radians();
    let mp_r = mp.to_radians();
    let f_r = f.to_radians();

    // Major longitude terms (×1e-6 deg)
    let sum_l = 6_288_774.0 * mp_r.sin()
        + 1_274_027.0 * (2.0 * d_r - mp_r).sin()
        + 658_314.0 * (2.0 * d_r).sin()
        + 213_618.0 * (2.0 * mp_r).sin()
        - 185_116.0 * m_r.sin()
        - 114_332.0 * (2.0 * f_r).sin();

    // Major latitude terms
    let sum_b = 5_128_122.0 * f_r.sin()
        + 280_602.0 * (mp_r + f_r).sin()
        + 277_693.0 * (mp_r - f_r).sin()
        + 173_237.0 * (2.0 * d_r - f_r).sin();

    // Major distance terms (km)
    let sum_r = -20_905_355.0 * mp_r.cos()
        - 3_699_111.0 * (2.0 * d_r - mp_r).cos()
        - 2_955_968.0 * (2.0 * d_r).cos()
        - 569_925.0 * (2.0 * mp_r).cos();

    let ecl_lon_rad = (lp + sum_l / 1_000_000.0)
        .to_radians()
        .rem_euclid(std::f64::consts::TAU);
    let ecl_lat_rad = (sum_b / 1_000_000.0).to_radians();
    let dist_km = 385_000.56 + sum_r / 1_000.0;

    // Ecliptic → equatorial with IAU 2006 mean obliquity
    let eps = precession_iau2006::mean_obliquity_iau2006(jd_tt).value();
    let ce = eps.cos();
    let se = eps.sin();
    let ra = (ecl_lon_rad.sin() * ce - ecl_lat_rad.tan() * se).atan2(ecl_lon_rad.cos());
    let ra_val = ra.rem_euclid(std::f64::consts::TAU);
    let dec_val = (ecl_lat_rad.sin() * ce + ecl_lat_rad.cos() * se * ecl_lon_rad.sin())
        .clamp(-1.0, 1.0)
        .asin();

    MoonMeeusCh47 {
        ra: Radians::new(ra_val),
        dec: Radians::new(dec_val),
        dist: Kilometers::new(dist_km),
        ecl_lon: Radians::new(ecl_lon_rad),
        ecl_lat: Radians::new(ecl_lat_rad),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn meeus_example_1992_apr_12() {
        // Meeus Example 47.a: 1992 April 12, 00:00 TD
        // JDE = 2448724.5
        let jd = JulianDate::new(2_448_724.5);
        let moon = moon_position_meeus_ch47(jd);

        // Expected: λ ≈ 133.17°, β ≈ -3.23°, Δ ≈ 368409 km (reduced-term)
        let lon_deg = moon.ecl_lon.value().to_degrees();
        let lat_deg = moon.ecl_lat.value().to_degrees();

        // Loose tolerances for the reduced 6-term model
        assert!(
            (lon_deg - 133.17).abs() < 1.0,
            "ecl lon = {lon_deg}°, expected ≈ 133.17°"
        );
        assert!(
            (lat_deg - (-3.23)).abs() < 0.5,
            "ecl lat = {lat_deg}°, expected ≈ -3.23°"
        );
        assert!(
            (moon.dist.value() - 368_409.0).abs() < 3000.0,
            "dist = {} km, expected ≈ 368409 km",
            moon.dist.value()
        );
    }
}
