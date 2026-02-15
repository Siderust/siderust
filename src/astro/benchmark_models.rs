// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Benchmark-facing astronomical models used by lab adapters.
//!
//! These helpers keep adapter crates thin by centralizing model formulas in
//! `siderust`. They expose deterministic, low-level `f64` APIs that mirror
//! benchmark protocol inputs/outputs.

use crate::astro::{nutation_iau2000b, precession_iau2006, sidereal};
use crate::time::JulianDate;

#[inline]
fn mat_mul(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

#[inline]
fn mat_transpose(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

#[inline]
fn wrap_0_2pi(a: f64) -> f64 {
    a.rem_euclid(std::f64::consts::TAU)
}

/// Convert equatorial (RA/Dec) to ecliptic-of-date (lon/lat).
///
/// The transform matches the adapter's benchmark path:
/// `R = Rx(-eps_a(t)) * P_iau2006(t)` where `P` is GCRS->mean-equator-of-date.
pub fn equatorial_to_ecliptic_of_date_iau2006(
    jd_tt: JulianDate,
    ra_rad: f64,
    dec_rad: f64,
) -> (f64, f64) {
    let prec = precession_iau2006::precession_matrix_iau2006(jd_tt);
    let eps_a = precession_iau2006::mean_obliquity_iau2006(jd_tt).value();

    let v_in = [
        dec_rad.cos() * ra_rad.cos(),
        dec_rad.cos() * ra_rad.sin(),
        dec_rad.sin(),
    ];
    let prec_mat = *prec.as_matrix();
    let v_prec = mat_mul(&prec_mat, v_in);

    let ce = eps_a.cos();
    let se = eps_a.sin();
    let v_ecl = [
        v_prec[0],
        v_prec[1] * ce + v_prec[2] * se,
        -v_prec[1] * se + v_prec[2] * ce,
    ];

    let lon = wrap_0_2pi(v_ecl[1].atan2(v_ecl[0]));
    let lat = v_ecl[2].clamp(-1.0, 1.0).asin();
    (lon, lat)
}

/// Convert ecliptic-of-date (lon/lat) to equatorial (RA/Dec).
///
/// Inverse of [`equatorial_to_ecliptic_of_date_iau2006`]:
/// `P^T(t) * Rx(+eps_a(t))`.
pub fn ecliptic_of_date_to_equatorial_iau2006(
    jd_tt: JulianDate,
    ecl_lon_rad: f64,
    ecl_lat_rad: f64,
) -> (f64, f64) {
    let prec = precession_iau2006::precession_matrix_iau2006(jd_tt);
    let eps_a = precession_iau2006::mean_obliquity_iau2006(jd_tt).value();
    let (sl, cl) = ecl_lon_rad.sin_cos();
    let (sb, cb) = ecl_lat_rad.sin_cos();

    let v_ecl = [cb * cl, cb * sl, sb];
    let ce = eps_a.cos();
    let se = eps_a.sin();
    let v_eq = [
        v_ecl[0],
        v_ecl[1] * ce - v_ecl[2] * se,
        v_ecl[1] * se + v_ecl[2] * ce,
    ];

    let prec_t = mat_transpose(prec.as_matrix());
    let v_back = mat_mul(&prec_t, v_eq);
    let ra = wrap_0_2pi(v_back[1].atan2(v_back[0]));
    let dec = v_back[2].clamp(-1.0, 1.0).asin();
    (ra, dec)
}

/// Convert equatorial (RA/Dec) to horizontal (Az/Alt), using GAST IAU 2006.
///
/// Azimuth is normalized to `[0, 2pi)`.
pub fn equatorial_to_horizontal_gast_iau2006(
    jd_ut1: JulianDate,
    jd_tt: JulianDate,
    ra_rad: f64,
    dec_rad: f64,
    obs_lon_rad: f64,
    obs_lat_rad: f64,
) -> (f64, f64) {
    let nut = nutation_iau2000b::nutation_iau2000b(jd_tt);
    let gast = sidereal::gast_iau2006(jd_ut1, jd_tt, nut.dpsi, nut.true_obliquity());
    let last_rad = gast.value() + obs_lon_rad;
    let ha_rad = last_rad - ra_rad;

    let sh = ha_rad.sin();
    let ch = ha_rad.cos();
    let sd = dec_rad.sin();
    let cd = dec_rad.cos();
    let sp = obs_lat_rad.sin();
    let cp = obs_lat_rad.cos();

    let x = -ch * cd * sp + sd * cp;
    let y = -sh * cd;
    let z = ch * cd * cp + sd * sp;
    let r = (x * x + y * y).sqrt();
    let az = wrap_0_2pi(if r != 0.0 { y.atan2(x) } else { 0.0 });
    let alt = z.atan2(r);
    (az, alt)
}

/// Convert horizontal (Az/Alt) to equatorial (RA/Dec), using GAST IAU 2006.
pub fn horizontal_to_equatorial_gast_iau2006(
    jd_ut1: JulianDate,
    jd_tt: JulianDate,
    az_rad: f64,
    alt_rad: f64,
    obs_lon_rad: f64,
    obs_lat_rad: f64,
) -> (f64, f64) {
    let nut = nutation_iau2000b::nutation_iau2000b(jd_tt);
    let gast = sidereal::gast_iau2006(jd_ut1, jd_tt, nut.dpsi, nut.true_obliquity());
    let last_rad = gast.value() + obs_lon_rad;

    let sp = obs_lat_rad.sin();
    let cp = obs_lat_rad.cos();
    let dec_rad = (sp * alt_rad.sin() + cp * alt_rad.cos() * az_rad.cos())
        .clamp(-1.0, 1.0)
        .asin();
    let ha_rad = (-az_rad.sin() * alt_rad.cos())
        .atan2(alt_rad.sin() * cp - alt_rad.cos() * az_rad.cos() * sp);
    let ra_rad = wrap_0_2pi(last_rad - ha_rad);
    (ra_rad, dec_rad)
}

/// Simplified Meeus Ch.47 lunar apparent geocentric equatorial position.
///
/// This is intentionally the reduced-term benchmark model (major terms only).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MoonMeeusCh47Simplified {
    pub ra_rad: f64,
    pub dec_rad: f64,
    pub dist_km: f64,
    pub ecl_lon_rad: f64,
    pub ecl_lat_rad: f64,
}

/// Compute simplified lunar position terms used in benchmark adapters.
pub fn moon_position_meeus_ch47_simplified(jd_tt: JulianDate) -> MoonMeeusCh47Simplified {
    let jd_tt_val = jd_tt.value();
    let date2 = jd_tt_val - 2_451_545.0;
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

    // Major longitude terms (x1e-6 deg)
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

    let ecl_lon_rad = (lp + sum_l / 1_000_000.0).to_radians();
    let ecl_lat_rad = (sum_b / 1_000_000.0).to_radians();
    let dist_km = 385_000.56 + sum_r / 1_000.0;

    // Ecliptic -> equatorial with IAU 2006 mean obliquity
    let eps = precession_iau2006::mean_obliquity_iau2006(jd_tt).value();
    let ce = eps.cos();
    let se = eps.sin();
    let ra = (ecl_lon_rad.sin() * ce - ecl_lat_rad.tan() * se).atan2(ecl_lon_rad.cos());
    let ra_rad = wrap_0_2pi(ra);
    let dec_rad = (ecl_lat_rad.sin() * ce + ecl_lat_rad.cos() * se * ecl_lon_rad.sin())
        .clamp(-1.0, 1.0)
        .asin();

    MoonMeeusCh47Simplified {
        ra_rad,
        dec_rad,
        dist_km,
        ecl_lon_rad,
        ecl_lat_rad,
    }
}
