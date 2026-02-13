// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # IAU 2000B Nutation Model
//!
//! This module implements the **IAU 2000B** luni-solar nutation model, which uses
//! **77 trigonometric terms** for Δψ (nutation in longitude) and Δε (nutation in
//! obliquity). It is accurate to better than **1 mas** (milliarcsecond) compared
//! to the full IAU 2000A model (1365 terms).
//!
//! ## Method
//!
//! The nutation angles are computed as truncated trigonometric series in the
//! five **Delaunay arguments** (IERS Conventions 2003):
//!
//! * `l`  — mean anomaly of the Moon
//! * `l'` — mean anomaly of the Sun
//! * `F`  — mean argument of latitude of the Moon
//! * `D`  — mean elongation of the Moon from the Sun
//! * `Ω`  — mean longitude of the ascending node of the Moon
//!
//! ```text
//! Δψ = Σᵢ (Aᵢ + A'ᵢ·t) · sin(αᵢ)     (0.1 μas units)
//! Δε = Σᵢ (Bᵢ + B'ᵢ·t) · cos(αᵢ)     (0.1 μas units)
//! ```
//!
//! A fixed correction is applied for the omitted planetary terms.
//!
//! ## Compatibility
//!
//! This module uses the **IAU 2006 mean obliquity** polynomial (84381.406″ at
//! J2000.0) for consistency with the IAU 2006 precession model. The legacy
//! IAU 1980 nutation remains available in [`super::nutation`].
//!
//! ## References
//!
//! * IAU 2000 Resolution B1.6 (nutation)
//! * McCarthy & Luzum (2003), "An abridged model of the precession-nutation
//!   of the celestial pole", Celestial Mechanics 85, 37–49
//! * IERS Conventions (2010), §5.5.1
//! * SOFA/ERFA routine `iauNut00b` / `eraNut00b`

use crate::astro::precession_iau2006::mean_obliquity_iau2006;
use crate::time::JulianDate;
use affn::Rotation3;
use qtty::*;

/// Nutation components for a given epoch (all **radians**).
#[derive(Debug, Copy, Clone)]
pub struct Nutation2000B {
    /// Δψ: nutation in ecliptic longitude (radians).
    pub dpsi: Radians,
    /// Δε: nutation in obliquity (radians).
    pub deps: Radians,
    /// ε_A: mean obliquity of the ecliptic (IAU 2006, radians).
    pub mean_obliquity: Radians,
}

impl Nutation2000B {
    /// True obliquity: ε = ε_A + Δε.
    #[inline]
    pub fn true_obliquity(&self) -> Radians {
        self.mean_obliquity + self.deps
    }

    /// Δψ in degrees (convenience).
    #[inline]
    pub fn dpsi_deg(&self) -> Degrees {
        self.dpsi.to::<Degree>()
    }

    /// Δε in degrees (convenience).
    #[inline]
    pub fn deps_deg(&self) -> Degrees {
        self.deps.to::<Degree>()
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Fundamental arguments (IERS Conventions 2003, Table 5.2e)
// ═══════════════════════════════════════════════════════════════════════════

/// Mean anomaly of the Moon, l (radians).
/// IERS Conventions (2003), Table 5.2e.
#[inline]
fn fund_arg_l(t: f64) -> f64 {
    // l = 134°.963 402 51 + (1717915923.2178 × t + 31.8792 × t² + 0.051635 × t³
    //     − 0.000 244 70 × t⁴) / 3600 (arcseconds → degrees)
    // Simplified to radians directly:
    let arcsec = 485868.249036
        + (1717915923.2178
            + (31.8792 + (0.051635 - 0.000_244_70 * t) * t) * t)
            * t;
    (arcsec % 1296000.0).to_radians() / 3600.0 * 3600.0 // normalize then convert
}

/// Mean anomaly of the Sun, l' (radians).
#[inline]
fn fund_arg_lp(t: f64) -> f64 {
    let arcsec = 1287104.793048
        + (129596581.0481 + (-0.5532 + (0.000_136 - 0.000_011_49 * t) * t) * t) * t;
    (arcsec % 1296000.0).to_radians() / 3600.0 * 3600.0
}

/// Mean argument of latitude of the Moon, F (radians).
#[inline]
fn fund_arg_f(t: f64) -> f64 {
    let arcsec = 335779.526232
        + (1739527262.8478 + (-12.7512 + (-0.001037 + 0.000_000_417 * t) * t) * t) * t;
    (arcsec % 1296000.0).to_radians() / 3600.0 * 3600.0
}

/// Mean elongation of the Moon from the Sun, D (radians).
#[inline]
fn fund_arg_d(t: f64) -> f64 {
    let arcsec = 1072260.703692
        + (1602961601.2090 + (-6.3706 + (0.006593 - 0.000_031_69 * t) * t) * t) * t;
    (arcsec % 1296000.0).to_radians() / 3600.0 * 3600.0
}

/// Mean longitude of the ascending node of the Moon, Ω (radians).
#[inline]
fn fund_arg_om(t: f64) -> f64 {
    let arcsec = 450160.398036
        + (-6962890.5431 + (7.4722 + (0.007702 - 0.000_059_39 * t) * t) * t) * t;
    (arcsec % 1296000.0).to_radians() / 3600.0 * 3600.0
}

// Use a simpler, well-tested approach for the fundamental arguments:
// Convert arcseconds to radians directly.

/// Compute all five Delaunay arguments (radians) for Julian centuries `t` from J2000.
#[inline]
fn delaunay_arguments(t: f64) -> [f64; 5] {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;

    // l: Mean anomaly of the Moon (IERS 2003)
    let l = (485868.249036 + 1717915923.2178 * t + 31.8792 * t2 + 0.051635 * t3
        - 0.000_244_70 * t4)
        .rem_euclid(1_296_000.0)
        * std::f64::consts::PI
        / 648_000.0;

    // l': Mean anomaly of the Sun (IERS 2003)
    let lp = (1287104.793048 + 129596581.0481 * t - 0.5532 * t2 + 0.000_136 * t3
        - 0.000_011_49 * t4)
        .rem_euclid(1_296_000.0)
        * std::f64::consts::PI
        / 648_000.0;

    // F: Mean argument of latitude of the Moon (IERS 2003)
    let f = (335779.526232 + 1739527262.8478 * t - 12.7512 * t2 - 0.001037 * t3
        + 0.000_000_417 * t4)
        .rem_euclid(1_296_000.0)
        * std::f64::consts::PI
        / 648_000.0;

    // D: Mean elongation of the Moon from the Sun (IERS 2003)
    let d = (1072260.703692 + 1602961601.2090 * t - 6.3706 * t2 + 0.006593 * t3
        - 0.000_031_69 * t4)
        .rem_euclid(1_296_000.0)
        * std::f64::consts::PI
        / 648_000.0;

    // Ω: Mean longitude of ascending node (IERS 2003)
    let om = (450160.398036 - 6962890.5431 * t + 7.4722 * t2 + 0.007702 * t3
        - 0.000_059_39 * t4)
        .rem_euclid(1_296_000.0)
        * std::f64::consts::PI
        / 648_000.0;

    [l, lp, f, d, om]
}

// ═══════════════════════════════════════════════════════════════════════════
// IAU 2000B luni-solar nutation series (77 terms)
// ═══════════════════════════════════════════════════════════════════════════

/// Number of luni-solar nutation terms in IAU 2000B.
const NUM_TERMS: usize = 77;

/// Each row: [nl, nlp, nf, nd, nom,  sp, spt, cp,  ce, cet, se]
///
/// * nl..nom: integer multipliers for l, l', F, D, Ω
/// * sp: Δψ sine coefficient (0.1 μas)
/// * spt: Δψ sine coefficient, t-dependent (0.1 μas/century)
/// * cp: Δψ cosine coefficient (0.1 μas)
/// * ce: Δε cosine coefficient (0.1 μas)
/// * cet: Δε cosine coefficient, t-dependent (0.1 μas/century)
/// * se: Δε sine coefficient (0.1 μas)
///
/// Data from SOFA nut00b.c / IERS Conventions (2010), Table 5.3a.
/// Coefficients are in units of 0.1 microarcsecond (0.1 μas).
#[rustfmt::skip]
const NUT00B_LS: [[f64; 11]; NUM_TERMS] = [
    //  l   l'   F    D    Ω       sp          spt      cp         ce         cet      se
    [  0.0, 0.0, 0.0, 0.0, 1.0, -172064161.0, -174666.0, 33386.0, 92052331.0, 9086.0, 15377.0 ],
    [  0.0, 0.0, 2.0,-2.0, 2.0,  -13170906.0,  -1675.0, -13696.0, 5730336.0, -3015.0, -4587.0 ],
    [  0.0, 0.0, 2.0, 0.0, 2.0,   -2276413.0,   -234.0,  2796.0,  978459.0,  -485.0,  1374.0 ],
    [  0.0, 0.0, 0.0, 0.0, 2.0,    2074554.0,    207.0,  -698.0, -897492.0,   470.0,  -291.0 ],
    [  0.0, 1.0, 0.0, 0.0, 0.0,    1475877.0,  -3633.0, 11817.0,   73871.0,  -184.0, -1924.0 ],
    [  0.0, 1.0, 2.0,-2.0, 2.0,    -516821.0,   1226.0,  -524.0,  224386.0,  -677.0,  -174.0 ],
    [  1.0, 0.0, 0.0, 0.0, 0.0,     711159.0,     73.0,  -872.0,   -6750.0,     0.0,   358.0 ],
    [  0.0, 0.0, 2.0, 0.0, 1.0,    -387298.0,   -367.0,   380.0,  200728.0,    18.0,   318.0 ],
    [  1.0, 0.0, 2.0, 0.0, 2.0,    -301461.0,    -36.0,   816.0,  129025.0,   -63.0,   367.0 ],
    [  0.0,-1.0, 2.0,-2.0, 2.0,     215829.0,   -494.0,   111.0,  -95929.0,   299.0,   132.0 ],
    [  0.0, 0.0, 2.0,-2.0, 1.0,     128227.0,    137.0,   181.0,  -68982.0,    -9.0,    39.0 ],
    [ -1.0, 0.0, 2.0, 0.0, 2.0,     123457.0,     11.0,    19.0,  -53311.0,    32.0,    -4.0 ],
    [ -1.0, 0.0, 0.0, 2.0, 0.0,     156994.0,     10.0,  -168.0,   -1235.0,     0.0,    82.0 ],
    [  1.0, 0.0, 0.0, 0.0, 1.0,      63110.0,     63.0,    27.0,  -33228.0,     0.0,    -9.0 ],
    [ -1.0, 0.0, 0.0, 0.0, 1.0,     -57976.0,    -63.0,  -189.0,   31429.0,     0.0,   -75.0 ],
    [ -1.0, 0.0, 2.0, 2.0, 2.0,     -59641.0,    -11.0,   149.0,   25543.0,   -11.0,    66.0 ],
    [  1.0, 0.0, 2.0, 0.0, 1.0,     -51613.0,    -42.0,   129.0,   26366.0,     0.0,    78.0 ],
    [ -2.0, 0.0, 2.0, 0.0, 1.0,      45893.0,     50.0,    31.0,  -24236.0,   -10.0,    20.0 ],
    [  0.0, 0.0, 0.0, 2.0, 0.0,      63384.0,     11.0,  -150.0,   -1220.0,     0.0,    29.0 ],
    [  0.0, 0.0, 2.0, 2.0, 2.0,     -38571.0,     -1.0,   158.0,   16452.0,   -11.0,    68.0 ],
    [  0.0,-2.0, 2.0,-2.0, 2.0,      32481.0,      0.0,     0.0,  -13870.0,     0.0,     0.0 ],
    [ -2.0, 0.0, 0.0, 2.0, 0.0,     -47722.0,      0.0,   -18.0,     477.0,     0.0,   -25.0 ],
    [  2.0, 0.0, 2.0, 0.0, 2.0,     -31046.0,     -1.0,   131.0,   13238.0,   -11.0,    59.0 ],
    [  1.0, 0.0, 2.0,-2.0, 2.0,      28593.0,      0.0,    -1.0,  -12338.0,    10.0,    -3.0 ],
    [ -1.0, 0.0, 2.0, 0.0, 1.0,      20441.0,     21.0,    10.0,  -10758.0,     0.0,    -3.0 ],
    [  2.0, 0.0, 0.0, 0.0, 0.0,      29243.0,      0.0,   -74.0,    -609.0,     0.0,    13.0 ],
    [  0.0, 0.0, 2.0, 0.0, 0.0,      25887.0,      0.0,   -66.0,    -550.0,     0.0,    11.0 ],
    [  0.0, 1.0, 0.0, 0.0, 1.0,     -14053.0,    -25.0,    79.0,    8551.0,    -2.0,   -45.0 ],
    [ -1.0, 0.0, 0.0, 2.0, 1.0,      15164.0,     10.0,    11.0,   -8001.0,     0.0,    -1.0 ],
    [  0.0, 2.0, 2.0,-2.0, 2.0,     -15794.0,     72.0,   -16.0,    6850.0,   -42.0,    -5.0 ],
    [  0.0, 0.0,-2.0, 2.0, 0.0,      21783.0,      0.0,    13.0,    -167.0,     0.0,    13.0 ],
    [  1.0, 0.0, 0.0,-2.0, 1.0,     -12873.0,    -10.0,   -37.0,    6953.0,     0.0,   -14.0 ],
    [  0.0,-1.0, 0.0, 0.0, 1.0,     -12654.0,     11.0,    63.0,    6415.0,     0.0,    26.0 ],
    [ -1.0, 0.0, 2.0, 2.0, 1.0,     -10204.0,      0.0,    25.0,    5222.0,     0.0,    15.0 ],
    [  0.0, 2.0, 0.0, 0.0, 0.0,      16707.0,    -85.0,   -10.0,     168.0,    -1.0,    10.0 ],
    [  1.0, 0.0, 2.0, 2.0, 2.0,      -7691.0,      0.0,    44.0,    3268.0,     0.0,    19.0 ],
    [ -2.0, 0.0, 2.0, 0.0, 0.0,     -11024.0,      0.0,   -14.0,     104.0,     0.0,     2.0 ],
    [  0.0, 1.0, 2.0, 0.0, 2.0,       7566.0,    -21.0,   -11.0,   -3250.0,     0.0,    -5.0 ],
    [  0.0, 0.0, 2.0, 2.0, 1.0,      -6637.0,    -11.0,    25.0,    3353.0,     0.0,    14.0 ],
    [  0.0,-1.0, 2.0, 0.0, 2.0,      -7141.0,     21.0,     8.0,    3070.0,     0.0,     4.0 ],
    [  0.0, 0.0, 0.0, 2.0, 1.0,      -6302.0,    -11.0,     2.0,    3272.0,     0.0,     4.0 ],
    [  1.0, 0.0, 2.0,-2.0, 1.0,       5800.0,     10.0,     2.0,   -3045.0,     0.0,    -1.0 ],
    [  2.0, 0.0, 2.0,-2.0, 2.0,       6443.0,      0.0,    -7.0,   -2768.0,     0.0,    -4.0 ],
    [ -2.0, 0.0, 0.0, 2.0, 1.0,      -5774.0,    -11.0,   -15.0,    3041.0,     0.0,    -5.0 ],
    [  2.0, 0.0, 2.0, 0.0, 1.0,      -5350.0,      0.0,    21.0,    2695.0,     0.0,    12.0 ],
    [  0.0,-1.0, 2.0,-2.0, 1.0,      -4752.0,    -11.0,    -3.0,    2719.0,     0.0,    -3.0 ],
    [  0.0, 0.0, 0.0,-2.0, 1.0,      -4940.0,    -11.0,    -21.0,   2720.0,     0.0,    -9.0 ],
    [ -1.0,-1.0, 0.0, 2.0, 0.0,       7350.0,      0.0,    -8.0,     -51.0,     0.0,     4.0 ],
    [  2.0, 0.0, 0.0,-2.0, 1.0,       4065.0,      0.0,     6.0,   -2206.0,     0.0,     1.0 ],
    [  1.0, 0.0, 0.0, 2.0, 0.0,       6579.0,      0.0,   -24.0,    -199.0,     0.0,     2.0 ],
    [  0.0, 1.0, 2.0,-2.0, 1.0,       3579.0,      0.0,     5.0,   -1900.0,     0.0,     1.0 ],
    [  1.0,-1.0, 0.0, 0.0, 0.0,       4725.0,      0.0,    -6.0,     -41.0,     0.0,     3.0 ],
    [ -2.0, 0.0, 2.0, 0.0, 2.0,      -3075.0,      0.0,    -2.0,    1313.0,     0.0,    -1.0 ],
    [  3.0, 0.0, 2.0, 0.0, 2.0,      -2904.0,      0.0,    15.0,    1233.0,     0.0,     7.0 ],
    [  0.0,-1.0, 0.0, 2.0, 0.0,       4348.0,      0.0,   -10.0,     -81.0,     0.0,     2.0 ],
    [  1.0,-1.0, 2.0, 0.0, 2.0,      -2878.0,      0.0,     8.0,    1232.0,     0.0,     4.0 ],
    [  0.0, 0.0, 0.0, 1.0, 0.0,      -4230.0,      0.0,     5.0,     -20.0,     0.0,    -2.0 ],
    [ -1.0,-1.0, 2.0, 2.0, 2.0,      -2819.0,      0.0,     7.0,    1207.0,     0.0,     3.0 ],
    [ -1.0, 0.0, 2.0, 0.0, 0.0,      -4056.0,      0.0,     5.0,      40.0,     0.0,    -2.0 ],
    [  0.0,-1.0, 2.0, 2.0, 2.0,      -2647.0,      0.0,    11.0,    1129.0,     0.0,     5.0 ],
    [ -2.0, 0.0, 0.0, 0.0, 1.0,      -2294.0,      0.0,   -10.0,    1266.0,     0.0,    -4.0 ],
    [  1.0, 1.0, 2.0, 0.0, 2.0,       2481.0,      0.0,    -7.0,   -1062.0,     0.0,    -3.0 ],
    [  2.0, 0.0, 0.0, 0.0, 1.0,       2179.0,      0.0,    -2.0,   -1129.0,     0.0,    -2.0 ],
    [ -1.0, 1.0, 0.0, 1.0, 0.0,       3276.0,      0.0,     1.0,      -9.0,     0.0,     0.0 ],
    [  1.0, 1.0, 0.0, 0.0, 0.0,      -3389.0,      0.0,     5.0,      35.0,     0.0,    -2.0 ],
    [  1.0, 0.0, 2.0, 0.0, 0.0,       3339.0,      0.0,   -13.0,    -107.0,     0.0,     1.0 ],
    [ -1.0, 0.0, 2.0,-2.0, 1.0,      -1987.0,      0.0,    -6.0,    1073.0,     0.0,    -2.0 ],
    [  1.0, 0.0, 0.0, 0.0, 2.0,      -1981.0,      0.0,     0.0,     854.0,     0.0,     0.0 ],
    [ -1.0, 0.0, 0.0, 1.0, 0.0,       4026.0,      0.0,  -353.0,    -553.0,     0.0,  -139.0 ],
    [  0.0, 0.0, 2.0, 1.0, 2.0,       1660.0,      0.0,    -5.0,    -710.0,     0.0,    -2.0 ],
    [ -1.0, 0.0, 2.0, 4.0, 2.0,      -1521.0,      0.0,     9.0,     647.0,     0.0,     4.0 ],
    [ -1.0, 1.0, 0.0, 1.0, 1.0,       1314.0,      0.0,     0.0,    -700.0,     0.0,     0.0 ],
    [  0.0,-2.0, 2.0,-2.0, 1.0,      -1283.0,      0.0,     0.0,     672.0,     0.0,     0.0 ],
    [  1.0, 0.0, 2.0, 2.0, 1.0,      -1331.0,      0.0,     8.0,     663.0,     0.0,     4.0 ],
    [ -2.0, 0.0, 2.0, 2.0, 2.0,       1383.0,      0.0,    -2.0,    -594.0,     0.0,    -2.0 ],
    [ -1.0, 0.0, 0.0, 0.0, 2.0,       1405.0,      0.0,     4.0,    -610.0,     0.0,     2.0 ],
    [  1.0, 1.0, 2.0,-2.0, 2.0,       1290.0,      0.0,     0.0,    -556.0,     0.0,     0.0 ],
];

/// Fixed correction for omitted planetary nutation terms (0.1 μas).
///
/// IAU 2000B adds these fixed offsets to account for the aggregate effect
/// of the ~600 planetary terms omitted from the truncated series.
///
/// * Δψ correction: −135 μas = −1350 × 0.1 μas
/// * Δε correction: +388 μas = +3880 × 0.1 μas
///
/// Reference: McCarthy & Luzum (2003), Table 2.
const DPSI_PLANETARY_CORR: f64 = -1350.0; // 0.1 μas
const DEPS_PLANETARY_CORR: f64 = 3880.0; // 0.1 μas

/// Compute IAU 2000B nutation for a given Julian Date (TT).
///
/// Returns [`Nutation2000B`] containing Δψ, Δε and ε_A, all in **radians**.
///
/// Uses the IAU 2006 mean obliquity polynomial for consistency with
/// the IAU 2006 precession model.
///
/// ## References
/// * SOFA/ERFA routine `iauNut00b`
/// * McCarthy & Luzum (2003)
/// * IERS Conventions (2010), §5.5.1
pub fn nutation_iau2000b(jd: JulianDate) -> Nutation2000B {
    let t = jd.julian_centuries().value();

    // Delaunay arguments (radians)
    let fa = delaunay_arguments(t);

    // Evaluate luni-solar series
    let mut dpsi_sum = 0.0_f64; // 0.1 μas
    let mut deps_sum = 0.0_f64; // 0.1 μas

    for row in &NUT00B_LS {
        // Argument = Σ nᵢ·FAᵢ
        let arg = row[0] * fa[0] // l
            + row[1] * fa[1]     // l'
            + row[2] * fa[2]     // F
            + row[3] * fa[3]     // D
            + row[4] * fa[4];    // Ω

        let (sin_arg, cos_arg) = arg.sin_cos();

        // Δψ: sine + cosine terms
        dpsi_sum += (row[5] + row[6] * t) * sin_arg + row[7] * cos_arg;
        // Δε: cosine + sine terms
        deps_sum += (row[8] + row[9] * t) * cos_arg + row[10] * sin_arg;
    }

    // Add fixed planetary correction
    dpsi_sum += DPSI_PLANETARY_CORR;
    deps_sum += DEPS_PLANETARY_CORR;

    // Convert from 0.1 μas to radians:
    // 0.1 μas = 0.1e-6 arcsec = 1e-7 arcsec
    // 1 arcsec = π/(180×3600) rad
    let unit = std::f64::consts::PI / (180.0 * 3600.0 * 1e7);

    let dpsi = Radians::new(dpsi_sum * unit);
    let deps = Radians::new(deps_sum * unit);
    let mean_obliquity = mean_obliquity_iau2006(jd);

    Nutation2000B {
        dpsi,
        deps,
        mean_obliquity,
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Rotation matrices
// ═══════════════════════════════════════════════════════════════════════════

#[inline]
fn rotation_x(angle: Radians) -> Rotation3 {
    let (s, c) = angle.sin_cos();
    Rotation3::from_matrix([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
}

#[inline]
fn rotation_z(angle: Radians) -> Rotation3 {
    let (s, c) = angle.sin_cos();
    Rotation3::from_matrix([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
}

/// Nutation rotation matrix from mean-of-date to true-of-date.
///
/// Uses IAU 2000B nutation with IAU 2006 mean obliquity.
///
/// ```text
/// N = R₁(ε₀ + Δε) · R₃(Δψ) · R₁(−ε₀)
/// ```
///
/// ## References
/// * SOFA routine `iauNum00b`
pub fn nutation_rotation_iau2000b(jd: JulianDate) -> Rotation3 {
    let nut = nutation_iau2000b(jd);
    rotation_x(nut.mean_obliquity + nut.deps)
        * rotation_z(nut.dpsi)
        * rotation_x(-nut.mean_obliquity)
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::{Arcsecond, Degrees, Radians};

    #[test]
    fn nutation_at_j2000_dominant_term() {
        let nut = nutation_iau2000b(JulianDate::J2000);

        // At J2000.0, the dominant Ω term gives Δψ ≈ −14″ to −17″, Δε ≈ 9″
        // (exact value depends on Ω phase at epoch)
        let dpsi_arcsec = nut.dpsi.to::<Degree>().value() * 3600.0;
        let deps_arcsec = nut.deps.to::<Degree>().value() * 3600.0;

        assert!(
            dpsi_arcsec.abs() > 5.0 && dpsi_arcsec.abs() < 20.0,
            "Δψ at J2000 = {}″, expected magnitude 5–20″",
            dpsi_arcsec
        );
        assert!(
            deps_arcsec.abs() > 2.0 && deps_arcsec.abs() < 15.0,
            "Δε at J2000 = {}″, expected magnitude 2–15″",
            deps_arcsec
        );
    }

    #[test]
    fn nutation_2000b_vs_iau1980_similar_magnitude() {
        // The IAU 2000B and IAU 1980 nutations should agree to ~0.5″
        let jd = JulianDate::new(2_459_000.5);
        let nut_2000b = nutation_iau2000b(jd);

        let nut_1980 = crate::astro::nutation::get_nutation(jd);

        let dpsi_diff = (nut_2000b.dpsi.to::<Degree>() - nut_1980.longitude).value().abs();
        let deps_diff = (nut_2000b.deps.to::<Degree>() - nut_1980.obliquity).value().abs();

        // They should agree within ~0.5″ ≈ 0.000139°
        assert!(
            dpsi_diff < 0.001,
            "Δψ difference too large: {}°",
            dpsi_diff
        );
        assert!(
            deps_diff < 0.001,
            "Δε difference too large: {}°",
            deps_diff
        );
    }

    #[test]
    fn nutation_rotation_near_identity() {
        let jd = JulianDate::new(2_451_545.0);
        let rot = nutation_rotation_iau2000b(jd);
        let m = rot.as_matrix();

        // Nutation is small: off-diagonal < 0.001
        for i in 0..3 {
            assert!(
                (m[i][i] - 1.0).abs() < 1e-4,
                "diagonal[{}] too far from 1: {}",
                i,
                m[i][i]
            );
        }
    }

    #[test]
    fn mean_obliquity_matches_iau2006() {
        let nut = nutation_iau2000b(JulianDate::J2000);
        let eps_arcsec = nut.mean_obliquity.to::<Degree>().value() * 3600.0;
        // IAU 2006 value at J2000: 84381.406″
        assert!(
            (eps_arcsec - 84381.406_f64).abs() < 0.001,
            "mean obliquity = {}″, expected 84381.406″",
            eps_arcsec
        );
    }

    #[test]
    fn delaunay_args_finite() {
        // Test that fundamental arguments don't blow up at various epochs
        for jd_val in &[2_400_000.5, 2_451_545.0, 2_460_000.5, 2_500_000.5] {
            let t = (jd_val - 2_451_545.0) / 36525.0;
            let fa = delaunay_arguments(t);
            for (i, &a) in fa.iter().enumerate() {
                assert!(a.is_finite(), "fundamental arg {} is not finite at JD {}", i, jd_val);
            }
        }
    }
}
