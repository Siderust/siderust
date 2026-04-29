// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # IAU 2000A / 2006A Nutation
//!
//! Full-precision nutation model with 678 luni-solar and 687 planetary
//! terms (MHB2000).  The public entry point [`nutation_iau2006a`] applies
//! the IAU 2006 precession adjustments (Wallace & Capitaine 2006, Eqs. 5)
//! on top of the raw IAU 2000A angles.
//!
//! ## References
//!
//! * Mathews, Herring & Buffett (2002), *J. Geophys. Res.* 107, B4
//! * Wallace & Capitaine (2006), *Astron. Astrophys.* 459, 981
//! * SOFA/ERFA routines `eraNut00a`, `eraNut06a`

use super::nut00a_tables::{NUT00A_LS, NUT00A_PL};
use super::NutationAngles;
use crate::astro::precession::mean_obliquity_iau2006;
use crate::time::JulianDate;
use crate::qtty::*;

// ═══════════════════════════════════════════════════════════════════════════
//  Fundamental arguments, planetary longitudes (IERS 2003)
// ═══════════════════════════════════════════════════════════════════════════

/// Mercury mean longitude (IERS 2003), radians.
#[inline]
fn fa_mercury(t: f64) -> f64 {
    (4.402608842 + 2608.7903141574 * t) % std::f64::consts::TAU
}

/// Venus mean longitude (IERS 2003), radians.
#[inline]
fn fa_venus(t: f64) -> f64 {
    (3.176146697 + 1021.3285546211 * t) % std::f64::consts::TAU
}

/// Earth mean longitude (IERS 2003), radians.
#[inline]
fn fa_earth(t: f64) -> f64 {
    (1.753470314 + 628.3075849991 * t) % std::f64::consts::TAU
}

/// Mars mean longitude (IERS 2003), radians.
#[inline]
fn fa_mars(t: f64) -> f64 {
    (6.203480913 + 334.0612426700 * t) % std::f64::consts::TAU
}

/// Jupiter mean longitude (IERS 2003), radians.
#[inline]
fn fa_jupiter(t: f64) -> f64 {
    (0.599546497 + 52.9690962641 * t) % std::f64::consts::TAU
}

/// Saturn mean longitude (IERS 2003), radians.
#[inline]
fn fa_saturn(t: f64) -> f64 {
    (0.874016757 + 21.3299104960 * t) % std::f64::consts::TAU
}

/// Uranus mean longitude (IERS 2003), radians.
#[inline]
fn fa_uranus(t: f64) -> f64 {
    (5.481293872 + 7.4781598567 * t) % std::f64::consts::TAU
}

/// General accumulated precession in longitude (IERS 2003), radians.
#[inline]
fn fa_pa(t: f64) -> f64 {
    (0.024381750 + 0.00000538691 * t) * t
}

// ═══════════════════════════════════════════════════════════════════════════
//  Fundamental arguments, luni-solar (IERS 2003 / MHB2000 mix)
// ═══════════════════════════════════════════════════════════════════════════

const AS2R: f64 = std::f64::consts::PI / (180.0 * 3600.0);
const TURNAS: f64 = 1_296_000.0;

/// Moon mean anomaly l (IERS 2003), radians.
#[inline]
fn fa_l_iers03(t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    ((485868.249036 + 1717915923.2178 * t + 31.8792 * t2 + 0.051635 * t3 - 0.000_244_70 * t4)
        % TURNAS)
        * AS2R
}

/// Sun mean anomaly l' (MHB2000), radians.
///
/// Uses the MHB2000 polynomial, NOT the IERS 2003 form, to match
/// the original `eraNut00a` implementation exactly.
#[inline]
fn fa_lp_mhb(t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    ((1287104.79305 + 129596581.0481 * t - 0.5532 * t2 + 0.000136 * t3 - 0.000_011_49 * t4)
        % TURNAS)
        * AS2R
}

/// Moon argument of latitude F (IERS 2003), radians.
#[inline]
fn fa_f_iers03(t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    ((335779.526232 + 1739527262.8478 * t - 12.7512 * t2 - 0.001037 * t3 + 0.000_000_417 * t4)
        % TURNAS)
        * AS2R
}

/// Moon elongation D (MHB2000), radians.
#[inline]
fn fa_d_mhb(t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    ((1072260.70369 + 1602961601.2090 * t - 6.3706 * t2 + 0.006593 * t3 - 0.000_031_69 * t4)
        % TURNAS)
        * AS2R
}

/// Moon ascending node Ω (IERS 2003), radians.
#[inline]
fn fa_om_iers03(t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    ((450160.398036 - 6962890.5431 * t + 7.4722 * t2 + 0.007702 * t3 - 0.000_059_39 * t4) % TURNAS)
        * AS2R
}

// ═══════════════════════════════════════════════════════════════════════════
//  MHB2000 Delaunay arguments for the planetary section
//  (slightly different from IERS 2003; reproduced faithfully)
// ═══════════════════════════════════════════════════════════════════════════

/// Moon mean anomaly (MHB2000), radians.
#[inline]
fn fa_l_mhb(t: f64) -> f64 {
    (2.35555598 + 8328.6914269554 * t) % std::f64::consts::TAU
}

/// Moon argument of latitude (MHB2000), radians.
#[inline]
fn fa_f_mhb(t: f64) -> f64 {
    (1.627905234 + 8433.466158131 * t) % std::f64::consts::TAU
}

/// Moon elongation (MHB2000), radians.
#[inline]
fn fa_d_mhb_pl(t: f64) -> f64 {
    (5.198466741 + 7771.3771468121 * t) % std::f64::consts::TAU
}

/// Moon ascending node (MHB2000), radians.
#[inline]
fn fa_om_mhb(t: f64) -> f64 {
    (2.18243920 - 33.757045 * t) % std::f64::consts::TAU
}

/// Neptune mean longitude (MHB2000), radians.
#[inline]
fn fa_neptune_mhb(t: f64) -> f64 {
    (5.321159000 + 3.8127774000 * t) % std::f64::consts::TAU
}

// ═══════════════════════════════════════════════════════════════════════════
//  IAU 2000A raw nutation (eraNut00a equivalent)
// ═══════════════════════════════════════════════════════════════════════════

/// 0.1 µas → radians.
const U2R: f64 = AS2R / 1e7;

/// Compute IAU 2000A nutation (MHB2000, 1365 terms).
///
/// Returns `(dpsi, deps)` in **radians**, referenced to the ecliptic of
/// date with the Lieske et al. (1977) obliquity (84381.448″), consistent
/// with the original MHB2000 model.
///
/// This is the raw IAU 2000A result; for use with IAU 2006 precession
/// call [`nutation_iau2006a`] instead, which applies the P03 corrections.
///
/// # Numerical convention: summation order
///
/// Both the luni-solar and the planetary series are summed in **reverse
/// table order** (`.iter().rev()`), so that the smallest-amplitude terms
/// are accumulated first and the dominant `Ω`-term (with `dpsi ≈ −17″`,
/// roughly seven orders of magnitude larger than the smallest planetary
/// terms at ~0.1 µas) is summed last. This matches the SOFA/ERFA
/// convention (`eraNut00a`) and preserves ~1 µas of precision in `Δψ` /
/// `Δε` even at extreme epochs (|t| ≳ 30 centuries) where catastrophic
/// cancellation in a "largest-first" summation would otherwise dominate
/// the IEEE-754 round-off budget.
///
/// Do not reorder the loop without updating the IAU compliance tests in
/// `tests/`.
pub(crate) fn nutation_iau2000a_raw(jd: JulianDate) -> (f64, f64) {
    let t = jd.julian_centuries();

    // ── Luni-solar fundamental arguments ──
    let el = fa_l_iers03(t);
    let elp = fa_lp_mhb(t);
    let f = fa_f_iers03(t);
    let d = fa_d_mhb(t);
    let om = fa_om_iers03(t);

    // ── Luni-solar nutation (678 terms, reverse order) ──
    let mut dp = 0.0_f64;
    let mut de = 0.0_f64;

    for row in NUT00A_LS.iter().rev() {
        let arg = (row[0] * el + row[1] * elp + row[2] * f + row[3] * d + row[4] * om)
            % std::f64::consts::TAU;
        let (sarg, carg) = arg.sin_cos();
        dp += (row[5] + row[6] * t) * sarg + row[7] * carg;
        de += (row[8] + row[9] * t) * carg + row[10] * sarg;
    }

    let dpsils = dp * U2R;
    let depsls = de * U2R;

    // ── Planetary fundamental arguments ──
    let al = fa_l_mhb(t);
    let af = fa_f_mhb(t);
    let ad = fa_d_mhb_pl(t);
    let aom = fa_om_mhb(t);
    let alme = fa_mercury(t);
    let alve = fa_venus(t);
    let alea = fa_earth(t);
    let alma = fa_mars(t);
    let alju = fa_jupiter(t);
    let alsa = fa_saturn(t);
    let alur = fa_uranus(t);
    let alne = fa_neptune_mhb(t);
    let apa = fa_pa(t);

    // ── Planetary nutation (687 terms, reverse order) ──
    dp = 0.0;
    de = 0.0;

    for row in NUT00A_PL.iter().rev() {
        let arg = (f64::from(row[0]) * al
            + f64::from(row[1]) * af
            + f64::from(row[2]) * ad
            + f64::from(row[3]) * aom
            + f64::from(row[4]) * alme
            + f64::from(row[5]) * alve
            + f64::from(row[6]) * alea
            + f64::from(row[7]) * alma
            + f64::from(row[8]) * alju
            + f64::from(row[9]) * alsa
            + f64::from(row[10]) * alur
            + f64::from(row[11]) * alne
            + f64::from(row[12]) * apa)
            % std::f64::consts::TAU;
        let (sarg, carg) = arg.sin_cos();
        dp += f64::from(row[13]) * sarg + f64::from(row[14]) * carg;
        de += f64::from(row[15]) * sarg + f64::from(row[16]) * carg;
    }

    let dpsipl = dp * U2R;
    let depspl = de * U2R;

    (dpsils + dpsipl, depsls + depspl)
}

// ═══════════════════════════════════════════════════════════════════════════
//  IAU 2006A nutation (eraNut06a equivalent)
// ═══════════════════════════════════════════════════════════════════════════

/// Compute IAU 2006/2000A nutation.
///
/// Applies IAU 2006 precession corrections to the raw IAU 2000A model:
///   * Ecliptic obliquity adjustment (P03 constant 0.4697×10⁻⁶)
///   * Secular variation in Earth's dynamical form factor J₂
///
/// Reference: Wallace & Capitaine (2006), Eqs. 5.
pub(crate) fn nutation_iau2006a(jd: JulianDate) -> NutationAngles {
    let t = jd.julian_centuries();
    let fj2 = -2.7774e-6 * t;

    let (dp, de) = nutation_iau2000a_raw(jd);

    let dpsi = dp + dp * (0.4697e-6 + fj2);
    let deps = de + de * fj2;

    NutationAngles {
        dpsi: Radians::new(dpsi),
        deps: Radians::new(deps),
        mean_obliquity: mean_obliquity_iau2006(jd),
    }
}

/// Compute pure IAU 2000A nutation (without IAU 2006 P03/J2 corrections).
///
/// This is useful when callers need an explicit "IAU 2000A" model selection
/// separate from the combined IAU 2006A convention.
pub(crate) fn nutation_iau2000a(jd: JulianDate) -> NutationAngles {
    let (dpsi, deps) = nutation_iau2000a_raw(jd);
    NutationAngles {
        dpsi: Radians::new(dpsi),
        deps: Radians::new(deps),
        mean_obliquity: mean_obliquity_iau2006(jd),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// At J2000.0, IAU 2000A and 2000B should agree to within 1 mas.
    #[test]
    fn nut2000a_agrees_with_2000b_at_j2000() {
        let jd = JulianDate::J2000;
        let n2a = nutation_iau2006a(jd);
        let n2b = crate::astro::nutation::nutation_iau2000b(jd);

        let dpsi_diff = (n2a.dpsi - n2b.dpsi).value().abs();
        let deps_diff = (n2a.deps - n2b.deps).value().abs();

        // IAU 2000B is accurate to ~1 mas relative to 2000A
        let one_mas = std::f64::consts::PI / (180.0 * 3600.0 * 1000.0);
        assert!(
            dpsi_diff < one_mas,
            "Δψ diff = {:.3e} rad > 1 mas ({:.3e})",
            dpsi_diff,
            one_mas
        );
        assert!(
            deps_diff < one_mas,
            "Δε diff = {:.3e} rad > 1 mas ({:.3e})",
            deps_diff,
            one_mas
        );
    }

    /// At a future epoch (J2020), 2000A and 2000B still agree to ~1 mas.
    #[test]
    fn nut2000a_agrees_with_2000b_at_j2020() {
        let jd = JulianDate::new(2_458_849.5); // 2020-01-01
        let n2a = nutation_iau2006a(jd);
        let n2b = crate::astro::nutation::nutation_iau2000b(jd);

        let dpsi_diff = (n2a.dpsi - n2b.dpsi).value().abs();
        let deps_diff = (n2a.deps - n2b.deps).value().abs();

        let one_mas = std::f64::consts::PI / (180.0 * 3600.0 * 1000.0);
        assert!(
            dpsi_diff < one_mas,
            "Δψ diff at J2020 = {:.3e} rad > 1 mas",
            dpsi_diff
        );
        assert!(
            deps_diff < one_mas,
            "Δε diff at J2020 = {:.3e} rad > 1 mas",
            deps_diff
        );
    }

    /// Cross-check against Python/ERFA nut06a at a reference epoch.
    /// Reference: erfa.nut06a(2458849.5, 0.0), 2020-01-01T12:00 TT
    #[test]
    fn nut2006a_matches_erfa_at_j2020() {
        let jd = JulianDate::new(2_458_849.5);
        let n = nutation_iau2006a(jd);

        let erfa_dpsi: f64 = -7.996558234083069e-05; // rad
        let erfa_deps: f64 = -8.251412879483328e-06; // rad

        let one_uas = std::f64::consts::PI / (180.0 * 3600.0 * 1e6); // 1 µas in rad
        let dpsi_diff = (n.dpsi.value() - erfa_dpsi).abs();
        let deps_diff = (n.deps.value() - erfa_deps).abs();

        assert!(
            dpsi_diff < one_uas,
            "Δψ diff vs ERFA = {:.3e} rad ({:.3} µas)",
            dpsi_diff,
            dpsi_diff / one_uas
        );
        assert!(
            deps_diff < one_uas,
            "Δε diff vs ERFA = {:.3e} rad ({:.3} µas)",
            deps_diff,
            deps_diff / one_uas
        );
    }

    /// The P03 J₂ correction should be tiny for near-current epochs.
    #[test]
    fn p03_correction_is_small() {
        let jd = JulianDate::J2000;
        let (dp_raw, de_raw) = nutation_iau2000a_raw(jd);
        let n06a = nutation_iau2006a(jd);

        // At J2000.0, t=0, fj2=0 so only the 0.4697e-6 factor on dpsi differs.
        let dpsi_corr = (n06a.dpsi.value() - dp_raw).abs();
        let deps_corr = (n06a.deps.value() - de_raw).abs();

        // Correction to dpsi should be dp_raw * 0.4697e-6 ≈ very small
        assert!(
            dpsi_corr < 1e-9,
            "dpsi P03 correction too large: {:.3e}",
            dpsi_corr
        );
        // deps correction at t=0 should be zero (fj2=0)
        assert!(
            deps_corr < 1e-15,
            "deps P03 correction should be zero at t=0: {:.3e}",
            deps_corr
        );
    }
}
