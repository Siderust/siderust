// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Moon Ephemeris Cache — Chebyshev Segment Interpolation
//!
//! Precomputes the Moon's geocentric ecliptic Cartesian coordinates (X, Y, Z)
//! on Chebyshev nodes within fixed-duration segments, then evaluates any
//! intermediate time via Clenshaw recurrence in **O(degree)** — replacing the
//! full ELP2000 series summation (~200 µs → ~1 µs per query).
//!
//! A companion [`NutationCache`] stores pre-evaluated nutation triplets
//! (Δψ, Δε, ε₀) at regular intervals and linearly interpolates, removing
//! the 63-term IAU 1980 series from the per-query hot path.
//!
//! ## Design
//!
//! | Parameter | Default | Notes |
//! |-----------|---------|-------|
//! | Segment duration | 4 days | Moon moves ≈52° in ecliptic longitude |
//! | Chebyshev degree | 8 | 9 nodes per segment, sub-arcsecond accuracy |
//! | Nutation step | 2 hours | Linear interpolation error < 0.001″ |
//!
//! ## Accuracy
//!
//! The Moon's shortest significant perturbation period is ~5 days (evection).
//! Chebyshev degree 8 over 4-day segments yields interpolation errors well
//! below 1 arcsecond in geocentric position — far smaller than the ~0.5°
//! atmospheric refraction uncertainty at the horizon.

use crate::astro::nutation::get_nutation;
use crate::astro::precession::precession_rotation_from_j2000;
use crate::astro::sidereal::{calculate_gst, calculate_lst, unmodded_gst};
use crate::bodies::solar_system::Moon;
use crate::coordinates::centers::ObserverSite;
use qtty::*;

// =============================================================================
// Constants
// =============================================================================

/// Chebyshev polynomial degree for position interpolation.
const CHEB_DEGREE: usize = 8;

/// Number of Chebyshev nodes per segment (degree + 1).
const CHEB_NODES: usize = CHEB_DEGREE + 1;

/// Duration of each Chebyshev segment in days.
const SEGMENT_DAYS: Days = Days::new(4.0);

/// J2000 mean obliquity ε₀ (IAU 2006): 84381.406″ converted to radians.
/// Used for ecliptic → equatorial rotation (constant for J2000 frame).
const J2000_OBLIQUITY_RAD: qtty::Quantity<Radian> =
    qtty::Quantity::<Radian>::new(84381.406 / 3600.0 * std::f64::consts::PI / 180.0);

/// Nutation cache step in days (2 hours).
const NUT_STEP_DAYS: Days = Hours::new(2.0).to_const::<Day>();

// =============================================================================
// Chebyshev nodes (precomputed for CHEB_NODES = 9)
// =============================================================================

/// Compute Chebyshev nodes on [-1, 1] at runtime.
#[inline]
fn cheb_nodes_on_unit() -> [f64; CHEB_NODES] {
    let mut nodes = [0.0_f64; CHEB_NODES];
    let n = CHEB_NODES as f64;
    for (k, node) in nodes.iter_mut().enumerate().take(CHEB_NODES) {
        let arg = std::f64::consts::PI * (2.0 * k as f64 + 1.0) / (2.0 * n);
        *node = arg.cos();
    }
    nodes
}

// =============================================================================
// Chebyshev coefficient computation & evaluation
// =============================================================================

/// Compute Chebyshev coefficients c_0..c_{n-1} from function values at
/// the n Chebyshev nodes.
///
/// Uses the DCT-like formula:
///   c_j = (2/n) Σ_{k=0}^{n-1} f(x_k) · cos(j·π·(2k+1)/(2n))
/// with c_0 scaled by 1/n instead of 2/n.
#[inline]
fn compute_cheb_coeffs(values: &[Kilometers; CHEB_NODES]) -> [Kilometers; CHEB_NODES] {
    let mut coeffs = [Kilometers::zero(); CHEB_NODES];
    let n = CHEB_NODES as f64;

    for (j, coeff) in coeffs.iter_mut().enumerate().take(CHEB_NODES) {
        let mut sum = Kilometers::zero();
        for (k, value) in values.iter().enumerate().take(CHEB_NODES) {
            let arg = std::f64::consts::PI * (j as f64) * (2.0 * k as f64 + 1.0) / (2.0 * n);
            sum += *value * arg.cos();
        }
        *coeff = if j == 0 { sum / n } else { 2.0 * sum / n };
    }
    coeffs
}

/// Evaluate a Chebyshev expansion at x ∈ [-1, 1] via Clenshaw recurrence.
///
/// f(x) ≈ c_0 + x·d_1 − d_2
/// where d_k = c_k + 2x·d_{k+1} − d_{k+2}, k = n-1, ..., 1.
#[inline]
fn clenshaw_eval(coeffs: &[Kilometers; CHEB_NODES], x: f64) -> Kilometers {
    let two_x = 2.0 * x;
    let mut d_k1 = Kilometers::zero(); // d_{k+1}
    let mut d_k2 = Kilometers::zero(); // d_{k+2}

    // k = CHEB_DEGREE down to 1
    let mut k = CHEB_DEGREE;
    while k >= 1 {
        let d_k = coeffs[k] + two_x * d_k1 - d_k2;
        d_k2 = d_k1;
        d_k1 = d_k;
        k -= 1;
    }

    coeffs[0] + x * d_k1 - d_k2
}

// =============================================================================
// MoonPositionCache
// =============================================================================

/// Chebyshev interpolation cache for the Moon's geocentric ecliptic
/// Cartesian coordinates (X, Y, Z) in kilometers.
///
/// The time domain is divided into segments of [`SEGMENT_DAYS`] days.
/// Within each segment the three coordinates are approximated by degree-8
/// Chebyshev polynomials fitted at the canonical nodes.
pub struct MoonPositionCache {
    /// Modified Julian Date of the first segment's start.
    mjd_start: ModifiedJulianDate,
    /// Number of segments.
    num_segments: usize,
    /// Chebyshev coefficients for X coordinate: [segment][CHEB_NODES].
    cx: Vec<[Kilometers; CHEB_NODES]>,
    /// Chebyshev coefficients for Y coordinate.
    cy: Vec<[Kilometers; CHEB_NODES]>,
    /// Chebyshev coefficients for Z coordinate.
    cz: Vec<[Kilometers; CHEB_NODES]>,
}

impl MoonPositionCache {
    /// Build the cache covering `[mjd_start, mjd_end]` (Modified Julian Dates).
    ///
    /// Adds a small padding on each side to accommodate Brent probes
    /// near the boundaries.
    pub fn new(mjd_start: ModifiedJulianDate, mjd_end: ModifiedJulianDate) -> Self {
        //let jd_start = JulianDate::from(mjd_start).value();
        //let jd_end = JulianDate::from(mjd_end).value();
        let pad = Days::new(1.0); // 1-day padding on each side
        let t0 = mjd_start - pad;
        let span = mjd_end + pad - t0;
        let num_segments = ((span / SEGMENT_DAYS).value().ceil() as usize).max(1);

        let nodes = cheb_nodes_on_unit();
        let mut cx = Vec::with_capacity(num_segments);
        let mut cy = Vec::with_capacity(num_segments);
        let mut cz = Vec::with_capacity(num_segments);

        for seg in 0..num_segments {
            let seg_start = t0 + seg as f64 * SEGMENT_DAYS;
            let seg_mid = seg_start + SEGMENT_DAYS * 0.5;
            let seg_half = SEGMENT_DAYS * 0.5;

            let mut vx = [Kilometers::zero(); CHEB_NODES];
            let mut vy = [Kilometers::zero(); CHEB_NODES];
            let mut vz = [Kilometers::zero(); CHEB_NODES];

            for k in 0..CHEB_NODES {
                let mjd_k = seg_mid + seg_half * nodes[k];
                let pos = Moon::get_geo_position::<Kilometer>(mjd_k.into());
                vx[k] = pos.x();
                vy[k] = pos.y();
                vz[k] = pos.z();
            }

            cx.push(compute_cheb_coeffs(&vx));
            cy.push(compute_cheb_coeffs(&vy));
            cz.push(compute_cheb_coeffs(&vz));
        }

        Self {
            mjd_start: t0,
            num_segments,
            cx,
            cy,
            cz,
        }
    }

    /// Evaluate the cached geocentric ecliptic (X, Y, Z) in km at `jd`.
    ///
    /// Falls back to full ELP2000 if `jd` is outside the cached range.
    #[inline]
    pub fn get_position_km(&self, mjd: ModifiedJulianDate) -> (Kilometers, Kilometers, Kilometers) {
        let offset = mjd - self.mjd_start;
        let seg_idx = (offset / SEGMENT_DAYS).value() as usize;

        if seg_idx >= self.num_segments {
            // Fallback: outside cache range
            let pos = Moon::get_geo_position::<Kilometer>(mjd.into());
            return (pos.x(), pos.y(), pos.z());
        }

        // Map jd into [-1, 1] within the segment
        let seg_start = self.mjd_start + seg_idx as f64 * SEGMENT_DAYS;
        let seg_mid = seg_start + SEGMENT_DAYS * 0.5;
        let x = (mjd - seg_mid) / (SEGMENT_DAYS * 0.5);
        let x = x.value(); // dimensionless
        let px = clenshaw_eval(&self.cx[seg_idx], x);
        let py = clenshaw_eval(&self.cy[seg_idx], x);
        let pz = clenshaw_eval(&self.cz[seg_idx], x);
        (px, py, pz)
    }
}

// =============================================================================
// NutationCache
// =============================================================================

/// Linear interpolation cache for nutation parameters (Δψ, Δε, ε₀),
/// stored as radians at regular 2-hour intervals.
pub struct NutationCache {
    /// Modified Julian Date of the first entry.
    mjd_start: ModifiedJulianDate,
    /// Number of entries.
    num_entries: usize,
    /// Pre-evaluated [dpsi_rad, deps_rad, eps0_rad] at each node.
    values: Vec<[Radians; 3]>,
}

impl NutationCache {
    /// Build the nutation cache covering `[jd_start, jd_end]`.
    pub fn new(mjd_start: ModifiedJulianDate, mjd_end: ModifiedJulianDate) -> Self {
        let pad = Days::new(1.0); // 1-day padding
        let t0 = mjd_start - pad;
        let t1 = mjd_end + pad;
        let span = t1 - t0;
        let num_entries = ((span / NUT_STEP_DAYS).value().ceil() as usize) + 1;

        let mut values = Vec::with_capacity(num_entries);
        for i in 0..num_entries {
            let mjd = t0 + i as f64 * NUT_STEP_DAYS;
            let nut = get_nutation(mjd.into());
            values.push([
                nut.longitude.to::<Radian>(),
                nut.obliquity.to::<Radian>(),
                nut.ecliptic.to::<Radian>(),
            ]);
        }

        Self {
            mjd_start: t0,
            num_entries,
            values,
        }
    }

    /// Interpolate (Δψ, Δε, ε₀) in radians at `jd`.
    #[inline]
    pub fn get_nutation_rad(&self, mjd: ModifiedJulianDate) -> (Radians, Radians, Radians) {
        let offset = mjd - self.mjd_start;
        let frac = (offset / NUT_STEP_DAYS).value();
        let idx = frac as usize;

        if idx + 1 >= self.num_entries {
            // Fallback: outside cache range — compute directly
            let nut = get_nutation(mjd.into());
            return (
                nut.longitude.to::<Radian>(),
                nut.obliquity.to::<Radian>(),
                nut.ecliptic.to::<Radian>(),
            );
        }

        let t = frac - idx as f64; // interpolation parameter [0, 1)
        let a = &self.values[idx];
        let b = &self.values[idx + 1];

        (
            a[0] + t * (b[0] - a[0]),
            a[1] + t * (b[1] - a[1]),
            a[2] + t * (b[2] - a[2]),
        )
    }

    /// Build the nutation rotation matrix from cached values at `jd`.
    ///
    /// Equivalent to [`crate::astro::nutation::nutation_rotation`] but uses
    /// interpolated nutation parameters instead of the full 63-term series.
    #[inline]
    pub fn nutation_rotation(&self, mjd: ModifiedJulianDate) -> affn::Rotation3 {
        let (dpsi, deps, eps0) = self.get_nutation_rad(mjd);
        let dpsi = dpsi.value();
        let deps = deps.value();
        let eps0 = eps0.value();

        // R1(ε0+Δε) · R3(Δψ) · R1(−ε0)
        affn::Rotation3::from_x_rotation(eps0 + deps)
            * affn::Rotation3::from_z_rotation(dpsi)
            * affn::Rotation3::from_x_rotation(-eps0)
    }
}

// =============================================================================
// MoonAltitudeContext — combines caches + site data
// =============================================================================

/// Pre-built context for fast repeated Moon altitude queries at a fixed site.
///
/// Holds:
/// * [`MoonPositionCache`] — Chebyshev interpolation of ELP2000 geocentric XYZ
/// * [`NutationCache`] — linear interpolation of IAU 1980 nutation values
/// * Precomputed observer ITRF position in km
///
/// The [`altitude_rad`] method reproduces the full transform chain
/// (ecliptic → equatorial → topocentric → precession → nutation → horizontal)
/// but replaces the two most expensive steps (ELP2000 and nutation) with
/// cache lookups.
pub struct MoonAltitudeContext {
    pos_cache: MoonPositionCache,
    nut_cache: NutationCache,
    /// Observer ITRF position in km: [x, y, z].
    site_itrf_km: [Kilometers; 3],
    /// Observer geodetic latitude in radians.
    lat_rad: f64,
    /// Observer geodetic longitude in degrees (for LST computation).
    lon: qtty::Degrees,
}

impl MoonAltitudeContext {
    /// Build an altitude context covering the MJD period for a given site.
    ///
    /// The caches are padded by 1 day on each side to accommodate Brent probes.
    pub fn new(
        mjd_start: ModifiedJulianDate,
        mjd_end: ModifiedJulianDate,
        site: ObserverSite,
    ) -> Self {
        // Convert ModifiedJulianDate to JulianDate for internal cache usage
        let pos_cache = MoonPositionCache::new(mjd_start, mjd_end);
        let nut_cache = NutationCache::new(mjd_start, mjd_end);

        // Precompute site ITRF position in km
        let site_ecef = site.geocentric_itrf::<Kilometer>();
        let site_itrf_km = [site_ecef.x(), site_ecef.y(), site_ecef.z()];

        Self {
            pos_cache,
            nut_cache,
            site_itrf_km,
            lat_rad: site.lat.to::<Radian>().value(),
            lon: site.lon,
        }
    }

    /// Compute the Moon's topocentric altitude in radians at `mjd`.
    ///
    /// Replicates the full transform chain of
    /// [`Moon::get_horizontal`] → [`Moon::get_apparent_topocentric_equ`]
    /// but replaces ELP2000 and nutation with cached interpolations.
    #[inline]
    pub fn altitude_rad(&self, mjd: ModifiedJulianDate) -> Quantity<Radian> {
        // ---------------------------------------------------------------
        // 1. Geocentric ecliptic Cartesian (km) — from Chebyshev cache
        // ---------------------------------------------------------------
        let (x_ecl, y_ecl, z_ecl) = self.pos_cache.get_position_km(mjd);

        // ---------------------------------------------------------------
        // 2. Ecliptic → EquatorialMeanJ2000 (constant rotation about +X)
        // ---------------------------------------------------------------
        let (sin_e, cos_e) = J2000_OBLIQUITY_RAD.sin_cos();
        let x_eq = x_ecl;
        let y_eq = cos_e * y_ecl - sin_e * z_ecl;
        let z_eq = sin_e * y_ecl + cos_e * z_ecl;

        // ---------------------------------------------------------------
        // 3. Topocentric correction: subtract observer position in J2000 eq
        // ---------------------------------------------------------------
        let gmst_rad = unmodded_gst(mjd.into()).to::<Radian>().value();
        let (sin_g, cos_g) = gmst_rad.sin_cos();

        let sx = self.site_itrf_km[0];
        let sy = self.site_itrf_km[1];
        let sz = self.site_itrf_km[2];

        // ITRF → EquatorialMeanJ2000 via R_z(-GMST)
        let site_eq_x = sx * cos_g - sy * sin_g;
        let site_eq_y = sx * sin_g + sy * cos_g;
        let site_eq_z = sz;

        let x_topo = x_eq - site_eq_x;
        let y_topo = y_eq - site_eq_y;
        let z_topo = z_eq - site_eq_z;

        // ---------------------------------------------------------------
        // 4. Precession: J2000 → mean-of-date
        // ---------------------------------------------------------------
        let rot_prec = precession_rotation_from_j2000(mjd.into());
        let [x_mod, y_mod, z_mod] =
            rot_prec.apply_array([x_topo.value(), y_topo.value(), z_topo.value()]);

        // ---------------------------------------------------------------
        // 5. Nutation: mean-of-date → true-of-date (from cache)
        // ---------------------------------------------------------------
        let rot_nut = self.nut_cache.nutation_rotation(mjd);
        let [x_tod, y_tod, z_tod] = rot_nut.apply_array([x_mod, y_mod, z_mod]);

        // ---------------------------------------------------------------
        // 6. Equatorial true-of-date → RA, Dec
        // ---------------------------------------------------------------
        let ra_rad = y_tod.atan2(x_tod);
        let r_xy = (x_tod * x_tod + y_tod * y_tod).sqrt();
        let dec_rad = z_tod.atan2(r_xy);

        // ---------------------------------------------------------------
        // 7. GAST → LST → HA → altitude
        // ---------------------------------------------------------------
        let gst = calculate_gst(mjd.into());
        let lst = calculate_lst(gst, self.lon);
        let ra_deg = qtty::Degrees::new(ra_rad.to_degrees());
        let ha_rad = (lst - ra_deg).normalize().to::<Radian>().value();

        let sin_alt =
            dec_rad.sin() * self.lat_rad.sin() + dec_rad.cos() * self.lat_rad.cos() * ha_rad.cos();

        Quantity::<Radian>::new(sin_alt.asin())
    }
}

// =============================================================================
// find_and_label_crossings — avoids probe evaluations
// =============================================================================

use crate::calculus::math_core::intervals::LabeledCrossing;
use crate::calculus::math_core::root_finding;
use crate::time::{ModifiedJulianDate, Period};

type Mjd = ModifiedJulianDate;
type Days = qtty::Quantity<qtty::Day>;

/// Tiny epsilon for deduplication (same as intervals.rs).
const DEDUPE_EPS: Days = Days::new(1e-8);

/// Combined scan + Brent root-finding + labelling in a single pass.
///
/// Unlike the generic `find_crossings` + `label_crossings` pipeline in
/// `intervals.rs`, this function records the crossing direction directly
/// from the sign change that triggered the Brent solve — **eliminating
/// the 2 extra probe evaluations per crossing**.
pub fn find_and_label_crossings<V, F>(
    period: Period<Mjd>,
    step: Days,
    f: &F,
    threshold: qtty::Quantity<V>,
) -> (Vec<LabeledCrossing>, bool)
where
    V: qtty::Unit,
    F: Fn(ModifiedJulianDate) -> qtty::Quantity<V>,
{
    let g = |t: Mjd| -> qtty::Quantity<V> { f(t) - threshold };

    let t_start = period.start;
    let t_end = period.end;

    let start_val = g(t_start);
    let start_above = start_val > 0.0;

    let mut labeled = Vec::new();
    let mut t = t_start;
    let mut prev = start_val;

    while t < t_end {
        let next_t = (t + step).min(t_end);
        let next_v = g(next_t);

        if prev.value() * next_v.value() < 0.0 {
            if let Some(root) =
                root_finding::brent_with_values(Period::new(t, next_t), prev, next_v, g)
            {
                if root >= t_start && root <= t_end {
                    // Direction from sign change: prev < 0 → next > 0 means entering (+1)
                    let direction = if prev < 0.0 { 1 } else { -1 };
                    labeled.push(LabeledCrossing { t: root, direction });
                }
            }
        }

        t = next_t;
        prev = next_v;
    }

    // Sort and deduplicate (should already be sorted from linear scan)
    labeled.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    labeled.dedup_by(|a, b| {
        let dup = (a.t - b.t).abs() < DEDUPE_EPS;
        if dup {
            // Keep the earlier one (b), discard a
        }
        dup
    });

    (labeled, start_above)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::calculus::lunar::moon_altitude_rad;
    use crate::coordinates::centers::ObserverSite;
    use crate::observatories::ROQUE_DE_LOS_MUCHACHOS;
    use crate::time::JulianDate;
    use qtty::Radians;

    #[test]
    fn chebyshev_position_accuracy() {
        // Compare cached vs. direct ELP2000 at random times within a segment
        let mjd_start: ModifiedJulianDate = JulianDate::J2000.into();
        let mjd_end = mjd_start + Days::new(30.0);
        let cache = MoonPositionCache::new(mjd_start, mjd_end);

        for i in 0..100 {
            let mjd = mjd_start + Days::new((i as f64) * 0.3 + 0.1); // sample every ~7 hours
            let (cx, cy, cz) = cache.get_position_km(mjd);
            let direct = Moon::get_geo_position::<Kilometer>(JulianDate::from(mjd));
            let (dx, dy, dz) = (direct.x(), direct.y(), direct.z());
            let err =
                ((cx - dx).value().powi(2) + (cy - dy).value().powi(2) + (cz - dz).value().powi(2))
                    .sqrt();
            // Error should be < 1 km (≈ 0.5 arcsecond at Moon distance)
            assert!(
                err < 1.0,
                "Chebyshev error at MJD {mjd}: {err:.6} km (x:{} vs {}, y:{} vs {}, z:{} vs {})",
                cx.value(),
                dx.value(),
                cy.value(),
                dy.value(),
                cz.value(),
                dz.value()
            );
        }
    }

    #[test]
    fn nutation_cache_accuracy() {
        let mjd_start: ModifiedJulianDate = JulianDate::J2000.into();
        let mjd_end: ModifiedJulianDate = mjd_start + Days::new(30.0);
        let cache = NutationCache::new(mjd_start, mjd_end);

        for i in 0..100 {
            let mjd = mjd_start + Days::new((i as f64) * 0.3 + 0.05);
            let (dpsi, deps, eps0) = cache.get_nutation_rad(mjd);
            let direct = get_nutation(mjd.into());
            let d_dpsi = direct.longitude.to::<Radian>().value();
            let d_deps = direct.obliquity.to::<Radian>().value();
            let d_eps0 = direct.ecliptic.to::<Radian>().value();

            let err_dpsi = (dpsi.value() - d_dpsi).abs();
            let err_deps = (deps.value() - d_deps).abs();
            let err_eps0 = (eps0.value() - d_eps0).abs();

            // Errors should be < 5e-10 radians (≈ 0.1 milliarcseconds)
            assert!(
                err_dpsi < 5e-10,
                "Nutation dpsi error at MJD {mjd}: {err_dpsi:.2e}"
            );
            assert!(
                err_deps < 5e-10,
                "Nutation deps error at MJD {mjd}: {err_deps:.2e}"
            );
            assert!(
                err_eps0 < 5e-10,
                "Nutation eps0 error at MJD {mjd}: {err_eps0:.2e}"
            );
        }
    }

    #[test]
    fn cached_altitude_matches_direct() {
        let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);
        let mjd_start: ModifiedJulianDate = JulianDate::J2000.into();
        let mjd_end = mjd_start + Days::new(7.0);
        let ctx = MoonAltitudeContext::new(mjd_start, mjd_end, site);

        for i in 0..50 {
            let mjd = ModifiedJulianDate::new(mjd_start.value() + (i as f64) * 0.14 + 0.01);
            let cached_alt = ctx.altitude_rad(mjd).value();
            let direct_alt = moon_altitude_rad(mjd, &site).value();

            let err_deg = (cached_alt - direct_alt).to_degrees().abs();
            // Should match within ~0.01° (limited by interpolation + nutation cache)
            assert!(
                err_deg < 0.02,
                "Altitude error at MJD {}: cached={:.4}° direct={:.4}° err={:.6}°",
                mjd.value(),
                cached_alt.to_degrees(),
                direct_alt.to_degrees(),
                err_deg
            );
        }
    }

    #[test]
    fn find_and_label_crossings_sine_wave() {
        // Test with a known sine wave: sin(2π(t+0.05)) crosses 0 at known times
        let f = |t: Mjd| Radians::new((2.0 * std::f64::consts::PI * (t.value() + 0.05)).sin());
        let period = Period::new(Mjd::new(0.0), Mjd::new(1.0));
        let step = Days::new(0.01);

        let (labeled, _start_above) = find_and_label_crossings(period, step, &f, Radians::new(0.0));

        // Should find 2 crossings (at t ≈ -0.05 + 0.5 = 0.45 and t ≈ -0.05 + 1.0 = 0.95)
        assert_eq!(labeled.len(), 2, "Expected 2 crossings, got {:?}", labeled);

        // First crossing should be exiting (sin going negative), second entering
        assert_eq!(labeled[0].direction, -1);
        assert_eq!(labeled[1].direction, 1);
    }
}
