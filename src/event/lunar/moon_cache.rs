// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Moon Ephemeris Cache, Chebyshev Segment Interpolation
//!
//! Interpolates geocentric lunar positions in four-day Chebyshev segments
//! for repeated altitude evaluations. Each evaluation rotates the cached
//! position into the terrestrial frame, subtracts the observer's fixed
//! position, and projects onto the local horizon. The configured default
//! Earth orientation model is evaluated at the observation epoch.

#![allow(unreachable_pub, missing_docs)]

use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::{EclipticMeanJ2000, ECEF};
use crate::coordinates::surface::EnuBasis;
use crate::coordinates::transform::{context::DefaultEphemeris, frame_rotation_with, AstroContext};
use crate::ephemeris::Ephemeris;
use crate::qtty::*;
use crate::time::JulianDate;
use alloc::vec::Vec;
use cheby;

// =============================================================================
// Constants
// =============================================================================

/// Chebyshev polynomial degree for position interpolation.
const CHEB_DEGREE: usize = 8;

/// Number of Chebyshev nodes per segment (degree + 1).
const CHEB_NODES: usize = CHEB_DEGREE + 1;

/// Duration of each Chebyshev segment in days.
const SEGMENT_DAYS: Days = Days::new(4.0);

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
    /// Chebyshev coefficients for X coordinate: [segment][CHEB_NODES], stored as raw km.
    cx: Vec<[f64; CHEB_NODES]>,
    /// Chebyshev coefficients for Y coordinate.
    cy: Vec<[f64; CHEB_NODES]>,
    /// Chebyshev coefficients for Z coordinate.
    cz: Vec<[f64; CHEB_NODES]>,
}

impl MoonPositionCache {
    /// Build the cache covering `[mjd_start, mjd_end]` (Modified Julian Dates).
    ///
    /// Adds a small padding on each side to accommodate Brent probes
    /// near the boundaries.
    ///
    /// # Arguments
    ///
    /// * `mjd_start`, inclusive lower bound of the cached window.
    /// * `mjd_end`, inclusive upper bound of the cached window.
    ///
    /// # Returns
    ///
    /// A populated [`MoonPositionCache`] holding Chebyshev coefficients
    /// for X, Y, Z over each segment.
    pub fn new(mjd_start: ModifiedJulianDate, mjd_end: ModifiedJulianDate) -> Self {
        let pad = Days::new(1.0); // 1-day padding on each side
        let t0 = crate::time::ModifiedJulianDate::new((mjd_start.raw() - pad).value());
        let span = mjd_end.raw() + pad - t0.raw();
        let num_segments = ((span / SEGMENT_DAYS).ceil() as usize).max(1);

        let nodes: [f64; CHEB_NODES] = cheby::nodes();
        let mut cx = Vec::with_capacity(num_segments);
        let mut cy = Vec::with_capacity(num_segments);
        let mut cz = Vec::with_capacity(num_segments);

        for seg in 0..num_segments {
            let seg_start = crate::time::ModifiedJulianDate::new(
                (t0.raw() + seg as f64 * SEGMENT_DAYS).value(),
            );
            let seg_mid = crate::time::ModifiedJulianDate::new(
                (seg_start.raw() + SEGMENT_DAYS * 0.5).value(),
            );
            let seg_half = SEGMENT_DAYS * 0.5;

            let mut vx = [0.0; CHEB_NODES];
            let mut vy = [0.0; CHEB_NODES];
            let mut vz = [0.0; CHEB_NODES];

            for k in 0..CHEB_NODES {
                let mjd_k = crate::time::ModifiedJulianDate::new(
                    (seg_mid.raw() + seg_half * nodes[k]).value(),
                );
                let pos = DefaultEphemeris::moon_geocentric(mjd_k.to::<crate::JD>());
                vx[k] = pos.x().value();
                vy[k] = pos.y().value();
                vz[k] = pos.z().value();
            }

            cx.push(cheby::fit_coeffs(&vx));
            cy.push(cheby::fit_coeffs(&vy));
            cz.push(cheby::fit_coeffs(&vz));
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
    ///
    /// # Arguments
    ///
    /// * `mjd`, evaluation instant.
    ///
    /// # Returns
    ///
    /// `(x, y, z)` in `Kilometers` in the EclipticMeanJ2000 frame.
    #[inline]
    pub fn get_position_km(&self, mjd: ModifiedJulianDate) -> (Kilometers, Kilometers, Kilometers) {
        let offset = mjd.raw() - self.mjd_start.raw();

        if offset < Days::zero() {
            return direct_position(mjd);
        }

        let seg_idx = (offset / SEGMENT_DAYS) as usize;

        if seg_idx >= self.num_segments {
            return direct_position(mjd);
        }

        // Map jd into [-1, 1] within the segment
        let seg_start = crate::time::ModifiedJulianDate::new(
            (self.mjd_start.raw() + seg_idx as f64 * SEGMENT_DAYS).value(),
        );
        let seg_mid =
            crate::time::ModifiedJulianDate::new((seg_start.raw() + SEGMENT_DAYS * 0.5).value());
        let x = (mjd.raw() - seg_mid.raw()) / (SEGMENT_DAYS * 0.5);
        let px = Kilometers::new(cheby::evaluate(&self.cx[seg_idx], x));
        let py = Kilometers::new(cheby::evaluate(&self.cy[seg_idx], x));
        let pz = Kilometers::new(cheby::evaluate(&self.cz[seg_idx], x));
        (px, py, pz)
    }
}

#[inline]
fn direct_position(mjd: ModifiedJulianDate) -> (Kilometers, Kilometers, Kilometers) {
    let pos = DefaultEphemeris::moon_geocentric(mjd.to::<crate::JD>());
    (pos.x(), pos.y(), pos.z())
}

/// Precomputed lunar ephemeris and observer geometry for altitude searches.
pub struct MoonAltitudeContext {
    pos_cache: MoonPositionCache,
    site_ecef_km: [Kilometers; 3],
    basis: EnuBasis,
}

impl MoonAltitudeContext {
    pub fn new(
        mjd_start: ModifiedJulianDate,
        mjd_end: ModifiedJulianDate,
        site: Geodetic<ECEF>,
    ) -> Self {
        Self {
            pos_cache: MoonPositionCache::new(mjd_start, mjd_end),
            site_ecef_km: *site.to_cartesian::<Kilometer>().as_array(),
            basis: EnuBasis::geodetic(site.lon, site.lat),
        }
    }

    /// Geometric topocentric altitude at the TT observation epoch.
    #[inline]
    pub fn altitude_rad(&self, mjd: ModifiedJulianDate) -> Radians {
        let jd: JulianDate = mjd.to::<crate::JD>();
        let ctx: AstroContext = AstroContext::default();
        let (x, y, z) = self.pos_cache.get_position_km(mjd);
        let fixed = frame_rotation_with::<EclipticMeanJ2000, ECEF, _>(jd, &ctx) * [x, y, z];
        let relative = core::array::from_fn(|i| (fixed[i] - self.site_ecef_km[i]).value());
        let [north, east, up] = self.basis.to_horizontal(relative);
        Radians::new(up.atan2((north * north + east * east).sqrt()))
    }
}

// =============================================================================
// find_and_label_crossings, avoids probe evaluations
// =============================================================================

use crate::event::search::intervals::LabeledCrossing;
use crate::event::search::root_finding;
use crate::time::{Interval, ModifiedJulianDate};

#[allow(dead_code)]
type Mjd = ModifiedJulianDate;
type Days = crate::qtty::Quantity<crate::qtty::Day>;

/// Tiny epsilon for deduplication (same as intervals.rs).
#[allow(dead_code)]
const DEDUPE_EPS: Days = Days::new(1e-8);

/// Combined scan + Brent root-finding + labelling in a single pass.
///
/// Unlike the generic `find_crossings` + `label_crossings` pipeline in
/// `intervals.rs`, this function records the crossing direction directly
/// from the sign change that triggered the Brent solve, **eliminating
/// the 2 extra probe evaluations per crossing**.
///
/// # Arguments
///
/// * `period`, MJD search window.
/// * `step`, scan step size used to bracket sign changes.
/// * `f`, scalar function `mjd → Quantity<V>` to root‑find against `threshold`.
/// * `threshold`, value subtracted from `f` before root‑finding.
///
/// # Returns
///
/// A tuple `(crossings, start_above)` where `crossings` is the sorted,
/// deduplicated `Vec<LabeledCrossing>` of root MJDs annotated with the
/// crossing direction (`+1` entering, `−1` leaving), and `start_above`
/// is `true` if `f(start) > threshold`.
#[allow(dead_code)]
pub fn find_and_label_crossings<V, F>(
    period: Interval<ModifiedJulianDate>,
    step: Days,
    f: &F,
    threshold: crate::qtty::Quantity<V>,
) -> (Vec<LabeledCrossing>, bool)
where
    V: crate::qtty::Unit,
    F: Fn(ModifiedJulianDate) -> crate::qtty::Quantity<V>,
{
    let g = |t: Mjd| -> crate::qtty::Quantity<V> { f(t) - threshold };

    let t_start = period.start;
    let t_end = period.end;

    let start_val = g(t_start);
    let start_above = start_val > crate::qtty::Quantity::<V>::new(0.0);

    let mut labeled = Vec::new();
    let mut t = t_start;
    let mut prev = start_val;

    while t < t_end {
        let next_t = {
            let t_nx = crate::time::ModifiedJulianDate::new((t.raw() + step).value());
            if t_nx.raw() <= t_end.raw() {
                t_nx
            } else {
                t_end
            }
        };
        let next_v = g(next_t);

        if prev.signum() * next_v.signum() < 0.0 {
            if let Some(root) =
                root_finding::brent_with_values(Interval::new(t, next_t), prev, next_v, g)
            {
                if root >= t_start && root <= t_end {
                    // Direction from sign change: prev < 0 → next > 0 means entering (+1)
                    let direction = if prev < crate::qtty::Quantity::<V>::new(0.0) {
                        1
                    } else {
                        -1
                    };
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
        let dup = (a.t.raw() - b.t.raw()).abs() < DEDUPE_EPS;
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
    use crate::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
    use crate::event::lunar::moon_altitude_rad;
    use crate::qtty::Radians;

    #[test]
    fn chebyshev_position_accuracy() {
        // Compare cached vs. direct ELP2000 at random times within a segment
        let mjd_start: ModifiedJulianDate = crate::J2000.to::<crate::MJD>();
        let mjd_end =
            crate::time::ModifiedJulianDate::new((mjd_start.raw() + Days::new(30.0)).value());
        let cache = MoonPositionCache::new(mjd_start, mjd_end);

        for i in 0..100 {
            let mjd = crate::time::ModifiedJulianDate::new(
                (mjd_start.raw() + Days::new((i as f64) * 0.3 + 0.1)).value(),
            ); // sample every ~7 hours
            let (cx, cy, cz) = cache.get_position_km(mjd);
            let direct = DefaultEphemeris::moon_geocentric(mjd.to::<crate::JD>());
            let (dx, dy, dz) = (direct.x(), direct.y(), direct.z());
            let err = (cx - dx).abs().max((cy - dy).abs()).max((cz - dz).abs());
            // Error should be < 1 km (≈ 0.5 arcsecond at Moon distance)
            assert!(
                err < Kilometers::new(1.0),
                "Chebyshev max-axis error at MJD {mjd}: {err} (x:{} vs {}, y:{} vs {}, z:{} vs {})",
                cx,
                dx,
                cy,
                dy,
                cz,
                dz
            );
        }
    }

    #[test]
    fn moon_position_cache_falls_back_before_start() {
        let start = ModifiedJulianDate::new(60000.0);
        let end = ModifiedJulianDate::new(60010.0);
        let cache = MoonPositionCache::new(start, end);

        let t = ModifiedJulianDate::new(59998.0);
        let cached = cache.get_position_km(t);
        let direct = DefaultEphemeris::moon_geocentric(t.to::<crate::JD>());

        assert!((cached.0.value() - direct.x().value()).abs() < 1e-12);
        assert!((cached.1.value() - direct.y().value()).abs() < 1e-12);
        assert!((cached.2.value() - direct.z().value()).abs() < 1e-12);
    }

    #[test]
    fn cached_altitude_matches_direct() {
        let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
        let mjd_start: ModifiedJulianDate = crate::J2000.to::<crate::MJD>();
        let mjd_end =
            crate::time::ModifiedJulianDate::new((mjd_start.raw() + Days::new(7.0)).value());
        let ctx = MoonAltitudeContext::new(mjd_start, mjd_end, site);

        for i in 0..50 {
            let mjd = crate::time::ModifiedJulianDate::new(
                (mjd_start.raw() + Days::new((i as f64) * 0.14 + 0.01)).value(),
            );
            let cached_alt = ctx.altitude_rad(mjd);
            let direct_alt = moon_altitude_rad(mjd, &site);

            let err_deg = (cached_alt - direct_alt).abs().to::<Degree>();
            // Should match within ~0.01° (limited by ephemeris interpolation)
            assert!(
                err_deg < Degrees::new(0.02),
                "Altitude error at MJD {}: cached={} direct={} err={}",
                mjd,
                cached_alt.to::<Degree>(),
                direct_alt.to::<Degree>(),
                err_deg
            );
        }
    }

    #[test]
    fn find_and_label_crossings_sine_wave() {
        // Test with a known sine wave: sin(2π(t+0.05)) crosses 0 at known times
        let f =
            |t: Mjd| Radians::new((2.0 * core::f64::consts::PI * (t.raw().value() + 0.05)).sin());
        let period = Interval::new(
            crate::time::ModifiedJulianDate::new((Days::new(0.0)).value()),
            crate::time::ModifiedJulianDate::new((Days::new(1.0)).value()),
        );
        let step = Days::new(0.01);

        let (labeled, _start_above) = find_and_label_crossings(period, step, &f, Radians::new(0.0));

        // Should find 2 crossings (at t ≈ -0.05 + 0.5 = 0.45 and t ≈ -0.05 + 1.0 = 0.95)
        assert_eq!(labeled.len(), 2, "Expected 2 crossings, got {:?}", labeled);

        // First crossing should be exiting (sin going negative), second entering
        assert_eq!(labeled[0].direction, -1);
        assert_eq!(labeled[1].direction, 1);
    }
}
