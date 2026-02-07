// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Star Altitude Analytical Model
//!
//! For a fixed star (static ICRS direction), the altitude as seen from a
//! terrestrial observer varies **sinusoidally** due to Earth's rotation:
//!
//! ```text
//! sin(alt) = A + B · cos(HA)
//! ```
//!
//! where:
//!
//! - `A = sin(δ) · sin(φ)`  — vertical offset (set by declination & latitude)
//! - `B = cos(δ) · cos(φ)`  — amplitude (always ≥ 0)
//! - `HA = LST − α`         — hour angle
//!
//! All stars share the same period (one sidereal day ≈ 23 h 56 m 4 s) since
//! the rotation rate is universal.  Differences between stars appear only in
//! A, B (set by δ and φ) and a phase shift (set by α).
//!
//! This module solves threshold crossings **analytically**:
//!
//! ```text
//! cos(H₀) = (sin(h) − A) / B
//! ```
//!
//! yielding O(1) bracket discovery per sidereal cycle — much faster than the
//! uniform scan used for bodies with non‑trivial orbital motion (Sun, Moon).

use crate::astro::nutation::corrected_ra_with_nutation;
use crate::astro::precession;
use crate::astro::sidereal::unmodded_gst;
use crate::astro::JulianDate;
use crate::coordinates::centers::{Geocentric, ObserverSite};
use crate::coordinates::frames::EquatorialMeanJ2000;
use crate::coordinates::spherical;
use crate::time::{ModifiedJulianDate, Period};
use qtty::*;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Earth's sidereal rotation rate in degrees per solar day.
///
/// This is the coefficient of (JD − J2000) in the IAU 2006 GST polynomial
/// used by [`crate::astro::sidereal::unmodded_gst`].
const OMEGA_SID: f64 = 360.985_647_366_29;

/// One sidereal day in solar days, derived from [`OMEGA_SID`].
///
/// ≈ 0.997 269 57 solar days ≈ 23 h 56 m 4.09 s.
const SIDEREAL_DAY_DAYS: f64 = 360.0 / OMEGA_SID;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Result of analytically solving for threshold crossing hour angles.
#[derive(Debug, Clone, Copy)]
pub(crate) enum ThresholdResult {
    /// Star is always above the threshold (circumpolar at this level).
    AlwaysAbove,
    /// Star never reaches the threshold.
    NeverAbove,
    /// Crossings exist at hour angles ±H₀  (0 < H₀ < π, radians).
    ///
    /// The star **rises** through the threshold at HA = 2π − H₀
    /// and **sets** through it at HA = H₀.
    Crossings { h0_rad: f64 },
}

/// Precomputed parameters for the sinusoidal star‑altitude model.
///
/// Constructed by precessing J2000 coordinates to a reference epoch and
/// caching the trigonometric constants that fully define the sinusoid.
///
/// # Model
///
/// ```text
/// sin(alt) = A + B · cos(HA)
///
/// max altitude = asin(A + B)   (upper culmination, HA = 0)
/// min altitude = asin(A − B)   (lower culmination, HA = π)
/// ```
#[derive(Debug, Clone, Copy)]
pub(crate) struct StarAltitudeParams {
    /// `sin(δ) · sin(φ)` — vertical offset of the sinusoid.
    a: f64,
    /// `cos(δ) · cos(φ)` — amplitude of the sinusoidal term (≥ 0).
    b: f64,
    /// RA corrected for precession + nutation (degrees, unwrapped).
    ra_corrected_deg: f64,
    /// Observer longitude (degrees, east positive).
    lon_deg: f64,
}

impl StarAltitudeParams {
    /// Build analytical parameters by precessing J2000 RA/Dec to
    /// `epoch_jd` and applying nutation correction.
    ///
    /// The epoch should be near the centre of the search window so that
    /// precession drift at the edges stays within the Brent bracket.
    /// Over a full year the drift is < 50″ of RA, i.e. < 3 s of time.
    pub fn from_j2000(
        ra_j2000: Degrees,
        dec_j2000: Degrees,
        site: &ObserverSite,
        epoch_jd: JulianDate,
    ) -> Self {
        // 1. Build a J2000 spherical position (unit distance, irrelevant)
        let pos = spherical::Position::<Geocentric, EquatorialMeanJ2000, LightYear>::new(
            ra_j2000,
            dec_j2000,
            LightYears::new(1.0),
        );

        // 2. Precess to mean‑of‑date
        let mod_pos = precession::precess_from_j2000(pos, epoch_jd);

        // 3. Apply nutation to RA
        let ra_corrected = corrected_ra_with_nutation(&mod_pos.direction(), epoch_jd);
        let dec = mod_pos.polar();

        // 4. Precompute trig constants
        let lat_rad = site.lat.to::<Radian>().value();
        let dec_rad = dec.to::<Radian>().value();

        Self {
            a: dec_rad.sin() * lat_rad.sin(),
            b: dec_rad.cos() * lat_rad.cos(),
            ra_corrected_deg: ra_corrected.value(),
            lon_deg: site.lon.value(),
        }
    }

    /// Solve for the threshold crossing hour angle.
    ///
    /// * `AlwaysAbove`  — the star's minimum altitude exceeds the threshold
    /// * `NeverAbove`   — the star's maximum altitude is below the threshold
    /// * `Crossings`    — two crossings per sidereal day at HA = ±H₀
    pub fn threshold_ha(&self, threshold_rad: f64) -> ThresholdResult {
        let sin_h = threshold_rad.sin();

        if self.b.abs() < 1e-15 {
            // Degenerate: observer at a pole or star at a pole — altitude
            // is constant.
            return if self.a > sin_h {
                ThresholdResult::AlwaysAbove
            } else {
                ThresholdResult::NeverAbove
            };
        }

        let cos_h0 = (sin_h - self.a) / self.b;

        if cos_h0 <= -1.0 {
            ThresholdResult::AlwaysAbove
        } else if cos_h0 >= 1.0 {
            ThresholdResult::NeverAbove
        } else {
            ThresholdResult::Crossings {
                h0_rad: cos_h0.acos(),
            }
        }
    }

    /// Unwrapped hour angle (degrees) at a given MJD.
    ///
    /// `HA = GST(t) + λ − α`
    #[inline]
    fn hour_angle_deg(&self, mjd: ModifiedJulianDate) -> f64 {
        let gst = unmodded_gst(mjd.to_julian_day());
        gst.value() + self.lon_deg - self.ra_corrected_deg
    }

    /// Predict all threshold crossing MJD times within `period` for a
    /// given threshold hour angle H₀.
    ///
    /// Returns `(mjd, direction)` pairs sorted chronologically, where
    /// direction = **+1** (rising above threshold) or **−1** (setting below).
    pub fn predict_crossings(
        &self,
        period: Period<ModifiedJulianDate>,
        h0_rad: f64,
    ) -> Vec<(ModifiedJulianDate, i32)> {
        let h0_deg = h0_rad.to_degrees();
        let t_start = period.start.value();
        let t_end = period.end.value();

        // HA at the start of the period
        let ha_start = self.hour_angle_deg(period.start);

        let mut crossings = Vec::new();

        // Rising: HA = 360° − H₀  (altitude increases through threshold)
        // Setting: HA = H₀         (altitude decreases through threshold)
        for &(ha_target_deg, dir) in &[(360.0 - h0_deg, 1_i32), (h0_deg, -1_i32)] {
            // Time offset from t_start to the first occurrence
            let dha = (ha_target_deg - ha_start).rem_euclid(360.0);
            let dt_first = dha / OMEGA_SID;

            let mut dt = dt_first;
            while t_start + dt <= t_end + 1e-12 {
                let t_cross = (t_start + dt).clamp(t_start, t_end);
                crossings.push((ModifiedJulianDate::new(t_cross), dir));
                dt += SIDEREAL_DAY_DAYS;
            }
        }

        crossings.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        crossings
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn greenwich() -> ObserverSite {
        ObserverSite::new(
            Degrees::new(0.0),
            Degrees::new(51.4769),
            Quantity::<Meter>::new(0.0),
        )
    }

    #[test]
    fn polaris_circumpolar_at_greenwich() {
        let params = StarAltitudeParams::from_j2000(
            Degrees::new(37.95),  // Polaris RA
            Degrees::new(89.26),  // Polaris Dec
            &greenwich(),
            JulianDate::J2000,
        );
        match params.threshold_ha(0.0) {
            ThresholdResult::AlwaysAbove => {} // expected
            other => panic!("Polaris should be circumpolar at 51°N, got {:?}", other),
        }
    }

    #[test]
    fn sirius_has_crossings_at_greenwich() {
        let params = StarAltitudeParams::from_j2000(
            Degrees::new(101.287), // Sirius RA
            Degrees::new(-16.716), // Sirius Dec
            &greenwich(),
            JulianDate::J2000,
        );
        match params.threshold_ha(0.0) {
            ThresholdResult::Crossings { h0_rad } => {
                assert!(h0_rad > 0.0 && h0_rad < std::f64::consts::PI);
            }
            other => panic!("Sirius should have horizon crossings at 51°N, got {:?}", other),
        }
    }

    #[test]
    fn never_visible_star_at_greenwich() {
        // Dec = −80° observed from 51°N: max alt = asin(sin(-80°)sin(51°) + cos(-80°)cos(51°))
        // ≈ asin(−0.766 × 0.629 + 0.174 × 0.777) ≈ asin(−0.347) ≈ −20.3°
        // → never above 0°
        let params = StarAltitudeParams::from_j2000(
            Degrees::new(0.0),
            Degrees::new(-80.0),
            &greenwich(),
            JulianDate::J2000,
        );
        match params.threshold_ha(0.0) {
            ThresholdResult::NeverAbove => {} // expected
            other => panic!("Star at Dec=−80° should never be visible at 51°N, got {:?}", other),
        }
    }

    #[test]
    fn predict_crossings_count_7_days() {
        let params = StarAltitudeParams::from_j2000(
            Degrees::new(101.287),
            Degrees::new(-16.716),
            &greenwich(),
            JulianDate::J2000,
        );
        let period = Period::new(
            ModifiedJulianDate::new(60000.0),
            ModifiedJulianDate::new(60007.0),
        );
        if let ThresholdResult::Crossings { h0_rad } = params.threshold_ha(0.0) {
            let crossings = params.predict_crossings(period, h0_rad);
            // 7 days ≈ 7.02 sidereal days → expect ~7 rises + ~7 sets = ~14 crossings
            assert!(
                crossings.len() >= 13 && crossings.len() <= 16,
                "expected ~14 crossings in 7 days, got {}",
                crossings.len()
            );

            // All crossing times should be within the period
            for (t, _) in &crossings {
                assert!(t.value() >= period.start.value() - 1e-10);
                assert!(t.value() <= period.end.value() + 1e-10);
            }

            // Crossings should alternate rise/set (roughly)
            let rises: usize = crossings.iter().filter(|(_, d)| *d > 0).count();
            let sets: usize = crossings.iter().filter(|(_, d)| *d < 0).count();
            assert!((rises as i32 - sets as i32).unsigned_abs() <= 1);
        } else {
            panic!("Sirius should have crossings at greenwich");
        }
    }
}
