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
use crate::astro::sidereal::{unmodded_gst, SIDEREAL_DAY};
use crate::coordinates::centers::{Geocentric, ObserverSite};
use crate::coordinates::spherical;
use crate::time::JulianDate;
use crate::time::{ModifiedJulianDate, Period};
use qtty::*;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// One sidereal day in solar-day units.
const SIDEREAL_DAY_DAYS: Days = SIDEREAL_DAY.to_const::<Day>();

/// Threshold for treating `B = cos(δ) cos(φ)` as numerically zero.
const DEGENERATE_B_EPS: Quantity<Unitless> = Quantity::new(1e-15);

/// Small temporal epsilon used to include crossings at the right boundary.
const CROSSING_EDGE_EPS: Days = Days::new(1e-12);

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
    Crossings { h0: Radians },
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
    a: Quantity<Unitless>,
    /// `cos(δ) · cos(φ)` — amplitude of the sinusoidal term (≥ 0).
    b: Quantity<Unitless>,
    /// RA corrected for precession + nutation (unwrapped degrees).
    ra_corrected: Degrees,
    /// Observer longitude (east positive).
    lon: Degrees,
}

impl StarAltitudeParams {
    /// Build analytical parameters by precessing J2000 RA/Dec to
    /// `epoch_jd` and applying nutation correction.
    ///
    /// The epoch should be near the centre of the search window so that
    /// precession drift at the edges stays within the Brent bracket.
    /// Over a full year the drift is < 50″ of RA, i.e. < 3 s of time.
    pub fn from_j2000(
        equatorial_j2000: spherical::direction::EquatorialMeanJ2000,
        site: &ObserverSite,
        epoch_jd: JulianDate,
    ) -> Self {
        // 1. Build a J2000 spherical position (unit distance is irrelevant)
        let pos = equatorial_j2000.position::<Geocentric, LightYear>(LightYears::new(1.0));

        // 2. Precess to mean‑of‑date
        let mod_pos = precession::precess_from_j2000(pos, epoch_jd);

        // 3. Apply nutation to RA
        let ra_corrected = corrected_ra_with_nutation(&mod_pos.direction(), epoch_jd);
        let dec = mod_pos.polar;

        // 4. Precompute trig constants
        let lat_rad = site.lat.to::<Radian>();
        let dec_rad = dec.to::<Radian>();

        Self {
            a: Quantity::new(dec_rad.sin() * lat_rad.sin()),
            b: Quantity::new(dec_rad.cos() * lat_rad.cos()),
            ra_corrected,
            lon: site.lon,
        }
    }

    /// Solve for the threshold crossing hour angle.
    ///
    /// * `AlwaysAbove`  — the star's minimum altitude exceeds the threshold
    /// * `NeverAbove`   — the star's maximum altitude is below the threshold
    /// * `Crossings`    — two crossings per sidereal day at HA = ±H₀
    pub fn threshold_ha(&self, threshold: Radians) -> ThresholdResult {
        let sin_h: Quantity<Unitless> = Quantity::new(threshold.sin());

        if self.b.abs() < DEGENERATE_B_EPS {
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
                h0: Radians::new(cos_h0.acos()),
            }
        }
    }

    /// Unwrapped hour angle (degrees) at a given Mjd.
    ///
    /// `HA = GST(t) + λ − α`
    #[inline]
    fn hour_angle(&self, mjd: ModifiedJulianDate) -> Degrees {
        let jd: JulianDate = mjd.into();
        unmodded_gst(jd) + self.lon - self.ra_corrected
    }

    /// Predict all threshold crossing Mjd times within `period` for a
    /// given threshold hour angle H₀.
    ///
    /// Returns `(mjd, direction)` pairs sorted chronologically, where
    /// direction = **+1** (rising above threshold) or **−1** (setting below).
    pub fn predict_crossings(
        &self,
        period: Period<ModifiedJulianDate>,
        h0: Radians,
    ) -> Vec<(ModifiedJulianDate, i32)> {
        let h0_deg = h0.to::<Degree>();

        // HA at the start of the period
        let ha_start = self.hour_angle(period.start);

        let mut crossings = Vec::new();

        // Rising: HA = 360° − H₀  (altitude increases through threshold)
        // Setting: HA = H₀         (altitude decreases through threshold)
        for &(ha_target, dir) in &[(Degrees::FULL_TURN - h0_deg, 1_i32), (h0_deg, -1_i32)] {
            // Time offset from t_start to the first occurrence
            let dha = (ha_target - ha_start).normalize();
            let phase = (dha / Degrees::FULL_TURN).simplify();
            let dt_first: Days = (SIDEREAL_DAY_DAYS * phase).to::<Day>();

            let mut dt = dt_first;
            while period.start + dt <= period.end + CROSSING_EDGE_EPS {
                let t_cross = (period.start + dt).max(period.start).min(period.end);
                crossings.push((t_cross, dir));
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

    fn equatorial_j2000(ra: Degrees, dec: Degrees) -> spherical::direction::EquatorialMeanJ2000 {
        spherical::direction::EquatorialMeanJ2000::new(ra, dec)
    }

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
            equatorial_j2000(
                Degrees::new(37.95), // Polaris RA
                Degrees::new(89.26), // Polaris Dec
            ),
            &greenwich(),
            JulianDate::J2000,
        );
        match params.threshold_ha(Radians::new(0.0)) {
            ThresholdResult::AlwaysAbove => {} // expected
            other => panic!("Polaris should be circumpolar at 51°N, got {:?}", other),
        }
    }

    #[test]
    fn sirius_has_crossings_at_greenwich() {
        let params = StarAltitudeParams::from_j2000(
            equatorial_j2000(
                Degrees::new(101.287), // Sirius RA
                Degrees::new(-16.716), // Sirius Dec
            ),
            &greenwich(),
            JulianDate::J2000,
        );
        match params.threshold_ha(Radians::new(0.0)) {
            ThresholdResult::Crossings { h0 } => {
                assert!(h0 > Radians::zero() && h0 < Radians::HALF_TURN);
            }
            other => panic!(
                "Sirius should have horizon crossings at 51°N, got {:?}",
                other
            ),
        }
    }

    #[test]
    fn never_visible_star_at_greenwich() {
        // Dec = −80° observed from 51°N: max alt = asin(sin(-80°)sin(51°) + cos(-80°)cos(51°))
        // ≈ asin(−0.766 × 0.629 + 0.174 × 0.777) ≈ asin(−0.347) ≈ −20.3°
        // → never above 0°
        let params = StarAltitudeParams::from_j2000(
            equatorial_j2000(Degrees::new(0.0), Degrees::new(-80.0)),
            &greenwich(),
            JulianDate::J2000,
        );
        match params.threshold_ha(Radians::new(0.0)) {
            ThresholdResult::NeverAbove => {} // expected
            other => panic!(
                "Star at Dec=−80° should never be visible at 51°N, got {:?}",
                other
            ),
        }
    }

    #[test]
    fn predict_crossings_count_7_days() {
        let params = StarAltitudeParams::from_j2000(
            equatorial_j2000(Degrees::new(101.287), Degrees::new(-16.716)),
            &greenwich(),
            JulianDate::J2000,
        );
        let period = Period::new(
            ModifiedJulianDate::new(60000.0),
            ModifiedJulianDate::new(60007.0),
        );
        if let ThresholdResult::Crossings { h0 } = params.threshold_ha(Radians::new(0.0)) {
            let crossings = params.predict_crossings(period, h0);
            // 7 days ≈ 7.02 sidereal days → expect ~7 rises + ~7 sets = ~14 crossings
            assert!(
                crossings.len() >= 13 && crossings.len() <= 16,
                "expected ~14 crossings in 7 days, got {}",
                crossings.len()
            );

            // All crossing times should be within the period
            for (t, _) in &crossings {
                assert!(*t >= period.start - Days::new(1e-10));
                assert!(*t <= period.end + Days::new(1e-10));
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
