// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # ΔT (Delta T) — UT↔TT Correction Layer
//!
//! This module implements the piecewise polynomial model for **ΔT = TT − UT**
//! from Chapter 9 of *Jean Meeus — Astronomical Algorithms (2nd ed. 1998)*.
//!
//! ## Integration with Time Scales
//!
//! The correction is applied **automatically** by the [`UT`](super::UT) time
//! scale marker.  When you convert from `Time<UT>` to any TT-based scale
//! (`.to::<JD>()`, `.to::<MJD>()`, etc.), `UT::to_jd_tt` adds ΔT.
//! The inverse (`UT::from_jd_tt`) uses a three-iteration fixed-point solver.
//!
//! [`Time::from_utc`](super::Time::from_utc) creates a `Time<UT>` internally
//! and then converts to the target scale, so external callers get the ΔT
//! correction without calling any function from this module.
//!
//! ## Quick Example
//! ```rust
//! use siderust::time::{UT, JD, Time};
//!
//! // UT-based Julian Day -> JD(TT) with ΔT applied
//! let ut = Time::<UT>::new(2_451_545.0);
//! let jd_tt = ut.to::<JD>();
//! println!("JD(TT) = {jd_tt}");
//!
//! // Query the raw ΔT value
//! let dt = ut.delta_t();
//! println!("ΔT = {dt}");
//! ```
//!
//! ## Scientific References
//! * Stephenson & Houlden (1986): *Atlas of Historical Eclipse Maps*.
//! * Morrison & Stephenson (2004): "Historical values of the Earth's clock error".
//! * IERS Conventions (2020): official ΔT data tables.
//!
//! ## Valid Time Range
//! The algorithm is valid from ancient times through approximately 2030, with
//! typical uncertainties ≤ ±2 s before 1800 CE and ≤ ±0.5 s since 1900.

use super::instant::Time;
use super::scales::UT;
use super::JulianDate;
use qtty::{Days, Seconds, Simplify};

/// Total number of tabulated terms (biennial 1620–1992).
const TERMS: usize = 187;

/// Biennial ΔT table from 1620 to 1992 (in seconds), compiled by J. Meeus.
#[rustfmt::skip]
const DELTA_T: [Seconds; TERMS] = qtty::qtty_vec!(
    Seconds;
    124.0,115.0,106.0, 98.0, 91.0, 85.0, 79.0, 74.0, 70.0, 65.0,
     62.0, 58.0, 55.0, 53.0, 50.0, 48.0, 46.0, 44.0, 42.0, 40.0,
     37.0, 35.0, 33.0, 31.0, 28.0, 26.0, 24.0, 22.0, 20.0, 18.0,
     16.0, 14.0, 13.0, 12.0, 11.0, 10.0,  9.0,  9.0,  9.0,  9.0,
      9.0,  9.0,  9.0,  9.0, 10.0, 10.0, 10.0, 10.0, 10.0, 11.0,
     11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 12.0, 12.0, 12.0, 12.0,
     12.0, 12.0, 13.0, 13.0, 13.0, 13.0, 14.0, 14.0, 14.0, 15.0,
     15.0, 15.0, 15.0, 16.0, 16.0, 16.0, 16.0, 16.0, 17.0, 17.0,
     17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 16.0, 16.0, 15.0, 14.0,
     13.7, 13.1, 12.7, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 12.3,
     12.0, 11.4, 10.6,  9.6,  8.6,  7.5,  6.6,  6.0,  5.7,  5.6,
      5.7,  5.9,  6.2,  6.5,  6.8,  7.1,  7.3,  7.5,  7.7,  7.8,
      7.9,  7.5,  6.4,  5.4,  2.9,  1.6, -1.0, -2.7, -3.6, -4.7,
     -5.4, -5.2, -5.5, -5.6, -5.8, -5.9, -6.2, -6.4, -6.1, -4.7,
     -2.7,  0.0,  2.6,  5.4,  7.7, 10.5, 13.4, 16.0, 18.2, 20.2,
     21.2, 22.4, 23.5, 23.9, 24.3, 24.0, 23.9, 23.9, 23.7, 24.0,
     24.3, 25.3, 26.2, 27.3, 28.2, 29.1, 30.0, 30.7, 31.4, 32.2,
     33.1, 34.0, 35.0, 36.5, 38.3, 40.2, 42.2, 44.5, 46.5, 48.5,
     50.5, 52.2, 53.8, 54.9, 55.8, 56.9, 58.3,
);

// ------------------------------------------------------------------------------------
// ΔT Approximation Sections by Time Interval
// ------------------------------------------------------------------------------------

/// **Years < 948 CE**
/// Quadratic formula from Stephenson & Houlden (1986).
#[inline]
fn delta_t_ancient(jd: JulianDate) -> Seconds {
    const DT_A0_S: Seconds = Seconds::new(1_830.0);
    const DT_A1_S: Seconds = Seconds::new(-405.0);
    const DT_A2_S: Seconds = Seconds::new(46.5);
    const JD_EPOCH_948_UT: JulianDate = JulianDate::new(2_067_314.5);
    let c = days_ratio(jd - JD_EPOCH_948_UT, JulianDate::JULIAN_CENTURY);
    DT_A0_S + DT_A1_S * c + DT_A2_S * c * c
}

/// **Years 948–1600 CE**
/// Second polynomial from Stephenson & Houlden (1986).
#[inline]
fn delta_t_medieval(jd: JulianDate) -> Seconds {
    const JD_EPOCH_1850_UT: JulianDate = JulianDate::new(2_396_758.5);
    const DT_A2_S: Seconds = Seconds::new(22.5);

    let c = days_ratio(jd - JD_EPOCH_1850_UT, JulianDate::JULIAN_CENTURY);
    DT_A2_S * c * c
}

/// **Years 1600–1992**
/// Bicubic interpolation from the biennial `DELTA_T` table.
#[inline]
fn delta_t_table(jd: JulianDate) -> Seconds {
    const JD_TABLE_START_1620: JulianDate = JulianDate::new(2_312_752.5);
    const BIENNIAL_STEP_D: Days = Days::new(730.5);

    let mut i = days_ratio(jd - JD_TABLE_START_1620, BIENNIAL_STEP_D) as usize;
    if i > TERMS - 3 {
        i = TERMS - 3;
    }
    let a: Seconds = DELTA_T[i + 1] - DELTA_T[i];
    let b: Seconds = DELTA_T[i + 2] - DELTA_T[i + 1];
    let c: Seconds = a - b;
    let n = days_ratio(
        jd - (JD_TABLE_START_1620 + BIENNIAL_STEP_D * i as f64),
        BIENNIAL_STEP_D,
    );
    DELTA_T[i + 1] + n / 2.0 * (a + b + n * c)
}

/// **Years 1992–2010**
/// Interpolation from Meeus's estimated ΔT for 1990, 2000, and 2010.
#[inline]
fn delta_t_recent(jd: JulianDate) -> Seconds {
    const DT: [Seconds; 3] = [Seconds::new(56.86), Seconds::new(63.83), Seconds::new(70.0)];
    const JD_YEAR_2000_UT: JulianDate = JulianDate::new(2_451_544.5);
    const DECADE_D: Days = Days::new(3_652.5);

    let a = DT[1] - DT[0];
    let b = DT[2] - DT[1];
    let c = b - a;
    let n = days_ratio(jd - JD_YEAR_2000_UT, DECADE_D);
    DT[1] + n / 2.0 * (a + b + n * c)
}

/// **Years > 2010**
/// Extrapolated via Equation (9.1) from Meeus.
#[inline]
fn delta_t_extrapolated(jd: JulianDate) -> Seconds {
    const JD_EPOCH_1810_UT: JulianDate = JulianDate::new(2_382_148.0);
    const DT_OFFSET_S: Seconds = Seconds::new(-15.0);
    const QUADRATIC_DIVISOR_D2_PER_S: f64 = 41_048_480.0;

    let t = days_ratio(jd - JD_EPOCH_1810_UT, Days::new(1.0));
    DT_OFFSET_S + Seconds::new((t * t) / QUADRATIC_DIVISOR_D2_PER_S)
}

#[inline]
fn days_ratio(num: Days, den: Days) -> f64 {
    (num / den).simplify().value()
}

/// Returns **ΔT** in seconds for a Julian Day on the **UT** axis.
#[inline]
pub(crate) fn delta_t_seconds_from_ut(jd_ut: JulianDate) -> Seconds {
    match jd_ut {
        jd if jd < JulianDate::new(2_067_314.5) => delta_t_ancient(jd),
        jd if jd < JulianDate::new(2_305_447.5) => delta_t_medieval(jd),
        jd if jd < JulianDate::new(2_448_622.5) => delta_t_table(jd),
        jd if jd <= JulianDate::new(2_455_197.5) => delta_t_recent(jd),
        _ => delta_t_extrapolated(jd_ut),
    }
}

// ── Time<UT> convenience method ───────────────────────────────────────────

impl Time<UT> {
    /// Returns **ΔT = TT − UT** in seconds for this UT epoch.
    ///
    /// This is a convenience accessor; the same correction is applied
    /// automatically when converting to any TT-based scale (`.to::<JD>()`).
    #[inline]
    pub fn delta_t(&self) -> Seconds {
        delta_t_seconds_from_ut(JulianDate::from_days(self.quantity()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::{Day, Days};

    #[test]
    fn delta_t_ancient_sample() {
        let dt = delta_t_seconds_from_ut(JulianDate::new(2_000_000.0));
        assert!((dt - Seconds::new(2_734.342_214_024_879_5)).abs() < Seconds::new(1e-6));
    }

    #[test]
    fn delta_t_medieval_sample() {
        let dt = delta_t_seconds_from_ut(JulianDate::new(2_100_000.0));
        assert!((dt - Seconds::new(1_485.280_240_204_242_3)).abs() < Seconds::new(1e-6));
    }

    #[test]
    fn delta_t_table_sample() {
        let dt = delta_t_seconds_from_ut(JulianDate::new(2_312_752.5));
        assert!((dt - Seconds::new(115.0)).abs() < Seconds::new(1e-6));
    }

    #[test]
    fn delta_t_table_upper_clip() {
        let dt = delta_t_table(JulianDate::new(2_449_356.0));
        assert!((dt - Seconds::new(59.3)).abs() < Seconds::new(1e-6));
    }

    #[test]
    fn delta_t_2000() {
        // IERS reference value: ~63.83 ±0.1 s
        let dt = delta_t_seconds_from_ut(JulianDate::J2000);
        assert!((dt - Seconds::new(63.83)).abs() < Seconds::new(0.5));
    }

    #[test]
    fn delta_t_recent_sample() {
        let dt = delta_t_seconds_from_ut(JulianDate::new(2_453_371.5));
        assert!((dt - Seconds::new(67.016_266_923_586_13)).abs() < Seconds::new(1e-6));
    }

    #[test]
    fn delta_t_extrapolated_sample() {
        let dt = delta_t_seconds_from_ut(JulianDate::new(2_457_000.0));
        assert!((dt - Seconds::new(121.492_798_369_147_89)).abs() < Seconds::new(1e-6));
    }

    #[test]
    fn ut_scale_applies_delta_t() {
        let ut = Time::<UT>::new(2_451_545.0);
        let jd_tt = ut.to::<crate::time::JD>();
        let offset = jd_tt - JulianDate::new(2_451_545.0);
        let expected = delta_t_seconds_from_ut(JulianDate::new(2_451_545.0)).to::<Day>();
        assert!((offset - expected).abs() < Days::new(1e-9));
    }

    #[test]
    fn ut_scale_roundtrip() {
        let jd_tt = JulianDate::new(2_451_545.0);
        let ut: Time<UT> = jd_tt.to::<UT>();
        let back: JulianDate = ut.to::<crate::time::JD>();
        assert!((back - jd_tt).abs() < Days::new(1e-12));
    }

    #[test]
    fn delta_t_convenience_method() {
        let ut = Time::<UT>::new(2_451_545.0);
        let dt = ut.delta_t();
        assert!((dt - Seconds::new(63.83)).abs() < Seconds::new(0.5));
    }
}
