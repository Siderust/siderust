// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # ΔT (Delta T) Calculation Module
//!
//! This module estimates the difference between **Universal Time (UT)** — the civil time based on the Earth's rotation —
//! and **Terrestrial Time (TT)** — a uniform time scale used in astronomical ephemerides. This difference is known as:
//! **ΔT = TT − UT**.
//!
//! ## Why ΔT Exists
//! The length of a solar day is not constant. The Earth's rotation slows and fluctuates due to gravitational forces
//! from the Moon and Sun (tidal friction), atmospheric changes, and internal geophysical processes.
//! As a result, uniform atomic time (TT) diverges from Earth's rotational time (UT).
//! Accurate astronomical modeling requires converting between these time standards by **adding ΔT**.
//!
//! ## Scientific References
//! * Stephenson & Houlden (1986): *Atlas of Historical Eclipse Maps* (polynomials for years < 948 and 948–1600).
//! * Morrison & Stephenson (2004): "Historical values of the Earth's clock error".
//! * IERS Conventions (2020): official ΔT data tables.
//!
//! This module re-implements the algorithm from Chapter 9 of **Jean Meeus – _Astronomical Algorithms_ (2nd ed. 1998)**.
//!
//! ## API
//! * [`julian_ephemeris_day`] — Converts a [`JulianDate`] into **Julian Ephemeris Day (JDE)** by adding ΔT / 86400.
//! * [`delta_t_seconds`] — Returns ΔT in seconds for a given Julian Day (JD).
//!
//! ## Quick Example
//! ```rust
//! use siderust::astro::JulianDate;
//! use siderust::astro::dynamical_time::julian_ephemeris_day;
//!
//! let jd = JulianDate::new(2_451_545.0); // January 1, 2000 (Epoch J2000)
//! let jde = julian_ephemeris_day(jd);
//! println!("JDE = {jde}");
//! ```
//!
//! ## Valid Time Range
//! The algorithm is valid from ancient times through approximately 2030, with typical uncertainties:
//! ≤ ±2 s before 1800 CE, and ≤ ±0.5 s since 1900.
//! For future predictions, more accurate IERS-provided curves are recommended.

use crate::astro::JulianDate;
use qtty::Days;

/// Total number of tabulated terms (biennial 1620–1992).
const TERMS: usize = 187;

/// Biennial ΔT table from 1620 to 1992 (in seconds), compiled by J. Meeus.
#[rustfmt::skip]
const DELTA_T: [f64; TERMS] = [
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
];

// ------------------------------------------------------------------------------------
// ΔT Approximation Sections by Time Interval
// ------------------------------------------------------------------------------------

/// **Years < 948 CE**
/// Quadratic formula from Stephenson & Houlden (1986).
#[inline]
fn delta_t_ancient(jd: f64) -> f64 {
    let c = (jd - 2_067_314.5) / 36_525.0;
    1_830.0 - 405.0 * c + 46.5 * c * c
}

/// **Years 948–1600 CE**
/// Second polynomial from Stephenson & Houlden (1986).
#[inline]
fn delta_t_medieval(jd: f64) -> f64 {
    let c = (jd - 2_396_758.5) / 36_525.0;
    22.5 * c * c
}

/// **Years 1600–1992**
/// Bicubic interpolation from the biennial `DELTA_T` table.
#[inline]
fn delta_t_table(jd: f64) -> f64 {
    let mut i = ((jd - 2_312_752.5) / 730.5) as usize;
    if i > TERMS - 3 {
        i = TERMS - 3;
    }
    let a = DELTA_T[i + 1] - DELTA_T[i];
    let b = DELTA_T[i + 2] - DELTA_T[i + 1];
    let c = a - b;
    let n = (jd - (2_312_752.5 + 730.5 * i as f64)) / 730.5;
    DELTA_T[i + 1] + n / 2.0 * (a + b + n * c)
}

/// **Years 1992–2010**
/// Interpolation from Meeus's estimated ΔT for 1990, 2000, and 2010.
#[inline]
fn delta_t_recent(jd: f64) -> f64 {
    const DT: [f64; 3] = [56.86, 63.83, 70.0];
    let a = DT[1] - DT[0];
    let b = DT[2] - DT[1];
    let c = b - a;
    let n = (jd - 2_451_544.5) / 3_652.5;
    DT[1] + n / 2.0 * (a + b + n * c)
}

/// **Years > 2010**
/// Extrapolated via Equation (9.1) from Meeus.
#[inline]
fn delta_t_extrapolated(jd: f64) -> f64 {
    let t = jd - 2_382_148.0;
    -15.0 + (t * t) / 41_048_480.0
}

/// Returns **ΔT** in seconds for a Julian Day (JD).
#[inline]
pub fn delta_t_seconds(jd: f64) -> f64 {
    match jd {
        jd if jd < 2_067_314.5 => delta_t_ancient(jd),
        jd if jd < 2_305_447.5 => delta_t_medieval(jd),
        jd if jd < 2_448_622.5 => delta_t_table(jd),
        jd if jd <= 2_455_197.5 => delta_t_recent(jd),
        _ => delta_t_extrapolated(jd),
    }
}

/// Converts a [`JulianDate`] to **Julian Ephemeris Day (JDE)** by adding ΔT/86400.
#[inline]
pub fn julian_ephemeris_day(jd: JulianDate) -> JulianDate {
    pub const SECONDS_PER_DAY: f64 = 86_400.0;
    jd + Days::new(delta_t_seconds(jd.value())) / SECONDS_PER_DAY
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn delta_t_ancient_sample() {
        let dt = delta_t_seconds(2_000_000.0);
        assert!((dt - 2_734.342_214_024_879_5).abs() < 1e-6);
    }

    #[test]
    fn delta_t_medieval_sample() {
        let dt = delta_t_seconds(2_100_000.0);
        assert!((dt - 1_485.280_240_204_242_3).abs() < 1e-6);
    }

    #[test]
    fn delta_t_table_sample() {
        let dt = delta_t_seconds(2_312_752.5);
        assert!((dt - 115.0).abs() < 1e-6);
    }

    #[test]
    fn delta_t_table_upper_clip() {
        let dt = delta_t_table(2_449_356.0);
        assert!((dt - 59.3).abs() < 1e-6);
    }

    #[test]
    fn delta_t_2000() {
        // IERS reference value: ~63.83 ±0.1 s
        let dt = delta_t_seconds(JulianDate::J2000.value());
        assert!((dt - 63.83).abs() < 0.5);
    }

    #[test]
    fn delta_t_recent_sample() {
        let dt = delta_t_seconds(2_453_371.5);
        assert!((dt - 67.016_266_923_586_13).abs() < 1e-6);
    }

    #[test]
    fn delta_t_extrapolated_sample() {
        let dt = delta_t_seconds(2_457_000.0);
        assert!((dt - 121.492_798_369_147_89).abs() < 1e-6);
    }

    #[test]
    fn jde_offset_matches_delta_t() {
        let jd = JulianDate::J2000;
        let jde = julian_ephemeris_day(jd);
        let offset = (jde - jd).value();
        let expected = delta_t_seconds(jd.value()) / 86_400.0;
        assert!((offset - expected).abs() < 1e-9);
    }
}
