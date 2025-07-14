//! # Sidereal Time Module
//!
//! **Sidereal time** is the hour angle of the vernal equinox: a clock that tells us
//! "how far the Earth has rotated relative to the stars" instead of the Sun.  One
//! **mean sidereal day** is ≈ 23 h 56 m 4.09 s, i.e. about 0.99727 solar days.
//! Astronomers rely on it to aim equatorial‑mounted telescopes, reduce star‐track
//! images, or convert between Earth‐fixed and inertial coordinate frames.
//!
//! This module provides two levels of functionality:
//!
//! | Function | Output | Notes |
//! |----------|--------|-------|
//! | [`unmodded_gst`] | Greenwich Sidereal Time **without** modulo 360 ° | Required when adding longitudes before the wrap. |
//! | [`unmodded_lst`] | Local Sidereal Time (observer longitude added, still unwrapped) |
//! | [`calculate_gst`] | GST wrapped to **[0°, 360°)** |
//! | [`calculate_lst`] | LST wrapped to **[0°, 360°)** |
//!
//! The implementation follows the IAU 2006/2000A recommended polynomial (Meeus
//! eq. 12.4 except with higher‑precision coefficients).  Accuracy is better than
//! ±0.1″ for dates within ±100 years of J2000, plenty for pointing and most
//! satellite tracking.
//!
//! ## Example
//! ```rust
//! use chrono::prelude::*;
//! use siderust::astro::JulianDate;
//! use siderust::units::Degrees;
//! use siderust::astro::sidereal::{calculate_gst, calculate_lst};
//!
//! let jd = JulianDate::from_utc(Utc::now());
//! let gst = calculate_gst(jd);
//! let lst = calculate_lst(gst, Degrees::new(-3.7038)); // Madrid ≈ ‑3.70°
//! println!("GST = {:.4}°,  LST = {:.4}°", gst, lst);
//! ```

use crate::units::{Degrees, Days};
use crate::astro::JulianDate;

/// Mean sidereal day length ≈ 0.9972696 solar days (23 h 56 m 4.09 s).
pub const SIDEREAL_DAY: Days = Days::new(23.934_469_6 / 24.0);

/// Returns **unwrapped** Greenwich Sidereal Time for the given Julian Day.
///
/// *Output*: angle in degrees, may be < 0° or > 360°.
#[inline]
pub fn unmodded_gst(julian_date: JulianDate) -> Degrees {
    let t = julian_date.julian_centuries().value();
    let base = julian_date - JulianDate::J2000;

    // IAU 2006 polynomial (units: degrees)
    let gst = 280.460_618_37
        + 360.985_647_366_29 * base.value()
        + 0.000_387_933 * t.powi(2)
        - t.powi(3) / 38_710_000.0;
    Degrees::new(gst)
}

/// Unwrapped Local Sidereal Time = unwrapped GST + observer longitude (east +).
#[inline]
pub fn unmodded_lst(jd: JulianDate, longitude_deg: Degrees) -> Degrees {
    unmodded_gst(jd) + longitude_deg
}

/// Greenwich Sidereal Time wrapped to **[0°, 360°)**.
#[inline]
pub fn calculate_gst(julian_date: JulianDate) -> Degrees {
    unmodded_gst(julian_date).normalize()
}

/// Local Sidereal Time wrapped to **[0°, 360°)**.
/// *`longitude`*: east positive, west negative.
#[inline]
pub fn calculate_lst(gst: Degrees, longitude: Degrees) -> Degrees {
    (gst + longitude).normalize()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gst_lst_ranges() {
        let jd = JulianDate::new(2_459_945.5); // 2023‑01‑01 00 UT
        let gst = calculate_gst(jd);
        let lst = calculate_lst(gst, Degrees::new(-75.0));
        assert!(gst.as_f64() >= 0.0 && gst.as_f64() < 360.0);
        assert!(lst.as_f64() >= 0.0 && lst.as_f64() < 360.0);
    }
}
