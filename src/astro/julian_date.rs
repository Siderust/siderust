//! Astro Date: Julian Date
//!
//! This module provides the `JulianDate` struct, representing the Julian Date,
//! a continuous count of days since the beginning of the Julian Period used in astronomy.

use crate::units::*;

use chrono::{DateTime, Utc};
use std::ops::{Add, Sub, AddAssign, SubAssign};

#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
/// Represents a Julian Day number (continuous count of days since the Julian Period).
pub struct JulianDate(Days);

impl JulianDate {
    const _J2000_: Days = Days::new(2_451_545.0);

    pub const J2000: JulianDate = JulianDate(Self::_J2000_); // Reference JD for J2000.0 epoch
    pub const JULIAN_YEAR: Days = Days::new(365.25);

    pub const fn new(jd: f64) -> Self {
        JulianDate(Days::new(jd))
    }

    #[inline]
    pub fn julian_millennias(&self) -> f64 {
        (self.0 - Self::_J2000_).value() / 365_250.0
    }

    #[inline]
    pub fn julian_centuries(&self) -> Centuries {
        Centuries::new((self.0 - Self::_J2000_).value() / 36_525.0)
    }

    #[inline]
    pub fn julian_years(&self) -> JulianYears {
        JulianYears::new((self.value() - Self::J2000.value()) / 365.25)
    }

    /// Converts JD(TT) to JD(TDB)
    pub fn tt_to_tdb(jd_tt: JulianDate) -> JulianDate {
        //let t = julian_centuries(jd_tt);
        let e = (357.53 + 0.98560028 * (jd_tt - Self::J2000).value()).to_radians();

        let delta_t = Days::new((1.658e-3 * e.sin() + 1.4e-6 * (2.0 * e).sin()) / 86400.0); // Convert to days
        jd_tt + delta_t
    }

    #[inline]
    pub fn value(&self) -> f64 {
        self.0.value()
    }

    pub fn to_utc(&self) -> Option<DateTime<Utc>> {
        let unix_epoch_jd = 2440587.5; // Julian Day for Unix epoch (1970-01-01T00:00:00Z)
        let seconds_since_epoch = (self.value() - unix_epoch_jd) * 86400.0;
        let secs = seconds_since_epoch.floor() as i64;
        let nanos = ((seconds_since_epoch - secs as f64) * 1e9) as u32;
        DateTime::<Utc>::from_timestamp(secs, nanos)
    }

    pub fn from_utc(datetime: DateTime<Utc>) -> Self {
        let unix_epoch_jd = 2440587.5; // Julian Day for Unix epoch (1970-01-01T00:00:00Z)
        let seconds_since_epoch = datetime.timestamp() as f64;
        let nanos = datetime.timestamp_subsec_nanos() as f64 / 1e9;
        let jd = unix_epoch_jd + (seconds_since_epoch + nanos) / 86400.0;
        JulianDate::new(jd)
    }

    pub const fn min(&self, other: JulianDate) -> JulianDate {
        JulianDate(self.0.min(other.0))
    }
}

impl std::fmt::Display for JulianDate {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Julian Day: {}", self.0)
    }
}

impl Add<Days> for JulianDate {
    type Output = JulianDate;

    fn add(self, days: Days) -> JulianDate {
        JulianDate(self.0 + days)
    }
}

impl AddAssign<Days> for JulianDate {
    fn add_assign(&mut self, rhs: Days) {
        self.0 += rhs;
    }
}

impl Add<Years> for JulianDate {
    type Output = JulianDate;

    fn add(self, years: Years) -> JulianDate {
        self + Self::JULIAN_YEAR* years.value()
    }
}

impl Sub for JulianDate {
    type Output = Days;

    fn sub(self, other: JulianDate) -> Days {
        self.0 - other.0
    }
}

impl Sub<Days> for JulianDate {
    type Output = JulianDate;

    fn sub(self, other: Days) -> JulianDate {
        JulianDate(self.0 - other)
    }
}

impl SubAssign<Days> for JulianDate {
    fn sub_assign(&mut self, rhs: Days) {
        self.0 -= rhs;
    }
}

impl std::ops::Div<Days> for JulianDate {
    type Output = f64;
    fn div(self, days: Days) -> f64 {
        self.value() / days.value()
    }
}

impl std::ops::Div<f64> for JulianDate {
    type Output = f64;
    fn div(self, val: f64) -> f64 {
        self.value() / val
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::DateTime;

    #[test]
    fn test_julian_day_creation() {
        let jd = JulianDate::new(2451545.0);
        assert_eq!(jd.value(), 2451545.0);
    }

    #[test]
    fn test_to_naive_datetime() {
        let jd = JulianDate::new(2451545.0);
        let datetime = jd.to_utc();
        assert_eq!(datetime, DateTime::from_timestamp(946728000, 0));
    }

    #[test]
    fn test_from_naive_datetime() {
        let datetime = DateTime::from_timestamp(946728000, 0).unwrap();
        let jd = JulianDate::from_utc(datetime);
        assert_eq!(jd.value(), 2451545.0);
    }
}
