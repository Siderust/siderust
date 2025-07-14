//! Time unit: Julian Day
//!
//! This module provides the `JulianDay` struct, representing the Julian Day number,
//! a continuous count of days since the beginning of the Julian Period used in astronomy.

use super::Days;

use chrono::{DateTime, Utc};
use std::ops::{Add, Sub, AddAssign, SubAssign};

#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
/// Represents a Julian Day number (continuous count of days since the Julian Period).
pub struct JulianDay(Days);

impl JulianDay {
    const _J2000_: Days = Days::new(2_451_545.0);

    pub const J2000: JulianDay = JulianDay(Self::_J2000_); // Reference JD for J2000.0 epoch
    pub const JULIAN_YEAR: Days = Days::new(365.25);

    pub const fn new(jd: f64) -> Self {
        JulianDay(Days::new(jd))
    }

    #[inline]
    pub fn julian_millennias(&self) -> f64 {
        (self.0 - Self::_J2000_).value() / 365_250.0
    }

    #[inline]
    pub fn julian_centuries(&self) -> super::Centuries {
        super::Centuries::new((self.0 - Self::_J2000_).value() / 36_525.0)
    }

    #[inline]
    pub fn julian_years(&self) -> super::JulianYear {
        super::JulianYear::new((self.value() - Self::J2000.value()) / 365.25)
    }

    /// Converts JD(TT) to JD(TDB)
    pub fn tt_to_tdb(jd_tt: JulianDay) -> JulianDay {
        //let t = julian_centuries(jd_tt);
        let e = (357.53 + 0.98560028 * (jd_tt - Self::J2000).value()).to_radians();

        let delta_t = Days::new((1.658e-3 * e.sin() + 1.4e-6 * (2.0 * e).sin()) / 86400.0); // Convert to days
        jd_tt + delta_t
    }

    #[inline]
    pub fn value(&self) -> f64 {
        self.0.value()
    }

    #[inline]
    pub fn days(&self) -> Days {
        self.0
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
        JulianDay::new(jd)
    }
}

impl std::fmt::Display for JulianDay {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Julian Day: {}", self.0)
    }
}

impl Add<Days> for JulianDay {
    type Output = JulianDay;

    fn add(self, days: Days) -> JulianDay {
        JulianDay(self.0 + days)
    }
}

impl AddAssign<Days> for JulianDay {
    fn add_assign(&mut self, rhs: Days) {
        self.0 += rhs;
    }
}

impl Add<super::Years> for JulianDay {
    type Output = JulianDay;

    fn add(self, years: super::Years) -> JulianDay {
        self + Self::JULIAN_YEAR* years.value()
    }
}

impl Sub for JulianDay {
    type Output = Days;

    fn sub(self, other: JulianDay) -> Days {
        self.0 - other.0
    }
}

impl Sub<Days> for JulianDay {
    type Output = JulianDay;

    fn sub(self, other: Days) -> JulianDay {
        JulianDay(self.0 - other)
    }
}

impl SubAssign<Days> for JulianDay {
    fn sub_assign(&mut self, rhs: Days) {
        self.0 -= rhs;
    }
}

impl std::ops::Div<Days> for JulianDay {
    type Output = f64;
    fn div(self, days: Days) -> f64 {
        self.value() / days.value()
    }
}

impl std::ops::Div<f64> for JulianDay {
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
        let jd = JulianDay::new(2451545.0);
        assert_eq!(jd.value(), 2451545.0);
    }

    #[test]
    fn test_to_naive_datetime() {
        let jd = JulianDay::new(2451545.0);
        let datetime = jd.to_utc();
        assert_eq!(datetime, DateTime::from_timestamp(946728000, 0));
    }

    #[test]
    fn test_from_naive_datetime() {
        let datetime = DateTime::from_timestamp(946728000, 0).unwrap();
        let jd = JulianDay::from_utc(datetime);
        assert_eq!(jd.value(), 2451545.0);
    }
}
