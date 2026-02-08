// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use chrono::{DateTime, Utc};
use qtty::Days;
use std::ops::{Add, Sub};

#[cfg(feature = "serde")]
use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// Represents Modified Julian Date (MJD), which is the Julian Date
/// minus 2400000.5, used in various scientific and technical applications.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct ModifiedJulianDate(Days);

impl ModifiedJulianDate {
    pub const fn new(mjd: f64) -> Self {
        ModifiedJulianDate(Days::new(mjd))
    }

    /// Returns the inner Modified Julian Day value.
    pub const fn value(&self) -> f64 {
        self.0.value()
    }

    pub const fn to_julian_day(&self) -> super::JulianDate {
        super::JulianDate::new(self.value() + 2400000.5)
    }

    pub fn to_utc(&self) -> Option<DateTime<Utc>> {
        let jd = self.to_julian_day().value();
        let unix_epoch_jd = 2440587.5; // Julian Day for Unix epoch (1970-01-01T00:00:00Z)
        let seconds_since_epoch = (jd - unix_epoch_jd) * 86400.0;
        let secs = seconds_since_epoch.floor() as i64;
        let nanos = ((seconds_since_epoch - secs as f64) * 1e9) as u32;
        DateTime::<Utc>::from_timestamp(secs, nanos)
    }

    pub fn from_utc(datetime: DateTime<Utc>) -> Self {
        let unix_epoch_jd = 2440587.5; // Julian Day for Unix epoch (1970-01-01T00:00:00Z)
        let seconds_since_epoch = datetime.timestamp() as f64;
        let nanos = datetime.timestamp_subsec_nanos() as f64 / 1e9;
        let jd = unix_epoch_jd + (seconds_since_epoch + nanos) / 86400.0;
        ModifiedJulianDate::new(jd - 2400000.5)
    }

    /// Returns the minimum of two MJD values.
    pub const fn min(self, other: ModifiedJulianDate) -> ModifiedJulianDate {
        ModifiedJulianDate(self.0.min_const(other.0))
    }

    /// Returns the maximum of two MJD values.
    pub const fn max(self, other: ModifiedJulianDate) -> ModifiedJulianDate {
        ModifiedJulianDate(self.0.max_const(other.0))
    }
}

#[cfg(feature = "serde")]
impl Serialize for ModifiedJulianDate {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_f64(self.value())
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for ModifiedJulianDate {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let v = f64::deserialize(deserializer)?;
        Ok(ModifiedJulianDate::new(v))
    }
}

impl From<ModifiedJulianDate> for super::JulianDate {
    fn from(mjd: ModifiedJulianDate) -> Self {
        mjd.to_julian_day()
    }
}

impl From<super::JulianDate> for ModifiedJulianDate {
    fn from(jd: super::JulianDate) -> Self {
        ModifiedJulianDate::new(jd.value() - 2400000.5)
    }
}

// Arithmetic operations
impl Add<Days> for ModifiedJulianDate {
    type Output = ModifiedJulianDate;

    fn add(self, rhs: Days) -> Self::Output {
        ModifiedJulianDate(self.0 + rhs)
    }
}

impl Sub<Days> for ModifiedJulianDate {
    type Output = ModifiedJulianDate;

    fn sub(self, rhs: Days) -> Self::Output {
        ModifiedJulianDate(self.0 - rhs)
    }
}

impl Sub<ModifiedJulianDate> for ModifiedJulianDate {
    type Output = Days;

    fn sub(self, rhs: ModifiedJulianDate) -> Self::Output {
        self.0 - rhs.0
    }
}

impl std::fmt::Display for ModifiedJulianDate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "MJD {}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::DateTime;

    #[test]
    fn test_modified_julian_day_creation() {
        let mjd = ModifiedJulianDate::new(51544.5);
        assert_eq!(mjd.value(), 51544.5);
    }

    #[test]
    fn test_to_julian_day() {
        let mjd = ModifiedJulianDate::new(51544.5);
        assert_eq!(mjd.to_julian_day().value(), 2451545.0);
    }

    #[test]
    fn test_to_utc() {
        let mjd = ModifiedJulianDate::new(51544.5);
        let datetime = mjd.to_utc();
        assert_eq!(datetime, DateTime::from_timestamp(946728000, 0));
    }

    #[test]
    fn test_from_utc() {
        let datetime = DateTime::from_timestamp(946728000, 0);
        let mjd = ModifiedJulianDate::from_utc(datetime.unwrap());
        assert_eq!(mjd.value(), 51544.5);
    }

    #[test]
    fn test_add_days() {
        let mjd = ModifiedJulianDate::new(59000.0);
        let result = mjd + Days::new(1.5);
        assert_eq!(result.value(), 59001.5);
    }

    #[test]
    fn test_sub_days() {
        let mjd = ModifiedJulianDate::new(59000.0);
        let result = mjd - Days::new(1.5);
        assert_eq!(result.value(), 58998.5);
    }

    #[test]
    fn test_sub_mjd() {
        let mjd1 = ModifiedJulianDate::new(59001.0);
        let mjd2 = ModifiedJulianDate::new(59000.0);
        let diff = mjd1 - mjd2;
        assert_eq!(diff.value(), 1.0);
    }

    #[test]
    fn test_comparison() {
        let mjd1 = ModifiedJulianDate::new(59000.0);
        let mjd2 = ModifiedJulianDate::new(59001.0);
        assert!(mjd1 < mjd2);
        assert!(mjd2 > mjd1);
    }
}
