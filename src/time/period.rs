// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Time Period implementation
//!
//! This module provides the `Period<T>` type, representing a time interval
//! between two time instants.

use super::TimeInstant;
use chrono::{DateTime, Utc};
use qtty::Days;
use serde::{ser::SerializeStruct, Deserialize, Deserializer, Serialize, Serializer};

/// Represents a time period between two instants.
///
/// A `Period` is defined by a start and end time instant of type `T`,
/// where `T` implements the `TimeInstant` trait. This allows for periods
/// defined in different time systems (Julian Date, Modified Julian Date, UTC, etc.).
///
/// # Examples
///
/// ```
/// use siderust::time::{Period, ModifiedJulianDate};
///
/// let start = ModifiedJulianDate::new(59000.0);
/// let end = ModifiedJulianDate::new(59001.0);
/// let period = Period::new(start, end);
///
/// // Duration in days
/// let duration = period.duration();
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Period<T: TimeInstant> {
    pub start: T,
    pub end: T,
}

impl<T: TimeInstant> Period<T> {
    /// Creates a new period between two time instants.
    ///
    /// # Arguments
    ///
    /// * `start` - The start time instant
    /// * `end` - The end time instant
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::time::{Period, JulianDate};
    ///
    /// let start = JulianDate::new(2451545.0);
    /// let end = JulianDate::new(2451546.0);
    /// let period = Period::new(start, end);
    /// ```
    pub fn new(start: T, end: T) -> Self {
        Period { start, end }
    }

    /// Returns the duration of the period as the difference between end and start.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::time::{Period, JulianDate};
    /// use qtty::Days;
    ///
    /// let start = JulianDate::new(2451545.0);
    /// let end = JulianDate::new(2451546.5);
    /// let period = Period::new(start, end);
    ///
    /// let duration = period.duration();
    /// assert_eq!(duration, Days::new(1.5));
    /// ```
    pub fn duration(&self) -> T::Duration {
        self.end.difference(&self.start)
    }
}

// Specific implementation for periods with Days duration (JD and MJD)
impl<T: TimeInstant<Duration = Days>> Period<T> {
    /// Returns the duration of the period in days as a floating-point value.
    ///
    /// This method is available for time instants with `Days` as their duration type
    /// (e.g., `JulianDate` and `ModifiedJulianDate`).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::time::{Period, ModifiedJulianDate};
    ///
    /// let start = ModifiedJulianDate::new(59000.0);
    /// let end = ModifiedJulianDate::new(59001.5);
    /// let period = Period::new(start, end);
    ///
    /// assert_eq!(period.duration_days(), 1.5);
    /// ```
    pub fn duration_days(&self) -> f64 {
        self.duration().value()
    }
}

// Specific implementation for UTC periods
impl Period<DateTime<Utc>> {
    /// Returns the duration in days as a floating-point value.
    ///
    /// This converts the chrono::Duration to days.
    pub fn duration_days(&self) -> f64 {
        self.duration().num_seconds() as f64 / 86400.0
    }

    /// Returns the duration in seconds.
    pub fn duration_seconds(&self) -> i64 {
        self.duration().num_seconds()
    }
}

// Serde support for Period<ModifiedJulianDate>
impl Serialize for Period<crate::time::ModifiedJulianDate> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut s = serializer.serialize_struct("Period", 2)?;
        s.serialize_field("start_mjd", &self.start.value())?;
        s.serialize_field("end_mjd", &self.end.value())?;
        s.end()
    }
}

impl<'de> Deserialize<'de> for Period<crate::time::ModifiedJulianDate> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct Raw {
            start_mjd: f64,
            end_mjd: f64,
        }

        let raw = Raw::deserialize(deserializer)?;
        Ok(Period::new(
            crate::time::ModifiedJulianDate::new(raw.start_mjd),
            crate::time::ModifiedJulianDate::new(raw.end_mjd),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time::{JulianDate, ModifiedJulianDate};

    #[test]
    fn test_period_creation_jd() {
        let start = JulianDate::new(2451545.0);
        let end = JulianDate::new(2451546.0);
        let period = Period::new(start, end);

        assert_eq!(period.start, start);
        assert_eq!(period.end, end);
    }

    #[test]
    fn test_period_creation_mjd() {
        let start = ModifiedJulianDate::new(59000.0);
        let end = ModifiedJulianDate::new(59001.0);
        let period = Period::new(start, end);

        assert_eq!(period.start, start);
        assert_eq!(period.end, end);
    }

    #[test]
    fn test_period_duration_jd() {
        let start = JulianDate::new(2451545.0);
        let end = JulianDate::new(2451546.5);
        let period = Period::new(start, end);

        assert_eq!(period.duration_days(), 1.5);
    }

    #[test]
    fn test_period_duration_mjd() {
        let start = ModifiedJulianDate::new(59000.0);
        let end = ModifiedJulianDate::new(59001.5);
        let period = Period::new(start, end);

        assert_eq!(period.duration_days(), 1.5);
    }

    #[test]
    fn test_period_duration_utc() {
        let start = DateTime::from_timestamp(0, 0).unwrap();
        let end = DateTime::from_timestamp(86400, 0).unwrap(); // 1 day later
        let period = Period::new(start, end);

        assert_eq!(period.duration_days(), 1.0);
        assert_eq!(period.duration_seconds(), 86400);
    }
}
