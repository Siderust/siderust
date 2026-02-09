// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Time Period implementation
//!
//! This module provides the `Period<T>` type, representing a time interval
//! between two time instants.

use super::TimeInstant;
use chrono::{DateTime, Utc};
use qtty::Days;

#[cfg(feature = "serde")]
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
    /// use qtty::Days;
    ///
    /// let start = ModifiedJulianDate::new(59000.0);
    /// let end = ModifiedJulianDate::new(59001.5);
    /// let period = Period::new(start, end);
    ///
    /// assert_eq!(period.duration_days(), Days::new(1.5));
    /// ```
    pub fn duration_days(&self) -> Days {
        self.duration()
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

// Serde support for Period<ModifiedJulianDate> (= Period<Time<MJD>>)
//
// Uses the historical field names `start_mjd` / `end_mjd` for backward
// compatibility with existing JSON reference data.
#[cfg(feature = "serde")]
impl Serialize for Period<crate::time::ModifiedJulianDate> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut s = serializer.serialize_struct("Period", 2)?;
        s.serialize_field("start_mjd", &self.start.quantity())?;
        s.serialize_field("end_mjd", &self.end.quantity())?;
        s.end()
    }
}

#[cfg(feature = "serde")]
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

// Serde support for Period<JulianDate> (= Period<Time<JD>>)
#[cfg(feature = "serde")]
impl Serialize for Period<crate::time::JulianDate> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut s = serializer.serialize_struct("Period", 2)?;
        s.serialize_field("start_jd", &self.start.quantity())?;
        s.serialize_field("end_jd", &self.end.quantity())?;
        s.end()
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for Period<crate::time::JulianDate> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct Raw {
            start_jd: f64,
            end_jd: f64,
        }

        let raw = Raw::deserialize(deserializer)?;
        Ok(Period::new(
            crate::time::JulianDate::new(raw.start_jd),
            crate::time::JulianDate::new(raw.end_jd),
        ))
    }
}

/// Returns the gaps (complement) of `periods` within the bounding `outer` period.
///
/// Given a sorted, non-overlapping list of sub-periods and a bounding period,
/// this returns the time intervals NOT covered by any sub-period.
///
/// Both `outer` and every element of `periods` must have `start <= end`.
/// The function runs in O(n) time with a single pass.
///
/// # Arguments
/// * `outer` - The bounding period
/// * `periods` - Sorted, non-overlapping sub-periods within `outer`
///
/// # Returns
/// The complement periods (gaps) in chronological order.
pub fn complement_within<T: TimeInstant>(
    outer: Period<T>,
    periods: &[Period<T>],
) -> Vec<Period<T>> {
    let mut gaps = Vec::new();
    let mut cursor = outer.start;
    for p in periods {
        if p.start > cursor {
            gaps.push(Period::new(cursor, p.start));
        }
        if p.end > cursor {
            cursor = p.end;
        }
    }
    if cursor < outer.end {
        gaps.push(Period::new(cursor, outer.end));
    }
    gaps
}

/// Returns the intersection of two sorted, non-overlapping period lists.
///
/// Uses an O(n+m) merge algorithm to find all overlapping spans.
///
/// # Arguments
/// * `a` - First sorted, non-overlapping period list
/// * `b` - Second sorted, non-overlapping period list
///
/// # Returns
/// Periods where both `a` and `b` overlap, in chronological order.
pub fn intersect_periods<T: TimeInstant>(a: &[Period<T>], b: &[Period<T>]) -> Vec<Period<T>> {
    let mut result = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        let start = if a[i].start >= b[j].start {
            a[i].start
        } else {
            b[j].start
        };
        let end = if a[i].end <= b[j].end {
            a[i].end
        } else {
            b[j].end
        };
        if start < end {
            result.push(Period::new(start, end));
        }
        if a[i].end <= b[j].end {
            i += 1;
        } else {
            j += 1;
        }
    }
    result
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

        assert_eq!(period.duration_days(), Days::new(1.5));
    }

    #[test]
    fn test_period_duration_mjd() {
        let start = ModifiedJulianDate::new(59000.0);
        let end = ModifiedJulianDate::new(59001.5);
        let period = Period::new(start, end);

        assert_eq!(period.duration_days(), Days::new(1.5));
    }

    #[test]
    fn test_period_duration_utc() {
        let start = DateTime::from_timestamp(0, 0).unwrap();
        let end = DateTime::from_timestamp(86400, 0).unwrap(); // 1 day later
        let period = Period::new(start, end);

        assert_eq!(period.duration_days(), 1.0);
        assert_eq!(period.duration_seconds(), 86400);
    }

    #[test]
    fn test_complement_within_gaps() {
        let outer = Period::new(ModifiedJulianDate::new(0.0), ModifiedJulianDate::new(10.0));
        let periods = vec![
            Period::new(ModifiedJulianDate::new(2.0), ModifiedJulianDate::new(4.0)),
            Period::new(ModifiedJulianDate::new(6.0), ModifiedJulianDate::new(8.0)),
        ];
        let gaps = complement_within(outer, &periods);
        assert_eq!(gaps.len(), 3);
        assert_eq!(gaps[0].start.value(), 0.0);
        assert_eq!(gaps[0].end.value(), 2.0);
        assert_eq!(gaps[1].start.value(), 4.0);
        assert_eq!(gaps[1].end.value(), 6.0);
        assert_eq!(gaps[2].start.value(), 8.0);
        assert_eq!(gaps[2].end.value(), 10.0);
    }

    #[test]
    fn test_complement_within_empty() {
        let outer = Period::new(ModifiedJulianDate::new(0.0), ModifiedJulianDate::new(10.0));
        let gaps = complement_within(outer, &[]);
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps[0].start.value(), 0.0);
        assert_eq!(gaps[0].end.value(), 10.0);
    }

    #[test]
    fn test_complement_within_full() {
        let outer = Period::new(ModifiedJulianDate::new(0.0), ModifiedJulianDate::new(10.0));
        let periods = vec![Period::new(
            ModifiedJulianDate::new(0.0),
            ModifiedJulianDate::new(10.0),
        )];
        let gaps = complement_within(outer, &periods);
        assert!(gaps.is_empty());
    }

    #[test]
    fn test_intersect_periods_overlap() {
        let a = vec![Period::new(
            ModifiedJulianDate::new(0.0),
            ModifiedJulianDate::new(5.0),
        )];
        let b = vec![Period::new(
            ModifiedJulianDate::new(3.0),
            ModifiedJulianDate::new(8.0),
        )];
        let overlap = intersect_periods(&a, &b);
        assert_eq!(overlap.len(), 1);
        assert_eq!(overlap[0].start.value(), 3.0);
        assert_eq!(overlap[0].end.value(), 5.0);
    }

    #[test]
    fn test_intersect_periods_no_overlap() {
        let a = vec![Period::new(
            ModifiedJulianDate::new(0.0),
            ModifiedJulianDate::new(3.0),
        )];
        let b = vec![Period::new(
            ModifiedJulianDate::new(5.0),
            ModifiedJulianDate::new(8.0),
        )];
        let overlap = intersect_periods(&a, &b);
        assert!(overlap.is_empty());
    }

    #[test]
    fn test_complement_intersect_roundtrip() {
        // above(min) ∩ complement(above(max)) = between(min, max)
        let outer = Period::new(ModifiedJulianDate::new(0.0), ModifiedJulianDate::new(10.0));
        let above_min = vec![
            Period::new(ModifiedJulianDate::new(1.0), ModifiedJulianDate::new(3.0)),
            Period::new(ModifiedJulianDate::new(5.0), ModifiedJulianDate::new(9.0)),
        ];
        let above_max = vec![
            Period::new(ModifiedJulianDate::new(2.0), ModifiedJulianDate::new(4.0)),
            Period::new(ModifiedJulianDate::new(7.0), ModifiedJulianDate::new(8.0)),
        ];
        let below_max = complement_within(outer, &above_max);
        let between = intersect_periods(&above_min, &below_max);
        // above_min: [1,3), [5,9)
        // above_max: [2,4), [7,8)
        // below_max (complement): [0,2), [4,7), [8,10)
        // intersection: [1,2), [5,7), [8,9)
        assert_eq!(between.len(), 3);
        assert_eq!(between[0].start.value(), 1.0);
        assert_eq!(between[0].end.value(), 2.0);
        assert_eq!(between[1].start.value(), 5.0);
        assert_eq!(between[1].end.value(), 7.0);
        assert_eq!(between[2].start.value(), 8.0);
        assert_eq!(between[2].end.value(), 9.0);
    }
}
