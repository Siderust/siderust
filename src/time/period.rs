// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Time period / interval implementation.
//!
//! This module provides:
//! - [`Interval<T>`]: generic interval over any [`TimeInstant`]
//! - [`Period<S>`]: scale-based alias for `Interval<Time<S>>`

use super::{Time, TimeInstant, TimeScale};
use chrono::{DateTime, Utc};
use qtty::Days;
use std::fmt;

#[cfg(feature = "serde")]
use serde::{ser::SerializeStruct, Deserialize, Deserializer, Serialize, Serializer};

/// Target type adapter for [`Interval<Time<S>>::to`].
///
/// This allows converting a period of `Time<S>` either to another time scale
/// marker (`MJD`, `JD`, `UT`, ...) or directly to `chrono::DateTime<Utc>`.
pub trait PeriodTimeTarget<S: TimeScale> {
    type Instant: TimeInstant;

    fn convert(value: Time<S>) -> Self::Instant;
}

impl<S: TimeScale, T: TimeScale> PeriodTimeTarget<S> for T {
    type Instant = Time<T>;

    #[inline]
    fn convert(value: Time<S>) -> Self::Instant {
        value.to::<T>()
    }
}

impl<S: TimeScale, T: TimeScale> PeriodTimeTarget<S> for Time<T> {
    type Instant = Time<T>;

    #[inline]
    fn convert(value: Time<S>) -> Self::Instant {
        value.to::<T>()
    }
}

impl<S: TimeScale> PeriodTimeTarget<S> for DateTime<Utc> {
    type Instant = DateTime<Utc>;

    #[inline]
    fn convert(value: Time<S>) -> Self::Instant {
        value
            .to_utc()
            .expect("time instant out of chrono::DateTime<Utc> representable range")
    }
}

/// Target type adapter for [`Interval<DateTime<Utc>>::to`].
pub trait PeriodUtcTarget {
    type Instant: TimeInstant;

    fn convert(value: DateTime<Utc>) -> Self::Instant;
}

impl<S: TimeScale> PeriodUtcTarget for S {
    type Instant = Time<S>;

    #[inline]
    fn convert(value: DateTime<Utc>) -> Self::Instant {
        Time::<S>::from_utc(value)
    }
}

impl<S: TimeScale> PeriodUtcTarget for Time<S> {
    type Instant = Time<S>;

    #[inline]
    fn convert(value: DateTime<Utc>) -> Self::Instant {
        Time::<S>::from_utc(value)
    }
}

impl PeriodUtcTarget for DateTime<Utc> {
    type Instant = DateTime<Utc>;

    #[inline]
    fn convert(value: DateTime<Utc>) -> Self::Instant {
        value
    }
}

/// Represents an interval between two instants.
///
/// An `Interval` is defined by a start and end time instant of type `T`,
/// where `T` implements the `TimeInstant` trait. This allows for periods
/// defined in different time systems (Julian Date, Modified Julian Date, UTC, etc.).
///
/// # Examples
///
/// ```
/// use siderust::time::{Interval, ModifiedJulianDate};
///
/// let start = ModifiedJulianDate::new(59000.0);
/// let end = ModifiedJulianDate::new(59001.0);
/// let period = Interval::new(start, end);
///
/// // Duration in days
/// let duration = period.duration();
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval<T: TimeInstant> {
    pub start: T,
    pub end: T,
}

/// Time-scale period alias.
///
/// This follows the same marker pattern as [`Time<S>`], so callers can write:
/// `Period<MJD>`, `Period<JD>`, etc.
pub type Period<S> = Interval<Time<S>>;

/// UTC interval alias.
pub type UtcPeriod = Interval<DateTime<Utc>>;

impl<T: TimeInstant> Interval<T> {
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
    /// use siderust::time::{Interval, JulianDate};
    ///
    /// let start = JulianDate::new(2451545.0);
    /// let end = JulianDate::new(2451546.0);
    /// let period = Interval::new(start, end);
    /// ```
    pub fn new(start: T, end: T) -> Self {
        Interval { start, end }
    }

    /// Returns the duration of the period as the difference between end and start.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::time::{Interval, JulianDate};
    /// use qtty::Days;
    ///
    /// let start = JulianDate::new(2451545.0);
    /// let end = JulianDate::new(2451546.5);
    /// let period = Interval::new(start, end);
    ///
    /// let duration = period.duration();
    /// assert_eq!(duration, Days::new(1.5));
    /// ```
    pub fn duration(&self) -> T::Duration {
        self.end.difference(&self.start)
    }

    /// Returns the overlapping sub-period between `self` and `other`.
    ///
    /// Periods are treated as half-open ranges `[start, end)`: if one period
    /// ends exactly when the other starts, the intersection is empty and `None`
    /// is returned.
    pub fn intersection(&self, other: &Self) -> Option<Self> {
        let start = if self.start >= other.start {
            self.start
        } else {
            other.start
        };
        let end = if self.end <= other.end {
            self.end
        } else {
            other.end
        };

        if start < end {
            Some(Self::new(start, end))
        } else {
            None
        }
    }
}

// Display implementation
impl<T: TimeInstant + fmt::Display> fmt::Display for Interval<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} to {}", self.start, self.end)
    }
}

impl<S: TimeScale> Interval<Time<S>> {
    /// Convert this period to another time scale.
    ///
    /// Each endpoint is converted preserving the represented absolute interval.
    ///
    /// Supported targets:
    /// - Any time-scale marker (`JD`, `MJD`, `UT`, ...)
    /// - `chrono::DateTime<Utc>`
    ///
    /// # Examples
    ///
    /// ```
    /// use chrono::{DateTime, Utc};
    /// use siderust::time::{Interval, JD, MJD, Period, Time};
    ///
    /// let period_jd = Period::new(Time::<JD>::new(2451545.0), Time::<JD>::new(2451546.0));
    /// let period_mjd = period_jd.to::<MJD>();
    /// let _period_utc: Interval<DateTime<Utc>> = period_jd.to::<DateTime<Utc>>();
    ///
    /// assert!((period_mjd.start.value() - 51544.5).abs() < 1e-12);
    /// assert!((period_mjd.end.value() - 51545.5).abs() < 1e-12);
    /// ```
    #[inline]
    pub fn to<Target>(&self) -> Interval<<Target as PeriodTimeTarget<S>>::Instant>
    where
        Target: PeriodTimeTarget<S>,
    {
        Interval::new(Target::convert(self.start), Target::convert(self.end))
    }
}

// Specific implementation for periods with Days duration (JD and MJD)
impl<T: TimeInstant<Duration = Days>> Interval<T> {
    /// Returns the duration of the period in days as a floating-point value.
    ///
    /// This method is available for time instants with `Days` as their duration type
    /// (e.g., `JulianDate` and `ModifiedJulianDate`).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::time::{Interval, ModifiedJulianDate};
    /// use qtty::Days;
    ///
    /// let start = ModifiedJulianDate::new(59000.0);
    /// let end = ModifiedJulianDate::new(59001.5);
    /// let period = Interval::new(start, end);
    ///
    /// assert_eq!(period.duration_days(), Days::new(1.5));
    /// ```
    pub fn duration_days(&self) -> Days {
        self.duration()
    }
}

// Specific implementation for UTC periods
impl Interval<DateTime<Utc>> {
    /// Convert this UTC interval to another target.
    ///
    /// Supported targets:
    /// - Any time-scale marker (`JD`, `MJD`, `UT`, ...)
    /// - Any `Time<...>` alias (`JulianDate`, `ModifiedJulianDate`, ...)
    /// - `chrono::DateTime<Utc>`
    #[inline]
    pub fn to<Target>(&self) -> Interval<<Target as PeriodUtcTarget>::Instant>
    where
        Target: PeriodUtcTarget,
    {
        Interval::new(Target::convert(self.start), Target::convert(self.end))
    }

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

// Serde support for Period<MJD> (= Interval<Time<MJD>>)
//
// Uses the historical field names `start_mjd` / `end_mjd` for backward
// compatibility with existing JSON reference data.
#[cfg(feature = "serde")]
impl Serialize for Interval<crate::time::ModifiedJulianDate> {
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
impl<'de> Deserialize<'de> for Interval<crate::time::ModifiedJulianDate> {
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
        Ok(Interval::new(
            crate::time::ModifiedJulianDate::new(raw.start_mjd),
            crate::time::ModifiedJulianDate::new(raw.end_mjd),
        ))
    }
}

// Serde support for Period<JD> (= Interval<Time<JD>>)
#[cfg(feature = "serde")]
impl Serialize for Interval<crate::time::JulianDate> {
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
impl<'de> Deserialize<'de> for Interval<crate::time::JulianDate> {
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
        Ok(Interval::new(
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
    outer: Interval<T>,
    periods: &[Interval<T>],
) -> Vec<Interval<T>> {
    let mut gaps = Vec::new();
    let mut cursor = outer.start;
    for p in periods {
        if p.start > cursor {
            gaps.push(Interval::new(cursor, p.start));
        }
        if p.end > cursor {
            cursor = p.end;
        }
    }
    if cursor < outer.end {
        gaps.push(Interval::new(cursor, outer.end));
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
pub fn intersect_periods<T: TimeInstant>(a: &[Interval<T>], b: &[Interval<T>]) -> Vec<Interval<T>> {
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
            result.push(Interval::new(start, end));
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
    use crate::time::{JulianDate, ModifiedJulianDate, JD, MJD};

    #[test]
    fn test_period_creation_jd() {
        let start = JulianDate::new(2451545.0);
        let end = JulianDate::new(2451546.0);
        let period = Period::new(start, end);

        assert_eq!(period.start, start);
        assert_eq!(period.end, end);
    }

    #[test]
    fn test_period_scale_conversion_jd_to_mjd() {
        let period_jd = Period::new(Time::<JD>::new(2_451_545.0), Time::<JD>::new(2_451_546.0));
        let period_mjd = period_jd.to::<MJD>();

        assert!((period_mjd.start.value() - 51_544.5).abs() < 1e-12);
        assert!((period_mjd.end.value() - 51_545.5).abs() < 1e-12);
    }

    #[test]
    fn test_period_scale_conversion_roundtrip() {
        let original = Period::new(Time::<MJD>::new(59_000.125), Time::<MJD>::new(59_001.75));
        let roundtrip = original.to::<JD>().to::<MJD>();

        assert!((roundtrip.start.value() - original.start.value()).abs() < 1e-12);
        assert!((roundtrip.end.value() - original.end.value()).abs() < 1e-12);
    }

    #[test]
    fn test_period_scale_conversion_to_utc() {
        let start_utc = DateTime::from_timestamp(1_700_000_000, 0).unwrap();
        let end_utc = DateTime::from_timestamp(1_700_000_600, 0).unwrap();
        let period_jd = Period::new(
            Time::<JD>::from_utc(start_utc),
            Time::<JD>::from_utc(end_utc),
        );

        let period_utc = period_jd.to::<DateTime<Utc>>();
        let start_delta_ns = period_utc.start.timestamp_nanos_opt().unwrap()
            - start_utc.timestamp_nanos_opt().unwrap();
        let end_delta_ns =
            period_utc.end.timestamp_nanos_opt().unwrap() - end_utc.timestamp_nanos_opt().unwrap();
        assert!(start_delta_ns.abs() < 10_000);
        assert!(end_delta_ns.abs() < 10_000);
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
        let period = Interval::new(start, end);

        assert_eq!(period.duration_days(), 1.0);
        assert_eq!(period.duration_seconds(), 86400);
    }

    #[test]
    fn test_period_to_conversion() {
        let mjd_start = ModifiedJulianDate::new(59000.0);
        let mjd_end = ModifiedJulianDate::new(59001.0);
        let mjd_period = Period::new(mjd_start, mjd_end);

        let utc_period = mjd_period.to::<DateTime<Utc>>();

        // The converted period should have approximately the same duration (within 1 second due to ΔT)
        let duration_secs = utc_period.duration().num_seconds();
        assert!(
            (duration_secs - 86400).abs() <= 1,
            "Duration was {} seconds",
            duration_secs
        );

        // Convert back and check that it's close (within small tolerance due to floating point)
        let back_to_mjd = utc_period.to::<ModifiedJulianDate>();
        let start_diff = (back_to_mjd.start.quantity() - mjd_start.quantity())
            .value()
            .abs();
        let end_diff = (back_to_mjd.end.quantity() - mjd_end.quantity())
            .value()
            .abs();
        assert!(start_diff < 1e-6, "Start difference: {}", start_diff);
        assert!(end_diff < 1e-6, "End difference: {}", end_diff);
    }

    #[test]
    fn test_period_display() {
        let start = ModifiedJulianDate::new(59000.0);
        let end = ModifiedJulianDate::new(59001.0);
        let period = Period::new(start, end);

        let display = format!("{}", period);
        assert!(display.contains("MJD 59000"));
        assert!(display.contains("MJD 59001"));
        assert!(display.contains("to"));
    }

    #[test]
    fn test_period_intersection_overlap() {
        let a = Period::new(ModifiedJulianDate::new(0.0), ModifiedJulianDate::new(5.0));
        let b = Period::new(ModifiedJulianDate::new(3.0), ModifiedJulianDate::new(8.0));

        let overlap = a.intersection(&b).expect("expected overlap");
        assert_eq!(overlap.start.quantity(), Days::new(3.0));
        assert_eq!(overlap.end.quantity(), Days::new(5.0));
    }

    #[test]
    fn test_period_intersection_disjoint() {
        let a = Period::new(ModifiedJulianDate::new(0.0), ModifiedJulianDate::new(3.0));
        let b = Period::new(ModifiedJulianDate::new(5.0), ModifiedJulianDate::new(8.0));

        assert_eq!(a.intersection(&b), None);
    }

    #[test]
    fn test_period_intersection_touching_edges() {
        let a = Period::new(ModifiedJulianDate::new(0.0), ModifiedJulianDate::new(3.0));
        let b = Period::new(ModifiedJulianDate::new(3.0), ModifiedJulianDate::new(8.0));

        assert_eq!(a.intersection(&b), None);
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
        assert_eq!(gaps[0].start.quantity(), Days::new(0.0));
        assert_eq!(gaps[0].end.quantity(), Days::new(2.0));
        assert_eq!(gaps[1].start.quantity(), Days::new(4.0));
        assert_eq!(gaps[1].end.quantity(), Days::new(6.0));
        assert_eq!(gaps[2].start.quantity(), Days::new(8.0));
        assert_eq!(gaps[2].end.quantity(), Days::new(10.0));
    }

    #[test]
    fn test_complement_within_empty() {
        let outer = Period::new(ModifiedJulianDate::new(0.0), ModifiedJulianDate::new(10.0));
        let gaps = complement_within(outer, &[]);
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps[0].start.quantity(), Days::new(0.0));
        assert_eq!(gaps[0].end.quantity(), Days::new(10.0));
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
        assert_eq!(overlap[0].start.quantity(), Days::new(3.0));
        assert_eq!(overlap[0].end.quantity(), Days::new(5.0));
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
        assert_eq!(between[0].start.quantity(), Days::new(1.0));
        assert_eq!(between[0].end.quantity(), Days::new(2.0));
        assert_eq!(between[1].start.quantity(), Days::new(5.0));
        assert_eq!(between[1].end.quantity(), Days::new(7.0));
        assert_eq!(between[2].start.quantity(), Days::new(8.0));
        assert_eq!(between[2].end.quantity(), Days::new(9.0));
    }
}
