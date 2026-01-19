//! Time Module
//!
//! This module provides time-related types and abstractions for astronomical calculations.
//! It includes:
//! - `JulianDate` (JD): Continuous count of days since the Julian Period
//! - `ModifiedJulianDate` (MJD): Julian Date minus 2400000.5
//! - `Period<T>`: A generic time period/interval between two time instants
//! - `TimeInstant`: Trait for types representing a point in time

use chrono::{DateTime, Utc};
use qtty::Days;

mod julian_date;
mod modified_julian_date;
mod period;

pub use julian_date::JulianDate;
pub use modified_julian_date::ModifiedJulianDate;
pub use period::Period;

/// Trait for types that represent a point in time.
///
/// Types implementing this trait can be used as time instants in `Period<T>`
/// and provide conversions to/from UTC and basic arithmetic operations.
pub trait TimeInstant: Copy + Clone + PartialEq + PartialOrd + Sized {
    /// The duration type used for arithmetic operations
    type Duration;

    /// Convert this time instant to UTC DateTime
    fn to_utc(&self) -> Option<DateTime<Utc>>;

    /// Create a time instant from UTC DateTime
    fn from_utc(datetime: DateTime<Utc>) -> Self;

    /// Compute the difference between two time instants
    fn difference(&self, other: &Self) -> Self::Duration;

    /// Add a duration to this time instant
    fn add_duration(&self, duration: Self::Duration) -> Self;

    /// Subtract a duration from this time instant
    fn sub_duration(&self, duration: Self::Duration) -> Self;
}

impl TimeInstant for JulianDate {
    type Duration = Days;

    fn to_utc(&self) -> Option<DateTime<Utc>> {
        JulianDate::to_utc(self)
    }

    fn from_utc(datetime: DateTime<Utc>) -> Self {
        JulianDate::from_utc(datetime)
    }

    fn difference(&self, other: &Self) -> Self::Duration {
        *self - *other
    }

    fn add_duration(&self, duration: Self::Duration) -> Self {
        *self + duration
    }

    fn sub_duration(&self, duration: Self::Duration) -> Self {
        *self - duration
    }
}

impl TimeInstant for ModifiedJulianDate {
    type Duration = Days;

    fn to_utc(&self) -> Option<DateTime<Utc>> {
        ModifiedJulianDate::to_utc(self)
    }

    fn from_utc(datetime: DateTime<Utc>) -> Self {
        ModifiedJulianDate::from_utc(datetime)
    }

    fn difference(&self, other: &Self) -> Self::Duration {
        *self - *other
    }

    fn add_duration(&self, duration: Self::Duration) -> Self {
        *self + duration
    }

    fn sub_duration(&self, duration: Self::Duration) -> Self {
        *self - duration
    }
}

impl TimeInstant for DateTime<Utc> {
    type Duration = chrono::Duration;

    fn to_utc(&self) -> Option<DateTime<Utc>> {
        Some(*self)
    }

    fn from_utc(datetime: DateTime<Utc>) -> Self {
        datetime
    }

    fn difference(&self, other: &Self) -> Self::Duration {
        *self - *other
    }

    fn add_duration(&self, duration: Self::Duration) -> Self {
        *self + duration
    }

    fn sub_duration(&self, duration: Self::Duration) -> Self {
        *self - duration
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::TimeZone;

    #[test]
    fn timeinstant_for_julian_date_handles_arithmetic() {
        let jd = JulianDate::new(2_451_545.0);
        let other = jd + Days::new(2.0);

        assert_eq!(jd.difference(&other), Days::new(-2.0));
        assert_eq!(jd.add_duration(Days::new(1.5)).value(), 2_451_546.5);
        assert_eq!(other.sub_duration(Days::new(0.5)).value(), 2_451_546.5);
    }

    #[test]
    fn timeinstant_for_modified_julian_date_roundtrips_utc() {
        let dt = DateTime::from_timestamp(946_684_800, 123_000_000).unwrap(); // 2000-01-01T00:00:00.123Z
        let mjd = ModifiedJulianDate::from_utc(dt);
        let back = mjd.to_utc().expect("mjd to utc");

        assert_eq!(mjd.difference(&mjd), Days::new(0.0));
        assert_eq!(mjd.add_duration(Days::new(1.0)).value(), mjd.value() + 1.0);
        assert_eq!(mjd.sub_duration(Days::new(0.5)).value(), mjd.value() - 0.5);
        let delta_ns = back.timestamp_nanos_opt().unwrap() - dt.timestamp_nanos_opt().unwrap();
        assert!(delta_ns.abs() < 10_000, "nanos differ by {}", delta_ns);
    }

    #[test]
    fn timeinstant_for_datetime_uses_chrono_durations() {
        let base = Utc.with_ymd_and_hms(2024, 1, 1, 0, 0, 0).unwrap();
        let later = Utc.with_ymd_and_hms(2024, 1, 2, 6, 0, 0).unwrap();
        let diff = later.difference(&base);

        assert_eq!(diff.num_hours(), 30);
        assert_eq!(
            base.add_duration(diff + chrono::Duration::hours(6)),
            later + chrono::Duration::hours(6)
        );
        assert_eq!(later.sub_duration(diff), base);
        assert_eq!(TimeInstant::to_utc(&later), Some(later));
    }
}
