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
