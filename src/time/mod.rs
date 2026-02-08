// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Time Module
//!
//! This module provides time-related types and abstractions for astronomical calculations.
//!
//! # Core types
//!
//! - [`Time<S>`] — generic instant parameterised by a [`TimeScale`] marker.
//! - [`TimeScale`] — trait that defines a time scale (epoch offset + conversions).
//! - [`JulianDate`] — type alias for `Time<JD>`.
//! - [`ModifiedJulianDate`] — type alias for `Time<MJD>`.
//! - [`Period<T>`] — a time interval between two [`TimeInstant`]s.
//! - [`TimeInstant`] — trait for points in time usable with [`Period`].
//!
//! # Time scales
//!
//! The following markers implement [`TimeScale`]:
//!
//! | Marker | Scale |
//! |--------|-------|
//! | [`JD`] | Julian Date |
//! | [`MJD`] | Modified Julian Date |
//! | [`TDB`] | Barycentric Dynamical Time |
//! | [`TT`] | Terrestrial Time |
//! | [`TAI`] | International Atomic Time |
//! | [`GPS`] | GPS Time |
//! | [`UnixTime`] | Unix / POSIX time |

use chrono::{DateTime, Utc};
use qtty::Days;

pub(crate) mod scales;
pub(crate) mod time;
mod period;

// ── Re-exports ────────────────────────────────────────────────────────────

pub use time::{Time, TimeScale};
pub use scales::{JD, MJD, TDB, TT, TAI, GPS, UnixTime};
pub use period::{complement_within, intersect_periods, Period};

// ── Backward-compatible type aliases ──────────────────────────────────────

/// Julian Date — continuous count of days since the Julian Period.
///
/// This is a type alias for [`Time<JD>`].  All historical call-sites
/// (`JulianDate::new(...)`, `JulianDate::J2000`, `.julian_centuries()`, …)
/// continue to work without modification.
pub type JulianDate = Time<JD>;

/// Modified Julian Date — `JD − 2 400 000.5`.
///
/// This is a type alias for [`Time<MJD>`].
pub type ModifiedJulianDate = Time<MJD>;

// ═══════════════════════════════════════════════════════════════════════════
// TimeInstant trait
// ═══════════════════════════════════════════════════════════════════════════

/// Trait for types that represent a point in time.
///
/// Types implementing this trait can be used as time instants in [`Period<T>`]
/// and provide conversions to/from UTC and basic arithmetic operations.
pub trait TimeInstant: Copy + Clone + PartialEq + PartialOrd + Sized {
    /// The duration type used for arithmetic operations.
    type Duration;

    /// Convert this time instant to UTC DateTime.
    fn to_utc(&self) -> Option<DateTime<Utc>>;

    /// Create a time instant from UTC DateTime.
    fn from_utc(datetime: DateTime<Utc>) -> Self;

    /// Compute the difference between two time instants.
    fn difference(&self, other: &Self) -> Self::Duration;

    /// Add a duration to this time instant.
    fn add_duration(&self, duration: Self::Duration) -> Self;

    /// Subtract a duration from this time instant.
    fn sub_duration(&self, duration: Self::Duration) -> Self;
}

// ── Blanket impl for all Time<S> ──────────────────────────────────────────

impl<S: TimeScale> TimeInstant for Time<S> {
    type Duration = Days;

    #[inline]
    fn to_utc(&self) -> Option<DateTime<Utc>> {
        Time::to_utc(self)
    }

    #[inline]
    fn from_utc(datetime: DateTime<Utc>) -> Self {
        Time::from_utc(datetime)
    }

    #[inline]
    fn difference(&self, other: &Self) -> Self::Duration {
        *self - *other
    }

    #[inline]
    fn add_duration(&self, duration: Self::Duration) -> Self {
        *self + duration
    }

    #[inline]
    fn sub_duration(&self, duration: Self::Duration) -> Self {
        *self - duration
    }
}

// ── DateTime<Utc> keeps its own implementation ────────────────────────────

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

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

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

    #[test]
    fn cross_scale_conversion_via_to() {
        let jd = JulianDate::new(2_451_545.0);
        let mjd: ModifiedJulianDate = jd.to::<MJD>();
        assert!((mjd.value() - 51_544.5).abs() < 1e-10);

        let back: JulianDate = mjd.to::<JD>();
        assert!((back.value() - 2_451_545.0).abs() < 1e-10);
    }

    #[test]
    fn cross_scale_conversion_via_from() {
        let jd = JulianDate::new(2_451_545.0);
        let mjd: ModifiedJulianDate = jd.into();
        assert!((mjd.value() - 51_544.5).abs() < 1e-10);

        let back: JulianDate = mjd.into();
        assert!((back.value() - 2_451_545.0).abs() < 1e-10);
    }

    #[test]
    fn backward_compat_to_julian_day() {
        let mjd = ModifiedJulianDate::new(51_544.5);
        let jd = mjd.to_julian_day();
        assert_eq!(jd.value(), 2_451_545.0);
    }
}
