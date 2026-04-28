// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Compatibility façade over `tempoch`.
//!
//! `siderust` historically exposed non-generic `JulianDate` /
//! `ModifiedJulianDate` wrappers plus `Period<T>` / `Interval<T>` over those
//! wrapper types. `tempoch 0.4.x` moved to generic encoded time aliases such as
//! `JulianDate<S> = EncodedTime<S, JD>`. This module keeps the older siderust
//! surface intact while delegating all actual time-scale logic to `tempoch`.

use chrono::{DateTime, Utc};

pub use tempoch::{
    constats, delta_t_seconds, delta_t_seconds_extrapolated, eop, ContinuousScale, ConversionError,
    ConversionTarget, CoordinateScale, EncodedTime, InfallibleConversionTarget, Interval,
    InvalidIntervalError, InvalidPeriodError, PeriodListError, Scale, ScaleKind, TAI, TCB, TCG,
    TDB, TT, Time, TimeContext, TimeDataError, UTC, UT1,
};

/// Backward-compatible alias: old siderust code used `UT` for the UT1 axis.
pub type UT = UT1;

/// Backward-compatible generic period alias over an instant type `T`.
pub type Period<T> = Interval<T>;

const JD_MINUS_MJD: f64 = 2_400_000.5;
const JD_J2000: f64 = 2_451_545.0;
const MJD_J2000: f64 = JD_J2000 - JD_MINUS_MJD;

/// TT-based Julian Date compatibility wrapper.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct JulianDate(f64);

impl JulianDate {
    pub const J2000: Self = Self(JD_J2000);
    pub const JULIAN_CENTURY: crate::qtty::Days = crate::qtty::Days::new(36_525.0);
    pub const JULIAN_YEAR: crate::qtty::Days = crate::qtty::Days::new(365.25);

    #[inline]
    pub const fn new(value: f64) -> Self {
        Self(value)
    }

    #[inline]
    pub const fn value(self) -> f64 {
        self.0
    }

    #[inline]
    pub fn from_utc(dt: DateTime<Utc>) -> Self {
        let tt = tempoch::Time::<UTC>::from_chrono(dt).to::<TT>();
        Self(tt.to::<tempoch::JD>().raw().value())
    }

    /// Construct a TT-aligned [`JulianDate`] directly from a [`tempoch::Time<UTC>`].
    ///
    /// Applies leap seconds via the same UTC→TT conversion path used by
    /// [`JulianDate::from_utc`], but without requiring a `chrono::DateTime`
    /// intermediary — downstream crates no longer need a direct `chrono` dependency
    /// just for this conversion.
    #[inline]
    pub fn from_tempoch_utc(time: tempoch::Time<UTC>) -> Self {
        let tt = time.to::<TT>();
        Self(tt.to::<tempoch::JD>().raw().value())
    }

    #[inline]
    pub fn to_utc(self) -> DateTime<Utc> {
        let tt: tempoch::JulianDate<TT> = self.into();
        tt.to_time()
            .to::<UTC>()
            .to_chrono()
            .expect("bundled TT->UTC conversions should stay within chrono's range")
    }

    #[inline]
    pub fn julian_centuries(self) -> f64 {
        (self.0 - JD_J2000) / 36_525.0
    }

    #[inline]
    pub fn julian_millennias(self) -> f64 {
        (self.0 - JD_J2000) / 365_250.0
    }

    #[inline]
    pub fn tt_to_tdb(jd_tt: Self) -> Self {
        let jd_tt: tempoch::JulianDate<TT> = jd_tt.into();
        let tdb = jd_tt.to_time().to::<TDB>();
        Self(tdb.to::<tempoch::JD>().raw().value())
    }

    #[inline]
    pub fn min(self, other: Self) -> Self {
        if self <= other { self } else { other }
    }

    #[inline]
    pub fn max(self, other: Self) -> Self {
        if self >= other { self } else { other }
    }

    #[inline]
    pub fn mean(self, other: Self) -> Self {
        Self((self.0 + other.0) * 0.5)
    }

    #[inline]
    pub fn quantity(self) -> crate::qtty::Days {
        crate::qtty::Days::new(self.0)
    }
}

impl core::fmt::Display for JulianDate {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Julian Day {:.9}", self.0)
    }
}

impl core::ops::Add<crate::qtty::Days> for JulianDate {
    type Output = Self;

    #[inline]
    fn add(self, rhs: crate::qtty::Days) -> Self {
        Self(self.0 + rhs.value())
    }
}

impl core::ops::Sub<crate::qtty::Days> for JulianDate {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: crate::qtty::Days) -> Self {
        Self(self.0 - rhs.value())
    }
}

impl core::ops::Sub for JulianDate {
    type Output = crate::qtty::Days;

    #[inline]
    fn sub(self, rhs: Self) -> crate::qtty::Days {
        crate::qtty::Days::new(self.0 - rhs.0)
    }
}

impl core::ops::Div<crate::qtty::Days> for JulianDate {
    type Output = f64;

    #[inline]
    fn div(self, rhs: crate::qtty::Days) -> Self::Output {
        self.0 / rhs.value()
    }
}

impl core::ops::AddAssign<crate::qtty::Days> for JulianDate {
    #[inline]
    fn add_assign(&mut self, rhs: crate::qtty::Days) {
        *self = *self + rhs;
    }
}

impl core::ops::SubAssign<crate::qtty::Days> for JulianDate {
    #[inline]
    fn sub_assign(&mut self, rhs: crate::qtty::Days) {
        *self = *self - rhs;
    }
}

/// TT-based Modified Julian Date compatibility wrapper.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct ModifiedJulianDate(f64);

impl ModifiedJulianDate {
    pub const J2000: Self = Self(MJD_J2000);

    #[inline]
    pub const fn new(value: f64) -> Self {
        Self(value)
    }

    #[inline]
    pub const fn value(self) -> f64 {
        self.0
    }

    #[inline]
    pub fn from_utc(dt: DateTime<Utc>) -> Self {
        let tt = tempoch::Time::<UTC>::from_chrono(dt).to::<TT>();
        Self(tt.to::<tempoch::MJD>().raw().value())
    }

    #[inline]
    pub fn to_utc(self) -> DateTime<Utc> {
        let tt: tempoch::ModifiedJulianDate<TT> = self.into();
        tt.to_time()
            .to::<UTC>()
            .to_chrono()
            .expect("bundled TT->UTC conversions should stay within chrono's range")
    }

    #[inline]
    pub fn min(self, other: Self) -> Self {
        if self <= other { self } else { other }
    }

    #[inline]
    pub fn max(self, other: Self) -> Self {
        if self >= other { self } else { other }
    }

    #[inline]
    pub fn mean(self, other: Self) -> Self {
        Self((self.0 + other.0) * 0.5)
    }

    #[inline]
    pub fn quantity(self) -> crate::qtty::Days {
        crate::qtty::Days::new(self.0)
    }
}

impl core::ops::Div<crate::qtty::Days> for ModifiedJulianDate {
    type Output = f64;

    #[inline]
    fn div(self, rhs: crate::qtty::Days) -> Self::Output {
        self.0 / rhs.value()
    }
}

impl core::fmt::Display for ModifiedJulianDate {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Modified Julian Day {:.9}", self.0)
    }
}

impl core::ops::Add<crate::qtty::Days> for ModifiedJulianDate {
    type Output = Self;

    #[inline]
    fn add(self, rhs: crate::qtty::Days) -> Self {
        Self(self.0 + rhs.value())
    }
}

impl core::ops::Sub<crate::qtty::Days> for ModifiedJulianDate {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: crate::qtty::Days) -> Self {
        Self(self.0 - rhs.value())
    }
}

impl core::ops::Sub for ModifiedJulianDate {
    type Output = crate::qtty::Days;

    #[inline]
    fn sub(self, rhs: Self) -> crate::qtty::Days {
        crate::qtty::Days::new(self.0 - rhs.0)
    }
}

impl core::ops::AddAssign<crate::qtty::Days> for ModifiedJulianDate {
    #[inline]
    fn add_assign(&mut self, rhs: crate::qtty::Days) {
        *self = *self + rhs;
    }
}

impl core::ops::SubAssign<crate::qtty::Days> for ModifiedJulianDate {
    #[inline]
    fn sub_assign(&mut self, rhs: crate::qtty::Days) {
        *self = *self - rhs;
    }
}

/// Backward-compatible aliases used throughout siderust.
pub type JD = JulianDate;
pub type MJD = ModifiedJulianDate;

impl From<ModifiedJulianDate> for JulianDate {
    #[inline]
    fn from(value: ModifiedJulianDate) -> Self {
        Self(value.0 + JD_MINUS_MJD)
    }
}

impl From<JulianDate> for ModifiedJulianDate {
    #[inline]
    fn from(value: JulianDate) -> Self {
        Self(value.0 - JD_MINUS_MJD)
    }
}

impl From<JulianDate> for tempoch::JulianDate<TT> {
    #[inline]
    fn from(value: JulianDate) -> Self {
        tempoch::JulianDate::<TT>::try_new(crate::qtty::Days::new(value.0))
            .expect("finite JulianDate should convert to tempoch::JulianDate<TT>")
    }
}

impl From<tempoch::JulianDate<TT>> for JulianDate {
    #[inline]
    fn from(value: tempoch::JulianDate<TT>) -> Self {
        Self(value.raw().value())
    }
}

impl From<ModifiedJulianDate> for tempoch::ModifiedJulianDate<TT> {
    #[inline]
    fn from(value: ModifiedJulianDate) -> Self {
        tempoch::ModifiedJulianDate::<TT>::try_new(crate::qtty::Days::new(value.0))
            .expect("finite ModifiedJulianDate should convert to tempoch::ModifiedJulianDate<TT>")
    }
}

impl From<tempoch::ModifiedJulianDate<TT>> for ModifiedJulianDate {
    #[inline]
    fn from(value: tempoch::ModifiedJulianDate<TT>) -> Self {
        Self(value.raw().value())
    }
}

/// Minimal instant trait kept for siderust's interval/root-finding helpers.
pub trait TimeInstant: Copy + PartialOrd {
    type Duration;

    fn difference(&self, other: &Self) -> Self::Duration;
    fn add_duration(&self, duration: Self::Duration) -> Self;
}

impl TimeInstant for JulianDate {
    type Duration = crate::qtty::Days;

    #[inline]
    fn difference(&self, other: &Self) -> Self::Duration {
        *self - *other
    }

    #[inline]
    fn add_duration(&self, duration: Self::Duration) -> Self {
        *self + duration
    }
}

impl TimeInstant for ModifiedJulianDate {
    type Duration = crate::qtty::Days;

    #[inline]
    fn difference(&self, other: &Self) -> Self::Duration {
        *self - *other
    }

    #[inline]
    fn add_duration(&self, duration: Self::Duration) -> Self {
        *self + duration
    }
}

#[inline]
pub fn complement_within<T: Copy + PartialOrd>(
    within: Interval<T>,
    periods: &[Interval<T>],
) -> Vec<Interval<T>> {
    within.complement(periods)
}

#[inline]
pub fn intersect_periods<T: Copy + PartialOrd>(a: &[Interval<T>], b: &[Interval<T>]) -> Vec<Interval<T>> {
    Interval::intersect_many(a, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 2024-01-01T00:00:00Z as a POSIX timestamp.
    const POSIX_2024_01_01: i64 = 1_704_067_200;

    #[test]
    fn from_tempoch_utc_matches_from_utc_for_known_instant() {
        let dt = chrono::DateTime::<chrono::Utc>::from_timestamp(POSIX_2024_01_01, 0)
            .expect("valid timestamp");
        let via_chrono = JulianDate::from_utc(dt);

        let utc_time = tempoch::Time::<UTC>::from_chrono(dt);
        let via_tempoch = JulianDate::from_tempoch_utc(utc_time);

        assert_eq!(
            via_chrono.value(),
            via_tempoch.value(),
            "from_tempoch_utc must produce identical JD to from_utc"
        );
    }

    /// Verify behaviour near a leap-second boundary (2016-12-31T23:59:59Z,
    /// one second before the 2017 leap second insertion).
    #[test]
    fn from_tempoch_utc_matches_from_utc_near_leap_second() {
        // 2016-12-31T23:59:59Z  →  POSIX 1483228799
        let posix_before_leap: i64 = 1_483_228_799;
        let dt = chrono::DateTime::<chrono::Utc>::from_timestamp(posix_before_leap, 0)
            .expect("valid timestamp");
        let via_chrono = JulianDate::from_utc(dt);

        let utc_time = tempoch::Time::<UTC>::from_chrono(dt);
        let via_tempoch = JulianDate::from_tempoch_utc(utc_time);

        assert_eq!(
            via_chrono.value(),
            via_tempoch.value(),
            "leap-second-adjacent instant must match from_utc exactly"
        );
    }
}
