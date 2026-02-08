// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Generic time–scale parameterised instant.
//!
//! [`Time<S>`] is the core type of the time module.  It stores a scalar
//! quantity in [`Days`] whose *meaning* is determined by the compile-time
//! marker `S: TimeScale`.  All arithmetic (addition/subtraction of
//! durations, difference between instants), UTC conversion, serialisation,
//! and display are implemented generically — no code duplication.
//!
//! Domain-specific methods that only make sense for a particular scale
//! (e.g. [`Time::<JD>::julian_centuries()`]) are placed in inherent `impl`
//! blocks gated on the concrete marker type.

use chrono::{DateTime, Utc};
use qtty::*;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Sub, SubAssign};

#[cfg(feature = "serde")]
use serde::{Deserialize, Deserializer, Serialize, Serializer};

// ═══════════════════════════════════════════════════════════════════════════
// TimeScale trait
// ═══════════════════════════════════════════════════════════════════════════

/// Marker trait for time scales.
///
/// A **time scale** defines:
///
/// 1. A human-readable **label** (e.g. `"JD"`, `"MJD"`, `"TAI"`).
/// 2. A pair of conversion functions between the scale's native quantity
///    (in [`Days`]) and **Julian Date in TT** (JD(TT)) — the canonical
///    internal representation used throughout the crate.
///
/// For pure *epoch counters* (JD, MJD, Unix Time, GPS) the conversions are
/// trivial constant offsets that the compiler will inline and fold away.
///
/// For *physical scales* (TT, TDB, TAI) the conversions may include
/// function-based corrections (e.g. the ≈1.7 ms TDB↔TT periodic term).
pub trait TimeScale: Copy + Clone + std::fmt::Debug + PartialEq + PartialOrd + 'static {
    /// Display label used by [`Time`] formatting.
    const LABEL: &'static str;

    /// Convert a quantity in this scale's native unit to an absolute JD(TT).
    fn to_jd_tt(value: Days) -> Days;

    /// Convert an absolute JD(TT) back to this scale's native quantity.
    fn from_jd_tt(jd_tt: Days) -> Days;
}

// ═══════════════════════════════════════════════════════════════════════════
// Time<S> — the generic instant
// ═══════════════════════════════════════════════════════════════════════════

/// A point on time scale `S`.
///
/// Internally stores a single `Days` quantity whose interpretation depends on
/// `S: TimeScale`.  The struct is `Copy` and zero-cost: `PhantomData` is
/// zero-sized, so `Time<S>` is layout-identical to `Days` (a single `f64`).
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct Time<S: TimeScale> {
    quantity: Days,
    _scale: PhantomData<S>,
}

impl<S: TimeScale> Time<S> {
    // ── constructors ──────────────────────────────────────────────────

    /// Create from a raw scalar (days since the scale's epoch).
    #[inline]
    pub const fn new(value: f64) -> Self {
        Self {
            quantity: Days::new(value),
            _scale: PhantomData,
        }
    }

    /// Create from a [`Days`] quantity.
    #[inline]
    pub const fn from_days(days: Days) -> Self {
        Self {
            quantity: days,
            _scale: PhantomData,
        }
    }

    // ── accessors ─────────────────────────────────────────────────────

    /// The underlying quantity in days.
    #[inline]
    pub const fn quantity(&self) -> Days {
        self.quantity
    }

    /// The underlying scalar value in days.
    #[inline]
    pub const fn value(&self) -> f64 {
        self.quantity.value()
    }

    /// Absolute Julian Day (TT) corresponding to this instant.
    #[inline]
    pub fn julian_day(&self) -> Days {
        S::to_jd_tt(self.quantity)
    }

    /// Absolute Julian Day (TT) as scalar.
    #[inline]
    pub fn julian_day_value(&self) -> f64 {
        self.julian_day().value()
    }

    /// Build an instant from an absolute Julian Day (TT).
    #[inline]
    pub fn from_julian_day(jd: Days) -> Self {
        Self::from_days(S::from_jd_tt(jd))
    }

    // ── cross-scale conversion (mirroring qtty's .to::<T>()) ─────────

    /// Convert this instant to another time scale.
    ///
    /// The conversion routes through the canonical JD(TT) intermediate:
    ///
    /// ```text
    /// self → JD(TT) → target
    /// ```
    ///
    /// For pure epoch-offset scales this compiles down to a single
    /// addition/subtraction.
    #[inline]
    pub fn to<T: TimeScale>(&self) -> Time<T> {
        Time::<T>::from_julian_day(S::to_jd_tt(self.quantity))
    }

    // ── UTC helpers ───────────────────────────────────────────────────

    /// Convert to a `chrono::DateTime<Utc>`.
    ///
    /// Returns `None` if the value falls outside chrono's representable range.
    pub fn to_utc(&self) -> Option<DateTime<Utc>> {
        const UNIX_EPOCH_JD: f64 = 2_440_587.5;
        let seconds_since_epoch = (self.julian_day_value() - UNIX_EPOCH_JD) * 86_400.0;
        let secs = seconds_since_epoch.floor() as i64;
        let nanos = ((seconds_since_epoch - secs as f64) * 1e9) as u32;
        DateTime::<Utc>::from_timestamp(secs, nanos)
    }

    /// Build an instant from a `chrono::DateTime<Utc>`.
    pub fn from_utc(datetime: DateTime<Utc>) -> Self {
        const UNIX_EPOCH_JD: f64 = 2_440_587.5;
        let seconds_since_epoch = datetime.timestamp() as f64;
        let nanos = datetime.timestamp_subsec_nanos() as f64 / 1e9;
        let jd = UNIX_EPOCH_JD + (seconds_since_epoch + nanos) / 86_400.0;
        Self::from_julian_day(Days::new(jd))
    }

    // ── min / max ─────────────────────────────────────────────────────

    /// Element-wise minimum.
    #[inline]
    pub const fn min(self, other: Self) -> Self {
        Self::from_days(self.quantity.min_const(other.quantity))
    }

    /// Element-wise maximum.
    #[inline]
    pub const fn max(self, other: Self) -> Self {
        Self::from_days(self.quantity.max_const(other.quantity))
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Generic trait implementations
// ═══════════════════════════════════════════════════════════════════════════

// ── Display ───────────────────────────────────────────────────────────────

impl<S: TimeScale> std::fmt::Display for Time<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {}", S::LABEL, self.quantity)
    }
}

// ── Serde ─────────────────────────────────────────────────────────────────

#[cfg(feature = "serde")]
impl<S: TimeScale> Serialize for Time<S> {
    fn serialize<Ser>(&self, serializer: Ser) -> Result<Ser::Ok, Ser::Error>
    where
        Ser: Serializer,
    {
        serializer.serialize_f64(self.value())
    }
}

#[cfg(feature = "serde")]
impl<'de, S: TimeScale> Deserialize<'de> for Time<S> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let v = f64::deserialize(deserializer)?;
        Ok(Self::new(v))
    }
}

// ── Arithmetic ────────────────────────────────────────────────────────────

impl<S: TimeScale> Add<Days> for Time<S> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Days) -> Self::Output {
        Self::from_days(self.quantity + rhs)
    }
}

impl<S: TimeScale> AddAssign<Days> for Time<S> {
    #[inline]
    fn add_assign(&mut self, rhs: Days) {
        self.quantity += rhs;
    }
}

impl<S: TimeScale> Sub<Days> for Time<S> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Days) -> Self::Output {
        Self::from_days(self.quantity - rhs)
    }
}

impl<S: TimeScale> SubAssign<Days> for Time<S> {
    #[inline]
    fn sub_assign(&mut self, rhs: Days) {
        self.quantity -= rhs;
    }
}

impl<S: TimeScale> Sub for Time<S> {
    type Output = Days;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        self.quantity - rhs.quantity
    }
}

impl<S: TimeScale> std::ops::Div<Days> for Time<S> {
    type Output = f64;
    #[inline]
    fn div(self, rhs: Days) -> Self::Output {
        self.value() / rhs.value()
    }
}

impl<S: TimeScale> std::ops::Div<f64> for Time<S> {
    type Output = f64;
    #[inline]
    fn div(self, rhs: f64) -> Self::Output {
        self.value() / rhs
    }
}

// ── From/Into Days ────────────────────────────────────────────────────────

impl<S: TimeScale> From<Days> for Time<S> {
    #[inline]
    fn from(days: Days) -> Self {
        Self::from_days(days)
    }
}

impl<S: TimeScale> From<Time<S>> for Days {
    #[inline]
    fn from(time: Time<S>) -> Self {
        time.quantity
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::scales::{JD, MJD};

    #[test]
    fn test_julian_day_creation() {
        let jd = Time::<JD>::new(2_451_545.0);
        assert_eq!(jd.value(), 2_451_545.0);
    }

    #[test]
    fn test_to_utc() {
        let jd = Time::<JD>::new(2_451_545.0);
        let datetime = jd.to_utc();
        assert_eq!(datetime, DateTime::from_timestamp(946_728_000, 0));
    }

    #[test]
    fn test_from_utc() {
        let datetime = DateTime::from_timestamp(946_728_000, 0).unwrap();
        let jd = Time::<JD>::from_utc(datetime);
        assert_eq!(jd.value(), 2_451_545.0);
    }

    #[test]
    fn test_julian_conversions() {
        let jd = Time::<JD>::J2000 + Days::new(365_250.0);
        assert!((jd.julian_millennias() - Millennia::new(1.0)).abs() < 1e-12);
        assert!((jd.julian_centuries().value() - 10.0).abs() < 1e-12);
        assert!((jd.julian_years().value() - 1000.0).abs() < 1e-9);
    }

    #[test]
    fn test_tt_to_tdb_and_min() {
        let jd_tdb = Time::<JD>::tt_to_tdb(Time::<JD>::J2000);
        assert!((jd_tdb - Time::<JD>::J2000).value().abs() < 1e-6);

        let later = Time::<JD>::J2000 + Days::new(1.0);
        assert_eq!(Time::<JD>::J2000.min(later), Time::<JD>::J2000);
    }

    #[test]
    fn test_into_days() {
        let jd = Time::<JD>::new(2_451_547.5);
        let days: Days = jd.into();
        assert_eq!(days.value(), 2_451_547.5);

        let roundtrip = Time::<JD>::from(days);
        assert_eq!(roundtrip, jd);
    }

    #[test]
    fn test_into_julian_years() {
        let jd = Time::<JD>::J2000 + Days::new(365.25 * 2.0);
        let years: JulianYears = jd.into();
        assert!((years.value() - 2.0).abs() < 1e-12);

        let roundtrip = Time::<JD>::from(years);
        assert!((roundtrip.value() - jd.value()).abs() < 1e-12);
    }

    #[test]
    fn test_into_centuries() {
        let jd = Time::<JD>::J2000 + Days::new(36_525.0 * 3.0);
        let centuries: Centuries = jd.into();
        assert!((centuries.value() - 3.0).abs() < 1e-12);

        let roundtrip = Time::<JD>::from(centuries);
        assert!((roundtrip.value() - jd.value()).abs() < 1e-12);
    }

    #[test]
    fn test_into_millennia() {
        let jd = Time::<JD>::J2000 + Days::new(365_250.0 * 1.5);
        let millennia: Millennia = jd.into();
        assert!((millennia.value() - 1.5).abs() < 1e-12);

        let roundtrip = Time::<JD>::from(millennia);
        assert!((roundtrip.value() - jd.value()).abs() < 1e-9);
    }

    #[test]
    fn test_mjd_creation() {
        let mjd = Time::<MJD>::new(51_544.5);
        assert_eq!(mjd.value(), 51_544.5);
    }

    #[test]
    fn test_mjd_into_jd() {
        let mjd = Time::<MJD>::new(51_544.5);
        let jd: Time<JD> = mjd.into();
        assert_eq!(jd.value(), 2_451_545.0);
    }

    #[test]
    fn test_mjd_to_utc() {
        let mjd = Time::<MJD>::new(51_544.5);
        let datetime = mjd.to_utc();
        assert_eq!(datetime, DateTime::from_timestamp(946_728_000, 0));
    }

    #[test]
    fn test_mjd_from_utc() {
        let datetime = DateTime::from_timestamp(946_728_000, 0).unwrap();
        let mjd = Time::<MJD>::from_utc(datetime);
        assert_eq!(mjd.value(), 51_544.5);
    }

    #[test]
    fn test_mjd_add_days() {
        let mjd = Time::<MJD>::new(59_000.0);
        let result = mjd + Days::new(1.5);
        assert_eq!(result.value(), 59_001.5);
    }

    #[test]
    fn test_mjd_sub_days() {
        let mjd = Time::<MJD>::new(59_000.0);
        let result = mjd - Days::new(1.5);
        assert_eq!(result.value(), 58_998.5);
    }

    #[test]
    fn test_mjd_sub_mjd() {
        let mjd1 = Time::<MJD>::new(59_001.0);
        let mjd2 = Time::<MJD>::new(59_000.0);
        let diff = mjd1 - mjd2;
        assert_eq!(diff.value(), 1.0);
    }

    #[test]
    fn test_mjd_comparison() {
        let mjd1 = Time::<MJD>::new(59_000.0);
        let mjd2 = Time::<MJD>::new(59_001.0);
        assert!(mjd1 < mjd2);
        assert!(mjd2 > mjd1);
    }

    #[test]
    fn test_display_jd() {
        let jd = Time::<JD>::new(2_451_545.0);
        let s = format!("{jd}");
        assert!(s.contains("Julian Day"));
    }

    #[test]
    fn test_display_mjd() {
        let mjd = Time::<MJD>::new(51_544.5);
        let s = format!("{mjd}");
        assert!(s.contains("MJD"));
    }

    #[test]
    fn test_add_assign_sub_assign() {
        let mut jd = Time::<JD>::new(2_451_545.0);
        jd += Days::new(1.0);
        assert_eq!(jd.value(), 2_451_546.0);
        jd -= Days::new(0.5);
        assert_eq!(jd.value(), 2_451_545.5);
    }

    #[test]
    fn test_add_years() {
        let jd = Time::<JD>::new(2_450_000.0);
        let with_years = jd + Years::new(1.0);
        let span: Days = with_years - jd;
        assert!((span.value() - Time::<JD>::JULIAN_YEAR.value()).abs() < 1e-9);
    }

    #[test]
    fn test_div_days_and_f64() {
        let jd = Time::<JD>::new(100.0);
        assert!((jd / Days::new(2.0) - 50.0).abs() < 1e-12);
        assert!((jd / 4.0 - 25.0).abs() < 1e-12);
    }

    #[test]
    fn test_to_method_jd_mjd() {
        let jd = Time::<JD>::new(2_451_545.0);
        let mjd = jd.to::<MJD>();
        assert!((mjd.value() - 51_544.5).abs() < 1e-10);
    }
}
