// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Time scales and instants
//!
//! ## Scientific scope
//!
//! Astronomical computations live on top of a small zoo of time scales —
//! TT (Terrestrial Time), TDB and TCB (the dynamical scales used by
//! barycentric ephemerides), TCG (geocentric coordinate time), TAI
//! (atomic time), UTC (civil time, with leap seconds), and UT1 (the
//! Earth-rotation scale). Mixing them silently produces errors of
//! seconds to tens of seconds — large enough to invalidate any
//! sub-arcsecond positional result. Following the IAU SOFA conventions
//! and the IERS Conventions (2010), every astronomical instant in this
//! crate carries its time scale at the type level so that conversions
//! between scales must be explicit.
//!
//! The reference epoch used throughout is **J2000.0 = JD 2 451 545.0
//! TT** (= 2000 January 1, 12:00:00 TT), as defined by IAU 2006
//! Resolution B1.
//!
//! ## Technical scope
//!
//! All time handling is delegated to the `tempoch` crate. This module
//! re-exports its full public surface (`Scale`, `EncodedTime`, `Time`,
//! the `TT`/`TDB`/`TCB`/`TCG`/`TAI`/`UTC`/`UT1` scale markers,
//! `Interval`, `delta_t_seconds`, …) and defines the v1 astronomy-facing
//! defaults as TT-based encoded dates plus the named same-scale helpers
//! (`to_jd`, `to_mjd`, `to_j2000_seconds`, `shifted_by`, `duration_since`)
//! provided by tempoch core:
//!
//! - [`JulianDate`] = `tempoch::JulianDate<TT>` = `EncodedTime<TT, JD>`.
//! - [`ModifiedJulianDate`] = `tempoch::ModifiedJulianDate<TT>` =
//!   `EncodedTime<TT, MJD>`.
//!
//! Advanced scale-aware work still uses `tempoch::JulianDate<S>` and
//! `tempoch::ModifiedJulianDate<S>` directly. The J2000.0 TT epoch is the
//! [`J2000`] constant (JD 2 451 545.0 TT).
//!
//! ## References
//!
//! - IAU SOFA Board (2021). *IAU SOFA Software Collection*. International
//!   Astronomical Union. <http://www.iausofa.org/>.
//! - Petit, G., Luzum, B. (Eds.) (2010). *IERS Conventions (2010)*. IERS
//!   Technical Note 36. Verlag des Bundesamts für Kartographie und
//!   Geodäsie, Frankfurt am Main.
//! - Seidelmann, P. K. (Ed.) (1992). *Explanatory Supplement to the
//!   Astronomical Almanac*. University Science Books. ISBN 0-935702-68-7.
//! - International Astronomical Union (2006). Resolution B3 ("Re-definition
//!   of Barycentric Dynamical Time, TDB"). IAU Transactions XXVIB.

pub use tempoch::{
    complement_within, constats, delta_t_seconds, delta_t_seconds_extrapolated, eop,
    ContinuousScale, ConversionError, ConversionTarget, CoordinateScale, EncodedTime,
    FormatForScale, InfallibleConversionTarget, InfallibleFormatForScale, Interval,
    InvalidIntervalError, PeriodListError, Scale, Time, TimeContext, TimeDataError, TimeInstant,
    JD, MJD, TAI, TCB, TCG, TDB, TT, UT1, UTC,
};

/// Julian year length in days, derived from [`qtty::time::JULIAN_YEAR`].
pub const JULIAN_YEAR_DAYS: qtty::Day = qtty::time::JULIAN_YEAR.to_const::<qtty::unit::Day>();

/// TT-scale Julian Date.  Equivalent to `tempoch::JulianDate<TT>`.
pub type JulianDate = tempoch::JulianDate<TT>;

/// TT-scale Modified Julian Date.  Equivalent to `tempoch::ModifiedJulianDate<TT>`.
pub type ModifiedJulianDate = tempoch::ModifiedJulianDate<TT>;

/// J2000.0 TT as a [`JulianDate`] (`JD 2 451 545.0 TT`); **`const`** everywhere.
pub const J2000: JulianDate = tempoch::JulianDate::<TT>::JD_EPOCH_J2000_0;

/// Runtime helper; prefer [`J2000`] in `const` contexts.
#[inline]
pub fn j2000_tt() -> JulianDate {
    J2000
}

/// Checked TT Julian Date constructor from a typed day quantity.
#[inline]
pub fn try_jd(raw: qtty::Day) -> Result<JulianDate, ConversionError> {
    if !raw.value().is_finite() {
        return Err(ConversionError::NonFinite);
    }
    JulianDate::try_new(raw)
}

/// Checked TT Julian Date constructor from a scalar JD value.
#[inline]
pub fn try_jd_f64(value: f64) -> Result<JulianDate, ConversionError> {
    try_jd(qtty::Day::new(value))
}

/// Checked TT Modified Julian Date constructor from a typed day quantity.
#[inline]
pub fn try_mjd(raw: qtty::Day) -> Result<ModifiedJulianDate, ConversionError> {
    if !raw.value().is_finite() {
        return Err(ConversionError::NonFinite);
    }
    ModifiedJulianDate::try_new(raw)
}

/// Checked TT Modified Julian Date constructor from a scalar MJD value.
#[inline]
pub fn try_mjd_f64(value: f64) -> Result<ModifiedJulianDate, ConversionError> {
    try_mjd(qtty::Day::new(value))
}

/// Helper: build a `ModifiedJulianDate` (TT scale) directly from a `chrono::DateTime<Utc>`.
///
/// This mirrors the common pattern `ModifiedJulianDate::from(dt)` but keeps a
/// crate-local helper for callers that prefer going through `siderust::time`.
#[inline]
pub fn modified_julian_date_from_chrono(dt: chrono::DateTime<chrono::Utc>) -> ModifiedJulianDate {
    ModifiedJulianDate::from(dt)
}

pub use crate::event::search::intervals::intersect as intersect_periods;
