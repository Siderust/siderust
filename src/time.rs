// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Time types for siderust, built directly on `tempoch`.
//!
//! All time handling is delegated to `tempoch`. This module re-exports the
//! full public surface of `tempoch` and adds two backward-compatible type
//! aliases so that unparameterised call sites continue to resolve to the
//! TT scale:
//!
//! - [`JulianDate`] = `tempoch::JulianDate<TT>` = `EncodedTime<TT, JD>`
//! - [`ModifiedJulianDate`] = `tempoch::ModifiedJulianDate<TT>` = `EncodedTime<TT, MJD>`
//!
//! A [`J2000`] constant is also re-exported for the J2000.0 TT epoch.

pub use tempoch::{
    constats, delta_t_seconds, delta_t_seconds_extrapolated, eop, ContinuousScale, ConversionError,
    ConversionTarget, CoordinateScale, EncodedTime, InfallibleConversionTarget,
    InfallibleRepresentationForScale, Interval, InvalidIntervalError, InvalidPeriodError, J2000_TT,
    JD, JulianTimeExt, MJD, PeriodListError, RepresentationForScale, Scale, ScaleKind, Time,
    TimeContext, TimeDataError, TimeInstant, TAI, TCB, TCG, TDB, TT, UT1, UTC,
};

/// Re-export `JulianDate` and `ModifiedJulianDate` under their generic names
/// so code that explicitly writes `JulianDate<TT>` or `JulianDate<TDB>` still compiles.
pub use tempoch::{
    JulianDate as JulianDateG, ModifiedJulianDate as ModifiedJulianDateG,
};

/// Backward-compatible alias: old siderust code used `UT` for the UT1 axis.
pub type UT = UT1;

/// Backward-compatible generic period alias over an instant type `T`.
pub type Period<T> = Interval<T>;

/// TT-scale Julian Date.  Equivalent to `tempoch::JulianDate<TT>`.
pub type JulianDate = tempoch::JulianDate<TT>;

/// TT-scale Modified Julian Date.  Equivalent to `tempoch::ModifiedJulianDate<TT>`.
pub type ModifiedJulianDate = tempoch::ModifiedJulianDate<TT>;

/// J2000.0 epoch as a [`JulianDate`] (TT, JD 2 451 545.0).
pub const J2000: JulianDate = J2000_TT;
