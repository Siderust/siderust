// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Event Search
//!
//! ## Scientific scope
//! Internal typed time-domain search utilities for Siderust event APIs.
//! Not a general numerical library.
//!
//! ## Technical scope
//! Root-finding, extremum detection, interval assembly, and bracketing
//! policies operating on `Interval<ModifiedJulianDate>` time windows and
//! typed closures `Fn(ModifiedJulianDate) → Quantity<V>`.
//!
//! Public users should call high-level event APIs such as
//! [`crate::event::altitude`], [`crate::event::azimuth`], or the solar/lunar
//! event modules.  The helpers in this module are internal implementation
//! details.
//!
//! ## References
//! - Brent, R.P. (1973). *Algorithms for Minimization without Derivatives*. Prentice-Hall.

pub(crate) mod bracketing;
pub(crate) mod extrema;
pub(crate) mod intervals;
pub(crate) mod root_finding;
