// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # DE440 Ephemeris Module
//!
//! This module provides a runtime evaluator for JPL DE440 ephemeris data,
//! which is embedded at compile time (Chebyshev polynomial coefficients
//! for Sun, Earth-Moon barycenter, and Moon).
//!
//! ## Architecture
//!
//! - **Data**: Coefficient tables are extracted from `de440.bsp` by the build
//!   script and embedded via `include_bytes!()`.
//! - **Evaluator**: [`chebyshev`] evaluates Chebyshev polynomials using the
//!   Clenshaw recurrence. [`eval`] handles segment lookup and time normalization.
//! - **Bodies**: [`bodies`] resolves composite body states (e.g., Earth from
//!   EMB and Moon, Earth heliocentric from barycentric).
//!
//! ## Reference Frame & Units
//!
//! DE440 data is in ICRF (≈ J2000 equatorial), units of km and km/s,
//! on the TDB timescale (seconds past J2000 epoch).
//!
//! The public API converts to ecliptic coordinates and AU to match the
//! `Ephemeris` trait signature used by the rest of siderust.

pub mod bodies;
pub mod chebyshev;
pub mod data;
pub mod eval;
