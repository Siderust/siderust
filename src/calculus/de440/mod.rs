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
//! - **[`data`]** — Coefficient tables extracted from `de440.bsp` by the build
//!   script and embedded via `include_bytes!()`.  Pre-built
//!   [`SegmentDescriptor`](eval::SegmentDescriptor)s (`SUN`, `EMB`, `MOON`)
//!   bundle per-body metadata for ergonomic evaluation.
//! - **[`chebyshev`]** — Clenshaw-recurrence Chebyshev polynomial evaluator
//!   (value, derivative, or both in one pass).
//! - **[`eval`]** — [`SegmentDescriptor`](eval::SegmentDescriptor): segment
//!   lookup, time-normalisation, and 3-D position / velocity evaluation.
//!   Physical constants (`SECONDS_PER_DAY`, `J2000_JD`) are sourced from
//!   `qtty` and `JulianDate`, not duplicated.
//! - **[`bodies`]** — Body-chain resolution (EMB → Earth, Moon geocentric,
//!   Earth heliocentric).  Time-scale conversion delegates to
//!   `JulianDate::tt_to_tdb()`. Unit conversion (`km → AU`) uses `qtty`
//!   unit ratios as the single source of truth.
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
