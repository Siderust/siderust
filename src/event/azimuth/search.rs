// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Search Options and Constants for Azimuth Calculus
//!
//! ## Scientific scope
//!
//! Numerical controls for the bracket‑and‑refine azimuth event search.
//! The default scan step (~10 min) is appropriate for diurnal motion of
//! solar‑system bodies and most stars at mid‑latitudes; bodies near
//! standstill or stationary points may need a finer step to avoid missed
//! brackets near a turning point.
//!
//! ## Technical scope
//!
//! Re‑exports [`SearchOpts`] from the altitude module (shared struct) and
//! provides crate‑internal scan‑step constants
//! ([`DEFAULT_SCAN_STEP`], [`EXTREMA_SCAN_STEP`]) consumed by
//! [`super::events`].
//!
//! ## References
//! None.

use crate::qtty::*;

/// Re-export the shared search options struct.
pub use crate::event::altitude::SearchOpts;

/// Default scan step for azimuth event detection: 10 minutes in days.
pub(crate) const DEFAULT_SCAN_STEP: Days = Minutes::new(10.0).to_const::<Day>();

/// Scan step used for azimuth extremum detection: 20 minutes in days.
pub(crate) const EXTREMA_SCAN_STEP: Days = Minutes::new(20.0).to_const::<Day>();
