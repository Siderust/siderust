// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Search Options and Constants for Azimuth Calculus
//!
//! Re-exports [`SearchOpts`] from the altitude module (shared struct) and
//! provides scan-step constants for azimuth event finding.

use qtty::*;

/// Re-export the shared search options struct.
pub use crate::calculus::altitude::SearchOpts;

/// Default scan step for azimuth event detection: 10 minutes in days.
pub(crate) const DEFAULT_SCAN_STEP: Days = Minutes::new(10.0).to_const::<Day>();

/// Scan step used for azimuth extremum detection: 20 minutes in days.
pub(crate) const EXTREMA_SCAN_STEP: Days = Minutes::new(20.0).to_const::<Day>();
