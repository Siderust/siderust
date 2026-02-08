// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Search Options and Constants
//!
//! Configuration for altitude search algorithms including tolerances,
//! scan steps, and target-specific defaults.

use qtty::*;

// ---------------------------------------------------------------------------
// Search Options
// ---------------------------------------------------------------------------

/// Options for controlling search precision and strategy.
#[derive(Debug, Clone, Copy)]
pub struct SearchOpts {
    /// Time tolerance for root/extremum refinement (days).
    /// Default: ~1 µs (1e-9 days).
    pub time_tolerance: Days,
    /// Scan step for coarse bracket detection (days).
    /// Default: 10 minutes. Override for slower-moving bodies like the Moon.
    pub scan_step_days: Option<Days>,
}

impl Default for SearchOpts {
    fn default() -> Self {
        Self {
            time_tolerance: Days::new(1e-9),
            scan_step_days: None,
        }
    }
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Default scan step: 10 minutes in days.
pub(crate) const DEFAULT_SCAN_STEP: Days = Minutes::new(10.0).to::<Day>();

/// Extrema scan step: 20 minutes in days (for culmination detection).
pub(crate) const EXTREMA_SCAN_STEP: Days = Minutes::new(20.0).to::<Day>();
