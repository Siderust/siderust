// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Search Options and Constants
//!
//! ## Scientific scope
//!
//! Numerical controls for the bracket‑and‑refine algorithm used to locate
//! altitude crossings and culminations. The default scan step (~10 min)
//! and refinement tolerance (~1 µs ≈ 86 µs in time) are sized for fast
//! diurnal motion of solar‑system bodies; very slow targets (Moon close
//! to standstill, deep‑sky tracking around culmination) may require a
//! coarser scan step or finer tolerance to avoid missed brackets or
//! spurious extrema.
//!
//! ## Technical scope
//!
//! Defines [`SearchOpts`], its [`Default`] impl, and the crate‑internal
//! step constants (`DEFAULT_SCAN_STEP`, `EXTREMA_SCAN_STEP`). No event
//! search is performed here; this module only carries configuration that
//! is consumed by [`super::events`] and [`super::provider`].
//!
//! ## References
//! None.

use crate::qtty::*;

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

/// Choose the best scan step given an optional provider hint and user overrides.
pub(crate) fn resolve_scan_step(hint: Option<Days>, opts: &SearchOpts, default_step: Days) -> Days {
    opts.scan_step_days.or(hint).unwrap_or(default_step)
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Default scan step: 10 minutes in days.
pub(crate) const DEFAULT_SCAN_STEP: Days = Minutes::new(10.0).to_const::<Day>();

/// Extrema scan step: 20 minutes in days (for culmination detection).
pub(crate) const EXTREMA_SCAN_STEP: Days = Minutes::new(20.0).to_const::<Day>();
