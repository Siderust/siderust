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

/// Crossing discovery algorithm.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossingAlgorithm {
    /// Choose Chebyshev roots for suitable default Sun/Moon windows and
    /// scan+Brent where scan is cheaper or more conservative for the target.
    Auto,
    /// Use the legacy uniform scan plus Brent refinement path.
    ScanBrent,
    /// Try Chebyshev root finding first, falling back per unsafe segment.
    ChebyshevRoots,
}

/// Options for Chebyshev-based crossing discovery.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChebyshevOptions {
    /// Segment length fitted by one Chebyshev polynomial, in days.
    pub segment_length: Days,
    /// Polynomial degree used for the fitted signal.
    pub degree: usize,
    /// Maximum allowed absolute tail coefficient sum before a segment is
    /// split or sent to scan+Brent fallback.
    pub max_tail_norm: f64,
    /// Maximum allowed residual `|sin_altitude(t) - sin(threshold)|` after
    /// validation/refinement.
    pub max_residual: f64,
    /// Whether candidate roots should be refined against the precise signal.
    pub refine: bool,
    /// Initial precise-refinement and classification margin, in days.
    pub refine_margin: Days,
    /// Minimum acceptable absolute polynomial slope at a root.
    pub min_slope: f64,
    /// Whether unsafe polynomial-tail segments should be split before falling
    /// back to scan+Brent.
    pub adaptive_split: bool,
    /// Maximum adaptive split depth per original segment.
    pub max_split_depth: usize,
}

impl Default for ChebyshevOptions {
    fn default() -> Self {
        Self {
            segment_length: Hours::new(12.0).to::<Day>(),
            degree: 10,
            max_tail_norm: 1e-6,
            max_residual: 1e-10,
            refine: true,
            refine_margin: Minutes::new(20.0).to::<Day>(),
            min_slope: 1e-8,
            adaptive_split: true,
            max_split_depth: 2,
        }
    }
}

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

/// Extended search options for callers that need explicit crossing-algorithm
/// control without changing the legacy [`SearchOpts`] layout.
#[derive(Debug, Clone, Copy)]
pub struct SearchOptsV2 {
    /// Time tolerance for root/extremum refinement (days).
    pub time_tolerance: Days,
    /// Scan step for scan+Brent fallback. `None` uses the target default.
    pub scan_step_days: Option<Days>,
    /// Crossing discovery algorithm.
    pub algorithm: CrossingAlgorithm,
    /// Chebyshev crossing-search controls.
    pub chebyshev: ChebyshevOptions,
}

impl SearchOptsV2 {
    /// Build extended options from the legacy option struct.
    #[inline]
    pub fn from_legacy(opts: SearchOpts) -> Self {
        Self {
            time_tolerance: opts.time_tolerance,
            scan_step_days: opts.scan_step_days,
            ..Self::default()
        }
    }

    /// Return the legacy subset of these options.
    #[inline]
    pub fn legacy(self) -> SearchOpts {
        SearchOpts {
            time_tolerance: self.time_tolerance,
            scan_step_days: self.scan_step_days,
        }
    }
}

impl Default for SearchOptsV2 {
    fn default() -> Self {
        let legacy = SearchOpts::default();
        Self {
            time_tolerance: legacy.time_tolerance,
            scan_step_days: legacy.scan_step_days,
            algorithm: CrossingAlgorithm::Auto,
            chebyshev: ChebyshevOptions::default(),
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
