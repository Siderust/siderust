// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Search Options and Constants
//!
//! ## Scientific scope
//!
//! Numerical controls for altitude crossing discovery and culmination
//! searches. The default refinement tolerance (~1 µs in time) is sized for
//! fast diurnal motion of solar-system bodies; very slow targets may require
//! a coarser scan step or finer tolerance.
//!
//! ## Technical scope
//!
//! Defines the stable public [`SearchOpts`] type and crate-internal crossing
//! search configuration consumed by [`super::events`] and [`super::provider`].
//!
//! ## References
//! None.

use crate::qtty::*;

// ---------------------------------------------------------------------------
// Stable public search options
// ---------------------------------------------------------------------------

/// Options for controlling search precision and strategy.
#[derive(Debug, Clone, Copy)]
pub struct SearchOpts {
    /// Time tolerance for root/extremum refinement (days).
    /// Default: ~1 µs (1e-9 days).
    pub time_tolerance: Days,
    /// Scan step for legacy scan+Brent compatibility (days).
    ///
    /// When set, crossing discovery uses the uniform scan plus Brent
    /// refinement path instead of the default Chebyshev-first engine.
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
// Internal crossing-search configuration
// ---------------------------------------------------------------------------

/// Internal crossing discovery mode for benchmarks and validation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CrossingAlgorithm {
    /// Chebyshev root finding first, with local scan+Brent fallback.
    ChebyshevFirst,
    /// Force the legacy uniform scan plus Brent refinement path.
    ScanBrent,
}

/// Internal options for Chebyshev-based crossing discovery.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ChebyshevOptions {
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

/// Internal extended search options used by event-search engines.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SearchOptsV2 {
    /// Time tolerance for root/extremum refinement (days).
    pub time_tolerance: Days,
    /// Scan step for scan+Brent fallback. `None` uses the target default.
    pub scan_step_days: Option<Days>,
    /// Crossing discovery mode.
    pub algorithm: CrossingAlgorithm,
    /// Chebyshev crossing-search controls.
    pub chebyshev: ChebyshevOptions,
}

impl SearchOptsV2 {
    /// Build internal options from the stable public option struct.
    #[inline]
    pub(crate) fn from_legacy(opts: SearchOpts) -> Self {
        Self {
            time_tolerance: opts.time_tolerance,
            scan_step_days: opts.scan_step_days,
            algorithm: if opts.scan_step_days.is_some() {
                CrossingAlgorithm::ScanBrent
            } else {
                CrossingAlgorithm::ChebyshevFirst
            },
            chebyshev: ChebyshevOptions::default(),
        }
    }

    /// Return the stable public subset of these options.
    #[inline]
    pub(crate) fn legacy(self) -> SearchOpts {
        SearchOpts {
            time_tolerance: self.time_tolerance,
            scan_step_days: self.scan_step_days,
        }
    }

    /// Benchmark/validation helper: force scan+Brent mode.
    #[inline]
    #[allow(dead_code)]
    pub(crate) fn scan_brent_legacy(opts: SearchOpts) -> Self {
        let mut out = Self::from_legacy(opts);
        out.algorithm = CrossingAlgorithm::ScanBrent;
        out
    }
}

impl Default for SearchOptsV2 {
    fn default() -> Self {
        let legacy = SearchOpts::default();
        Self {
            time_tolerance: legacy.time_tolerance,
            scan_step_days: legacy.scan_step_days,
            algorithm: CrossingAlgorithm::ChebyshevFirst,
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
