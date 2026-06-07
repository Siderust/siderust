// SPDX-License-Identifier: AGPL-3.0-only
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
//! search configuration consumed by [`super::events`].
//!
//! ## References
//! None.

use crate::qtty::*;

// ---------------------------------------------------------------------------
// Stable public search options
// ---------------------------------------------------------------------------

/// Options for controlling search precision.
#[derive(Debug, Clone, Copy)]
pub struct SearchOpts {
    /// Time tolerance for root/extremum refinement (days).
    /// Default: ~1 µs (1e-9 days).
    pub time_tolerance: Days,
    /// Optional scan-step override for validation or baseline searches.
    ///
    /// Normal users should leave this as `None`; the default engine is
    /// Chebyshev-first with internal scan+Brent fallback.
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

/// Internal options for Chebyshev-based crossing discovery.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ChebyshevSearchConfig {
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

impl Default for ChebyshevSearchConfig {
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

/// Internal search configuration used by event-search engines.
#[derive(Debug, Clone, Copy)]
pub(crate) struct InternalSearchConfig {
    /// Time tolerance for root/extremum refinement (days).
    pub time_tolerance: Days,
    /// Scan step for scan+Brent baseline. `None` uses the target default.
    pub scan_step_days: Option<Days>,
    /// Chebyshev crossing-search controls.
    pub chebyshev: ChebyshevSearchConfig,
}

impl InternalSearchConfig {
    /// Build internal configuration from the stable public option struct.
    #[inline]
    pub(crate) fn from_public_opts(opts: SearchOpts) -> Self {
        Self {
            time_tolerance: opts.time_tolerance,
            scan_step_days: opts.scan_step_days,
            chebyshev: ChebyshevSearchConfig::default(),
        }
    }

    /// Return the stable public subset of these options.
    #[inline]
    pub(crate) fn public_opts(self) -> SearchOpts {
        SearchOpts {
            time_tolerance: self.time_tolerance,
            scan_step_days: self.scan_step_days,
        }
    }

    /// Whether the caller requested a uniform scan+Brent baseline path.
    #[inline]
    pub(crate) fn uses_scan_baseline(self) -> bool {
        self.scan_step_days.is_some()
    }

    /// Test/bench helper: force scan+Brent baseline from public options.
    #[allow(dead_code)]
    pub(crate) fn scan_brent_baseline(opts: SearchOpts) -> Self {
        let mut config = Self::from_public_opts(opts);
        if config.scan_step_days.is_none() {
            config.scan_step_days = Some(DEFAULT_SCAN_STEP);
        }
        config
    }
}

impl Default for InternalSearchConfig {
    fn default() -> Self {
        Self::from_public_opts(SearchOpts::default())
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
