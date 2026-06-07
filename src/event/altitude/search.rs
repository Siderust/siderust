// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Search Options and Constants
//!
//! ## Scientific scope
//!
//! Numerical controls for altitude crossing discovery and culmination
//! searches. The default refinement tolerance (~1 µs in time) is sized for
//! fast diurnal motion of solar-system bodies.
//!
//! ## Technical scope
//!
//! Defines the stable public [`SearchOpts`] type and crate-internal crossing
//! search configuration consumed by [`super::events`].

use crate::qtty::*;

// ---------------------------------------------------------------------------
// Stable public search options (Option A: semantic tolerance only)
// ---------------------------------------------------------------------------

/// Options for controlling search precision.
#[derive(Debug, Clone, Copy)]
pub struct SearchOpts {
    /// Time tolerance for root/extremum refinement (days).
    /// Default: ~1 µs (1e-9 days).
    pub time_tolerance: Days,
}

impl Default for SearchOpts {
    fn default() -> Self {
        Self {
            time_tolerance: DEFAULT_TIME_TOLERANCE,
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

/// Internal scan+Brent fallback controls.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub(crate) struct FallbackSearchConfig {
    /// Scan step for scan+Brent fallback. `None` uses the target default hint.
    pub scan_step: Option<Days>,
    /// When true, the entire search window uses scan+Brent instead of the
    /// primary engine (tests/benches only).
    pub force_scan_baseline: bool,
}

/// Internal search configuration used by event-search engines.
#[derive(Debug, Clone, Copy)]
pub(crate) struct InternalSearchConfig {
    /// Time tolerance for root/extremum refinement (days).
    pub time_tolerance: Days,
    /// Chebyshev crossing-search controls.
    pub chebyshev: ChebyshevSearchConfig,
    /// Scan+Brent fallback controls.
    pub fallback: FallbackSearchConfig,
    /// Skip the solar daily predictor and use the generic Chebyshev engine
    /// (tests/benches only).
    pub disable_solar_daily_predictor: bool,
}

impl InternalSearchConfig {
    /// Build internal configuration from the stable public option struct.
    #[inline]
    pub(crate) fn from_public_opts(opts: SearchOpts) -> Self {
        Self {
            time_tolerance: opts.time_tolerance,
            chebyshev: ChebyshevSearchConfig::default(),
            fallback: FallbackSearchConfig::default(),
            disable_solar_daily_predictor: false,
        }
    }

    /// Return the stable public subset of these options.
    #[inline]
    pub(crate) fn public_opts(self) -> SearchOpts {
        SearchOpts {
            time_tolerance: self.time_tolerance,
        }
    }

    /// Whether the caller requested a uniform scan+Brent baseline path.
    #[inline]
    pub(crate) fn uses_scan_baseline(self) -> bool {
        self.fallback.force_scan_baseline
    }

    /// Resolve the scan step for fallback paths.
    #[inline]
    pub(crate) fn fallback_scan_step(self, default_step: Days) -> Days {
        self.fallback.scan_step.unwrap_or(default_step)
    }

    /// Test/bench helper: force scan+Brent baseline from public options.
    #[cfg(any(test, feature = "bench-internals"))]
    pub(crate) fn scan_brent_baseline_config(opts: SearchOpts) -> Self {
        Self {
            time_tolerance: opts.time_tolerance,
            chebyshev: ChebyshevSearchConfig::default(),
            fallback: FallbackSearchConfig {
                scan_step: Some(DEFAULT_SCAN_STEP),
                force_scan_baseline: true,
            },
            disable_solar_daily_predictor: true,
        }
    }

    /// Test/bench helper: generic Chebyshev-first baseline (no solar daily predictor).
    #[cfg(any(test, feature = "bench-internals"))]
    pub(crate) fn chebyshev_baseline_config(opts: SearchOpts) -> Self {
        Self {
            time_tolerance: opts.time_tolerance,
            chebyshev: ChebyshevSearchConfig::default(),
            fallback: FallbackSearchConfig::default(),
            disable_solar_daily_predictor: true,
        }
    }
}

impl Default for InternalSearchConfig {
    fn default() -> Self {
        Self::from_public_opts(SearchOpts::default())
    }
}

/// Choose the best scan step given an optional provider hint and internal config.
pub(crate) fn resolve_scan_step(
    hint: Option<Days>,
    opts: &InternalSearchConfig,
    default_step: Days,
) -> Days {
    opts.fallback.scan_step.or(hint).unwrap_or(default_step)
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// One mean solar day.
pub(crate) const ONE_DAY: Days = Days::new(1.0);

/// Default root/extremum refinement tolerance (~1 µs).
pub(crate) const DEFAULT_TIME_TOLERANCE: Days = Days::new(1e-9);

/// Minimum bracket width before treating endpoints as the root.
pub(crate) const ROOT_INTERVAL_EPS: Days = Days::new(1e-12);

/// Deduplicate crossing times closer than this interval.
pub(crate) const CROSSING_DEDUPE_EPS: Days = Days::new(1e-8);

/// Chebyshev segment scan hint for slowly varying diurnal bodies (2 hours).
pub(crate) const DIURNAL_CHEBY_SCAN_STEP: Days = Hours::new(2.0).to_const::<Day>();

/// Default scan step: 10 minutes in days.
pub(crate) const DEFAULT_SCAN_STEP: Days = Minutes::new(10.0).to_const::<Day>();

/// Extrema scan step: 20 minutes in days (for culmination detection).
pub(crate) const EXTREMA_SCAN_STEP: Days = Minutes::new(20.0).to_const::<Day>();
