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

// ---------------------------------------------------------------------------
// Unstable exports for the separate FFI crate (gated; not semver-stable).
// ---------------------------------------------------------------------------

#[cfg(feature = "unstable-event-search")]
/// Crossing algorithm selector for experimental callers.
pub enum CrossingAlgorithmFfi {
    /// Chebyshev-first with per-segment scan+Brent fallback.
    Auto,
    /// Force uniform scan plus Brent refinement.
    ScanBrent,
    /// Chebyshev polynomial roots with per-segment fallback.
    ChebyshevRoots,
}

#[cfg(feature = "unstable-event-search")]
#[allow(private_interfaces)]
/// Chebyshev tuning retained for experimental callers.
pub type ChebyshevOptionsFfi = ChebyshevOptions;

#[cfg(feature = "unstable-event-search")]
#[allow(private_interfaces)]
/// Extended search options retained for experimental callers.
pub struct SearchOptsV2Ffi {
    /// Time tolerance for refinement (days).
    pub time_tolerance: Days,
    /// Optional scan step override (days).
    pub scan_step_days: Option<Days>,
    /// Crossing discovery mode.
    pub algorithm: CrossingAlgorithmFfi,
    /// Chebyshev segment controls.
    pub chebyshev: ChebyshevOptions,
}

#[cfg(feature = "unstable-event-search")]
#[allow(private_interfaces)]
impl SearchOptsV2Ffi {
    /// Build experimental options from the stable public struct.
    pub fn from_legacy(opts: SearchOpts) -> Self {
        Self {
            time_tolerance: opts.time_tolerance,
            scan_step_days: opts.scan_step_days,
            algorithm: if opts.scan_step_days.is_some() {
                CrossingAlgorithmFfi::ScanBrent
            } else {
                CrossingAlgorithmFfi::Auto
            },
            chebyshev: ChebyshevOptions::default(),
        }
    }

    /// Convert to the internal engine configuration.
    pub fn to_internal(self) -> SearchOptsV2 {
        let mut opts = SearchOptsV2 {
            time_tolerance: self.time_tolerance,
            scan_step_days: self.scan_step_days,
            chebyshev: self.chebyshev,
            ..SearchOptsV2::default()
        };
        opts.algorithm = match self.algorithm {
            CrossingAlgorithmFfi::ScanBrent => CrossingAlgorithm::ScanBrent,
            CrossingAlgorithmFfi::Auto | CrossingAlgorithmFfi::ChebyshevRoots => {
                CrossingAlgorithm::ChebyshevFirst
            }
        };
        if opts.scan_step_days.is_some() {
            opts.algorithm = CrossingAlgorithm::ScanBrent;
        }
        opts
    }
}

#[cfg(feature = "unstable-event-search")]
impl Default for SearchOptsV2Ffi {
    fn default() -> Self {
        Self {
            time_tolerance: SearchOpts::default().time_tolerance,
            scan_step_days: None,
            algorithm: CrossingAlgorithmFfi::Auto,
            chebyshev: ChebyshevOptions::default(),
        }
    }
}

#[cfg(feature = "unstable-event-search")]
#[allow(private_interfaces, dead_code)]
pub(crate) fn compose_search_opts_v2_ffi(
    time_tolerance: Days,
    scan_step_days: Option<Days>,
    algorithm: CrossingAlgorithmFfi,
    chebyshev: ChebyshevOptionsFfi,
) -> SearchOptsV2Ffi {
    SearchOptsV2Ffi {
        time_tolerance,
        scan_step_days,
        algorithm,
        chebyshev,
    }
}

#[cfg(feature = "unstable-event-search")]
#[allow(private_interfaces, clippy::too_many_arguments)]
/// Build experimental search options from scalar fields (for the FFI crate).
pub fn search_opts_v2_ffi_from_scalars(
    time_tolerance_days: f64,
    has_scan_step: bool,
    scan_step_days: f64,
    algorithm: CrossingAlgorithmFfi,
    has_chebyshev: bool,
    segment_length_days: f64,
    degree: u32,
    max_tail_norm: f64,
    max_residual: f64,
    refine: bool,
    refine_margin_days: f64,
    min_slope: f64,
    adaptive_split: bool,
    max_split_depth: u32,
) -> SearchOptsV2Ffi {
    let defaults = SearchOptsV2Ffi::default();
    let time_tolerance = if time_tolerance_days > 0.0 {
        Days::new(time_tolerance_days)
    } else {
        defaults.time_tolerance
    };
    let scan_step = if has_scan_step && scan_step_days > 0.0 {
        Some(Days::new(scan_step_days))
    } else {
        None
    };
    let chebyshev = if has_chebyshev {
        compose_chebyshev_options_ffi(
            if segment_length_days > 0.0 {
                Days::new(segment_length_days)
            } else {
                defaults.chebyshev.segment_length
            },
            if degree > 0 {
                degree as usize
            } else {
                defaults.chebyshev.degree
            },
            if max_tail_norm > 0.0 {
                max_tail_norm
            } else {
                defaults.chebyshev.max_tail_norm
            },
            if max_residual > 0.0 {
                max_residual
            } else {
                defaults.chebyshev.max_residual
            },
            refine,
            if refine_margin_days > 0.0 {
                Days::new(refine_margin_days)
            } else {
                defaults.chebyshev.refine_margin
            },
            if min_slope > 0.0 {
                min_slope
            } else {
                defaults.chebyshev.min_slope
            },
            adaptive_split,
            if max_split_depth > 0 {
                max_split_depth as usize
            } else {
                defaults.chebyshev.max_split_depth
            },
        )
    } else {
        defaults.chebyshev
    };
    let mut opts = SearchOptsV2Ffi {
        time_tolerance,
        scan_step_days: scan_step,
        algorithm,
        chebyshev,
    };
    if opts.scan_step_days.is_some() {
        opts.algorithm = CrossingAlgorithmFfi::ScanBrent;
    }
    opts
}

#[cfg(feature = "unstable-event-search")]
#[allow(private_interfaces, clippy::too_many_arguments, dead_code)]
pub(crate) fn compose_chebyshev_options_ffi(
    segment_length: Days,
    degree: usize,
    max_tail_norm: f64,
    max_residual: f64,
    refine: bool,
    refine_margin: Days,
    min_slope: f64,
    adaptive_split: bool,
    max_split_depth: usize,
) -> ChebyshevOptionsFfi {
    ChebyshevOptions {
        segment_length,
        degree,
        max_tail_norm,
        max_residual,
        refine,
        refine_margin,
        min_slope,
        adaptive_split,
        max_split_depth,
    }
}
