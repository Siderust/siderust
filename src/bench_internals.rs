// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Internal baseline helpers for Criterion benchmarks.
//!
//! Enabled with the `bench-internals` feature; not part of the stable API.
//! Run gated benches with:
//!
//! ```bash
//! cargo bench --features bench-internals --bench solar_altitude
//! cargo bench --features bench-internals --bench moon_altitude
//! ```

use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::event::altitude::search::InternalSearchConfig;
use crate::event::altitude::{CrossingEvent, SearchOpts};
use crate::qtty::*;
use crate::time::{Interval, ModifiedJulianDate};

/// Run a solar below-threshold search using the internal scan+Brent baseline.
pub fn solar_below_threshold_scan_baseline(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<Interval<ModifiedJulianDate>> {
    crate::event::solar::solar_below_threshold_impl(
        site,
        window,
        threshold,
        InternalSearchConfig::scan_brent_baseline_config(opts),
    )
}

/// Run a solar below-threshold search using the generic Chebyshev engine only.
pub fn solar_below_threshold_chebyshev_baseline(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<Interval<ModifiedJulianDate>> {
    crate::event::solar::solar_below_threshold_impl(
        site,
        window,
        threshold,
        InternalSearchConfig::chebyshev_baseline_config(opts),
    )
}

/// Run a solar altitude-range search using the internal scan+Brent baseline.
pub fn solar_altitude_ranges_scan_baseline(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    h_min: Degrees,
    h_max: Degrees,
    opts: SearchOpts,
) -> Vec<Interval<ModifiedJulianDate>> {
    crate::event::solar::solar_altitude_ranges_impl(
        site,
        window,
        (h_min, h_max),
        InternalSearchConfig::scan_brent_baseline_config(opts),
    )
}

/// Run a solar crossing search using the internal scan+Brent baseline.
pub fn solar_crossings_scan_baseline(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<CrossingEvent> {
    crate::event::solar::solar_crossings_impl(
        site,
        window,
        threshold,
        InternalSearchConfig::scan_brent_baseline_config(opts),
    )
}

/// Run a lunar above-threshold search using the internal scan+Brent baseline.
pub fn lunar_above_threshold_scan_baseline(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<Interval<ModifiedJulianDate>> {
    crate::event::lunar::lunar_above_threshold_impl(
        site,
        window,
        threshold,
        InternalSearchConfig::scan_brent_baseline_config(opts),
    )
}

/// Run a lunar below-threshold search using the internal scan+Brent baseline.
pub fn lunar_below_threshold_scan_baseline(
    site: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<Interval<ModifiedJulianDate>> {
    crate::event::lunar::lunar_below_threshold_impl(
        site,
        window,
        threshold,
        InternalSearchConfig::scan_brent_baseline_config(opts),
    )
}
