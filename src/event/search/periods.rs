// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Generic threshold-period assembly.
//!
//! These helpers are body-agnostic. Callers own crossing discovery; this module
//! only converts labelled crossings into above, below, or in-band intervals
//! within the same MJD/TT query window.

use crate::time::{complement_within, Interval, ModifiedJulianDate};

use super::intervals::{self, LabeledCrossing};

/// Assemble periods where a signal is above a threshold.
pub(crate) fn assemble_above_threshold_periods(
    labelled_crossings: &[LabeledCrossing],
    window: Interval<ModifiedJulianDate>,
    start_above: bool,
) -> Vec<Interval<ModifiedJulianDate>> {
    intervals::build_above_periods_directed(labelled_crossings, window, start_above)
}

/// Assemble periods where a signal is below a threshold.
///
/// This is intentionally defined as the complement of the above-threshold
/// periods, making the invariant auditable at every call site.
pub(crate) fn complement_threshold_periods(
    window: Interval<ModifiedJulianDate>,
    above_threshold_periods: &[Interval<ModifiedJulianDate>],
) -> Vec<Interval<ModifiedJulianDate>> {
    complement_within(window, above_threshold_periods)
}

/// Assemble periods where `min_threshold <= signal <= max_threshold`.
pub(crate) fn assemble_in_range_periods(
    min_crossings: &[LabeledCrossing],
    start_above_min: bool,
    max_crossings: &[LabeledCrossing],
    start_above_max: bool,
    window: Interval<ModifiedJulianDate>,
) -> Vec<Interval<ModifiedJulianDate>> {
    let above_min = assemble_above_threshold_periods(min_crossings, window, start_above_min);
    let above_max = assemble_above_threshold_periods(max_crossings, window, start_above_max);
    let below_max = complement_threshold_periods(window, &above_max);
    intervals::intersect(&above_min, &below_max)
}
