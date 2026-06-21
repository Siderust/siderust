// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::starlight::ApparentMagnitude;

/// Convert an apparent magnitude to tenth-magnitude-star flux units.
#[must_use]
pub fn flux_10mag_units(magnitude: ApparentMagnitude) -> f64 {
    10.0_f64.powf(-0.4 * (magnitude.value() - 10.0))
}
