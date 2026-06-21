// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Photometric conversions used by the v1 starlight model.
//!
//! S10 denotes the flux of a tenth-magnitude star, following the traditional
//! diffuse night-sky-brightness convention used by Leinert et al. (1998).

use crate::starlight::ApparentMagnitude;

/// Convert an apparent magnitude to tenth-magnitude-star flux units.
///
/// A magnitude of 10 maps to one S10 flux unit. Brighter sources have larger
/// values according to the astronomical magnitude relation.
#[must_use]
pub fn flux_10mag_units(magnitude: ApparentMagnitude) -> f64 {
    10.0_f64.powf(-0.4 * (magnitude.value() - 10.0))
}
