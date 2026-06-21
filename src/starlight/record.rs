// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::EquatorialMeanJ2000;
use crate::starlight::ApparentMagnitude;

/// Generic stellar catalogue record consumed by the map builder.
#[derive(Debug, Clone)]
pub struct StellarCatalogueRecord {
    /// Optional source identifier from the upstream catalogue.
    pub source_id: Option<String>,
    /// Source direction in the EquatorialMeanJ2000 frame.
    pub direction: Direction<EquatorialMeanJ2000>,
    /// Optional B-band apparent magnitude.
    pub b_mag: Option<ApparentMagnitude>,
    /// Optional V-band apparent magnitude.
    pub v_mag: Option<ApparentMagnitude>,
    /// Multiplicative catalogue weight; must be finite and non-negative.
    pub weight: f64,
}
