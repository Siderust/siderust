// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Observation site and ground-station metadata
//!
//! ## Scientific scope
//!
//! An observation site aggregates a geodetic position (encoded as a
//! [`LocalFrame`] carrying typed geodetic latitude and longitude) with a
//! free-text site identifier and an optional terrain-mask elevation
//! cut-off. These are the primitives needed to transform topocentric
//! azimuth/elevation/range solutions into visibility windows and to apply
//! site-specific horizon masks.
//!
//! ## Technical scope
//!
//! - [`Location`] — observation-site / ground-station record.
//! - Build with [`Location::new`] and optionally attach a mask with
//!   [`Location::with_mask`].
//!
//! ## References
//!
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and
//!   Applications*, 4th ed. §3.4.
//! - Wertz, J. R. (2011). *Mission Geometry; Orbit and Constellation
//!   Design and Management*. §10.

use crate::mission::geometry::{LocalFrame, TerrainMask};

/// Observation site or ground-station location.
#[derive(Debug, Clone, PartialEq)]
pub struct Location {
    /// Free-text site identifier.
    pub id: String,
    /// Local east-north-up frame at the site (geodetic lat/lon).
    pub frame: LocalFrame,
    /// Optional terrain-mask elevation cut-off.
    pub mask: Option<TerrainMask>,
}

impl Location {
    /// Build a [`Location`] from a free-text identifier and a
    /// [`LocalFrame`].
    pub fn new(id: impl Into<String>, frame: LocalFrame) -> Self {
        Self {
            id: id.into(),
            frame,
            mask: None,
        }
    }

    /// Attach a terrain mask to this location.
    pub fn with_mask(mut self, mask: TerrainMask) -> Self {
        self.mask = Some(mask);
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::Quantity;

    #[test]
    fn location_with_mask_round_trip() {
        let frame = LocalFrame::new(Quantity::new(0.5), Quantity::new(0.3));
        let mask = TerrainMask::new(vec![(Quantity::new(0.0), Quantity::new(0.1))]);
        let loc = Location::new("VLBI-A", frame).with_mask(mask);
        assert_eq!(loc.id, "VLBI-A");
        assert!(loc.mask.is_some());
    }
}
