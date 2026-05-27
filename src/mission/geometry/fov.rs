// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Field-of-view geometry and instrument metadata.

use qtty::angular::Radians;

/// Field-of-view geometry attached to an instrument.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Fov {
    /// Right circular cone, parameterised by half-angle from the
    /// boresight to the cone edge.
    Conical {
        /// Half-angle of the cone.
        half_angle: Radians,
    },
    /// Rectangular pyramid with full-angle extents about the boresight.
    Rectangular {
        /// Full angle across the X axis.
        across_track: Radians,
        /// Full angle across the Y axis.
        along_track: Radians,
    },
}

impl Fov {
    /// Test whether an off-axis angle (measured from the boresight in the
    /// instrument frame, boresight = +Z) lies inside the field of view.
    ///
    /// For conical FoVs, the point is inside when `|off_axis| ≤ half_angle`.
    ///
    /// For rectangular FoVs a single off-axis angle is insufficient to fully
    /// describe containment in both axes; this method conservatively tests
    /// `|off_axis| ≤ min(across_track, along_track) / 2`. Callers that need
    /// a proper 2-axis rectangular test should decompose the off-axis vector
    /// into its x and y components and test each axis independently.
    pub fn contains_boresight(&self, off_axis: Radians) -> bool {
        let abs_off = off_axis.abs();
        match self {
            Fov::Conical { half_angle } => abs_off <= *half_angle,
            Fov::Rectangular {
                across_track,
                along_track,
            } => {
                let half_x = *across_track / 2.0;
                let half_y = *along_track / 2.0;
                // Use the smaller half-angle as the conservative single-axis bound.
                let bound = if half_x <= half_y { half_x } else { half_y };
                abs_off <= bound
            }
        }
    }
}

/// Instrument metadata: a typed FoV plus a free-text identifier.
#[derive(Debug, Clone, PartialEq)]
pub struct Instrument {
    /// Free-text identifier for the instrument.
    pub id: String,
    /// Field-of-view geometry.
    pub fov: Fov,
}

impl Instrument {
    /// Build an [`Instrument`] with the supplied identifier and FoV.
    pub fn new(id: impl Into<String>, fov: Fov) -> Self {
        Self { id: id.into(), fov }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::Quantity;

    #[test]
    fn fov_conical_contains_axis() {
        let fov = Fov::Conical {
            half_angle: Quantity::new(0.1),
        };
        assert!(fov.contains_boresight(Quantity::new(0.0)));
        assert!(fov.contains_boresight(Quantity::new(0.099)));
        assert!(!fov.contains_boresight(Quantity::new(0.2)));
    }

    #[test]
    fn fov_rectangular_uses_smaller_half_angle() {
        // across_track = 0.2 rad full → half = 0.1; along_track = 0.4 rad full → half = 0.2.
        // Conservative bound = min(0.1, 0.2) = 0.1.
        let fov = Fov::Rectangular {
            across_track: Quantity::new(0.2),
            along_track: Quantity::new(0.4),
        };
        assert!(fov.contains_boresight(Quantity::new(0.0)));
        assert!(fov.contains_boresight(Quantity::new(0.09)));
        assert!(!fov.contains_boresight(Quantity::new(0.15)));
    }

    #[test]
    fn instrument_construction_preserves_id() {
        let inst = Instrument::new(
            "demo",
            Fov::Conical {
                half_angle: Quantity::new(0.1),
            },
        );
        assert_eq!(inst.id, "demo");
    }
}
