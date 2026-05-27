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
    /// Test whether a unit vector expressed in the instrument frame
    /// (boresight = +Z) lies inside the field of view.
    ///
    /// The input is the cosine of the angle off boresight and (for
    /// rectangular FoVs) the off-axis offsets `(x_angle, y_angle)` in
    /// radians.
    pub fn contains_boresight(&self, off_axis: Radians) -> bool {
        match self {
            Fov::Conical { half_angle } => off_axis.value().abs() <= half_angle.value(),
            Fov::Rectangular {
                across_track,
                along_track,
            } => {
                off_axis.value().abs() <= across_track.value() / 2.0
                    || off_axis.value().abs() <= along_track.value() / 2.0
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
