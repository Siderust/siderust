// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Instruments and ground-station locations
//!
//! ## Scientific scope
//!
//! Typed primitives for instrument-attached field-of-view geometry and
//! observatory / ground-station location metadata. The library exposes
//! the pure geometric and metadata layer; mission-planning workflows,
//! scheduling, and database persistence belong in **SatOps**.
//!
//! ## Technical scope
//!
//! - [`Instrument`] couples a free-text identifier with a typed
//!   [`Fov`] field of view.
//! - [`Fov`] is a flat enum covering [`Fov::Conical`] (half-angle)
//!   and [`Fov::Rectangular`] (full angles in two orthogonal axes).
//! - [`Location`] aggregates a typed `affn::ellipsoidal::Position`,
//!   a free-text site identifier, and an optional [`TerrainMask`].
//! - [`TerrainMask`] is a piece-wise linear elevation cut-off keyed
//!   by azimuth in `[0, 2π)`. All inputs are typed radians.
//!
//! ## References
//!
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and
//!   Applications*, 4th ed. §3.4.
//! - Wertz, J. R. (2011). *Mission Geometry; Orbit and Constellation
//!   Design and Management*. §10.

#![forbid(unsafe_code)]

use qtty::angular::Radians;

use crate::mission_geometry::LocalFrame;

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

/// Piece-wise linear terrain-mask elevation profile keyed by azimuth.
///
/// The supplied samples are sorted by azimuth on construction; mask
/// elevation at an arbitrary azimuth is linearly interpolated between
/// adjacent samples and the profile wraps around at `2π`.
#[derive(Debug, Clone, PartialEq)]
pub struct TerrainMask {
    samples: Vec<(f64, f64)>,
}

impl TerrainMask {
    /// Build a terrain mask from `(azimuth, elevation)` samples. Empty
    /// inputs produce a flat zero mask.
    pub fn new(mut samples: Vec<(Radians, Radians)>) -> Self {
        samples.sort_by(|a, b| a.0.value().partial_cmp(&b.0.value()).unwrap());
        let s = samples
            .into_iter()
            .map(|(a, e)| (a.value().rem_euclid(core::f64::consts::TAU), e.value()))
            .collect();
        Self { samples: s }
    }

    /// Elevation of the mask at the supplied azimuth.
    pub fn elevation_at(&self, azimuth: Radians) -> Radians {
        use core::f64::consts::TAU;
        if self.samples.is_empty() {
            return qtty::Quantity::new(0.0);
        }
        let az = azimuth.value().rem_euclid(TAU);
        // Linear interpolation with wrap-around.
        let n = self.samples.len();
        for i in 0..n {
            let (a0, e0) = self.samples[i];
            let (a1, e1) = if i + 1 < n {
                self.samples[i + 1]
            } else {
                let (a, e) = self.samples[0];
                (a + TAU, e)
            };
            if az >= a0 && az <= a1 {
                let t = if a1 == a0 { 0.0 } else { (az - a0) / (a1 - a0) };
                return qtty::Quantity::new(e0 + t * (e1 - e0));
            }
        }
        // azimuth before first sample → wrap.
        let (a0, e0) = self.samples[n - 1];
        let (a1, e1) = self.samples[0];
        let a0w = a0 - TAU;
        let t = if a1 == a0w { 0.0 } else { (az - a0w) / (a1 - a0w) };
        qtty::Quantity::new(e0 + t * (e1 - e0))
    }

    /// True iff the supplied elevation clears the mask at the given
    /// azimuth.
    pub fn clears(&self, azimuth: Radians, elevation: Radians) -> bool {
        elevation.value() >= self.elevation_at(azimuth).value()
    }
}

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
    use approx::assert_abs_diff_eq;
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

    #[test]
    fn terrain_mask_interpolates_linearly() {
        use core::f64::consts::PI;
        let mask = TerrainMask::new(vec![
            (Quantity::new(0.0), Quantity::new(0.1)),
            (Quantity::new(PI), Quantity::new(0.2)),
        ]);
        // Halfway between 0 and π → midpoint elevation.
        let e = mask.elevation_at(Quantity::new(PI / 2.0));
        assert_abs_diff_eq!(e.value(), 0.15, epsilon = 1e-12);
    }

    #[test]
    fn empty_mask_is_zero() {
        let mask = TerrainMask::new(vec![]);
        assert_eq!(mask.elevation_at(Quantity::new(1.0)).value(), 0.0);
        assert!(mask.clears(Quantity::new(1.0), Quantity::new(0.0)));
    }

    #[test]
    fn location_with_mask_round_trip() {
        let frame = LocalFrame::new(Quantity::new(0.5), Quantity::new(0.3));
        let mask = TerrainMask::new(vec![(Quantity::new(0.0), Quantity::new(0.1))]);
        let loc = Location::new("VLBI-A", frame).with_mask(mask);
        assert_eq!(loc.id, "VLBI-A");
        assert!(loc.mask.is_some());
    }
}
