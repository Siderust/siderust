// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Piece-wise linear terrain-mask elevation profile.

use qtty::angular::Radians;
use qtty::Quantity;

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
            return Quantity::new(0.0);
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
                return Quantity::new(e0 + t * (e1 - e0));
            }
        }
        // azimuth before first sample → wrap.
        let (a0, e0) = self.samples[n - 1];
        let (a1, e1) = self.samples[0];
        let a0w = a0 - TAU;
        let t = if a1 == a0w {
            0.0
        } else {
            (az - a0w) / (a1 - a0w)
        };
        Quantity::new(e0 + t * (e1 - e0))
    }

    /// True iff the supplied elevation clears the mask at the given
    /// azimuth.
    pub fn clears(&self, azimuth: Radians, elevation: Radians) -> bool {
        elevation.value() >= self.elevation_at(azimuth).value()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use core::f64::consts::PI;

    #[test]
    fn terrain_mask_interpolates_linearly() {
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
}
