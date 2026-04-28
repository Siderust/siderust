// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed alt/az hemispherical sky grid sampler.
//!
//! Provides [`SkyGrid`], a configurable iterator over sky directions covering
//! (part of) the visible hemisphere in regular altitude–azimuth cells. The
//! grid yields frame-typed [`Direction<Horizontal>`](crate::coordinates::spherical::direction::Horizontal)
//! values; per-cell solid angles are exposed as typed [`Steradians`] through
//! [`SkyGrid::with_solid_angle`] / [`SkyGrid::iter_cells`].
//!
//! ## Module placement
//!
//! Hemispherical samplers are geometric primitives independent of any
//! particular astronomical model — they don't need ephemerides, sidereal
//! time, or even an observer. They live in [`crate::geometry`] (alongside
//! other sky-coverage helpers) rather than under `astro::` or
//! `calculus::horizontal`, both of which are reserved for time-dependent or
//! body-specific calculations.
//!
//! ## Cell layout
//!
//! - **Cell centres** are sampled (not edges):
//!   `alt = alt_min + (i + 0.5)·Δalt`, `az = (j + 0.5)·Δaz`.
//! - `n_alt = round((alt_max − alt_min) / Δalt)`,
//!   `n_az = round(360° / Δaz)` for the uniform grid.
//! - Iteration order is **altitude-outer ascending, azimuth-inner ascending**.
//! - The default altitude range is `[0°, 90°)` (upper hemisphere).
//!
//! ## Equal-area mode
//!
//! When constructed with [`SkyGrid::equal_area`], the azimuth count per
//! altitude band scales with `cos(alt)`:
//!
//! ```text
//! n_az(alt) = max(1, round(2π · cos(alt) / Δaz_horizon))
//! ```
//!
//! so cells subtend approximately the same solid angle. The horizon row uses
//! the user-supplied `az_step_at_horizon`; near the zenith, rows collapse to
//! a single cell. This is the iso-latitude ring partition commonly used by
//! HEALPix-like quick samplers.
//!
//! ## Solid angle
//!
//! Per-cell solid angles are returned as typed [`Steradians`]:
//!
//! ```text
//! dΩ = cos(alt) · Δalt_rad · Δaz_rad
//! ```
//!
//! The hemisphere integrates to `≈ 2π sr` (see the tests for tolerances).
//!
//! ## Example
//!
//! ```rust
//! use siderust::geometry::SkyGrid;
//! use siderust::qtty::DEG;
//!
//! let grid = SkyGrid::uniform(5.0 * DEG);
//! for dir in &grid {
//!     // dir: Direction<Horizontal>
//!     let _ = dir.polar;     // altitude
//!     let _ = dir.azimuth;   // azimuth (north-clockwise)
//! }
//! ```

use std::f64::consts::PI;

use crate::coordinates::frames;
use crate::coordinates::spherical;
use crate::qtty::{Degrees, Steradians};

// ─────────────────────────────────────────────────────────────────────────────
// Public types
// ─────────────────────────────────────────────────────────────────────────────

/// A single cell of a [`SkyGrid`]: a sky direction together with its
/// approximate solid angle.
#[derive(Debug, Clone, Copy)]
pub struct SkyGridCell {
    /// Sky direction (altitude = polar, azimuth) in the siderust Horizontal frame.
    ///
    /// - `direction.polar`   — altitude above the horizon, in `[-90°, +90°]`
    /// - `direction.azimuth` — azimuth from North (clockwise), in `[0°, 360°)`
    pub direction: spherical::direction::Horizontal,

    /// Approximate solid angle subtended by this cell.
    ///
    /// - Uniform grid: `dΩ = cos(alt) · Δalt_rad · Δaz_rad`.
    /// - Equal-area grid: `dΩ = cos(alt) · Δalt_rad · (2π / n_az(alt))`.
    pub solid_angle: Steradians,
}

/// Typed hemispherical alt/az grid sampler.
///
/// Constructors:
/// - [`SkyGrid::uniform`] — fixed `Δalt = Δaz = step`.
/// - [`SkyGrid::with_steps`] — independent altitude / azimuth steps.
/// - [`SkyGrid::equal_area`] — azimuth count scales with `cos(alt)`.
///
/// Builders:
/// - [`SkyGrid::with_alt_range`] — restrict the altitude range
///   (e.g. `10°..=90°` to apply a horizon mask).
///
/// The grid is *immutable*; `IntoIterator` and the explicit accessors
/// ([`SkyGrid::iter`], [`SkyGrid::iter_cells`], [`SkyGrid::with_solid_angle`])
/// each materialise a fresh traversal.
///
/// # Examples
///
/// Plain direction iteration via `IntoIterator`:
///
/// ```rust
/// use siderust::geometry::SkyGrid;
/// use siderust::qtty::DEG;
///
/// let grid = SkyGrid::uniform(10.0 * DEG);
/// let n: usize = (&grid).into_iter().count();
/// assert_eq!(n, 9 * 36);
/// ```
///
/// Direction + per-cell solid angle:
///
/// ```rust
/// use siderust::geometry::SkyGrid;
/// use siderust::qtty::{DEG, Steradians};
///
/// let grid = SkyGrid::uniform(5.0 * DEG);
/// let total: f64 = grid
///     .with_solid_angle()
///     .map(|(_, sr): (_, Steradians)| sr.value())
///     .sum();
/// assert!((total - 2.0 * std::f64::consts::PI).abs() / (2.0 * std::f64::consts::PI) < 0.01);
/// ```
#[derive(Debug, Clone)]
pub struct SkyGrid {
    alt_min: Degrees,
    alt_max: Degrees,
    alt_step: Degrees,
    az_step: Degrees,
    /// If true, az count per row scales with `cos(alt)`.
    equal_solid_angle: bool,
}

// ─────────────────────────────────────────────────────────────────────────────
// Constructors & builder
// ─────────────────────────────────────────────────────────────────────────────

impl SkyGrid {
    /// Uniform grid with equal altitude and azimuth steps, covering the upper
    /// hemisphere `[0°, 90°)`.
    pub fn uniform(step: Degrees) -> Self {
        Self::with_steps(step, step)
    }

    /// Uniform grid with independent altitude and azimuth steps, covering the
    /// upper hemisphere `[0°, 90°)`.
    pub fn with_steps(alt_step: Degrees, az_step: Degrees) -> Self {
        Self {
            alt_min: Degrees::new(0.0),
            alt_max: Degrees::new(90.0),
            alt_step,
            az_step,
            equal_solid_angle: false,
        }
    }

    /// Equal-area grid covering `[0°, 90°)`.
    ///
    /// `az_step_at_horizon` controls the azimuth granularity at the horizon
    /// row; higher altitude rings use proportionally fewer cells:
    ///
    /// ```text
    /// n_az(alt) = max(1, round(2π · cos(alt) / Δaz_horizon))
    /// ```
    pub fn equal_area(alt_step: Degrees, az_step_at_horizon: Degrees) -> Self {
        Self {
            alt_min: Degrees::new(0.0),
            alt_max: Degrees::new(90.0),
            alt_step,
            az_step: az_step_at_horizon,
            equal_solid_angle: true,
        }
    }

    /// Override the altitude range (builder pattern). Useful for applying a
    /// horizon mask such as `10°..=90°`.
    pub fn with_alt_range(mut self, lo: Degrees, hi: Degrees) -> Self {
        self.alt_min = lo;
        self.alt_max = hi;
        self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Accessors
// ─────────────────────────────────────────────────────────────────────────────

impl SkyGrid {
    fn n_alt(&self) -> usize {
        let span = self.alt_max.value() - self.alt_min.value();
        if span <= 0.0 || self.alt_step.value() <= 0.0 {
            return 0;
        }
        (span / self.alt_step.value()).round() as usize
    }

    fn n_az_uniform(&self) -> usize {
        if self.az_step.value() <= 0.0 {
            return 0;
        }
        (360.0_f64 / self.az_step.value()).round() as usize
    }

    fn n_az_equal_area(&self, alt_center_deg: f64) -> usize {
        let cos_alt = alt_center_deg.to_radians().cos().max(0.0);
        let az_horizon_rad = self.az_step.value().to_radians();
        if az_horizon_rad <= 0.0 {
            return 1;
        }
        let n = ((2.0 * PI * cos_alt) / az_horizon_rad).round() as usize;
        n.max(1)
    }

    /// Exact number of cells yielded by [`iter`](Self::iter).
    pub fn len(&self) -> usize {
        let n_alt = self.n_alt();
        if self.equal_solid_angle {
            (0..n_alt)
                .map(|i| {
                    let alt = self.alt_min.value() + (i as f64 + 0.5) * self.alt_step.value();
                    self.n_az_equal_area(alt)
                })
                .sum()
        } else {
            n_alt * self.n_az_uniform()
        }
    }

    /// Returns `true` if the grid contains no cells.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Iterate over directions only.
    ///
    /// Order: **altitude-outer ascending, azimuth-inner ascending**. The first
    /// azimuth in each row is `0.5·Δaz`.
    pub fn iter(&self) -> impl Iterator<Item = spherical::direction::Horizontal> + '_ {
        self.iter_cells().map(|c| c.direction)
    }

    /// Iterate over `(direction, solid_angle)` pairs.
    pub fn with_solid_angle(
        &self,
    ) -> impl Iterator<Item = (spherical::direction::Horizontal, Steradians)> + '_ {
        self.iter_cells().map(|c| (c.direction, c.solid_angle))
    }

    /// Iterate over full [`SkyGridCell`] records.
    pub fn iter_cells(&self) -> impl Iterator<Item = SkyGridCell> + '_ {
        let n_alt = self.n_alt();
        let alt_step_rad = self.alt_step.value().to_radians();
        let equal_area = self.equal_solid_angle;
        let alt_min = self.alt_min.value();
        let alt_step = self.alt_step.value();
        let az_step_uniform = self.az_step.value();
        let n_az_uniform = self.n_az_uniform();

        (0..n_alt).flat_map(move |i| {
            let alt_deg = alt_min + (i as f64 + 0.5) * alt_step;
            let alt_rad = alt_deg.to_radians();
            let cos_alt = alt_rad.cos();

            let n_az: usize;
            let az_step_deg: f64;

            if equal_area {
                n_az = {
                    let cos_clamped = cos_alt.max(0.0);
                    let az_horizon_rad = az_step_uniform.to_radians();
                    let raw = if az_horizon_rad > 0.0 {
                        ((2.0 * PI * cos_clamped) / az_horizon_rad).round() as usize
                    } else {
                        1
                    };
                    raw.max(1)
                };
                az_step_deg = 360.0 / n_az as f64;
            } else {
                n_az = n_az_uniform;
                az_step_deg = az_step_uniform;
            }

            let az_step_rad = az_step_deg.to_radians();
            let solid_angle = Steradians::new(cos_alt * alt_step_rad * az_step_rad);

            (0..n_az).map(move |j| {
                let az_deg = (j as f64 + 0.5) * az_step_deg;
                SkyGridCell {
                    direction: spherical::Direction::<frames::Horizontal>::new_unchecked(
                        Degrees::new(alt_deg),
                        Degrees::new(az_deg),
                    ),
                    solid_angle,
                }
            })
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// IntoIterator: lets `for dir in &grid {}` yield bare `Direction<Horizontal>`.
// ─────────────────────────────────────────────────────────────────────────────

impl<'a> IntoIterator for &'a SkyGrid {
    type Item = spherical::direction::Horizontal;
    type IntoIter = Box<dyn Iterator<Item = spherical::direction::Horizontal> + 'a>;

    fn into_iter(self) -> Self::IntoIter {
        Box::new(self.iter())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::DEG;

    /// 5° regular hemispherical grid: 18 alt rows × 72 az = 1296 cells.
    #[test]
    fn regular_5deg_hemisphere_count() {
        let grid = SkyGrid::uniform(5.0 * DEG);
        assert_eq!(grid.len(), 18 * 72);
        assert_eq!(grid.iter().count(), 18 * 72);
    }

    /// 10° regular grid: 9 × 36 = 324 cells.
    #[test]
    fn regular_10deg_count() {
        let grid = SkyGrid::uniform(10.0 * DEG);
        assert_eq!(grid.len(), 324);
    }

    /// `len()` matches counted iteration across many step sizes and modes.
    #[test]
    fn len_matches_iter_count() {
        for step in [1.0_f64, 5.0, 10.0, 15.0, 30.0, 45.0] {
            let g = SkyGrid::uniform(step * DEG);
            assert_eq!(g.len(), g.iter().count(), "uniform step={step}");
        }
        for step in [5.0_f64, 10.0, 15.0, 30.0] {
            let g = SkyGrid::equal_area(step * DEG, step * DEG);
            assert_eq!(g.len(), g.iter().count(), "equal_area step={step}");
        }
    }

    /// Uniform hemisphere: ∑ dΩ ≈ 2π sr (within 1%).
    #[test]
    fn uniform_solid_angle_sum_approx_two_pi() {
        let grid = SkyGrid::uniform(1.0 * DEG);
        let total: f64 = grid.with_solid_angle().map(|(_, sr)| sr.value()).sum();
        let expected = 2.0 * PI;
        let rel_err = (total - expected).abs() / expected;
        assert!(rel_err < 0.01, "uniform sum={total}, rel_err={rel_err}");
    }

    /// Equal-area hemisphere: ∑ dΩ ≈ 2π sr (within 0.1%).
    #[test]
    fn equal_area_solid_angle_sum_approx_two_pi() {
        let grid = SkyGrid::equal_area(1.0 * DEG, 1.0 * DEG);
        let total: f64 = grid.with_solid_angle().map(|(_, sr)| sr.value()).sum();
        let expected = 2.0 * PI;
        let rel_err = (total - expected).abs() / expected;
        assert!(
            rel_err < 0.001,
            "equal-area sum={total}, rel_err={rel_err}"
        );
    }

    /// All cells lie strictly inside the requested alt/az ranges (no
    /// duplicates, no boundary leakage).
    #[test]
    fn all_directions_inside_ranges() {
        let grid = SkyGrid::uniform(5.0 * DEG);
        for dir in grid.iter() {
            let alt = dir.polar.value();
            let az = dir.azimuth.value();
            assert!(
                (0.0..=90.0).contains(&alt),
                "alt {alt} out of [0, 90]"
            );
            assert!((0.0..360.0).contains(&az), "az {az} out of [0, 360)");
        }
    }

    /// First cell sits at `(alt_min + 0.5·Δalt, 0.5·Δaz)`; no zenith/horizon
    /// duplicates are emitted.
    #[test]
    fn boundary_cells_no_duplicates() {
        let grid = SkyGrid::uniform(10.0 * DEG);
        let cells: Vec<_> = grid.iter().collect();

        let first = cells.first().expect("non-empty");
        assert!((first.polar.value() - 5.0).abs() < 1e-9);
        assert!((first.azimuth.value() - 5.0).abs() < 1e-9);

        let last = cells.last().expect("non-empty");
        // Last alt row centre is 85°, last azimuth centre is 355°.
        assert!((last.polar.value() - 85.0).abs() < 1e-9);
        assert!((last.azimuth.value() - 355.0).abs() < 1e-9);

        // No two cells share both coordinates.
        for (i, a) in cells.iter().enumerate() {
            for b in &cells[i + 1..] {
                let dup = (a.polar.value() - b.polar.value()).abs() < 1e-12
                    && (a.azimuth.value() - b.azimuth.value()).abs() < 1e-12;
                assert!(!dup, "duplicate cell at ({}, {})", a.polar.value(), a.azimuth.value());
            }
        }
    }

    /// Equal-area: zenith ring collapses to a single cell, horizon ring is the
    /// densest, and all cells stay in range.
    #[test]
    fn equal_area_boundary_rings() {
        let grid = SkyGrid::equal_area(10.0 * DEG, 10.0 * DEG);
        let cells: Vec<_> = grid.iter().collect();

        // Zenith-most centre is 85°: round(2π·cos(85°) / 10°_rad) cells.
        let cos_top = 85.0_f64.to_radians().cos();
        let expected_top = (((2.0 * PI * cos_top) / 10.0_f64.to_radians()).round() as usize).max(1);
        let top_count = cells
            .iter()
            .filter(|c| (c.polar.value() - 85.0).abs() < 1e-9)
            .count();
        assert_eq!(top_count, expected_top);

        // All cells in canonical ranges.
        for c in &cells {
            assert!((0.0..=90.0).contains(&c.polar.value()));
            assert!((0.0..360.0).contains(&c.azimuth.value()));
        }
    }

    /// `with_alt_range` applies a horizon mask correctly.
    #[test]
    fn horizon_mask_via_with_alt_range() {
        let grid = SkyGrid::uniform(10.0 * DEG)
            .with_alt_range(30.0 * DEG, 60.0 * DEG);
        assert_eq!(grid.len(), 3 * 36);
        for dir in grid.iter() {
            assert!(dir.polar.value() >= 30.0 && dir.polar.value() <= 60.0);
        }
    }

    /// `IntoIterator` over `&SkyGrid` yields bare directions.
    #[test]
    fn into_iterator_yields_directions() {
        use crate::coordinates::frames::Horizontal;
        use crate::coordinates::spherical::Direction;

        let grid = SkyGrid::uniform(15.0 * DEG);
        let mut count = 0usize;
        for dir in &grid {
            // Compile-time frame check.
            let _: Direction<Horizontal> = dir;
            count += 1;
        }
        assert_eq!(count, grid.len());
    }

    /// Solid-angle is typed as `Steradians`.
    #[test]
    fn solid_angle_is_typed_steradians() {
        let grid = SkyGrid::uniform(30.0 * DEG);
        let cell = grid.iter_cells().next().unwrap();
        let _: Steradians = cell.solid_angle;
        assert!(cell.solid_angle.value() > 0.0);
    }
}
