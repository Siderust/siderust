// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Scattering phase functions
//!
//! ## Scientific scope
//!
//! A scattering phase function `P(λ, θ)` describes the angular
//! redistribution of light scattered by a medium. It is normalised so
//! that its integral over the sphere equals 1 (or 4π, depending on
//! convention) and feeds directly into single-scattering radiative
//! transfer, sky-brightness models, and aerosol retrieval.
//!
//! This module provides:
//!
//! - the dimensionless [`ScatteringFactor`] unit shared across the
//!   atmosphere submodules,
//! - a generic [`PhaseFunction`] trait,
//! - the closed-form [`RayleighPhaseFunction`], and
//! - a tabulated implementation [`TabulatedPhaseFunction`] for empirical
//!   or Mie phase functions, backed by [`optica::grid::Grid2D`].
//!
//! ## Technical scope
//!
//! All phase-function evaluators take typed [`Nanometers`] wavelengths
//! and typed [`crate::qtty::Radians`] scattering angles, and return typed
//! [`Quantity<ScatteringFactor>`] values.
//!
//! ## References
//!
//! - Hansen, J. E., & Travis, L. D. (1974). "Light scattering in
//!   planetary atmospheres". *Space Science Reviews* 16, 527.

use crate::atmosphere::rayleigh::rayleigh_phase;
use crate::ext_qtty::{Dimensionless, Quantity, Unit};
use crate::qtty::unit::{Degree, Nanometer};
use crate::qtty::{Nanometers, Radians};

/// Dimensionless value marker for phase functions and angular correction
/// factors.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ScatteringFactor;

impl Unit for ScatteringFactor {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "";
}

/// A wavelength-dependent scattering phase function.
pub trait PhaseFunction {
    /// Evaluate the phase function at `scattering_angle` for `wavelength`.
    fn phase(
        &self,
        wavelength: Nanometers,
        scattering_angle: Radians,
    ) -> Quantity<ScatteringFactor>;
}

/// Rayleigh phase function `P(θ) = (3 / 16π)·(1 + cos²θ)`.
#[derive(Debug, Clone, Copy, Default)]
pub struct RayleighPhaseFunction;

impl PhaseFunction for RayleighPhaseFunction {
    #[inline]
    fn phase(
        &self,
        _wavelength: Nanometers,
        scattering_angle: Radians,
    ) -> Quantity<ScatteringFactor> {
        rayleigh_phase(scattering_angle)
    }
}

/// Generic tabulated phase or correction function over wavelength and
/// scattering angle, backed by a [`optica::grid::Grid2D`].
///
/// Evaluation is infalible: out-of-range queries are clamped to the nearest
/// endpoint. Construct with [`from_raw_row_major`](Self::from_raw_row_major).
#[derive(Debug, Clone)]
pub struct TabulatedPhaseFunction {
    grid: optica::grid::Grid2D<Nanometer, Degree, ScatteringFactor>,
}

impl TabulatedPhaseFunction {
    /// Construct from row-major data.
    ///
    /// Rows are indexed by `angles_deg` (y-axis) and columns by
    /// `wavelengths_nm` (x-axis). Values are stored `[ny][nx]`.
    ///
    /// # Errors
    ///
    /// Returns [`optica::grid::GridError`] when axes are invalid or the
    /// value count does not match `wavelengths_nm.len() * angles_deg.len()`.
    pub fn from_raw_row_major(
        wavelengths_nm: &[f64],
        angles_deg: &[f64],
        row_major: &[f64],
    ) -> Result<Self, optica::grid::GridError> {
        Ok(Self {
            grid: optica::grid::Grid2D::from_raw_row_major(
                wavelengths_nm,
                angles_deg,
                row_major,
                optica::grid::OutOfRange::ClampToEndpoints,
            )?,
        })
    }

    /// Attach provenance to the underlying grid.
    pub fn with_provenance(mut self, provenance: optica::data::Provenance) -> Self {
        self.grid = self.grid.with_provenance(provenance);
        self
    }

    /// Interpolate at a typed wavelength and scattering angle.
    ///
    /// Out-of-range queries are clamped to the nearest endpoint.
    pub fn interp_at(
        &self,
        wavelength: Nanometers,
        scattering_angle: Radians,
    ) -> Quantity<ScatteringFactor> {
        self.grid
            .interp_at(wavelength, scattering_angle.to::<Degree>())
    }
}

impl PhaseFunction for TabulatedPhaseFunction {
    fn phase(
        &self,
        wavelength: Nanometers,
        scattering_angle: Radians,
    ) -> Quantity<ScatteringFactor> {
        self.interp_at(wavelength, scattering_angle)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::{unit::Radian, Degrees};

    #[test]
    fn rayleigh_phase_is_symmetric() {
        let r = RayleighPhaseFunction;
        let p0 = r
            .phase(Nanometers::new(550.0), Degrees::new(0.0).to::<Radian>())
            .value();
        let p180 = r
            .phase(Nanometers::new(550.0), Degrees::new(180.0).to::<Radian>())
            .value();
        let p90 = r
            .phase(Nanometers::new(550.0), Degrees::new(90.0).to::<Radian>())
            .value();
        assert!((p0 - p180).abs() < 1.0e-15);
        assert!(p0 > p90);
    }

    #[test]
    fn tabulated_phase_interpolates() {
        let t = TabulatedPhaseFunction::from_raw_row_major(
            &[500.0, 600.0],
            &[0.0, 90.0],
            &[1.0, 2.0, 3.0, 4.0],
        )
        .unwrap();
        let v = t
            .phase(Nanometers::new(550.0), Degrees::new(45.0).to::<Radian>())
            .value();
        assert!((v - 2.5).abs() < 1.0e-12);
    }
}
