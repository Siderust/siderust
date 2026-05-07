// SPDX-License-Identifier: AGPL-3.0-or-later
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
//! - an optional tabulated implementation [`TabulatedPhaseFunction`]
//!   (behind the `tables` feature) for empirical / Mie phase functions.
//!
//! ## Technical scope
//!
//! All phase-function evaluators take typed [`Nanometers`] wavelengths
//! and typed [`Radians`] scattering angles, and return typed
//! [`Quantity<ScatteringFactor>`] values.
//!
//! ## References
//!
//! - Hansen, J. E., & Travis, L. D. (1974). "Light scattering in
//!   planetary atmospheres". *Space Science Reviews* 16, 527.

use crate::atmosphere::rayleigh::rayleigh_phase;
use crate::ext_qtty::{Dimensionless, Quantity, Unit};
#[cfg(feature = "tables")]
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
/// scattering angle.
#[cfg(feature = "tables")]
#[derive(Debug, Clone)]
pub struct TabulatedPhaseFunction {
    grid: crate::tables::Grid2D<Nanometer, Degree, ScatteringFactor>,
}

#[cfg(feature = "tables")]
impl TabulatedPhaseFunction {
    /// Construct from row-major data with rows indexed by `angles_deg` and
    /// columns indexed by `wavelengths_nm`.
    pub fn from_raw_row_major(
        wavelengths_nm: Vec<f64>,
        angles_deg: Vec<f64>,
        row_major: Vec<f64>,
    ) -> Result<Self, crate::tables::TableError> {
        Ok(Self {
            grid: crate::tables::Grid2D::from_raw_row_major(wavelengths_nm, angles_deg, row_major)?,
        })
    }

    /// Attach provenance to the underlying grid.
    pub fn with_provenance(mut self, provenance: crate::provenance::Provenance) -> Self {
        self.grid = self.grid.with_provenance(provenance);
        self
    }

    /// Interpolate at a typed wavelength and scattering angle.
    pub fn interp_at(
        &self,
        wavelength: Nanometers,
        scattering_angle: Radians,
    ) -> Result<Quantity<ScatteringFactor>, crate::tables::TableError> {
        self.grid.interp_at(
            wavelength,
            scattering_angle.to::<Degree>(),
            crate::tables::OutOfRange::ClampToEndpoints,
            crate::tables::OutOfRange::ClampToEndpoints,
        )
    }
}

#[cfg(feature = "tables")]
impl PhaseFunction for TabulatedPhaseFunction {
    fn phase(
        &self,
        wavelength: Nanometers,
        scattering_angle: Radians,
    ) -> Quantity<ScatteringFactor> {
        self.interp_at(wavelength, scattering_angle)
            .expect("tabulated phase-function interpolation failed")
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

    #[cfg(feature = "tables")]
    #[test]
    fn tabulated_phase_interpolates() {
        let t = TabulatedPhaseFunction::from_raw_row_major(
            vec![500.0, 600.0],
            vec![0.0, 90.0],
            vec![1.0, 2.0, 3.0, 4.0],
        )
        .unwrap();
        let v = t
            .phase(Nanometers::new(550.0), Degrees::new(45.0).to::<Radian>())
            .value();
        assert!((v - 2.5).abs() < 1.0e-12);
    }
}
