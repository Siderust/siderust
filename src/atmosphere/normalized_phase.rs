// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Strict, unity-normalized atmospheric tables using Optica interpolation.

use crate::qtty::{Nanometers, Quantity, Radians};
use alloc::vec::Vec;
use core::f64::consts::{FRAC_PI_2, PI, TAU};
use optica::grid::{GridError, OutOfRange};
use optica::phase::{PhaseTable, ScatteringFactor};
#[cfg(not(feature = "std"))]
use qtty::{Scalar, Transcendental};

/// Invalid atmospheric phase data or a failed strict lookup.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum NormalizedPhaseError {
    /// Optica axis, shape, or query-domain error (x: nanometers; y: radians).
    #[error(transparent)]
    Grid(#[from] GridError),
    /// Wavelengths must be positive and strictly ascending.
    #[error("wavelengths must be positive and strictly ascending")]
    WavelengthDomain,
    /// Angles must be strictly ascending from exactly zero to exactly π.
    #[error("angles must be strictly ascending from 0 to pi radians")]
    AngularDomain,
    /// A row-major value was negative or nonfinite.
    #[error("phase value at row-major index {index} must be finite and nonnegative")]
    InvalidValue {
        /// Index in the supplied row-major data.
        index: usize,
    },
    /// The actual linear angular interpolant failed unity normalization.
    #[error("phase integral at wavelength index {index} differs from unity by more than 1e-6")]
    NotNormalized {
        /// Index in the supplied wavelength axis.
        index: usize,
    },
    /// A query coordinate was NaN or infinite.
    #[error("phase query coordinates must be finite")]
    NonFiniteQuery,
}

/// A validated atmospheric phase table with exclusively fallible evaluation.
///
/// Wavelengths are positive, finite and ascending; angles cover exactly
/// `[0, π]` radians in ascending order. Values are finite and nonnegative.
/// Each wavelength slice satisfies `2π ∫ P(θ) sin(θ) dθ = 1` to absolute
/// tolerance **1e-6**, checked on the actual piecewise-linear interpolant.
/// Wavelength interpolation is a convex combination of these slices and
/// preserves that tolerance, apart from floating-point rounding.
///
/// Segment integrals are analytic, with a small-interval series and compensated
/// summation. For half-widths below 0.01 rad, the sinc series' first omitted
/// term is below 3e-22; above that cutoff, cancellation in the closed form
/// is of order 1e-12 or less. The tolerance allows discretization error in
/// finely sampled analytic models; it is much larger than ordinary f64 integration roundoff.
/// It does not certify agreement with an unresolved analytic peak. Data are
/// never renormalized. External 4π-normalized data must first be divided by 4π.
///
/// The output uses **Optica's** dimensionless `ScatteringFactor` marker, distinct
/// from Siderust's marker. Physically the density is per steradian (sr⁻¹);
/// converting to reciprocal-solid-angle quantities requires an explicit marker
/// conversion at the consumer boundary.
///
/// Unlike [`super::TabulatedPhaseFunction`], this type rejects out-of-domain
/// queries. It deliberately implements neither infallible `PhaseFunction`
/// trait and does not expose its grid for mutation or infallible evaluation.
/// Construction and storage reuse Optica's [`PhaseTable`]; queries allocate
/// nothing and use its fallible interpolation with [`OutOfRange::Error`].
///
/// ```
/// use siderust::atmosphere::NormalizedPhaseTable;
/// use siderust::optica::phase::ScatteringFactor;
/// use siderust::qtty::{Nanometers, Quantity, Radians};
/// use core::f64::consts::PI;
/// let p = Quantity::<ScatteringFactor>::new(1.0 / (4.0 * PI));
/// let table = NormalizedPhaseTable::new(
///     &[Nanometers::new(400.0), Nanometers::new(700.0)],
///     &[Radians::new(0.0), Radians::new(PI)],
///     &[p, p, p, p],
/// )?;
/// let value = table.try_phase(Nanometers::new(550.0), Radians::new(1.0))?;
/// assert!((value.value() - p.value()).abs() < 1e-15);
/// assert!(table.try_phase(Nanometers::new(800.0), Radians::new(1.0)).is_err());
/// # Ok::<(), siderust::atmosphere::NormalizedPhaseError>(())
/// ```
#[derive(Debug, Clone)]
pub struct NormalizedPhaseTable {
    grid: PhaseTable,
}

impl NormalizedPhaseTable {
    /// Validate typed axes and row-major values (`values[iy * nx + ix]`).
    ///
    /// At least two samples on each axis are required by Optica. Axis and
    /// shape validation is delegated to Optica; the additional atmospheric
    /// invariants and normalization are checked once here.
    ///
    /// # Errors
    /// Returns [`NormalizedPhaseError`] for invalid axes, shape, values or
    /// angular normalization. Endpoint angles must equal `0.0` and `PI`.
    pub fn new(
        wavelengths: &[Nanometers],
        angles: &[Radians],
        values: &[Quantity<ScatteringFactor>],
    ) -> Result<Self, NormalizedPhaseError> {
        let xs: Vec<_> = wavelengths.iter().map(|v| v.value()).collect();
        let ys: Vec<_> = angles.iter().map(|v| v.value()).collect();
        let data: Vec<_> = values.iter().map(|v| v.value()).collect();
        let grid = PhaseTable::from_raw_row_major(&xs, &ys, &data, OutOfRange::Error)?;
        // Optica already guarantees finite, strictly monotone axes and shape.
        if xs[0] <= 0.0 || xs[0] >= xs[xs.len() - 1] {
            return Err(NormalizedPhaseError::WavelengthDomain);
        }
        if ys[0] != 0.0 || ys[ys.len() - 1] != PI {
            return Err(NormalizedPhaseError::AngularDomain);
        }
        for (index, &value) in data.iter().enumerate() {
            if !value.is_finite() || value < 0.0 {
                return Err(NormalizedPhaseError::InvalidValue { index });
            }
        }
        for ix in 0..xs.len() {
            let mut sum = 0.0;
            let mut correction = 0.0;
            for (iy, pair) in ys.windows(2).enumerate() {
                let integral = segment_integral(
                    pair[0],
                    pair[1],
                    data[iy * xs.len() + ix],
                    data[(iy + 1) * xs.len() + ix],
                );
                let adjusted = integral - correction;
                let next = sum + adjusted;
                correction = (next - sum) - adjusted;
                sum = next;
            }
            if !sum.is_finite() || (TAU * sum - 1.0).abs() > 1e-6 {
                return Err(NormalizedPhaseError::NotNormalized { index: ix });
            }
        }
        Ok(Self { grid })
    }

    /// Interpolate strictly inside the inclusive tabulated domain.
    ///
    /// # Errors
    /// Nonfinite coordinates return [`NormalizedPhaseError::NonFiniteQuery`].
    /// Finite coordinates outside the domain return an Optica [`GridError`].
    pub fn try_phase(
        &self,
        wavelength: Nanometers,
        scattering_angle: Radians,
    ) -> Result<Quantity<ScatteringFactor>, NormalizedPhaseError> {
        if !wavelength.value().is_finite() || !scattering_angle.value().is_finite() {
            return Err(NormalizedPhaseError::NonFiniteQuery);
        }
        Ok(self.grid.try_interp_at(wavelength, scattering_angle)?)
    }
}

// Exact integral of the linear segment times sin(theta). With half-width h
// and midpoint m, endpoint weights are sin(m)*sin(h) ± cos(m)*(sinc(h)-cos(h)).
// Factor out h before multiplying large phase values, and use a series for
// (sinc(h)-cos(h))/h to avoid cancellation. Reflect m near π to retain precision.
fn segment_integral(a: f64, b: f64, p: f64, q: f64) -> f64 {
    let h = (b - a) * 0.5;
    let m = a + h;
    let (sin_m, cos_m) = if m > FRAC_PI_2 {
        let distance = (PI - b) + h;
        (distance.sin(), -distance.cos())
    } else {
        (m.sin(), m.cos())
    };
    let (sinc, slope) = if h < 0.01 {
        let h2 = h * h;
        (
            1.0 + h2 * (-1.0 / 6.0 + h2 * (1.0 / 120.0 - h2 / 5040.0)),
            h * (1.0 / 3.0 + h2 * (-1.0 / 30.0 + h2 * (1.0 / 840.0 - h2 / 45360.0))),
        )
    } else {
        let sinc = h.sin() / h;
        (sinc, (sinc - h.cos()) / h)
    };
    (p * h) * (sin_m * sinc - cos_m * slope) + (q * h) * (sin_m * sinc + cos_m * slope)
}
