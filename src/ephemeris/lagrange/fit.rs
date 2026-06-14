// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Sun-Earth Lagrange Chebyshev fitting
//!
//! ## Scientific scope
//!
//! Fits short or long Sun-Earth L1-L5 ephemeris archives from the instantaneous
//! N-body equilibrium solver. The official production target is 32-day blocks
//! across 1900-2100 validated every six hours with a 100 metre absolute error
//! ceiling.
//!
//! ## Technical scope
//!
//! The fitter samples the solver at Chebyshev root nodes, writes flat records
//! `[mid_seconds, radius_seconds, x_coeffs..., y_coeffs..., z_coeffs...]`, and
//! validates against a regular time grid. Stored coefficients are kilometres;
//! public evaluators convert to astronomical units.
//!
//! ## References
//!
//! - Mason, J. C., Handscomb, D. C. (2003). *Chebyshev Polynomials*.
//! - Park, R. S., et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440
//!   and DE441". *The Astronomical Journal* 161, 105.

use super::solver::{solve_sun_earth_lagrange_dyn_with_config, SolverConfig, SolverError};
use super::{evaluate_records, SunEarthLagrangePoint};
use crate::coordinates::cartesian::Position;
use crate::coordinates::centers::Barycentric;
use crate::coordinates::frames::EclipticMeanJ2000;
use crate::ephemeris::DynEphemeris;
use crate::qtty::{AstronomicalUnit, Kilometer, Meters, Second, Seconds};
use crate::time::JulianDate;
use std::fmt;

const SECONDS_PER_DAY: f64 = crate::qtty::time::SECONDS_PER_DAY;
const J2000_JD: f64 = tempoch::J2000_JD_TT_DAY.value();
const DEFAULT_FROM_JD: f64 = 2_415_020.5;
const DEFAULT_TO_JD: f64 = 2_488_070.5;
const DEFAULT_BLOCK_SECONDS: f64 = 32.0 * SECONDS_PER_DAY;
const DEFAULT_VALIDATION_SECONDS: f64 = 6.0 * 3_600.0;
const DEFAULT_TOLERANCE_METERS: f64 = 100.0;

/// Chebyshev fitting configuration.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct FitConfig {
    /// First epoch to fit.
    pub from: JulianDate,
    /// End epoch to fit.
    pub to: JulianDate,
    /// Chebyshev block size.
    pub block: Second,
    /// Validation grid step.
    pub validation_step: Second,
    /// Maximum allowed absolute validation error.
    pub tolerance: Meters,
    /// Solver controls for each sample.
    pub solver: SolverConfig,
}

impl Default for FitConfig {
    fn default() -> Self {
        Self {
            from: crate::time::try_jd_f64(DEFAULT_FROM_JD).expect("valid default start JD"),
            to: crate::time::try_jd_f64(DEFAULT_TO_JD).expect("valid default end JD"),
            block: Seconds::new(DEFAULT_BLOCK_SECONDS),
            validation_step: Seconds::new(DEFAULT_VALIDATION_SECONDS),
            tolerance: Meters::new(DEFAULT_TOLERANCE_METERS),
            solver: SolverConfig::default(),
        }
    }
}

/// Summary of validation errors for one fitted point.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct FitStats {
    /// Maximum absolute position error over validation samples.
    pub max_abs_error: Meters,
    /// Root-mean-square position error over validation samples.
    pub rms_error: Meters,
}

/// Fitted records and validation statistics for one Lagrange point.
#[derive(Debug, Clone, PartialEq)]
pub struct FittedPoint {
    /// Lagrange point represented by these records.
    pub point: SunEarthLagrangePoint,
    /// Number of Chebyshev coefficients per coordinate.
    pub ncoeff: usize,
    /// Flat record payload.
    pub records: Vec<f64>,
    /// Validation statistics.
    pub stats: FitStats,
}

/// Errors returned by Chebyshev fitting.
#[derive(Debug, Clone, PartialEq)]
pub enum FitError {
    /// The supplied fitting configuration is invalid.
    InvalidConfig,
    /// A typed Julian Date could not be constructed for a sample.
    InvalidEpoch,
    /// The N-body solver failed for a sample.
    Solver(SolverError),
    /// No candidate order up to 48 coefficients met the requested tolerance.
    ToleranceNotMet {
        /// Best maximum absolute error reached.
        max_abs_error: Meters,
    },
}

impl fmt::Display for FitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfig => write!(f, "invalid Lagrange fit configuration"),
            Self::InvalidEpoch => write!(f, "invalid Julian Date while fitting Lagrange records"),
            Self::Solver(err) => write!(f, "solver failed while fitting Lagrange records: {err}"),
            Self::ToleranceNotMet { max_abs_error } => write!(
                f,
                "Lagrange Chebyshev fit did not meet tolerance; best max error={} m",
                max_abs_error.value()
            ),
        }
    }
}

impl std::error::Error for FitError {}

impl From<SolverError> for FitError {
    fn from(value: SolverError) -> Self {
        Self::Solver(value)
    }
}

/// Fits one Sun-Earth Lagrange point using the standard coefficient candidates.
///
/// # Arguments
///
/// - `ephemeris`: Runtime ephemeris provider.
/// - `point`: L1-L5 selector.
/// - `config`: Fit range, block size, validation step, and tolerance.
///
/// # Returns
///
/// Flat Chebyshev records and validation statistics.
///
/// # Errors
///
/// Returns [`FitError`] if the configuration is invalid, a solver sample fails,
/// or no candidate order satisfies the tolerance.
pub fn fit_sun_earth_lagrange(
    ephemeris: &dyn DynEphemeris,
    point: SunEarthLagrangePoint,
    config: FitConfig,
) -> Result<FittedPoint, FitError> {
    let mut best = None;
    for ncoeff in [16_usize, 24, 32, 40, 48] {
        let fitted = match ncoeff {
            16 => fit_with_order::<16>(ephemeris, point, config)?,
            24 => fit_with_order::<24>(ephemeris, point, config)?,
            32 => fit_with_order::<32>(ephemeris, point, config)?,
            40 => fit_with_order::<40>(ephemeris, point, config)?,
            48 => fit_with_order::<48>(ephemeris, point, config)?,
            _ => unreachable!(),
        };
        let max_error = fitted.stats.max_abs_error.value();
        best = Some(fitted.clone());
        if max_error <= config.tolerance.value() {
            return Ok(fitted);
        }
    }
    let max_abs_error = best.map_or(Meters::new(f64::INFINITY), |fitted| {
        fitted.stats.max_abs_error
    });
    Err(FitError::ToleranceNotMet { max_abs_error })
}

/// Evaluates fitted records at an epoch.
///
/// # Arguments
///
/// - `fitted`: Fitted Chebyshev payload.
/// - `jd`: Epoch as a typed Julian Date.
///
/// # Returns
///
/// Barycentric ecliptic J2000 position in astronomical units.
///
/// # Errors
///
/// Returns an ephemeris out-of-range error if `jd` lies outside the fitted records.
pub fn evaluate_fitted(
    fitted: &FittedPoint,
    jd: JulianDate,
) -> Result<
    Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>,
    crate::ephemeris::EphemerisError,
> {
    evaluate_records(&fitted.records, fitted.ncoeff, jd)
}

fn fit_with_order<const N: usize>(
    ephemeris: &dyn DynEphemeris,
    point: SunEarthLagrangePoint,
    config: FitConfig,
) -> Result<FittedPoint, FitError> {
    validate_config(config)?;
    let from_s = seconds_since_j2000(config.from);
    let to_s = seconds_since_j2000(config.to);
    let block_s = config.block.value();
    let mut records = Vec::new();
    let mut start = from_s;

    while start < to_s {
        let end = (start + block_s).min(to_s);
        let mid = 0.5 * (start + end);
        let radius = 0.5 * (end - start);
        let nodes = cheby::nodes::<N>();
        let mut xs = [0.0_f64; N];
        let mut ys = [0.0_f64; N];
        let mut zs = [0.0_f64; N];
        for (idx, node) in nodes.iter().enumerate() {
            let sample_s = mid + radius * node;
            let jd = jd_from_seconds(sample_s)?;
            // Each node gets an independent CR3BP seed. Propagating seeds backward
            // through Chebyshev nodes (ordered high-to-low within the block) causes
            // Newton overshoots to spurious near-zero-residual points.
            let solution = solve_sun_earth_lagrange_dyn_with_config(
                ephemeris,
                point,
                jd,
                config.solver,
                None,
            )?;
            xs[idx] = solution.position.x().value();
            ys[idx] = solution.position.y().value();
            zs[idx] = solution.position.z().value();
        }
        let cx = cheby::fit_coeffs(&xs);
        let cy = cheby::fit_coeffs(&ys);
        let cz = cheby::fit_coeffs(&zs);
        records.push(mid);
        records.push(radius);
        records.extend_from_slice(&cx);
        records.extend_from_slice(&cy);
        records.extend_from_slice(&cz);
        start = end;
    }

    let stats = validate_fit(ephemeris, point, config, N, &records)?;
    Ok(FittedPoint {
        point,
        ncoeff: N,
        records,
        stats,
    })
}

fn validate_fit(
    ephemeris: &dyn DynEphemeris,
    point: SunEarthLagrangePoint,
    config: FitConfig,
    ncoeff: usize,
    records: &[f64],
) -> Result<FitStats, FitError> {
    let mut sample_s = seconds_since_j2000(config.from);
    let end_s = seconds_since_j2000(config.to);
    let step_s = config.validation_step.value();
    let mut max_error = 0.0_f64;
    let mut sum_sq = 0.0_f64;
    let mut count = 0_usize;
    while sample_s <= end_s + 1.0e-9 {
        let jd = jd_from_seconds(sample_s)?;
        // Each validation sample uses an independent CR3BP seed. Forward seed
        // propagation also causes Newton overshoots across frame rotations.
        let truth =
            solve_sun_earth_lagrange_dyn_with_config(ephemeris, point, jd, config.solver, None)?;
        let fit = evaluate_records(records, ncoeff, jd).map_err(|_| FitError::InvalidConfig)?;
        let fit_km = fit.to_unit::<Kilometer>();
        let err_km = distance_km(truth.position, fit_km);
        let err_m = err_km * 1_000.0;
        max_error = max_error.max(err_m);
        sum_sq += err_m * err_m;
        count += 1;
        sample_s += step_s;
    }
    let rms = if count == 0 {
        0.0
    } else {
        (sum_sq / count as f64).sqrt()
    };
    Ok(FitStats {
        max_abs_error: Meters::new(max_error),
        rms_error: Meters::new(rms),
    })
}

fn validate_config(config: FitConfig) -> Result<(), FitError> {
    let from = seconds_since_j2000(config.from);
    let to = seconds_since_j2000(config.to);
    let block = config.block.value();
    let step = config.validation_step.value();
    if from.is_finite()
        && to.is_finite()
        && to > from
        && block.is_finite()
        && block > 0.0
        && step.is_finite()
        && step > 0.0
    {
        Ok(())
    } else {
        Err(FitError::InvalidConfig)
    }
}

fn seconds_since_j2000(jd: JulianDate) -> f64 {
    (jd.raw().value() - J2000_JD) * SECONDS_PER_DAY
}

fn jd_from_seconds(seconds: f64) -> Result<JulianDate, FitError> {
    crate::time::try_jd_f64(J2000_JD + seconds / SECONDS_PER_DAY)
        .map_err(|_| FitError::InvalidEpoch)
}

fn distance_km(
    lhs: Position<Barycentric, EclipticMeanJ2000, Kilometer>,
    rhs: Position<Barycentric, EclipticMeanJ2000, Kilometer>,
) -> f64 {
    let dx = lhs.x().value() - rhs.x().value();
    let dy = lhs.y().value() - rhs.y().value();
    let dz = lhs.z().value() - rhs.z().value();
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ephemeris::Vsop87Ephemeris;

    #[test]
    fn tiny_fixture_fits_and_evaluates() {
        let config = FitConfig {
            from: crate::time::try_jd_f64(J2000_JD).expect("valid JD"),
            to: crate::time::try_jd_f64(J2000_JD + 0.25).expect("valid JD"),
            block: Seconds::new(0.125 * SECONDS_PER_DAY),
            validation_step: Seconds::new(0.125 * SECONDS_PER_DAY),
            tolerance: Meters::new(5_000_000.0),
            solver: SolverConfig::default(),
        };
        let eph = Vsop87Ephemeris;
        let fitted = fit_sun_earth_lagrange(&eph, SunEarthLagrangePoint::L1, config)
            .expect("fixture should fit");
        assert!(!fitted.records.is_empty());
        let evaluated = evaluate_fitted(&fitted, config.from).expect("fit should evaluate");
        assert!(evaluated.x().value().is_finite());
    }
}
