// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # N-body Sun-Earth Lagrange solver
//!
//! ## Scientific scope
//!
//! Solves instantaneous Sun-Earth L1-L5 equilibrium positions in a rotating
//! Sun-Earth frame with lunar gravity included as a third-body perturbation.
//! The result is an epoch-local equilibrium point, not a halo orbit or a
//! propagated spacecraft trajectory.
//!
//! ## Technical scope
//!
//! The solver accepts either a compile-time [`Ephemeris`] backend or a runtime
//! [`DynEphemeris`] provider. It works internally in kilometres and seconds,
//! using Newton iteration with a finite-difference 3 × 3 Jacobian and CR3BP
//! analytic seeds.
//!
//! ## References
//!
//! - Szebehely, V. (1967). *Theory of Orbits: The Restricted Problem of Three Bodies*.
//! - Koon, W. S., Lo, M. W., Marsden, J. E., Ross, S. D. (2011). *Dynamical Systems,
//!   the Three-Body Problem and Space Mission Design*.

use super::SunEarthLagrangePoint;
use crate::coordinates::cartesian::Position;
use crate::coordinates::centers::Barycentric;
use crate::coordinates::frames::EclipticMeanJ2000;
use crate::ephemeris::{DynEphemeris, Ephemeris, EphemerisError};
use crate::qtty::{Kilometer, Kilometers, GM_EARTH, GM_MOON, GM_SUN};
use crate::time::JulianDate;
use std::fmt;

const MAX_ITERS: usize = 30;
const STEP_TOL_KM: f64 = 1.0e-6;
const FD_STEP_KM: f64 = 100.0;
const RESIDUAL_TOL_KM_S2: f64 = 1.0e-9;
const DERIVATIVE_STEP_DAYS: f64 = 0.01;
const SECONDS_PER_DAY: f64 = crate::qtty::time::SECONDS_PER_DAY;

/// Configuration for Newton iteration.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct SolverConfig {
    /// Maximum number of Newton iterations.
    pub max_iterations: usize,
    /// Euclidean position-step convergence tolerance in kilometres.
    pub step_tolerance_km: f64,
    /// Finite-difference Jacobian perturbation in kilometres.
    pub finite_difference_step_km: f64,
    /// Residual acceleration convergence tolerance in km/s².
    pub residual_tolerance_km_s2: f64,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            max_iterations: MAX_ITERS,
            step_tolerance_km: STEP_TOL_KM,
            finite_difference_step_km: FD_STEP_KM,
            residual_tolerance_km_s2: RESIDUAL_TOL_KM_S2,
        }
    }
}

/// Converged Sun-Earth Lagrange solution.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct SolverSolution {
    /// Barycentric ecliptic J2000 position in kilometres.
    pub position: Position<Barycentric, EclipticMeanJ2000, Kilometer>,
    /// Number of Newton iterations performed.
    pub iterations: usize,
    /// Final Newton step norm in kilometres.
    pub step_norm_km: f64,
    /// Final residual acceleration norm in km/s².
    pub residual_norm_km_s2: f64,
}

/// Errors returned by the Sun-Earth Lagrange solver.
#[derive(Debug, Clone, PartialEq)]
pub enum SolverError {
    /// The backing ephemeris could not provide a required body state.
    Ephemeris(EphemerisError),
    /// The Sun-Earth rotating frame could not be constructed.
    DegenerateFrame,
    /// Newton iteration failed to converge within the configured iteration cap.
    DidNotConverge {
        /// Number of iterations attempted.
        iterations: usize,
        /// Last position-step norm in kilometres.
        step_norm_km: f64,
        /// Last residual norm in km/s².
        residual_norm_km_s2: f64,
    },
    /// The finite-difference Jacobian was singular.
    SingularJacobian,
}

impl fmt::Display for SolverError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Ephemeris(err) => write!(f, "ephemeris error while solving Lagrange point: {err}"),
            Self::DegenerateFrame => write!(f, "degenerate Sun-Earth rotating frame"),
            Self::DidNotConverge {
                iterations,
                step_norm_km,
                residual_norm_km_s2,
            } => write!(
                f,
                "Lagrange solver did not converge after {iterations} iterations; step={step_norm_km} km residual={residual_norm_km_s2} km/s²"
            ),
            Self::SingularJacobian => write!(f, "singular finite-difference Jacobian"),
        }
    }
}

impl std::error::Error for SolverError {}

impl From<EphemerisError> for SolverError {
    fn from(value: EphemerisError) -> Self {
        Self::Ephemeris(value)
    }
}

/// Solves a Sun-Earth Lagrange point with a compile-time ephemeris backend.
///
/// # Arguments
///
/// - `point`: L1-L5 selector.
/// - `jd`: Epoch as a typed Julian Date.
///
/// # Returns
///
/// Converged barycentric position and diagnostics.
///
/// # Errors
///
/// Returns [`SolverError`] if ephemeris sampling, frame construction, Jacobian
/// inversion, or Newton convergence fails.
pub fn solve_sun_earth_lagrange<Eph: Ephemeris>(
    point: SunEarthLagrangePoint,
    jd: JulianDate,
) -> Result<SolverSolution, SolverError> {
    solve_sun_earth_lagrange_with_config::<Eph>(point, jd, SolverConfig::default())
}

/// Solves a Sun-Earth Lagrange point with an explicit solver configuration.
///
/// # Arguments
///
/// - `point`: L1-L5 selector.
/// - `jd`: Epoch as a typed Julian Date.
/// - `config`: Newton iteration controls.
///
/// # Returns
///
/// Converged barycentric position and diagnostics.
///
/// # Errors
///
/// Returns [`SolverError`] if ephemeris sampling, frame construction, Jacobian
/// inversion, or Newton convergence fails.
pub fn solve_sun_earth_lagrange_with_config<Eph: Ephemeris>(
    point: SunEarthLagrangePoint,
    jd: JulianDate,
    config: SolverConfig,
) -> Result<SolverSolution, SolverError> {
    let sampler = StaticSampler::<Eph>(std::marker::PhantomData);
    solve_with_sampler(&sampler, point, jd, config, None)
}

/// Solves a Sun-Earth Lagrange point with a runtime ephemeris backend.
///
/// # Arguments
///
/// - `ephemeris`: Runtime ephemeris provider.
/// - `point`: L1-L5 selector.
/// - `jd`: Epoch as a typed Julian Date.
///
/// # Returns
///
/// Converged barycentric position and diagnostics.
///
/// # Errors
///
/// Returns [`SolverError`] if ephemeris sampling, frame construction, Jacobian
/// inversion, or Newton convergence fails.
pub fn solve_sun_earth_lagrange_dyn(
    ephemeris: &dyn DynEphemeris,
    point: SunEarthLagrangePoint,
    jd: JulianDate,
) -> Result<SolverSolution, SolverError> {
    solve_sun_earth_lagrange_dyn_with_config(ephemeris, point, jd, SolverConfig::default(), None)
}

/// Solves a runtime-ephemeris point with explicit configuration and optional seed.
///
/// # Arguments
///
/// - `ephemeris`: Runtime ephemeris provider.
/// - `point`: L1-L5 selector.
/// - `jd`: Epoch as a typed Julian Date.
/// - `config`: Newton iteration controls.
/// - `seed`: Optional barycentric kilometre seed from a previous sample.
///
/// # Returns
///
/// Converged barycentric position and diagnostics.
///
/// # Errors
///
/// Returns [`SolverError`] if ephemeris sampling, frame construction, Jacobian
/// inversion, or Newton convergence fails.
pub fn solve_sun_earth_lagrange_dyn_with_config(
    ephemeris: &dyn DynEphemeris,
    point: SunEarthLagrangePoint,
    jd: JulianDate,
    config: SolverConfig,
    seed: Option<Position<Barycentric, EclipticMeanJ2000, Kilometer>>,
) -> Result<SolverSolution, SolverError> {
    let sampler = DynSampler { ephemeris };
    solve_with_sampler(&sampler, point, jd, config, seed.map(|pos| [pos.x().value(), pos.y().value(), pos.z().value()]))
}

trait Sampler {
    fn sun(&self, jd: JulianDate) -> Result<Vec3, EphemerisError>;
    fn earth(&self, jd: JulianDate) -> Result<Vec3, EphemerisError>;
    fn moon_geo(&self, jd: JulianDate) -> Result<Vec3, EphemerisError>;
}

struct StaticSampler<Eph>(std::marker::PhantomData<Eph>);

impl<Eph: Ephemeris> Sampler for StaticSampler<Eph> {
    fn sun(&self, jd: JulianDate) -> Result<Vec3, EphemerisError> {
        Ok(pos_au_to_km(Eph::try_sun_barycentric(jd)?))
    }

    fn earth(&self, jd: JulianDate) -> Result<Vec3, EphemerisError> {
        Ok(pos_au_to_km(Eph::try_earth_barycentric(jd)?))
    }

    fn moon_geo(&self, jd: JulianDate) -> Result<Vec3, EphemerisError> {
        Ok(pos_km_to_vec(Eph::try_moon_geocentric(jd)?))
    }
}

struct DynSampler<'a> {
    ephemeris: &'a dyn DynEphemeris,
}

impl Sampler for DynSampler<'_> {
    fn sun(&self, jd: JulianDate) -> Result<Vec3, EphemerisError> {
        Ok(pos_au_to_km(self.ephemeris.try_sun_barycentric(jd)?))
    }

    fn earth(&self, jd: JulianDate) -> Result<Vec3, EphemerisError> {
        Ok(pos_au_to_km(self.ephemeris.try_earth_barycentric(jd)?))
    }

    fn moon_geo(&self, jd: JulianDate) -> Result<Vec3, EphemerisError> {
        Ok(pos_km_to_vec(self.ephemeris.try_moon_geocentric(jd)?))
    }
}

fn solve_with_sampler(
    sampler: &dyn Sampler,
    point: SunEarthLagrangePoint,
    jd: JulianDate,
    config: SolverConfig,
    seed: Option<Vec3>,
) -> Result<SolverSolution, SolverError> {
    let frame = RotatingFrame::new(sampler, jd)?;
    let mut pos = seed.unwrap_or_else(|| frame.seed(point));
    let mut last_step = f64::INFINITY;
    let mut residual_norm = f64::INFINITY;

    for iter in 0..config.max_iterations {
        let residual = frame.residual(pos);
        residual_norm = norm(residual);
        if residual_norm <= config.residual_tolerance_km_s2 {
            return Ok(SolverSolution {
                position: Position::new(
                    Kilometers::new(pos[0]),
                    Kilometers::new(pos[1]),
                    Kilometers::new(pos[2]),
                ),
                iterations: iter,
                step_norm_km: last_step,
                residual_norm_km_s2: residual_norm,
            });
        }
        let jac = jacobian(&frame, pos, config.finite_difference_step_km);
        let mut delta = solve_linear(jac, neg(residual)).ok_or(SolverError::SingularJacobian)?;
        let max_step = 0.05 * frame.distance;
        let delta_norm = norm(delta);
        if delta_norm > max_step {
            delta = scale(delta, max_step / delta_norm);
        }
        let mut accepted = delta;
        let mut alpha = 1.0;
        for _ in 0..12 {
            let trial = add(pos, scale(delta, alpha));
            if norm(frame.residual(trial)) <= residual_norm || alpha <= 1.0 / 4096.0 {
                accepted = scale(delta, alpha);
                break;
            }
            alpha *= 0.5;
        }
        pos = add(pos, accepted);
        last_step = norm(accepted);
        if last_step <= config.step_tolerance_km {
            return Ok(SolverSolution {
                position: Position::new(
                    Kilometers::new(pos[0]),
                    Kilometers::new(pos[1]),
                    Kilometers::new(pos[2]),
                ),
                iterations: iter + 1,
                step_norm_km: last_step,
                residual_norm_km_s2: residual_norm,
            });
        }
    }

    Err(SolverError::DidNotConverge {
        iterations: config.max_iterations,
        step_norm_km: last_step,
        residual_norm_km_s2: residual_norm,
    })
}

type Vec3 = [f64; 3];
type Mat3 = [[f64; 3]; 3];

struct RotatingFrame {
    sun: Vec3,
    earth: Vec3,
    moon: Vec3,
    origin: Vec3,
    x_axis: Vec3,
    y_axis: Vec3,
    z_axis: Vec3,
    omega: f64,
    distance: f64,
    origin_accel: Vec3,
}

impl RotatingFrame {
    fn new(sampler: &dyn Sampler, jd: JulianDate) -> Result<Self, SolverError> {
        let sun = sampler.sun(jd)?;
        let earth = sampler.earth(jd)?;
        let moon = add(earth, sampler.moon_geo(jd)?);
        let rel = sub(earth, sun);
        let d = norm(rel);
        if d <= 0.0 || !d.is_finite() {
            return Err(SolverError::DegenerateFrame);
        }
        let x_axis = scale(rel, 1.0 / d);
        let velocity = relative_velocity(sampler, jd)?;
        let mut z_axis = cross(rel, velocity);
        let z_norm = norm(z_axis);
        if z_norm <= 0.0 || !z_norm.is_finite() {
            z_axis = [0.0, 0.0, 1.0];
        } else {
            z_axis = scale(z_axis, 1.0 / z_norm);
        }
        let y_axis = normalize(cross(z_axis, x_axis)).ok_or(SolverError::DegenerateFrame)?;
        let total_mu = GM_SUN.value() + GM_EARTH.value();
        let origin = scale(add(scale(sun, GM_SUN.value()), scale(earth, GM_EARTH.value())), 1.0 / total_mu);
        let omega = (total_mu / d.powi(3)).sqrt();
        let origin_accel = sun_earth_barycenter_external_accel(sun, earth, moon);
        Ok(Self { sun, earth, moon, origin, x_axis, y_axis, z_axis, omega, distance: d, origin_accel })
    }

    fn residual(&self, pos: Vec3) -> Vec3 {
        let bodies = [(self.sun, GM_SUN.value()), (self.earth, GM_EARTH.value()), (self.moon, GM_MOON.value())];
        let rho = sub(pos, self.origin);
        let omega_vec = scale(self.z_axis, self.omega);
        sub(sub(grav_accel(pos, &bodies), self.origin_accel), cross(omega_vec, cross(omega_vec, rho)))
    }

    fn seed(&self, point: SunEarthLagrangePoint) -> Vec3 {
        let rel = sub(self.earth, self.sun);
        let d = norm(rel);
        let mu = GM_EARTH.value() / (GM_SUN.value() + GM_EARTH.value());
        let gamma = (mu / 3.0).cbrt();
        let x = match point {
            SunEarthLagrangePoint::L1 => 1.0 - gamma,
            SunEarthLagrangePoint::L2 => 1.0 + gamma,
            SunEarthLagrangePoint::L3 => -1.0 - 5.0 * mu / 12.0,
            SunEarthLagrangePoint::L4 | SunEarthLagrangePoint::L5 => 0.5,
        };
        let y = match point {
            SunEarthLagrangePoint::L4 => 3.0_f64.sqrt() / 2.0,
            SunEarthLagrangePoint::L5 => -3.0_f64.sqrt() / 2.0,
            _ => 0.0,
        };
        let bary_x = (x - mu) * d;
        add(self.origin, add(scale(self.x_axis, bary_x), scale(self.y_axis, y * d)))
    }
}

fn relative_velocity(sampler: &dyn Sampler, jd: JulianDate) -> Result<Vec3, SolverError> {
    let jd_value = jd.raw().value();
    let before = crate::time::try_jd_f64(jd_value - DERIVATIVE_STEP_DAYS).map_err(|_| SolverError::DegenerateFrame)?;
    let after = crate::time::try_jd_f64(jd_value + DERIVATIVE_STEP_DAYS).map_err(|_| SolverError::DegenerateFrame)?;
    let rel_before = sub(sampler.earth(before)?, sampler.sun(before)?);
    let rel_after = sub(sampler.earth(after)?, sampler.sun(after)?);
    Ok(scale(sub(rel_after, rel_before), 1.0 / (2.0 * DERIVATIVE_STEP_DAYS * SECONDS_PER_DAY)))
}

fn pos_au_to_km<C>(pos: Position<C, EclipticMeanJ2000, crate::qtty::AstronomicalUnit>) -> Vec3
where
    C: crate::coordinates::centers::ReferenceCenter,
    C::Params: Clone,
{
    pos_km_to_vec(pos.to_unit::<Kilometer>())
}

fn pos_km_to_vec<C>(pos: Position<C, EclipticMeanJ2000, Kilometer>) -> Vec3
where
    C: crate::coordinates::centers::ReferenceCenter,
{
    [pos.x().value(), pos.y().value(), pos.z().value()]
}


fn sun_earth_barycenter_external_accel(sun: Vec3, earth: Vec3, moon: Vec3) -> Vec3 {
    let sun_due_moon = point_mass_accel(sun, moon, GM_MOON.value());
    let earth_due_moon = point_mass_accel(earth, moon, GM_MOON.value());
    let total_mu = GM_SUN.value() + GM_EARTH.value();
    scale(
        add(scale(sun_due_moon, GM_SUN.value()), scale(earth_due_moon, GM_EARTH.value())),
        1.0 / total_mu,
    )
}

fn point_mass_accel(pos: Vec3, body: Vec3, mu: f64) -> Vec3 {
    let r = sub(body, pos);
    let rnorm = norm(r);
    if rnorm > 0.0 {
        scale(r, mu / rnorm.powi(3))
    } else {
        [0.0; 3]
    }
}

fn grav_accel(pos: Vec3, bodies: &[(Vec3, f64); 3]) -> Vec3 {
    let mut acc = [0.0; 3];
    for (body, mu) in bodies {
        let r = sub(*body, pos);
        let rnorm = norm(r);
        if rnorm > 0.0 {
            acc = add(acc, scale(r, *mu / rnorm.powi(3)));
        }
    }
    acc
}

fn jacobian(frame: &RotatingFrame, pos: Vec3, h: f64) -> Mat3 {
    let mut jac = [[0.0; 3]; 3];
    for axis in 0..3 {
        let mut plus = pos;
        plus[axis] += h;
        let mut minus = pos;
        minus[axis] -= h;
        let diff = scale(sub(frame.residual(plus), frame.residual(minus)), 1.0 / (2.0 * h));
        for row in 0..3 {
            jac[row][axis] = diff[row];
        }
    }
    jac
}

fn solve_linear(mut a: Mat3, mut b: Vec3) -> Option<Vec3> {
    for col in 0..3 {
        let mut pivot = col;
        for row in (col + 1)..3 {
            if a[row][col].abs() > a[pivot][col].abs() {
                pivot = row;
            }
        }
        if a[pivot][col].abs() < 1.0e-30 {
            return None;
        }
        if pivot != col {
            a.swap(pivot, col);
            b.swap(pivot, col);
        }
        let div = a[col][col];
        for value in a[col].iter_mut().skip(col) {
            *value /= div;
        }
        b[col] /= div;
        let pivot_row = a[col];
        for (row_idx, row_values) in a.iter_mut().enumerate() {
            if row_idx != col {
                let factor = row_values[col];
                for (value, pivot_value) in row_values.iter_mut().zip(pivot_row).skip(col) {
                    *value -= factor * pivot_value;
                }
                b[row_idx] -= factor * b[col];
            }
        }
    }
    Some(b)
}

fn add(a: Vec3, b: Vec3) -> Vec3 { [a[0] + b[0], a[1] + b[1], a[2] + b[2]] }
fn sub(a: Vec3, b: Vec3) -> Vec3 { [a[0] - b[0], a[1] - b[1], a[2] - b[2]] }
fn neg(a: Vec3) -> Vec3 { [-a[0], -a[1], -a[2]] }
fn scale(a: Vec3, s: f64) -> Vec3 { [a[0] * s, a[1] * s, a[2] * s] }
fn dot(a: Vec3, b: Vec3) -> f64 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
fn norm(a: Vec3) -> f64 { dot(a, a).sqrt() }
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn normalize(a: Vec3) -> Option<Vec3> {
    let n = norm(a);
    (n > 0.0 && n.is_finite()).then(|| scale(a, 1.0 / n))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ephemeris::Vsop87Ephemeris;
    use crate::time::J2000;

    #[test]
    fn solver_converges_at_j2000_for_all_points() {
        for point in [
            SunEarthLagrangePoint::L1,
            SunEarthLagrangePoint::L2,
            SunEarthLagrangePoint::L3,
            SunEarthLagrangePoint::L4,
            SunEarthLagrangePoint::L5,
        ] {
            let solution = solve_sun_earth_lagrange::<Vsop87Ephemeris>(point, J2000)
                .expect("solver should converge");
            assert!(solution.step_norm_km <= STEP_TOL_KM * 10.0 || solution.residual_norm_km_s2 < 1.0e-6);
        }
    }
}
