// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Finite-burn / thrust-arc physical model.
//!
//! This module provides a minimal but unit-safe thrust model usable both
//! from precise-orbit-determination workflows (`siderust-pod-dynamics`,
//! where the constant scale factor / delta-v components are estimable)
//! and from future mission-design workflows.
//!
//! The model is intentionally simple — constant thrust over a closed
//! interval, expressed in a 3-component unit vector — and avoids
//! dragging in the full mission-design machinery (Sims-Flanagan,
//! Pontryagin, low-thrust optimisation), which lives in future
//! `siderust-mission-design` / `siderust-low-thrust` crates.

use qtty::{force::Newtons, Kilogram, Second};
use tempoch::{Time, UTC};

/// Standard sea-level gravitational acceleration `g₀` used to relate
/// specific impulse to mass flow rate, in m / s².
///
/// IUPAC / ISO 80000 standard value.
pub const G0_M_PER_S2: f64 = 9.80665;

/// A constant-thrust arc applied over a closed UTC interval.
///
/// `direction` is a unit 3-vector expressed in the spacecraft frame
/// chosen by the consuming pipeline (typically VNC, RTN, or body).
/// The pipeline is responsible for projecting the resulting acceleration
/// into the propagation frame.
#[derive(Clone, Debug)]
pub struct ThrustArc {
    /// Start time of the thrust arc.
    pub start: Time<UTC>,
    /// Stop time of the thrust arc.
    pub stop: Time<UTC>,
    /// Thrust magnitude in Newtons.
    pub thrust: Newtons,
    /// Specific impulse.
    pub isp: Second,
    /// Unit thrust direction `[x, y, z]` in the consuming pipeline's frame.
    /// Magnitude is *not* validated by `ThrustArc` itself — see
    /// [`ThrustArc::with_normalised_direction`].
    pub direction: [f64; 3],
}

/// Errors produced when evaluating a thrust arc.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum ManeuverError {
    /// Thrust direction has non-finite or zero magnitude.
    #[error("thrust direction has non-finite or zero magnitude (got {0:?})")]
    InvalidDirection([f64; 3]),
    /// Specific impulse must be strictly positive.
    #[error("specific impulse must be strictly positive (got {0} s)")]
    NonPositiveIsp(f64),
    /// Spacecraft mass must be strictly positive.
    #[error("spacecraft mass must be strictly positive (got {0} kg)")]
    NonPositiveMass(f64),
    /// Thrust arc is empty: start >= stop.
    #[error("thrust arc is empty: start {start:?} >= stop {stop:?}")]
    EmptyInterval {
        /// Start time of the arc.
        start: Time<UTC>,
        /// Stop time of the arc.
        stop: Time<UTC>,
    },
}

impl ThrustArc {
    /// Construct a thrust arc, normalising the direction vector and
    /// rejecting zero / non-finite directions, non-positive Isp, and
    /// empty intervals.
    pub fn with_normalised_direction(
        start: Time<UTC>,
        stop: Time<UTC>,
        thrust: Newtons,
        isp: Second,
        direction: [f64; 3],
    ) -> Result<Self, ManeuverError> {
        let isp_val = isp.value();
        if !isp_val.is_finite() || isp_val <= 0.0 {
            return Err(ManeuverError::NonPositiveIsp(isp_val));
        }
        // `!(stop > start)` deliberately rejects NaN-bearing instants too;
        // a literal `stop <= start` would silently accept them.
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(stop > start) {
            return Err(ManeuverError::EmptyInterval { start, stop });
        }
        let [x, y, z] = direction;
        if !x.is_finite() || !y.is_finite() || !z.is_finite() {
            return Err(ManeuverError::InvalidDirection(direction));
        }
        let mag = (x * x + y * y + z * z).sqrt();
        // Reject zero, negative, and NaN magnitudes in one expression.
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(mag > 0.0) {
            return Err(ManeuverError::InvalidDirection(direction));
        }
        Ok(Self {
            start,
            stop,
            thrust,
            isp,
            direction: [x / mag, y / mag, z / mag],
        })
    }

    /// True iff `epoch ∈ [start, stop)`.
    pub fn is_active(&self, epoch: Time<UTC>) -> bool {
        epoch >= self.start && epoch < self.stop
    }
}

/// Thrust acceleration in m / s², or `[0; 3]` outside the active interval.
///
/// Mass must be strictly positive; otherwise [`ManeuverError::NonPositiveMass`]
/// is returned.
pub fn thrust_acceleration(
    arc: &ThrustArc,
    spacecraft_mass: Kilogram,
    epoch: Time<UTC>,
) -> Result<[f64; 3], ManeuverError> {
    let mass_val = spacecraft_mass.value();
    if !mass_val.is_finite() || mass_val <= 0.0 {
        return Err(ManeuverError::NonPositiveMass(mass_val));
    }
    if !arc.is_active(epoch) {
        return Ok([0.0, 0.0, 0.0]);
    }
    let f_n = arc.thrust.value();
    let a = f_n / mass_val;
    let [ux, uy, uz] = arc.direction;
    Ok([a * ux, a * uy, a * uz])
}

/// Mass-flow rate `ṁ = F / (Isp · g₀)` in kg / s.
///
/// Returns [`ManeuverError::NonPositiveIsp`] for a non-positive Isp.
pub fn mass_flow_rate(thrust: Newtons, isp: Second) -> Result<f64, ManeuverError> {
    let isp_val = isp.value();
    if !isp_val.is_finite() || isp_val <= 0.0 {
        return Err(ManeuverError::NonPositiveIsp(isp_val));
    }
    Ok(thrust.value() / (isp_val * G0_M_PER_S2))
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::{TimeZone, Utc};

    fn epoch(secs_from_j2000: i64) -> Time<UTC> {
        let dt = Utc
            .timestamp_opt(946_728_000 + secs_from_j2000, 0)
            .single()
            .unwrap();
        Time::<UTC>::try_from_chrono(dt).unwrap()
    }

    fn arc() -> ThrustArc {
        ThrustArc::with_normalised_direction(
            epoch(0),
            epoch(60),
            Newtons::new(0.5),
            Second::new(1500.0),
            [1.0, 0.0, 0.0],
        )
        .unwrap()
    }

    #[test]
    fn zero_thrust_outside_interval() {
        let a = thrust_acceleration(&arc(), Kilogram::new(100.0), epoch(120)).unwrap();
        assert_eq!(a, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn acceleration_is_force_over_mass_in_direction() {
        let a = thrust_acceleration(&arc(), Kilogram::new(100.0), epoch(10)).unwrap();
        // 0.5 N / 100 kg = 0.005 m/s² along +X
        assert!((a[0] - 0.005).abs() < 1e-15);
        assert_eq!(a[1], 0.0);
        assert_eq!(a[2], 0.0);
    }

    #[test]
    fn mass_flow_rate_uses_g0() {
        let mdot = mass_flow_rate(Newtons::new(0.5), Second::new(1500.0)).unwrap();
        let expected = 0.5 / (1500.0 * G0_M_PER_S2);
        assert!(
            (mdot - expected).abs() < 1e-18,
            "mdot={mdot} expected={expected}"
        );
    }

    #[test]
    fn rejects_zero_direction() {
        let r = ThrustArc::with_normalised_direction(
            epoch(0),
            epoch(10),
            Newtons::new(1.0),
            Second::new(300.0),
            [0.0, 0.0, 0.0],
        );
        assert!(matches!(r, Err(ManeuverError::InvalidDirection(_))));
    }

    #[test]
    fn rejects_non_positive_isp() {
        let r = mass_flow_rate(Newtons::new(1.0), Second::new(0.0));
        assert!(matches!(r, Err(ManeuverError::NonPositiveIsp(_))));
    }

    #[test]
    fn rejects_non_positive_mass() {
        let r = thrust_acceleration(&arc(), Kilogram::new(0.0), epoch(5));
        assert!(matches!(r, Err(ManeuverError::NonPositiveMass(_))));
    }

    #[test]
    fn direction_is_normalised() {
        let a = ThrustArc::with_normalised_direction(
            epoch(0),
            epoch(1),
            Newtons::new(1.0),
            Second::new(300.0),
            [3.0, 0.0, 4.0],
        )
        .unwrap();
        let mag = (a.direction[0].powi(2) + a.direction[1].powi(2) + a.direction[2].powi(2)).sqrt();
        assert!((mag - 1.0).abs() < 1e-15);
        assert!((a.direction[0] - 0.6).abs() < 1e-15);
        assert!((a.direction[2] - 0.8).abs() < 1e-15);
    }
}
