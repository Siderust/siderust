// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Constant empirical acceleration in the RTN frame.

use principia::{AccelerationModel, PrincipiaError};

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::GCRS;
use crate::qtty::KmPerSecondsSquared;
use crate::time::TT;

/// Constant empirical acceleration in the Radial-Transverse-Normal (RTN) frame.
#[derive(Debug, Clone, Copy)]
pub struct EmpiricalAcceleration {
    /// Radial (along-position-vector) component of the empirical acceleration (km s⁻²).
    pub radial: KmPerSecondsSquared,
    /// Transverse (along-velocity-vector projected) component (km s⁻²).
    pub transverse: KmPerSecondsSquared,
    /// Normal (orbit-plane perpendicular) component (km s⁻²).
    pub normal: KmPerSecondsSquared,
}

impl EmpiricalAcceleration {
    /// Construct an empirical acceleration from its Radial, Transverse, Normal components.
    pub fn rtn(
        radial: KmPerSecondsSquared,
        transverse: KmPerSecondsSquared,
        normal: KmPerSecondsSquared,
    ) -> Self {
        Self {
            radial,
            transverse,
            normal,
        }
    }
}

impl AccelerationModel<DynamicsContext, TT, Geocentric, GCRS> for EmpiricalAcceleration {
    fn name(&self) -> &'static str {
        "empirical_rtn"
    }

    fn acceleration(
        &self,
        s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, PrincipiaError> {
        let rx = s.position.x().value();
        let ry = s.position.y().value();
        let rz = s.position.z().value();
        let vx = s.velocity.x().value();
        let vy = s.velocity.y().value();
        let vz = s.velocity.z().value();
        let r_norm = (rx * rx + ry * ry + rz * rz).sqrt();
        if r_norm == 0.0 {
            return Err(PrincipiaError::DegenerateGeometry {
                reason: "zero position magnitude — RTN frame undefined",
            });
        }
        let n_raw = [ry * vz - rz * vy, rz * vx - rx * vz, rx * vy - ry * vx];
        let n_norm = (n_raw[0] * n_raw[0] + n_raw[1] * n_raw[1] + n_raw[2] * n_raw[2]).sqrt();
        if n_norm == 0.0 {
            return Err(PrincipiaError::DegenerateGeometry {
                reason: "position and velocity are parallel — RTN frame undefined",
            });
        }
        let r_hat = [rx / r_norm, ry / r_norm, rz / r_norm];
        let n_hat = [n_raw[0] / n_norm, n_raw[1] / n_norm, n_raw[2] / n_norm];
        let t_hat = [
            n_hat[1] * r_hat[2] - n_hat[2] * r_hat[1],
            n_hat[2] * r_hat[0] - n_hat[0] * r_hat[2],
            n_hat[0] * r_hat[1] - n_hat[1] * r_hat[0],
        ];
        let a_rtn = [
            self.radial.value(),
            self.transverse.value(),
            self.normal.value(),
        ];
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            r_hat[0] * a_rtn[0] + t_hat[0] * a_rtn[1] + n_hat[0] * a_rtn[2],
            r_hat[1] * a_rtn[0] + t_hat[1] * a_rtn[1] + n_hat[1] * a_rtn[2],
            r_hat[2] * a_rtn[0] + t_hat[2] * a_rtn[1] + n_hat[2] * a_rtn[2],
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::time::JulianDate;

    fn leo() -> OrbitState {
        OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn pure_radial_empirical_is_collinear_with_position() {
        let a = EmpiricalAcceleration::rtn(
            KmPerSecondsSquared::new(1e-9),
            KmPerSecondsSquared::new(0.0),
            KmPerSecondsSquared::new(0.0),
        )
        .acceleration(&leo(), &DynamicsContext::empty())
        .unwrap();
        assert!(a.y().value().abs() < 1e-15);
        assert!(a.z().value().abs() < 1e-15);
    }
}
