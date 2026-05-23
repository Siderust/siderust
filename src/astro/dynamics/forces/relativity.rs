// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! First-order post-Newtonian Schwarzschild correction.

use principia::{AccelerationModel, PrincipiaError};

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::GCRS;
use crate::qtty::{KmPerSeconds, SPEED_OF_LIGHT_KM_S};
use crate::time::TT;

use super::{GravitationalParameter, DEGENERATE_RADIUS_KM, GM_EARTH};

/// First-order post-Newtonian (Schwarzschild) relativistic correction for a central body.
///
/// Models the GR perihelion/perigee precession term. The acceleration is:
/// `a = (GM / c²r³) · [(4GM/r − v²) r + 4(r·v) v]`
#[derive(Debug, Clone, Copy)]
pub struct CentralBodyRelativity1Pn {
    /// Standard gravitational parameter µ = GM (km³ s⁻²).
    pub gm: GravitationalParameter,
    /// Speed of light (km s⁻¹).
    pub c: KmPerSeconds,
}

impl CentralBodyRelativity1Pn {
    /// Pre-configured for the Earth–satellite system (GM_EARTH, c = 299792.458 km/s).
    pub fn earth() -> Self {
        Self {
            gm: GM_EARTH,
            c: SPEED_OF_LIGHT_KM_S,
        }
    }
}

impl AccelerationModel<DynamicsContext, TT, Geocentric, GCRS> for CentralBodyRelativity1Pn {
    fn name(&self) -> &'static str {
        "relativity_1pn"
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
        let r2 = rx * rx + ry * ry + rz * rz;
        let r = r2.sqrt();
        if r < DEGENERATE_RADIUS_KM {
            return Err(PrincipiaError::DegenerateGeometry {
                reason: "CentralBodyRelativity1Pn: radius near zero",
            });
        }
        let r3 = r2 * r;
        let v2 = vx * vx + vy * vy + vz * vz;
        let rdotv = rx * vx + ry * vy + rz * vz;
        let gm = self.gm.value();
        let c2 = self.c.value().powi(2);
        let pre = gm / (c2 * r3);
        let coeff_r = 4.0 * gm / r - v2;
        let coeff_v = 4.0 * rdotv;
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            pre * (coeff_r * rx + coeff_v * vx),
            pre * (coeff_r * ry + coeff_v * vy),
            pre * (coeff_r * rz + coeff_v * vz),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::time::JulianDate;

    fn leo_circular() -> OrbitState {
        OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn relativity_magnitude_order_of_magnitude() {
        let mag = CentralBodyRelativity1Pn::earth()
            .acceleration(&leo_circular(), &DynamicsContext::empty())
            .unwrap()
            .magnitude()
            .value();
        assert!(mag > 1e-12 && mag < 1e-9);
    }
}
