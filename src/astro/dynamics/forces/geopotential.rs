// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Spherical-harmonic geopotential acceleration model.

use principia::{spherical_harmonic_acceleration, AccelerationModel, PrincipiaError};

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::astro::era::earth_rotation_angle;
use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::GCRS;
use crate::time::{JD, TT};

#[derive(Debug, Clone, Copy)]
pub struct Geopotential {
    pub degree: usize,
    pub order: usize,
}

impl Geopotential {
    pub fn new(degree: usize, order: usize) -> Self {
        Self { degree, order }
    }

    pub fn full(degree: usize) -> Self {
        Self {
            degree,
            order: degree,
        }
    }
}

impl AccelerationModel<DynamicsContext, TT, Geocentric, GCRS> for Geopotential {
    fn name(&self) -> &'static str {
        "geopotential"
    }

    fn acceleration(
        &self,
        s: &OrbitState,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, PrincipiaError> {
        let provider = ctx
            .require_gravity_field()
            .map_err(DynamicsError::into_principia)?;
        let (pos_kernel, era_sin_cos) = if self.order > 0 {
            let era = earth_rotation_angle(s.epoch.to::<JD>());
            let (sin_t, cos_t) = era.value().sin_cos();
            let x = s.position.x().value();
            let y = s.position.y().value();
            let z = s.position.z().value();
            (
                [cos_t * x + sin_t * y, -sin_t * x + cos_t * y, z],
                Some((sin_t, cos_t)),
            )
        } else {
            (
                [
                    s.position.x().value(),
                    s.position.y().value(),
                    s.position.z().value(),
                ],
                None,
            )
        };
        let a_kernel = spherical_harmonic_acceleration(
            provider.as_ref(),
            pos_kernel,
            self.degree,
            self.order,
        )?;
        let [ax, ay, az] = if let Some((sin_t, cos_t)) = era_sin_cos {
            [
                cos_t * a_kernel[0] - sin_t * a_kernel[1],
                sin_t * a_kernel[0] + cos_t * a_kernel[1],
                a_kernel[2],
            ]
        } else {
            a_kernel
        };
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(ax, ay, az))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::forces::{TwoBody, EARTH_J2, GM_EARTH, J2, R_EARTH};
    use crate::astro::dynamics::gravity::LowDegreeEarth;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::time::JulianDate;

    fn make_state(x: f64, y: f64, z: f64) -> OrbitState {
        OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::<GCRS>::new(x, y, z),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    fn ctx_with_gravity() -> DynamicsContext {
        DynamicsContext::builder()
            .with_gravity(Arc::new(LowDegreeEarth))
            .build()
    }

    #[test]
    fn geopotential_n2m0_matches_two_body_plus_j2() {
        let ctx = ctx_with_gravity();
        let geopot = Geopotential::new(2, 0);
        let tb = TwoBody::new(GM_EARTH);
        let j2 = J2::new(GM_EARTH, R_EARTH, EARTH_J2);
        let s = make_state(7_000.0, 0.0, 0.0);
        let a_geo = geopot.acceleration(&s, &ctx).unwrap();
        let a_tb = tb.acceleration(&s, &DynamicsContext::empty()).unwrap();
        let a_j2 = j2.acceleration(&s, &DynamicsContext::empty()).unwrap();
        let rel = ((a_geo.x().value() - (a_tb.x().value() + a_j2.x().value())).abs())
            / a_geo.x().value().abs();
        assert!(rel < 5e-9);
    }

    #[test]
    fn geopotential_no_provider_returns_error() {
        let err = Geopotential::new(2, 2)
            .acceleration(&make_state(7000.0, 0.0, 0.0), &DynamicsContext::empty());
        assert!(matches!(
            err,
            Err(PrincipiaError::ContextDataUnavailable {
                what: "gravity field"
            })
        ));
    }
}
