// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Astronomy-specific perturbation models plus `principia`'s generic accelerations.

pub mod drag;
pub mod empirical;
pub mod geopotential;
pub mod relativity;
pub mod srp;
pub mod third_body;

use crate::qtty::{AstronomicalUnit, InverseSeconds, Kilometer, Kilometers, Pascals, Unit};

pub use crate::astro::dynamics::density::{ExponentialAtmosphere, Nrlmsise00LiteApprox};
pub use crate::astro::dynamics::units::{GravitationalParameter, GM_EARTH, GM_MOON, GM_SUN};
pub use drag::{DragForce, ExponentialDrag};
pub use empirical::EmpiricalAcceleration;
pub use geopotential::Geopotential;
pub use principia::{AccelerationModel, AccelerationPartials, CompositeModel, TwoBody, J2};
pub use relativity::CentralBodyRelativity1Pn;
pub use srp::{CannonballSrp, Conical, Cylindrical, EclipseModel, NoEclipse, ShadowModel};
pub use third_body::{MoonPerturbation, SunPerturbation, ThirdBody, ThirdBodyProvider};

/// Earth mean equatorial radius (GRS-80 / WGS-84), km.
pub const R_EARTH: Kilometers = Kilometers::new(6_378.137);
/// Earth mean rotation rate (sidereal), rad/s.
pub const OMEGA_EARTH_RAD_S: InverseSeconds = InverseSeconds::new(7.292_115_146_706_979e-5);
/// Minimum geocentric radius below which force-model computations are considered degenerate.
pub const DEGENERATE_RADIUS_KM: f64 = 100.0;
/// Solar radiation pressure at 1 AU, N/m².
pub const P0: Pascals = Pascals::new(4.560e-6);
/// Canonical Earth J2 coefficient used by the built-in low-degree field.
pub const EARTH_J2: f64 = 1.082_626_68e-3;
/// Astronomical unit in km.
pub const AU_IN_KM: f64 = AstronomicalUnit::RATIO / Kilometer::RATIO;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::{OrbitState, Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::time::JulianDate;

    fn leo() -> OrbitState {
        OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn composite_model_sums_components() {
        let ctx = DynamicsContext::empty();
        let model = CompositeModel::empty()
            .push(Box::new(TwoBody::new(GM_EARTH)))
            .push(Box::new(J2::new(GM_EARTH, R_EARTH, EARTH_J2)));
        let a = model.acceleration(&leo(), &ctx).unwrap();
        assert!(a.x().value() < 0.0);
    }
}
