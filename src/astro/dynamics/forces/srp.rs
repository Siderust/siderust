// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Cannonball solar radiation pressure acceleration.

use std::marker::PhantomData;

use affn::cartesian::Displacement;
use principia::{AccelerationModel, PrincipiaError};

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::GCRS;
use crate::qtty::{AreaToMass, Kilometer, Ratios, SrpCoefficient};
use crate::time::{JD, TT};

use super::third_body::sun_geocentric;
use super::{AU_IN_KM, P0, R_EARTH};

const R_SUN_KM: f64 = 695_700.0;

/// Shadow model to apply when computing the SRP eclipse fraction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ShadowModel {
    /// No shadowing — satellite is always in sunlight.
    None,
    /// Cylindrical shadow model (sharp umbra, no penumbra).
    Cylindrical,
    /// Conical (dual-cone) model with partial penumbra (default).
    #[default]
    Conical,
}

mod sealed {
    pub trait Sealed {}
}

/// Object-safe trait for eclipse models used by [`CannonballSrp`].
///
/// Implementors compute the fraction of solar flux seen by the satellite
/// (1.0 = full sunlight, 0.0 = full shadow).
pub trait EclipseModel: sealed::Sealed + Send + Sync + 'static {
    /// Return the shadow factor ν ∈ [0, 1] for the given satellite and Sun vectors.
    fn eclipse_factor(
        r_sat: Displacement<GCRS, Kilometer>,
        r_sun: Displacement<GCRS, Kilometer>,
    ) -> Ratios;
}

/// Eclipse model that always returns full illumination (ν = 1).
#[derive(Debug, Clone, Copy, Default)]
pub struct NoEclipse;
/// Cylindrical (hard) shadow model — sharp umbra boundary.
#[derive(Debug, Clone, Copy, Default)]
pub struct Cylindrical;
/// Dual-cone (penumbra + umbra) shadow model.
#[derive(Debug, Clone, Copy, Default)]
pub struct Conical;

impl sealed::Sealed for NoEclipse {}
impl sealed::Sealed for Cylindrical {}
impl sealed::Sealed for Conical {}

impl EclipseModel for NoEclipse {
    fn eclipse_factor(
        _r_sat: Displacement<GCRS, Kilometer>,
        _r_sun: Displacement<GCRS, Kilometer>,
    ) -> Ratios {
        Ratios::new(1.0)
    }
}
impl EclipseModel for Cylindrical {
    fn eclipse_factor(
        r_sat: Displacement<GCRS, Kilometer>,
        r_sun: Displacement<GCRS, Kilometer>,
    ) -> Ratios {
        let sat = [r_sat.x().value(), r_sat.y().value(), r_sat.z().value()];
        let sun = [r_sun.x().value(), r_sun.y().value(), r_sun.z().value()];
        Ratios::new(cylindrical_shadow_factor(sat, sun, R_EARTH.value()))
    }
}
impl EclipseModel for Conical {
    fn eclipse_factor(
        r_sat: Displacement<GCRS, Kilometer>,
        r_sun: Displacement<GCRS, Kilometer>,
    ) -> Ratios {
        let sat = [r_sat.x().value(), r_sat.y().value(), r_sat.z().value()];
        let sun = [r_sun.x().value(), r_sun.y().value(), r_sun.z().value()];
        Ratios::new(conical_shadow_factor(sat, sun, R_EARTH.value(), R_SUN_KM))
    }
}

/// Cannonball solar radiation pressure model with a configurable eclipse model.
#[derive(Debug, Clone, Copy)]
pub struct CannonballSrp<S: EclipseModel = Conical> {
    /// SRP reflectivity coefficient Cr (typically 1.0–2.0).
    pub cr: SrpCoefficient,
    /// Area-to-mass ratio (m² kg⁻¹).
    pub area_to_mass: AreaToMass,
    _shadow: PhantomData<fn() -> S>,
}

impl<S: EclipseModel> CannonballSrp<S> {
    /// Construct a cannonball SRP model with the given Cr and A/m ratio.
    pub fn new(cr: SrpCoefficient, area_to_mass: AreaToMass) -> Self {
        Self {
            cr,
            area_to_mass,
            _shadow: PhantomData,
        }
    }
}

fn cylindrical_shadow_factor(r_sat: [f64; 3], r_sun: [f64; 3], r_earth_km: f64) -> f64 {
    let r_sun_mag = (r_sun[0] * r_sun[0] + r_sun[1] * r_sun[1] + r_sun[2] * r_sun[2]).sqrt();
    if r_sun_mag == 0.0 {
        return 1.0;
    }
    let sun_hat = [
        r_sun[0] / r_sun_mag,
        r_sun[1] / r_sun_mag,
        r_sun[2] / r_sun_mag,
    ];
    let proj = r_sat[0] * sun_hat[0] + r_sat[1] * sun_hat[1] + r_sat[2] * sun_hat[2];
    if proj >= 0.0 {
        return 1.0;
    }
    let perp_sq = (r_sat[0] - proj * sun_hat[0]).powi(2)
        + (r_sat[1] - proj * sun_hat[1]).powi(2)
        + (r_sat[2] - proj * sun_hat[2]).powi(2);
    if perp_sq < r_earth_km * r_earth_km {
        0.0
    } else {
        1.0
    }
}

fn conical_shadow_factor(
    r_sat: [f64; 3],
    r_sun: [f64; 3],
    r_earth_km: f64,
    r_sun_radius_km: f64,
) -> f64 {
    let d_earth = (r_sat[0] * r_sat[0] + r_sat[1] * r_sat[1] + r_sat[2] * r_sat[2]).sqrt();
    if d_earth == 0.0 {
        return 1.0;
    }
    let sat_to_sun = [
        r_sun[0] - r_sat[0],
        r_sun[1] - r_sat[1],
        r_sun[2] - r_sat[2],
    ];
    let d_sun = (sat_to_sun[0] * sat_to_sun[0]
        + sat_to_sun[1] * sat_to_sun[1]
        + sat_to_sun[2] * sat_to_sun[2])
        .sqrt();
    if d_sun == 0.0 {
        return 1.0;
    }
    let rho_e = (r_earth_km / d_earth).min(1.0).asin();
    let rho_s = (r_sun_radius_km / d_sun).min(1.0).asin();
    let cos_theta =
        (-r_sat[0] * sat_to_sun[0] - r_sat[1] * sat_to_sun[1] - r_sat[2] * sat_to_sun[2])
            / (d_earth * d_sun);
    let theta = cos_theta.clamp(-1.0, 1.0).acos();
    if theta >= rho_s + rho_e {
        return 1.0;
    }
    if rho_e >= rho_s && theta <= rho_e - rho_s {
        return 0.0;
    }
    if rho_s > rho_e && theta <= rho_s - rho_e {
        return 1.0 - (rho_e / rho_s).powi(2);
    }
    let x = (theta * theta + rho_s * rho_s - rho_e * rho_e) / (2.0 * theta);
    let y = (theta * theta + rho_e * rho_e - rho_s * rho_s) / (2.0 * theta);
    let alpha = (x / rho_s).clamp(-1.0, 1.0).acos();
    let beta = (y / rho_e).clamp(-1.0, 1.0).acos();
    let semiperimeter_sq = (rho_s + rho_e + theta)
        * (rho_s + rho_e - theta)
        * (theta + rho_s - rho_e)
        * (theta - rho_s + rho_e);
    let triangle_term = semiperimeter_sq.max(0.0).sqrt() / 2.0;
    let a_overlap = rho_s * rho_s * alpha + rho_e * rho_e * beta - triangle_term;
    let a_sun = std::f64::consts::PI * rho_s * rho_s;
    (1.0 - a_overlap / a_sun).clamp(0.0, 1.0)
}

impl<S: EclipseModel> AccelerationModel<DynamicsContext, TT, Geocentric, GCRS>
    for CannonballSrp<S>
{
    fn name(&self) -> &'static str {
        "srp"
    }

    fn acceleration(
        &self,
        s: &OrbitState<Geocentric, GCRS>,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, PrincipiaError> {
        let eph = ctx
            .require_ephemeris()
            .map_err(DynamicsError::into_principia)?;
        let sun = sun_geocentric(eph, s.epoch.to::<JD>()).map_err(DynamicsError::into_principia)?;
        let r_sun_sat = Displacement::<GCRS, Kilometer>::new(
            s.position.x().value() - sun.x().value(),
            s.position.y().value() - sun.y().value(),
            s.position.z().value() - sun.z().value(),
        );
        let r = r_sun_sat.magnitude().value();
        if r == 0.0 {
            return Ok(Acceleration::<GCRS, AccelerationUnit>::new(0.0, 0.0, 0.0));
        }
        let r_sat_disp =
            Displacement::<GCRS, Kilometer>::new(s.position.x(), s.position.y(), s.position.z());
        let nu = S::eclipse_factor(r_sat_disp, sun).value();
        let mag_km_s2 = nu
            * self.cr.value()
            * P0.value()
            * (AU_IN_KM * AU_IN_KM / (r * r))
            * self.area_to_mass.value()
            / 1_000.0;
        let inv_r = 1.0 / r;
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            mag_km_s2 * r_sun_sat.x().value() * inv_r,
            mag_km_s2 * r_sun_sat.y().value() * inv_r,
            mag_km_s2 * r_sun_sat.z().value() * inv_r,
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::astro::dynamics::context::DynamicsContextBuilder;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::{
        cartesian::{Position as CartPosition, Velocity as CartVelocity},
        centers::{Barycentric, Geocentric as GeoCenter, Heliocentric},
        frames::{EclipticMeanJ2000, GCRS},
    };
    use crate::ephemeris::{AuPerDay, DynEphemeris, EphemerisError};
    use crate::qtty::{AstronomicalUnit, Kilometer};
    use crate::time::JulianDate;

    struct FixedSunEphemeris;

    impl DynEphemeris for FixedSunEphemeris {
        fn try_sun_barycentric(
            &self,
            _jd: JulianDate,
        ) -> Result<CartPosition<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError>
        {
            Ok(CartPosition::new(1.0, 0.0, 0.0))
        }
        fn sun_barycentric(
            &self,
            _jd: JulianDate,
        ) -> CartPosition<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
            CartPosition::new(1.0, 0.0, 0.0)
        }
        fn try_earth_barycentric(
            &self,
            _jd: JulianDate,
        ) -> Result<CartPosition<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError>
        {
            Ok(CartPosition::new(0.0, 0.0, 0.0))
        }
        fn earth_barycentric(
            &self,
            _jd: JulianDate,
        ) -> CartPosition<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
            CartPosition::new(0.0, 0.0, 0.0)
        }
        fn try_earth_heliocentric(
            &self,
            _jd: JulianDate,
        ) -> Result<CartPosition<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError>
        {
            Ok(CartPosition::new(0.0, 0.0, 0.0))
        }
        fn earth_heliocentric(
            &self,
            _jd: JulianDate,
        ) -> CartPosition<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> {
            CartPosition::new(0.0, 0.0, 0.0)
        }
        fn try_earth_barycentric_velocity(
            &self,
            _jd: JulianDate,
        ) -> Result<CartVelocity<EclipticMeanJ2000, AuPerDay>, EphemerisError> {
            Ok(CartVelocity::new(0.0, 0.0, 0.0))
        }
        fn earth_barycentric_velocity(
            &self,
            _jd: JulianDate,
        ) -> CartVelocity<EclipticMeanJ2000, AuPerDay> {
            CartVelocity::new(0.0, 0.0, 0.0)
        }
        fn try_moon_geocentric(
            &self,
            _jd: JulianDate,
        ) -> Result<CartPosition<GeoCenter, EclipticMeanJ2000, Kilometer>, EphemerisError> {
            Ok(CartPosition::new(0.0, 0.0, 0.0))
        }
        fn moon_geocentric(
            &self,
            _jd: JulianDate,
        ) -> CartPosition<GeoCenter, EclipticMeanJ2000, Kilometer> {
            CartPosition::new(0.0, 0.0, 0.0)
        }
    }

    fn ctx_stub() -> DynamicsContext {
        DynamicsContextBuilder::new()
            .with_ephemeris(Arc::new(FixedSunEphemeris))
            .build()
    }

    fn leo() -> OrbitState {
        OrbitState::new(
            JulianDate::JD_EPOCH_J2000_0.to_j2000s(),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn srp_order_of_magnitude_at_leo() {
        let srp = CannonballSrp::<NoEclipse>::new(SrpCoefficient::new(1.5), AreaToMass::new(0.01));
        let mag = srp
            .acceleration(&leo(), &ctx_stub())
            .unwrap()
            .magnitude()
            .value();
        assert!((1e-11..5e-10).contains(&mag));
    }

    #[test]
    fn shadow_model_none_ignores_eclipse() {
        let state = OrbitState::new(
            JulianDate::JD_EPOCH_J2000_0.to_j2000s(),
            Position::<GCRS>::new(-(R_EARTH.value() + 500.0), 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        );
        let srp = CannonballSrp::<NoEclipse>::new(SrpCoefficient::new(1.5), AreaToMass::new(0.01));
        assert!(
            srp.acceleration(&state, &ctx_stub())
                .unwrap()
                .magnitude()
                .value()
                > 0.0
        );
    }
}
