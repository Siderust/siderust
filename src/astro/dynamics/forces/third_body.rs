// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Third-body perturbation force model.

use std::sync::Arc;

use affn::cartesian::Displacement;
use principia::{AccelerationModel, PrincipiaError};

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::astro::precession::ecliptic_of_date_to_mean_equatorial_matrix;
use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::{EclipticMeanJ2000, GCRS};
use crate::ephemeris::DynEphemeris;
use crate::qtty::{AstronomicalUnit, Kilometer};
use crate::time::{JulianDate, JD, TT};

use super::{GravitationalParameter, GM_MOON, GM_SUN};

fn battin_third_body_accel(r_s: [f64; 3], r_b: [f64; 3], mu: f64) -> [f64; 3] {
    let r_b_sq = r_b[0] * r_b[0] + r_b[1] * r_b[1] + r_b[2] * r_b[2];
    let r_b_norm = r_b_sq.sqrt();
    if r_b_norm == 0.0 {
        return [0.0; 3];
    }
    let q = (r_s[0] * (r_s[0] - 2.0 * r_b[0])
        + r_s[1] * (r_s[1] - 2.0 * r_b[1])
        + r_s[2] * (r_s[2] - 2.0 * r_b[2]))
        / r_b_sq;
    let fq = q * (3.0 + q * (3.0 + q)) / (1.0 + (1.0 + q).powf(1.5));
    let denom = r_b_norm * r_b_sq * (1.0 + fq);
    if denom == 0.0 {
        return [0.0; 3];
    }
    let scale = -mu / denom;
    [
        scale * (r_b[0] * fq + r_s[0]),
        scale * (r_b[1] * fq + r_s[1]),
        scale * (r_b[2] * fq + r_s[2]),
    ]
}

pub(super) fn sun_geocentric(
    eph: &Arc<dyn DynEphemeris + Send + Sync>,
    jd: JulianDate,
) -> Result<Displacement<GCRS, Kilometer>, DynamicsError> {
    let sun_b = eph
        .try_sun_barycentric(jd)
        .map_err(|e| DynamicsError::EphemerisUnavailable {
            body: "Sun",
            source: Some(Box::new(e)),
        })?;
    let earth_b =
        eph.try_earth_barycentric(jd)
            .map_err(|e| DynamicsError::EphemerisUnavailable {
                body: "Earth",
                source: Some(Box::new(e)),
            })?;
    let d_ecl: Displacement<EclipticMeanJ2000, AstronomicalUnit> = sun_b - earth_b;
    let rot = ecliptic_of_date_to_mean_equatorial_matrix(JulianDate::JD_EPOCH_J2000_0);
    Ok(rot.apply_vec(d_ecl).to_unit::<Kilometer>())
}

pub(super) fn moon_geocentric(
    eph: &Arc<dyn DynEphemeris + Send + Sync>,
    jd: JulianDate,
) -> Result<Displacement<GCRS, Kilometer>, DynamicsError> {
    let m = eph
        .try_moon_geocentric(jd)
        .map_err(|e| DynamicsError::EphemerisUnavailable {
            body: "Moon",
            source: Some(Box::new(e)),
        })?;
    let m_ecl = Displacement::<EclipticMeanJ2000, Kilometer>::new(
        m.x().value(),
        m.y().value(),
        m.z().value(),
    );
    let rot = ecliptic_of_date_to_mean_equatorial_matrix(JulianDate::JD_EPOCH_J2000_0);
    Ok(rot.apply_vec(m_ecl))
}

/// Trait for objects that provide the position and gravitational parameter of a
/// perturbing third body relative to the reference frame origin.
pub trait ThirdBodyProvider: Send + Sync {
    /// Short human-readable name of the perturbing body (e.g. `"Sun"`).
    fn name(&self) -> &'static str;
    /// Standard gravitational parameter µ = GM (km³ s⁻²) of the body.
    fn gm(&self) -> GravitationalParameter;
    /// Geocentric position of the body in GCRS at `epoch`.
    fn position_relative_to_origin(
        &self,
        eph: &Arc<dyn DynEphemeris + Send + Sync>,
        epoch: JulianDate,
    ) -> Result<Displacement<GCRS, Kilometer>, DynamicsError>;
}

/// Third-body perturbation due to solar gravity.
#[derive(Debug, Clone, Copy, Default)]
pub struct SunPerturbation;
impl ThirdBodyProvider for SunPerturbation {
    fn name(&self) -> &'static str {
        "Sun"
    }
    fn gm(&self) -> GravitationalParameter {
        GM_SUN
    }
    fn position_relative_to_origin(
        &self,
        eph: &Arc<dyn DynEphemeris + Send + Sync>,
        epoch: JulianDate,
    ) -> Result<Displacement<GCRS, Kilometer>, DynamicsError> {
        sun_geocentric(eph, epoch)
    }
}

/// Third-body perturbation due to lunar gravity.
#[derive(Debug, Clone, Copy, Default)]
pub struct MoonPerturbation;
impl ThirdBodyProvider for MoonPerturbation {
    fn name(&self) -> &'static str {
        "Moon"
    }
    fn gm(&self) -> GravitationalParameter {
        GM_MOON
    }
    fn position_relative_to_origin(
        &self,
        eph: &Arc<dyn DynEphemeris + Send + Sync>,
        epoch: JulianDate,
    ) -> Result<Displacement<GCRS, Kilometer>, DynamicsError> {
        moon_geocentric(eph, epoch)
    }
}

/// Composite third-body force model: accumulates contributions from any number
/// of [`ThirdBodyProvider`]s using the Battin formulation for numerical stability.
pub struct ThirdBody {
    bodies: Vec<Box<dyn ThirdBodyProvider>>,
}

impl ThirdBody {
    /// Create an empty model with no perturbing bodies.
    pub fn new() -> Self {
        Self { bodies: Vec::new() }
    }
    /// Add solar gravity perturbation.
    pub fn with_sun(mut self) -> Self {
        self.bodies.push(Box::new(SunPerturbation));
        self
    }
    /// Add lunar gravity perturbation.
    pub fn with_moon(mut self) -> Self {
        self.bodies.push(Box::new(MoonPerturbation));
        self
    }
    /// Add an arbitrary third-body provider.
    pub fn with(mut self, body: Box<dyn ThirdBodyProvider>) -> Self {
        self.bodies.push(body);
        self
    }
    /// Convenience constructor: Sun + Moon perturbations.
    pub fn sun_and_moon() -> Self {
        Self::new().with_sun().with_moon()
    }
}

impl Default for ThirdBody {
    fn default() -> Self {
        Self::new()
    }
}

impl AccelerationModel<DynamicsContext, TT, Geocentric, GCRS> for ThirdBody {
    fn name(&self) -> &'static str {
        "third_body"
    }

    fn acceleration(
        &self,
        s: &OrbitState<Geocentric, GCRS>,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, PrincipiaError> {
        if self.bodies.is_empty() {
            return Ok(Acceleration::new(0.0, 0.0, 0.0));
        }
        let eph = ctx
            .require_ephemeris()
            .map_err(DynamicsError::into_principia)?;
        let jd = s.epoch.to::<JD>();
        let r_s = [
            s.position.x().value(),
            s.position.y().value(),
            s.position.z().value(),
        ];
        let mut total = [0.0; 3];
        for body in &self.bodies {
            let r_b = body
                .position_relative_to_origin(eph, jd)
                .map_err(DynamicsError::into_principia)?;
            let a = battin_third_body_accel(
                r_s,
                [r_b.x().value(), r_b.y().value(), r_b.z().value()],
                body.gm().value(),
            );
            total[0] += a[0];
            total[1] += a[1];
            total[2] += a[2];
        }
        Ok(Acceleration::new(total[0], total[1], total[2]))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::astro::dynamics::context::DynamicsContextBuilder;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::ephemeris::Vsop87Ephemeris;

    fn direct_third_body_accel(r_s: [f64; 3], r_b: [f64; 3], mu: f64) -> [f64; 3] {
        let d = [
            r_b[0] - r_s[0],
            r_b[1] - r_s[1],
            r_b[2] - r_s[2],
        ];
        let d3 = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).powf(1.5);
        let rb3 = (r_b[0] * r_b[0] + r_b[1] * r_b[1] + r_b[2] * r_b[2]).powf(1.5);
        [
            mu * (d[0] / d3 - r_b[0] / rb3),
            mu * (d[1] / d3 - r_b[1] / rb3),
            mu * (d[2] / d3 - r_b[2] / rb3),
        ]
    }

    fn leo() -> OrbitState {
        OrbitState::new(
            JulianDate::JD_EPOCH_J2000_0.to_j2000s(),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    fn ctx_with_eph() -> DynamicsContext {
        DynamicsContextBuilder::new()
            .with_ephemeris(Arc::new(Vsop87Ephemeris))
            .build()
    }

    #[test]
    fn empty_body_list_returns_zero_without_ephemeris() {
        let a = ThirdBody::new()
            .acceleration(&leo(), &DynamicsContext::empty())
            .unwrap();
        assert_eq!(a.magnitude().value(), 0.0);
    }

    #[test]
    fn sun_and_moon_nonzero_acceleration() {
        let mag = ThirdBody::sun_and_moon()
            .acceleration(&leo(), &ctx_with_eph())
            .unwrap()
            .magnitude()
            .value();
        assert!(mag > 0.0 && mag < 1e-4);
    }

    /// Verify the Battin formula matches the direct (numerically naive) formula
    /// for well-separated geometries (LEO, MEO, GEO spacecraft vs Moon/Sun).
    #[test]
    fn battin_matches_direct_formula_well_separated() {
        const MU_MOON: f64 = 4902.8;
        const MU_SUN: f64 = 1.327_124_4e11;

        let cases: &[(&str, [f64; 3], [f64; 3], f64)] = &[
            // LEO spacecraft, Moon roughly in +Y
            ("leo_moon", [7_000.0, 0.0, 0.0], [0.0, 384_400.0, 0.0], MU_MOON),
            // MEO spacecraft, Moon at an angle
            ("meo_moon", [20_000.0, 5_000.0, 0.0], [200_000.0, 300_000.0, 50_000.0], MU_MOON),
            // GEO spacecraft, Sun in +X
            ("geo_sun", [42_164.0, 0.0, 0.0], [1.496e8, 0.0, 0.0], MU_SUN),
            // Out-of-plane spacecraft, Moon at angle
            ("inclined_moon", [5_000.0, 3_000.0, 4_000.0], [-150_000.0, 330_000.0, 80_000.0], MU_MOON),
        ];

        for (name, r_s, r_b, mu) in cases {
            let battin = battin_third_body_accel(*r_s, *r_b, *mu);
            let direct = direct_third_body_accel(*r_s, *r_b, *mu);

            let mag = (direct[0] * direct[0] + direct[1] * direct[1] + direct[2] * direct[2]).sqrt();
            assert!(mag > 0.0, "{name}: reference magnitude is zero");

            for i in 0..3 {
                let rel_err = ((battin[i] - direct[i]) / mag).abs();
                assert!(
                    rel_err < 1e-10,
                    "{name}[{i}]: Battin={:.6e}, direct={:.6e}, rel_err={:.3e}",
                    battin[i], direct[i], rel_err
                );
            }
        }
    }
}
