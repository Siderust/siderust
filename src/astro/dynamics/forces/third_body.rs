// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Third-body perturbations from Sun and Moon — `ThirdBodySunMoon` force model.
//!
//! ## Scope
//!
//! Provides [`ThirdBodySunMoon`], fetching Sun and Moon geocentric positions
//! from the [`DynEphemeris`] in the [`DynamicsContext`].
//!
//! ## Equations
//!
//! The Battin formulation for third-body acceleration (Vallado §8.4):
//!
//! ```text
//! a = μ_b · [ (d − r) / |d − r|³ − d / |d|³ ]
//! ```
//!
//! where `d` is the geocentric position of the third body (km) and `r` is the
//! satellite geocentric position (km).
//!
//! ## Units
//!
//! km, km/s, km/s².
//!
//! ## Frame/center assumptions
//!
//! `<Geocentric, GCRS>`.  Body positions are transformed from ecliptic to GCRS
//! using [`ecliptic_of_date_to_mean_equatorial_matrix`] at J2000.
//!
//! ## References
//!
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, §8.4.
//! * Montenbruck & Gill, *Satellite Orbits*, §3.2.

use std::sync::Arc;

use affn::cartesian::Displacement;

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::astro::precession::ecliptic_of_date_to_mean_equatorial_matrix;
use crate::calculus::ephemeris::DynEphemeris;
use crate::coordinates::frames::{EclipticMeanJ2000, GCRS};
use crate::qtty::{AstronomicalUnit, Kilometer};
use crate::time::JulianDate;

use super::traits::{ForceModel, MU_MOON, MU_SUN};

// ---- Private helpers --------------------------------------------------------

/// Compute the geocentric Sun displacement (GCRS, km) from an ephemeris.
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
    let rot = ecliptic_of_date_to_mean_equatorial_matrix(JulianDate::J2000);
    let d_gcrs: Displacement<GCRS, AstronomicalUnit> = rot.apply_vec(d_ecl);
    Ok(d_gcrs.to_unit::<Kilometer>())
}

/// Compute the geocentric Moon displacement (GCRS, km) from an ephemeris.
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
    let rot = ecliptic_of_date_to_mean_equatorial_matrix(JulianDate::J2000);
    Ok(rot.apply_vec(m_ecl))
}

// ---- ThirdBodySunMoon -------------------------------------------------------

/// Third-body point-mass perturbations from Sun and Moon.
///
/// Body positions are fetched from the [`DynEphemeris`] provider in the
/// [`DynamicsContext`].  The Battin formula gives the acceleration:
///
/// ```text
/// a = μ_b · ( (d − r) / |d − r|³ − d / |d|³ )
/// ```
///
/// where `d` is the geocentric position of the third body (km) and `r` is the
/// satellite geocentric position (km).
#[derive(Debug, Clone, Copy, Default)]
pub struct ThirdBodySunMoon;

impl ForceModel for ThirdBodySunMoon {
    #[inline]
    fn acceleration(
        &self,
        s: &OrbitState,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        let eph = ctx.require_ephemeris()?;
        let jd = s.epoch_jd();
        type BodyAccel = (f64, Result<Displacement<GCRS, Kilometer>, DynamicsError>);
        let bodies: [BodyAccel; 2] = [
            (MU_SUN.value(), sun_geocentric(eph, jd)),
            (MU_MOON.value(), moon_geocentric(eph, jd)),
        ];
        let mut ax = 0.0_f64;
        let mut ay = 0.0_f64;
        let mut az = 0.0_f64;
        for (mu, d_res) in bodies {
            let d = d_res?;
            let drx = d.x().value() - s.position.x().value();
            let dry = d.y().value() - s.position.y().value();
            let drz = d.z().value() - s.position.z().value();
            let dr_n = (drx * drx + dry * dry + drz * drz).sqrt();
            let d_n = d.magnitude().value();
            if dr_n == 0.0 || d_n == 0.0 {
                continue;
            }
            let dr3 = dr_n * dr_n * dr_n;
            let d3 = d_n * d_n * d_n;
            ax += mu * (drx / dr3 - d.x().value() / d3);
            ay += mu * (dry / dr3 - d.y().value() / d3);
            az += mu * (drz / dr3 - d.z().value() / d3);
        }
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(ax, ay, az))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::astro::dynamics::context::DynamicsContextBuilder;
    use crate::astro::dynamics::state::OrbitState;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::calculus::ephemeris::Vsop87Ephemeris;
    use crate::coordinates::frames::GCRS;
    use crate::time::JulianDate;

    fn leo_at(epoch: JulianDate) -> OrbitState {
        OrbitState::new_at_jd(
            epoch,
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn third_body_nonzero_acceleration() {
        let ctx = DynamicsContextBuilder::new()
            .with_ephemeris(Arc::new(Vsop87Ephemeris))
            .build();
        let s = leo_at(JulianDate::new(2_451_545.0));
        let a = ThirdBodySunMoon.acceleration(&s, &ctx).unwrap();
        let mag = a.magnitude().value();
        assert!(mag > 0.0, "third-body acceleration should be nonzero");
        assert!(
            mag < 1e-4,
            "third-body acceleration unrealistically large: {mag}"
        );
    }
}
