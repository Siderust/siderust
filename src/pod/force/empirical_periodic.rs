// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Periodic empirical accelerations (1-CPR, 2-CPR) in the RTN frame.
//!
//! Upstream `siderust` ships only [`crate::astro::dynamics::forces::EmpiricalAcceleration`],
//! a *constant* RTN acceleration. Sequential-filter POD typically also
//! estimates once-per-revolution (1-CPR) and twice-per-revolution (2-CPR)
//! sinusoidal coefficients. Those are POD-specific composition; they live
//! here, not in the generic `siderust-dynamics` crate.
//!
//! The phase of the periodic term is driven by an externally-supplied
//! orbital period and a reference epoch, so the model is force-model-pure
//! (no implicit two-body propagation inside).
//!
//! ## Equations
//!
//! Define the phase angle
//! ```text
//! θ(t) = n · ω · (t - t_ref),     ω = 2π / T_orbit
//! ```
//! where `n ∈ {1, 2}` selects the harmonic.
//!
//! ```text
//! a_R(t) = R_cos · cos θ + R_sin · sin θ
//! a_T(t) = T_cos · cos θ + T_sin · sin θ
//! a_N(t) = N_cos · cos θ + N_sin · sin θ
//! a_GCRS(t) = R_{GCRS←RTN}(state) · [a_R, a_T, a_N]
//! ```

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::forces::AccelerationModel;
use crate::astro::dynamics::frames::{LocalOrbitalFrame, RTN};
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::GCRS;
use crate::time::JulianDate;
use principia::PrincipiaError;
use qtty::{KmPerSecondsSquared, Second};
use tempoch::{JD, TT};

/// Harmonic order of the periodic empirical acceleration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PeriodicHarmonic {
    /// Once-per-revolution.
    OncePerRev,
    /// Twice-per-revolution.
    TwicePerRev,
}

impl PeriodicHarmonic {
    /// Numerical multiplier (1 or 2).
    pub fn multiplier(self) -> f64 {
        match self {
            PeriodicHarmonic::OncePerRev => 1.0,
            PeriodicHarmonic::TwicePerRev => 2.0,
        }
    }
}

/// Periodic empirical acceleration in the RTN frame, parameterised by the
/// orbital period and a reference epoch.
///
/// # Example
///
/// ```
/// use siderust::pod::force::{EmpiricalPeriodicAcceleration, PeriodicHarmonic};
/// use qtty::{KmPerSecondsSquared, Second};
/// use siderust::time::JulianDate;
///
/// let f = EmpiricalPeriodicAcceleration::new(
///     PeriodicHarmonic::OncePerRev,
///     JulianDate::new(2_451_545.0),
///     Second::new(5400.0),
///     KmPerSecondsSquared::new(0.0), KmPerSecondsSquared::new(0.0),
///     KmPerSecondsSquared::new(1e-12), KmPerSecondsSquared::new(0.0),
///     KmPerSecondsSquared::new(0.0), KmPerSecondsSquared::new(0.0),
/// );
/// assert_eq!(f.harmonic, PeriodicHarmonic::OncePerRev);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct EmpiricalPeriodicAcceleration {
    /// Harmonic of the orbit period.
    pub harmonic: PeriodicHarmonic,
    /// Reference epoch defining `θ = 0`.
    pub epoch_ref: JulianDate,
    /// Orbital period `T_orbit` used to compute `ω = 2π / T_orbit`.
    pub period: Second,
    /// Radial cosine coefficient.
    pub r_cos: KmPerSecondsSquared,
    /// Radial sine coefficient.
    pub r_sin: KmPerSecondsSquared,
    /// Transverse cosine coefficient.
    pub t_cos: KmPerSecondsSquared,
    /// Transverse sine coefficient.
    pub t_sin: KmPerSecondsSquared,
    /// Normal cosine coefficient.
    pub n_cos: KmPerSecondsSquared,
    /// Normal sine coefficient.
    pub n_sin: KmPerSecondsSquared,
}

impl EmpiricalPeriodicAcceleration {
    /// Construct a periodic empirical acceleration from typed RTN coefficients.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        harmonic: PeriodicHarmonic,
        epoch_ref: JulianDate,
        period: Second,
        r_cos: KmPerSecondsSquared,
        r_sin: KmPerSecondsSquared,
        t_cos: KmPerSecondsSquared,
        t_sin: KmPerSecondsSquared,
        n_cos: KmPerSecondsSquared,
        n_sin: KmPerSecondsSquared,
    ) -> Self {
        Self {
            harmonic,
            epoch_ref,
            period,
            r_cos,
            r_sin,
            t_cos,
            t_sin,
            n_cos,
            n_sin,
        }
    }

    fn phase(&self, epoch: JulianDate) -> f64 {
        let dt_days = epoch.value() - self.epoch_ref.value();
        let dt_s = dt_days * 86_400.0;
        let t = self.period.value();
        if t <= 0.0 || !t.is_finite() {
            return 0.0;
        }
        let omega = 2.0 * core::f64::consts::PI / t;
        self.harmonic.multiplier() * omega * dt_s
    }
}

impl AccelerationModel<DynamicsContext, TT, Geocentric, GCRS> for EmpiricalPeriodicAcceleration {
    fn name(&self) -> &'static str {
        "empirical_periodic"
    }

    fn acceleration(
        &self,
        s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, PrincipiaError> {
        let theta = self.phase(s.epoch.to::<JD>());
        let (sin_t, cos_t) = theta.sin_cos();
        let a_rtn = [
            self.r_cos.value() * cos_t + self.r_sin.value() * sin_t,
            self.t_cos.value() * cos_t + self.t_sin.value() * sin_t,
            self.n_cos.value() * cos_t + self.n_sin.value() * sin_t,
        ];
        let rtn_frame = LocalOrbitalFrame::<RTN>::try_from_state(s)?;
        let a_gcrs = rtn_frame.rotation_inverse().apply_array(a_rtn);
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            a_gcrs[0], a_gcrs[1], a_gcrs[2],
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;

    fn s0() -> OrbitState {
        OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::<GCRS>::new(7_000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5450, 0.0),
        )
    }

    #[test]
    fn one_cpr_at_reference_epoch_returns_cosine_amplitude() {
        let f = EmpiricalPeriodicAcceleration::new(
            PeriodicHarmonic::OncePerRev,
            JulianDate::new(2_451_545.0),
            Second::new(5400.0),
            KmPerSecondsSquared::new(1e-9),
            KmPerSecondsSquared::new(0.0),
            KmPerSecondsSquared::new(0.0),
            KmPerSecondsSquared::new(0.0),
            KmPerSecondsSquared::new(0.0),
            KmPerSecondsSquared::new(0.0),
        );
        let a = f.acceleration(&s0(), &DynamicsContext::empty()).unwrap();
        // r̂ at epoch is +X in GCRS, so a_GCRS_x ≈ 1e-9 km/s².
        assert!((a.x().value() - 1e-9).abs() < 1e-15);
    }

    #[test]
    fn two_cpr_doubles_phase() {
        let f1 = EmpiricalPeriodicAcceleration::new(
            PeriodicHarmonic::OncePerRev,
            JulianDate::new(2_451_545.0),
            Second::new(5400.0),
            KmPerSecondsSquared::new(0.0),
            KmPerSecondsSquared::new(0.0),
            KmPerSecondsSquared::new(1e-9),
            KmPerSecondsSquared::new(0.0),
            KmPerSecondsSquared::new(0.0),
            KmPerSecondsSquared::new(0.0),
        );
        let mut f2 = f1;
        f2.harmonic = PeriodicHarmonic::TwicePerRev;
        // After T_orbit/4 the 1cpr has cos(π/2)=0, the 2cpr has cos(π)=-1.
        let dt_quarter = 5400.0 / 4.0;
        let epoch = JulianDate::new(2_451_545.0 + dt_quarter / 86_400.0);
        let mut s = s0();
        s.epoch = epoch.into();
        let a1 = f1.acceleration(&s, &DynamicsContext::empty()).unwrap();
        let a2 = f2.acceleration(&s, &DynamicsContext::empty()).unwrap();
        let mag1 = (a1.x().value().powi(2) + a1.y().value().powi(2)).sqrt();
        let mag2 = (a2.x().value().powi(2) + a2.y().value().powi(2)).sqrt();
        assert!(mag1 < 1e-13, "1cpr mag at T/4 should be ≈0, got {mag1:e}");
        assert!(
            (mag2 - 1e-9).abs() < 1e-12,
            "2cpr mag at T/4 should be ≈1e-9, got {mag2:e}"
        );
    }
}
