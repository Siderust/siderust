// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Boolean / declarative force-model selection for a POD run.
//!
//! [`ForceModelConfig`] is the high-level "which force families are on?"
//! description that the POD service / CLI consumes. It maps to the
//! lower-level [`siderust::pod::dynamics::registry::ForceModelSpec`] vector via
//! [`ForceModelConfig::to_specs`], which is then materialised by
//! [`siderust::pod::dynamics::registry::ForceModelRegistry`].
//!
//! The split keeps the config layer free of heavy numerical types while
//! still guaranteeing that the registry-level translation is exhaustive
//! and unit-typed.

use crate::astro::dynamics::forces::ShadowModel;
use qtty::{AreaToMass, DragCoefficient, KmPerSecondsSquared, Second, SrpCoefficient};
use crate::time::JulianDate;

use super::empirical_periodic::PeriodicHarmonic;
use super::registry::{ForceModelParams, ForceModelSpec};

/// Declarative description of which force-model contributions are active.
///
/// # Example
///
/// ```
/// use siderust::pod::dynamics::ForceModelConfig;
/// let cfg = ForceModelConfig::default();
/// assert!(cfg.two_body && cfg.j2);
/// let specs = cfg.to_specs();
/// assert_eq!(specs.len(), 2);
/// assert_eq!(specs[0].name, "two_body");
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct ForceModelConfig {
    /// Enable two-body central gravity.
    pub two_body: bool,
    /// Enable J2 oblateness.
    pub j2: bool,
    /// Maximum spherical-harmonics degree/order. `None` disables non-J2 harmonics.
    pub harmonics: Option<(usize, usize)>,
    /// Enable Sun third-body perturbation.
    pub third_body_sun: bool,
    /// Enable Moon third-body perturbation.
    pub third_body_moon: bool,
    /// Optional drag (typed `C_D`, `A/m`).
    pub drag: Option<(DragCoefficient, AreaToMass)>,
    /// Optional cannonball SRP (typed `C_R`, `A/m`, shadow model).
    pub srp: Option<(SrpCoefficient, AreaToMass, ShadowModel)>,
    /// Enable central-body Schwarzschild relativistic correction.
    pub relativity: bool,
    /// Optional constant RTN empirical acceleration.
    pub empirical_constant: Option<(
        KmPerSecondsSquared,
        KmPerSecondsSquared,
        KmPerSecondsSquared,
    )>,
    /// Optional 1-CPR empirical acceleration:
    /// `(epoch_ref, period, [r_cos, r_sin, t_cos, t_sin, n_cos, n_sin])`.
    pub empirical_1cpr: Option<(JulianDate, Second, [KmPerSecondsSquared; 6])>,
    /// Optional 2-CPR empirical acceleration with the same shape.
    pub empirical_2cpr: Option<(JulianDate, Second, [KmPerSecondsSquared; 6])>,
}

impl Default for ForceModelConfig {
    /// Conservative default: two-body central gravity plus J2 only.
    /// All other perturbations are off and must be enabled
    /// explicitly by the caller.
    fn default() -> Self {
        Self {
            two_body: true,
            j2: true,
            harmonics: None,
            third_body_sun: false,
            third_body_moon: false,
            drag: None,
            srp: None,
            relativity: false,
            empirical_constant: None,
            empirical_1cpr: None,
            empirical_2cpr: None,
        }
    }
}

impl ForceModelConfig {
    /// Translate the boolean / option flags into an ordered list of
    /// registry specs. The order is deterministic and matches the field
    /// order on this struct, so it can be hashed into the run manifest.
    ///
    /// # Example
    ///
    /// ```
    /// use siderust::pod::dynamics::ForceModelConfig;
    /// let mut cfg = ForceModelConfig::default();
    /// cfg.relativity = true;
    /// let specs = cfg.to_specs();
    /// assert_eq!(specs.last().unwrap().name, "relativity");
    /// ```
    pub fn to_specs(&self) -> Vec<ForceModelSpec> {
        let mut out = Vec::new();
        if self.two_body {
            out.push(ForceModelSpec::named("two_body"));
        }
        if self.j2 {
            out.push(ForceModelSpec::named("j2"));
        }
        if let Some((d, o)) = self.harmonics {
            out.push(ForceModelSpec::with_params(
                "geopotential",
                ForceModelParams::Geopotential {
                    degree: d,
                    order: o,
                },
            ));
        }
        if self.third_body_sun {
            out.push(ForceModelSpec::named("third_body_sun"));
        }
        if self.third_body_moon {
            out.push(ForceModelSpec::named("third_body_moon"));
        }
        if let Some((cd, am)) = self.drag {
            out.push(ForceModelSpec::with_params(
                "drag",
                ForceModelParams::Drag {
                    cd,
                    area_to_mass: am,
                },
            ));
        }
        if let Some((cr, am, sh)) = self.srp {
            out.push(ForceModelSpec::with_params(
                "srp_cannonball",
                ForceModelParams::SrpCannonball {
                    cr,
                    area_to_mass: am,
                    shadow: sh,
                },
            ));
        }
        if self.relativity {
            out.push(ForceModelSpec::named("relativity"));
        }
        if let Some((r, t, n)) = self.empirical_constant {
            out.push(ForceModelSpec::with_params(
                "empirical_constant",
                ForceModelParams::EmpiricalConstant {
                    radial: r,
                    transverse: t,
                    normal: n,
                },
            ));
        }
        if let Some((epoch, per, c)) = self.empirical_1cpr {
            out.push(ForceModelSpec::with_params(
                "empirical_1cpr",
                ForceModelParams::EmpiricalPeriodic {
                    harmonic: PeriodicHarmonic::OncePerRev,
                    epoch_ref: epoch,
                    period: per,
                    coeffs: c,
                },
            ));
        }
        if let Some((epoch, per, c)) = self.empirical_2cpr {
            out.push(ForceModelSpec::with_params(
                "empirical_2cpr",
                ForceModelParams::EmpiricalPeriodic {
                    harmonic: PeriodicHarmonic::TwicePerRev,
                    epoch_ref: epoch,
                    period: per,
                    coeffs: c,
                },
            ));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_emits_two_specs() {
        let s = ForceModelConfig::default().to_specs();
        assert_eq!(s.len(), 2);
        assert_eq!(s[0].name, "two_body");
        assert_eq!(s[1].name, "j2");
    }

    #[test]
    fn full_config_round_trips_through_registry() {
        let cfg = ForceModelConfig {
            two_body: true,
            j2: true,
            harmonics: Some((4, 4)),
            third_body_sun: true,
            third_body_moon: true,
            drag: Some((DragCoefficient::new(2.2), AreaToMass::new(0.01))),
            srp: Some((
                SrpCoefficient::new(1.5),
                AreaToMass::new(0.02),
                ShadowModel::Conical,
            )),
            relativity: true,
            empirical_constant: Some((
                KmPerSecondsSquared::new(0.0),
                KmPerSecondsSquared::new(1e-12),
                KmPerSecondsSquared::new(0.0),
            )),
            empirical_1cpr: None,
            empirical_2cpr: None,
        };
        let specs = cfg.to_specs();
        assert_eq!(specs.len(), 9);
        let reg = crate::pod::dynamics::registry::ForceModelRegistry::with_builtins();
        let composite = reg.build(&specs).unwrap();
        assert_eq!(composite.len(), 9);
    }
}
