// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Named force-model registry.
//!
//! The POD service layer ingests a YAML/JSON config and needs to translate
//! string keys (`"two_body"`, `"j2"`, `"drag"`, …) into concrete
//! [`ForceModel`] instances. The registry centralises that mapping so that:
//!
//! * Built-in models (two-body, J2, geopotential, Sun/Moon third body,
//!   cannonball SRP, drag, central-body relativity, constant/1-CPR/2-CPR
//!   empirical accelerations) are always available out of the box.
//! * Downstream crates can plug in additional factories without forking
//!   this crate.
//!
//! Each factory implements the [`ForceModelFactory`] trait. A
//! [`ForceModelSpec`] (`name + ForceModelParams`) is consumed by the
//! registry to produce a heap-allocated `Box<DynSiderustForceModel>` ready for
//! insertion into [`siderust::astro::dynamics::forces::CompositeForce`].
//!
//! # Example
//!
//! ```
//! use siderust::pod::dynamics::registry::{ForceModelRegistry, ForceModelSpec};
//! let mut reg = ForceModelRegistry::with_builtins();
//! let composite = reg
//!     .build(&[ForceModelSpec::named("two_body"), ForceModelSpec::named("j2")])
//!     .unwrap();
//! assert_eq!(composite.len(), 2);
//! ```

use std::collections::BTreeMap;

use crate::astro::dynamics::atmosphere::DensityProvider;
use crate::astro::dynamics::forces::{
    CannonballSrp, CentralBodyRelativity1Pn, Conical, Cylindrical, DragForce,
    EmpiricalAcceleration, Geopotential, J2, NoEclipse, ShadowModel, ThirdBody, TwoBody,
};
use crate::astro::dynamics::{EARTH_J2, GM_EARTH, R_EARTH};
use qtty::{AreaToMass, DragCoefficient, KmPerSecondsSquared, Second, SrpCoefficient};
use crate::time::JulianDate;

use super::empirical_periodic::{EmpiricalPeriodicAcceleration, PeriodicHarmonic};
use super::pod_error::PodDynamicsError;
use crate::pod::dynamics::{DynSiderustForceModel, SiderustCompositeModel};

/// Concrete parameter payload supplied with a [`ForceModelSpec`].
///
/// Each variant matches the parameter shape expected by the corresponding
/// built-in factory. Custom factories typically use [`ForceModelParams::None`]
/// or [`ForceModelParams::Custom`].
#[non_exhaustive]
#[derive(Debug, Clone)]
pub enum ForceModelParams {
    /// Factory needs no parameters (e.g. two-body, J2, relativity).
    None,
    /// Spherical-harmonic geopotential truncation `(degree, order)`.
    Geopotential {
        /// Maximum degree N.
        degree: usize,
        /// Maximum order M (≤ N).
        order: usize,
    },
    /// Cannonball drag parameters.
    Drag {
        /// Drag coefficient `C_D`.
        cd: DragCoefficient,
        /// Effective area-to-mass ratio (m²/kg).
        area_to_mass: AreaToMass,
    },
    /// Cannonball SRP parameters.
    SrpCannonball {
        /// SRP reflectivity coefficient `C_R`.
        cr: SrpCoefficient,
        /// Effective area-to-mass ratio (m²/kg).
        area_to_mass: AreaToMass,
        /// Earth shadow model.
        shadow: ShadowModel,
    },
    /// Constant RTN empirical acceleration.
    EmpiricalConstant {
        /// Radial component.
        radial: KmPerSecondsSquared,
        /// Transverse (along-track) component.
        transverse: KmPerSecondsSquared,
        /// Normal (cross-track) component.
        normal: KmPerSecondsSquared,
    },
    /// Periodic empirical acceleration (1-CPR or 2-CPR).
    EmpiricalPeriodic {
        /// Harmonic order.
        harmonic: PeriodicHarmonic,
        /// Reference epoch defining `θ = 0`.
        epoch_ref: JulianDate,
        /// Orbital period `T_orbit`.
        period: Second,
        /// `[r_cos, r_sin, t_cos, t_sin, n_cos, n_sin]` typed coefficients.
        coeffs: [KmPerSecondsSquared; 6],
    },
    /// Free-form custom payload, opaque to the built-in registry.
    Custom(String),
}

/// Declarative spec consumed by [`ForceModelRegistry::build`].
///
/// # Example
///
/// ```
/// use siderust::pod::dynamics::registry::ForceModelSpec;
/// let s = ForceModelSpec::named("j2");
/// assert_eq!(s.name, "j2");
/// ```
#[derive(Debug, Clone)]
pub struct ForceModelSpec {
    /// Registry key (e.g. `"two_body"`).
    pub name: String,
    /// Parameters consumed by the matching factory.
    pub params: ForceModelParams,
}

impl ForceModelSpec {
    /// Build a no-parameter spec (factory must accept [`ForceModelParams::None`]).
    pub fn named(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            params: ForceModelParams::None,
        }
    }

    /// Build a spec with explicit parameters.
    pub fn with_params(name: impl Into<String>, params: ForceModelParams) -> Self {
        Self {
            name: name.into(),
            params,
        }
    }
}

/// Trait implemented by registry factories.
pub trait ForceModelFactory: Send + Sync {
    /// Stable string key under which this factory is registered.
    fn name(&self) -> &'static str;

    /// Construct an instance from the supplied parameters.
    fn build(&self, params: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError>;
}

/// String-keyed force-model registry.
///
/// Use [`ForceModelRegistry::with_builtins`] to get a registry pre-populated
/// with the built-in factories listed in the crate-level docs, then
/// [`ForceModelRegistry::register`] to add custom factories.
pub struct ForceModelRegistry {
    factories: BTreeMap<String, Box<dyn ForceModelFactory>>,
}

impl Default for ForceModelRegistry {
    fn default() -> Self {
        Self::with_builtins()
    }
}

impl ForceModelRegistry {
    /// Create an empty registry.
    pub fn new() -> Self {
        Self {
            factories: BTreeMap::new(),
        }
    }

    /// Create a registry pre-populated with the built-in factories:
    /// `two_body`, `j2`, `geopotential`, `third_body_sun`, `third_body_moon`,
    /// `third_body_sun_moon`, `drag`, `srp_cannonball`, `srp_boxwing` (stub),
    /// `relativity`, `empirical_constant`, `empirical_1cpr`, `empirical_2cpr`.
    ///
    /// # Example
    ///
    /// ```
    /// use siderust::pod::dynamics::registry::ForceModelRegistry;
    /// let reg = ForceModelRegistry::with_builtins();
    /// assert!(reg.is_registered("two_body"));
    /// assert!(reg.is_registered("empirical_2cpr"));
    /// ```
    pub fn with_builtins() -> Self {
        let mut r = Self::new();
        r.register(Box::new(TwoBodyFactory));
        r.register(Box::new(J2Factory));
        r.register(Box::new(GeopotentialFactory));
        r.register(Box::new(ThirdBodySunFactory));
        r.register(Box::new(ThirdBodyMoonFactory));
        r.register(Box::new(ThirdBodySunMoonFactory));
        r.register(Box::new(DragFactory));
        r.register(Box::new(SrpCannonballFactory));
        r.register(Box::new(SrpBoxwingFactory));
        r.register(Box::new(RelativityFactory));
        r.register(Box::new(EmpiricalConstantFactory));
        r.register(Box::new(Empirical1CprFactory));
        r.register(Box::new(Empirical2CprFactory));
        r
    }

    /// Register or replace a factory.
    pub fn register(&mut self, f: Box<dyn ForceModelFactory>) {
        self.factories.insert(f.name().to_string(), f);
    }

    /// Return `true` iff a factory with `name` is currently registered.
    pub fn is_registered(&self, name: &str) -> bool {
        self.factories.contains_key(name)
    }

    /// All currently-registered names, sorted lexicographically.
    pub fn names(&self) -> Vec<&str> {
        self.factories.keys().map(|s| s.as_str()).collect()
    }

    /// Resolve a single spec into a boxed [`ForceModel`].
    ///
    /// # Errors
    ///
    /// * [`PodDynamicsError::UnknownModel`] if `spec.name` is not registered.
    /// * Any error returned by the factory.
    pub fn build_one(
        &self,
        spec: &ForceModelSpec,
    ) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        let f = self
            .factories
            .get(&spec.name)
            .ok_or_else(|| PodDynamicsError::UnknownModel(spec.name.clone()))?;
        f.build(&spec.params)
    }

    /// Resolve a batch of specs into a [`CompositeForce`] in declared order.
    ///
    /// # Example
    ///
    /// ```
    /// use siderust::pod::dynamics::registry::{ForceModelRegistry, ForceModelSpec, ForceModelParams};
    /// let reg = ForceModelRegistry::with_builtins();
    /// let composite = reg.build(&[
    ///     ForceModelSpec::named("two_body"),
    ///     ForceModelSpec::with_params(
    ///         "geopotential",
    ///         ForceModelParams::Geopotential { degree: 4, order: 4 },
    ///     ),
    /// ]).unwrap();
    /// assert_eq!(composite.len(), 2);
    /// ```
    pub fn build(&self, specs: &[ForceModelSpec]) -> Result<SiderustCompositeModel, PodDynamicsError> {
        let mut out = SiderustCompositeModel::empty();
        for s in specs {
            out = out.push(self.build_one(s)?);
        }
        Ok(out)
    }
}

// ───────────────────────────── built-in factories ──────────────────────────

struct TwoBodyFactory;
impl ForceModelFactory for TwoBodyFactory {
    fn name(&self) -> &'static str {
        "two_body"
    }
    fn build(&self, _p: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        Ok(Box::new(TwoBody::new(GM_EARTH)))
    }
}

struct J2Factory;
impl ForceModelFactory for J2Factory {
    fn name(&self) -> &'static str {
        "j2"
    }
    fn build(&self, _p: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        Ok(Box::new(J2::new(GM_EARTH, R_EARTH, EARTH_J2)))
    }
}

struct GeopotentialFactory;
impl ForceModelFactory for GeopotentialFactory {
    fn name(&self) -> &'static str {
        "geopotential"
    }
    fn build(&self, p: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        match p {
            ForceModelParams::Geopotential { degree, order } => {
                Ok(Box::new(Geopotential::new(*degree, *order)))
            }
            _ => Err(PodDynamicsError::InvalidParameters {
                name: "geopotential".into(),
                reason: "expected ForceModelParams::Geopotential { degree, order }",
            }),
        }
    }
}

struct ThirdBodySunFactory;
impl ForceModelFactory for ThirdBodySunFactory {
    fn name(&self) -> &'static str {
        "third_body_sun"
    }
    fn build(&self, _p: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        Ok(Box::new(ThirdBody::new().with_sun()))
    }
}

struct ThirdBodyMoonFactory;
impl ForceModelFactory for ThirdBodyMoonFactory {
    fn name(&self) -> &'static str {
        "third_body_moon"
    }
    fn build(&self, _p: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        Ok(Box::new(ThirdBody::new().with_moon()))
    }
}

struct ThirdBodySunMoonFactory;
impl ForceModelFactory for ThirdBodySunMoonFactory {
    fn name(&self) -> &'static str {
        "third_body_sun_moon"
    }
    fn build(&self, _p: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        Ok(Box::new(ThirdBody::sun_and_moon()))
    }
}

struct DragFactory;
impl ForceModelFactory for DragFactory {
    fn name(&self) -> &'static str {
        "drag"
    }
    fn build(&self, p: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        match p {
            ForceModelParams::Drag { cd, area_to_mass } => {
                Ok(Box::new(DragForce::new(*cd, *area_to_mass)))
            }
            _ => Err(PodDynamicsError::InvalidParameters {
                name: "drag".into(),
                reason: "expected ForceModelParams::Drag { cd, area_to_mass }",
            }),
        }
    }
}

struct SrpCannonballFactory;
impl ForceModelFactory for SrpCannonballFactory {
    fn name(&self) -> &'static str {
        "srp_cannonball"
    }
    fn build(&self, p: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        match p {
            ForceModelParams::SrpCannonball {
                cr,
                area_to_mass,
                shadow,
            } => {
                let force: Box<DynSiderustForceModel> = match shadow {
                    ShadowModel::None => {
                        Box::new(CannonballSrp::<NoEclipse>::new(*cr, *area_to_mass))
                    }
                    ShadowModel::Cylindrical => {
                        Box::new(CannonballSrp::<Cylindrical>::new(*cr, *area_to_mass))
                    }
                    ShadowModel::Conical => {
                        Box::new(CannonballSrp::<Conical>::new(*cr, *area_to_mass))
                    }
                };
                Ok(force)
            }
            _ => Err(PodDynamicsError::InvalidParameters {
                name: "srp_cannonball".into(),
                reason: "expected ForceModelParams::SrpCannonball { cr, area_to_mass, shadow }",
            }),
        }
    }
}

struct SrpBoxwingFactory;
impl ForceModelFactory for SrpBoxwingFactory {
    fn name(&self) -> &'static str {
        "srp_boxwing"
    }
    fn build(&self, _p: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        Err(PodDynamicsError::FeatureNotImplemented("srp_boxwing"))
    }
}

struct RelativityFactory;
impl ForceModelFactory for RelativityFactory {
    fn name(&self) -> &'static str {
        "relativity"
    }
    fn build(&self, _p: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        Ok(Box::new(CentralBodyRelativity1Pn::earth()))
    }
}

struct EmpiricalConstantFactory;
impl ForceModelFactory for EmpiricalConstantFactory {
    fn name(&self) -> &'static str {
        "empirical_constant"
    }
    fn build(&self, p: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        match p {
            ForceModelParams::EmpiricalConstant {
                radial,
                transverse,
                normal,
            } => Ok(Box::new(EmpiricalAcceleration::rtn(
                *radial,
                *transverse,
                *normal,
            ))),
            _ => Err(PodDynamicsError::InvalidParameters {
                name: "empirical_constant".into(),
                reason:
                    "expected ForceModelParams::EmpiricalConstant { radial, transverse, normal }",
            }),
        }
    }
}

fn build_periodic(
    expected: PeriodicHarmonic,
    label: &'static str,
    p: &ForceModelParams,
) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
    match p {
        ForceModelParams::EmpiricalPeriodic {
            harmonic,
            epoch_ref,
            period,
            coeffs,
        } if *harmonic == expected => Ok(Box::new(EmpiricalPeriodicAcceleration::new(
            *harmonic, *epoch_ref, *period, coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4],
            coeffs[5],
        ))),
        _ => Err(PodDynamicsError::InvalidParameters {
            name: label.into(),
            reason: "expected ForceModelParams::EmpiricalPeriodic with matching harmonic",
        }),
    }
}

struct Empirical1CprFactory;
impl ForceModelFactory for Empirical1CprFactory {
    fn name(&self) -> &'static str {
        "empirical_1cpr"
    }
    fn build(&self, p: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        build_periodic(PeriodicHarmonic::OncePerRev, "empirical_1cpr", p)
    }
}

struct Empirical2CprFactory;
impl ForceModelFactory for Empirical2CprFactory {
    fn name(&self) -> &'static str {
        "empirical_2cpr"
    }
    fn build(&self, p: &ForceModelParams) -> Result<Box<DynSiderustForceModel>, PodDynamicsError> {
        build_periodic(PeriodicHarmonic::TwicePerRev, "empirical_2cpr", p)
    }
}

/// Lightweight type alias for atmosphere density providers consumable by the
/// drag factory indirectly through [`siderust::astro::dynamics::context::DynamicsContext`].
///
/// Built-in models such as
/// [`siderust::astro::dynamics::atmosphere::ExponentialAtmosphere`] and
/// [`siderust::astro::dynamics::atmosphere::Nrlmsise00LiteApprox`] implement
/// this trait directly, and downstream crates can plug in MSIS-86, JB2008, …
/// by implementing [`DensityProvider`].
pub type AtmosphereDensityProvider = dyn DensityProvider + Send + Sync;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builtins_have_expected_names() {
        let r = ForceModelRegistry::with_builtins();
        for name in [
            "two_body",
            "j2",
            "geopotential",
            "third_body_sun",
            "third_body_moon",
            "third_body_sun_moon",
            "drag",
            "srp_cannonball",
            "srp_boxwing",
            "relativity",
            "empirical_constant",
            "empirical_1cpr",
            "empirical_2cpr",
        ] {
            assert!(r.is_registered(name), "missing builtin: {name}");
        }
    }

    #[test]
    fn unknown_model_errors() {
        let r = ForceModelRegistry::with_builtins();
        let e = r.build_one(&ForceModelSpec::named("nope")).err().unwrap();
        assert!(matches!(e, PodDynamicsError::UnknownModel(ref s) if s == "nope"));
    }

    #[test]
    fn srp_boxwing_returns_feature_not_implemented() {
        let r = ForceModelRegistry::with_builtins();
        let e = r
            .build_one(&ForceModelSpec::named("srp_boxwing"))
            .err()
            .unwrap();
        assert!(matches!(
            e,
            PodDynamicsError::FeatureNotImplemented("srp_boxwing")
        ));
    }

    #[test]
    fn drag_requires_typed_params() {
        let r = ForceModelRegistry::with_builtins();
        let e = r.build_one(&ForceModelSpec::named("drag")).err().unwrap();
        assert!(matches!(e, PodDynamicsError::InvalidParameters { .. }));
        let ok = r.build_one(&ForceModelSpec::with_params(
            "drag",
            ForceModelParams::Drag {
                cd: DragCoefficient::new(2.2),
                area_to_mass: AreaToMass::new(0.01),
            },
        ));
        assert!(ok.is_ok());
    }
}
