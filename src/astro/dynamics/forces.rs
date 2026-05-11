// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Force models for spacecraft dynamics.
//!
//! ## Trait design
//!
//! [`ForceModel<C, F>`] is the central abstraction.  Each call to
//! [`acceleration`][ForceModel::acceleration] or [`partials`][ForceModel::partials]
//! receives a [`DynamicsContext`] reference so force models can query
//! ephemeris, atmosphere, and gravity-field providers **without holding them
//! as fields**.  Force models store only their tunable physical parameters
//! (C_D, C_R, A/m, J2 coefficient, truncation degree, …).
//!
//! ## Partial derivatives
//!
//! [`ForcePartials<F>`] is the linearization block
//! `A(t) = ∂a/∂[r, v]`, the lower half of the variational matrix
//! `F(t) = [[0, I], [A_r, A_v]]`:
//!
//! * `d_acc_d_pos` (`A_r = ∂a/∂r`): units km/s² per km = s⁻² (stored as raw f64 in frame `F`)
//! * `d_acc_d_vel` (`A_v = ∂a/∂v`): units km/s² per km/s = s⁻¹ (stored as raw f64 in frame `F`)
//!
//! Only [`TwoBody`] supplies analytic partials in this release.  All other
//! models return `Err(DynamicsError::Provider(_))` from the default impl.
//!
//! ## Provided models
//!
//! | Model | Type parameters | Description |
//! |-------|-----------------|-------------|
//! | [`TwoBody`] | `<Geocentric, GCRS>` | Central Newtonian gravity |
//! | [`J2`] | `<Geocentric, GCRS>` | Zonal oblateness perturbation |
//! | [`DragForce<D>`] | `<Geocentric, GCRS>` | Cannonball atmospheric drag |
//! | [`ThirdBodySunMoon`] | `<Geocentric, GCRS>` | Sun + Moon point-mass perturbation |
//! | [`CannonballSrp`] | `<Geocentric, GCRS>` | Cannonball solar radiation pressure |
//! | [`CompositeForce<C, F>`] | generic | Linear sum of any force models |

use std::sync::Arc;

use affn::cartesian::Displacement;
use affn::matrix3::FrameMatrix3;

use crate::astro::dynamics::atmosphere::{DensityProvider, ExponentialAtmosphere};
use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::astro::dynamics::units::{GravitationalParameter, GM_EARTH, GM_MOON, GM_SUN};
use crate::astro::precession::ecliptic_of_date_to_mean_equatorial_matrix;
use crate::calculus::ephemeris::DynEphemeris;
use crate::coordinates::centers::{Geocentric, ReferenceCenter};
use crate::coordinates::frames::{EclipticMeanJ2000, ReferenceFrame, GCRS};
use crate::qtty::{
    AreaToMass, AstronomicalUnit, DragCoefficient, J2Coefficient, Kilometer, Kilometers, Pascals,
    SrpCoefficient, Unit,
};
use crate::time::JulianDate;

// =============================================================================
// Earth constants (used by multiple force models)
// =============================================================================

/// Earth mean equatorial radius (GRS-80 / WGS-84).
pub const R_EARTH: Kilometers = Kilometers::new(6_378.137);

/// Earth mean rotation rate (sidereal), rad/s (IAU 2000).
pub const OMEGA_EARTH_RAD_S: f64 = 7.292_115_146_706_979e-5;

/// Solar radiation pressure at 1 AU.
pub const P0: Pascals = Pascals::new(4.560e-6);

/// Standard gravitational parameter of the Sun.
pub const MU_SUN: GravitationalParameter = GM_SUN;
/// Standard gravitational parameter of the Moon.
pub const MU_MOON: GravitationalParameter = GM_MOON;

/// Astronomical unit in km (derived from unit ratios, not a magic number).
const AU_IN_KM: f64 = AstronomicalUnit::RATIO / Kilometer::RATIO;

// =============================================================================
// ForcePartials
// =============================================================================

/// Frame-tagged Jacobian blocks of the acceleration: `∂a/∂[r, v]`.
///
/// These are the `A_r` and `A_v` matrices that form the lower half of the
/// variational dynamics matrix `F(t)`:
///
/// ```text
/// F(t) = [ 0   I  ]
///        [ A_r A_v ]
/// ```
///
/// where:
/// * `A_r = ∂a/∂r` (units: km/s² per km = s⁻²)
/// * `A_v = ∂a/∂v` (units: km/s² per km/s = s⁻¹)
///
/// Both matrices are stored as raw `f64` values in the tagged frame `F`.
/// For conservative forces such as gravity, `A_v = 0`.
#[derive(Debug, Clone, Copy)]
pub struct ForcePartials<F = GCRS> {
    /// `∂a/∂r` in frame `F` (units: s⁻²).
    pub d_acc_d_pos: FrameMatrix3<F>,
    /// `∂a/∂v` in frame `F` (units: s⁻¹). Zero for conservative forces.
    pub d_acc_d_vel: FrameMatrix3<F>,
}

impl<F> ForcePartials<F> {
    /// Zero partials — both Jacobian blocks are zero matrices.
    ///
    /// Used as a neutral element when compositing forces whose partials are
    /// not implemented.
    pub fn zero() -> Self {
        Self {
            d_acc_d_pos: FrameMatrix3::zero(),
            d_acc_d_vel: FrameMatrix3::zero(),
        }
    }

    /// Analytic `∂a/∂r` for the Newtonian two-body acceleration `a = -μ r / |r|³`.
    ///
    /// The formula is:
    ///
    /// ```text
    /// A_r = -μ/r³ I + 3μ/r⁵ r rᵀ
    /// ```
    ///
    /// `A_v = 0` because central gravity has no velocity dependence.
    ///
    /// # Arguments
    ///
    /// * `mu`  — standard gravitational parameter μ = GM (km³/s²).
    /// * `r`   — Cartesian position vector `[x, y, z]` (km) in frame `F`.
    pub fn two_body(mu: GravitationalParameter, r: [f64; 3]) -> Self {
        let r_norm = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
        if r_norm == 0.0 {
            return Self::zero();
        }
        let r3 = r_norm.powi(3);
        let r5 = r_norm.powi(5);
        let mu_val = mu.value();
        let diag = -mu_val / r3;
        let scale = 3.0 * mu_val / r5;
        let mut data = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                data[i][j] = scale * r[i] * r[j] + if i == j { diag } else { 0.0 };
            }
        }
        Self {
            d_acc_d_pos: FrameMatrix3::from_array(data),
            d_acc_d_vel: FrameMatrix3::zero(),
        }
    }

    /// Element-wise sum of `self` and `other`, returning a new [`ForcePartials`].
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        let a_pos = self.d_acc_d_pos.as_array();
        let b_pos = other.d_acc_d_pos.as_array();
        let mut out_pos = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                out_pos[i][j] = a_pos[i][j] + b_pos[i][j];
            }
        }
        let a_vel = self.d_acc_d_vel.as_array();
        let b_vel = other.d_acc_d_vel.as_array();
        let mut out_vel = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                out_vel[i][j] = a_vel[i][j] + b_vel[i][j];
            }
        }
        Self {
            d_acc_d_pos: FrameMatrix3::from_array(out_pos),
            d_acc_d_vel: FrameMatrix3::from_array(out_vel),
        }
    }

    /// Element-wise add `other` into `self` in place.
    pub fn add_in_place(&mut self, other: &Self) {
        let b_pos = *other.d_acc_d_pos.as_array();
        let cur_pos = *self.d_acc_d_pos.as_array();
        let mut out_pos = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                out_pos[i][j] = cur_pos[i][j] + b_pos[i][j];
            }
        }
        self.d_acc_d_pos = FrameMatrix3::from_array(out_pos);

        let b_vel = *other.d_acc_d_vel.as_array();
        let cur_vel = *self.d_acc_d_vel.as_array();
        let mut out_vel = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                out_vel[i][j] = cur_vel[i][j] + b_vel[i][j];
            }
        }
        self.d_acc_d_vel = FrameMatrix3::from_array(out_vel);
    }
}

// =============================================================================
// ForceModel trait
// =============================================================================

/// A force model evaluated on an inertial [`OrbitState<C, F>`].
///
/// All providers (ephemeris, atmosphere, gravity field) are accessed via
/// `ctx`.  Implementors store only their tunable physical parameters.
///
/// # Type parameters
///
/// * `C` — reference center (default [`Geocentric`]).
/// * `F` — reference frame (default [`GCRS`]).
///
/// # Units
///
/// The returned acceleration is in km/s² in frame `F`.
///
/// # Errors
///
/// Returns a [`DynamicsError`] when:
/// * a required provider is absent from `ctx` (e.g. ephemeris for third-body),
/// * the spacecraft is below the surface,
/// * geometry is degenerate.
///
/// # Partial derivatives
///
/// The default `partials` implementation returns
/// `Err(DynamicsError::Provider(_))` ("analytic partials not implemented").
/// Override it in force models that have a closed-form Jacobian.
///
/// The `partials` convention is the lower half of `F(t) = [[0, I], [A_r, A_v]]`:
///
/// * `A_r = ∂a/∂r` (3×3, units: s⁻²)
/// * `A_v = ∂a/∂v` (3×3, units: s⁻¹; zero for conservative forces)
pub trait ForceModel<C = Geocentric, F = GCRS>: Send + Sync
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    /// Compute the inertial acceleration acting on `state`.
    ///
    /// # Frame
    ///
    /// The result is expressed in frame `F` with units km/s².
    fn acceleration(
        &self,
        state: &OrbitState<C, F>,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<F, AccelerationUnit>, DynamicsError>;

    /// Analytic partial derivatives of the acceleration with respect to
    /// position and velocity: `(∂a/∂r, ∂a/∂v)` in frame `F`.
    ///
    /// The default implementation returns an error indicating that analytic
    /// partials are not available for this model.  Override this method in
    /// force models that have a closed-form Jacobian.
    fn partials(
        &self,
        _state: &OrbitState<C, F>,
        _ctx: &DynamicsContext,
    ) -> Result<ForcePartials<F>, DynamicsError> {
        Err(DynamicsError::Provider(Box::new(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "analytic partials not implemented for this force model",
        ))))
    }

    /// Human-readable name identifying this force model.
    ///
    /// The default returns the fully-qualified Rust type name.
    fn name(&self) -> &'static str {
        std::any::type_name::<Self>()
    }
}

// =============================================================================
// CompositeForce
// =============================================================================

/// Sum of N component force models.
///
/// `acceleration` sums child accelerations, short-circuiting on the first
/// error.  `partials` sums child [`ForcePartials`], also short-circuiting on
/// the first error.
///
/// # Type parameters
///
/// * `C` — reference center (default [`Geocentric`]).
/// * `F` — reference frame (default [`GCRS`]).
pub struct CompositeForce<C = Geocentric, F = GCRS>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    components: Vec<Box<dyn ForceModel<C, F>>>,
}

impl<C, F> Default for CompositeForce<C, F>
where
    C: ReferenceCenter + 'static,
    F: ReferenceFrame + 'static,
{
    fn default() -> Self {
        Self::empty()
    }
}

impl<C, F> CompositeForce<C, F>
where
    C: ReferenceCenter + 'static,
    F: ReferenceFrame + 'static,
{
    /// Empty composite force (returns zero acceleration).
    pub fn empty() -> Self {
        Self {
            components: Vec::new(),
        }
    }

    /// Append a force component.
    pub fn push(mut self, f: Box<dyn ForceModel<C, F>>) -> Self {
        self.components.push(f);
        self
    }

    /// Number of components.
    pub fn len(&self) -> usize {
        self.components.len()
    }

    /// True when no components have been added.
    pub fn is_empty(&self) -> bool {
        self.components.is_empty()
    }
}

impl<C, F> ForceModel<C, F> for CompositeForce<C, F>
where
    C: ReferenceCenter + 'static,
    F: ReferenceFrame + 'static,
{
    fn acceleration(
        &self,
        state: &OrbitState<C, F>,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<F, AccelerationUnit>, DynamicsError> {
        let mut acc = Acceleration::<F, AccelerationUnit>::new(0.0, 0.0, 0.0);
        for f in &self.components {
            acc = acc + f.acceleration(state, ctx)?;
        }
        Ok(acc)
    }

    fn partials(
        &self,
        state: &OrbitState<C, F>,
        ctx: &DynamicsContext,
    ) -> Result<ForcePartials<F>, DynamicsError> {
        let mut total = ForcePartials::zero();
        for f in &self.components {
            total.add_in_place(&f.partials(state, ctx)?);
        }
        Ok(total)
    }
}

// =============================================================================
// Two-body central gravity
// =============================================================================

/// Newtonian central-gravity acceleration `−μ r / |r|³`.
///
/// This model provides analytic partial derivatives via the
/// [`ForceModel::partials`] method.  All other force models in this module
/// use the default error-returning implementation.
#[derive(Debug, Clone, Copy)]
pub struct TwoBody {
    /// Gravitational parameter `GM` (km³/s²).
    pub gm: GravitationalParameter,
}

impl TwoBody {
    /// Earth two-body field with EGM2008 GM.
    pub fn earth() -> Self {
        Self { gm: GM_EARTH }
    }
}

impl ForceModel for TwoBody {
    #[inline]
    fn acceleration(
        &self,
        s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        let r = s.position.distance().value();
        let r2 = r * r;
        let k = -self.gm.value() / (r2 * r);
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            k * s.position.x().value(),
            k * s.position.y().value(),
            k * s.position.z().value(),
        ))
    }

    fn partials(
        &self,
        s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<ForcePartials<GCRS>, DynamicsError> {
        let r = [
            s.position.x().value(),
            s.position.y().value(),
            s.position.z().value(),
        ];
        Ok(ForcePartials::two_body(self.gm, r))
    }
}

// =============================================================================
// J2 perturbation
// =============================================================================

/// J2 (Earth oblateness) perturbation acceleration in the inertial frame.
///
/// `J2` is the unnormalised zonal harmonic coefficient (Earth: `1.082_626_68e-3`).
#[derive(Debug, Clone, Copy)]
pub struct J2 {
    /// `GM` (km³/s²).
    pub gm: GravitationalParameter,
    /// Equatorial radius (km).
    pub req: Kilometers,
    /// Unnormalised J₂ coefficient (dimensionless).
    pub j2: J2Coefficient,
}

impl J2 {
    /// Standard Earth values (EGM2008 / WGS-84).
    pub fn earth() -> Self {
        Self {
            gm: GM_EARTH,
            req: R_EARTH,
            j2: J2Coefficient::new(1.082_626_68e-3),
        }
    }
}

impl ForceModel for J2 {
    #[inline]
    fn acceleration(
        &self,
        s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        let r = s.position.distance().value();
        let r2 = r * r;
        let rx = s.position.x().value();
        let ry = s.position.y().value();
        let rz = s.position.z().value();
        let z2_over_r2 = (rz * rz) / r2;
        let req = self.req.value();
        let factor =
            1.5 * self.j2.value() * self.gm.value() * req * req / (r2 * r2 * r);
        let cx = 5.0 * z2_over_r2 - 1.0;
        let cz = 5.0 * z2_over_r2 - 3.0;
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            factor * rx * cx,
            factor * ry * cx,
            factor * rz * cz,
        ))
    }
}

// =============================================================================
// Drag force
// =============================================================================

/// Cannonball atmospheric drag acceleration.
///
/// ```text
/// a_drag = − ½ · Cd · (A/m) · ρ(h) · |v_rel| · v_rel
/// ```
///
/// where:
///
/// * `Cd` is the drag coefficient (typed [`DragCoefficient`]; typical LEO value ≈ 2.2);
/// * `A/m` is the area-to-mass ratio (typed [`AreaToMass`], m²/kg);
/// * `ρ(h)` is the atmospheric mass density from the embedded [`DensityProvider`];
/// * `h = |r| − R_⊕` is the geocentric altitude (km);
/// * `v_rel = v − ω_⊕ × r` is the velocity in the co-rotating atmosphere frame.
///
/// The embedded atmosphere `D` holds the density-model parameters (e.g.
/// exponential-atmosphere scale heights).  For a context-driven density
/// model see [`DynamicsContext::require_atmosphere`].
#[derive(Debug, Clone)]
pub struct DragForce<D: DensityProvider> {
    /// Drag coefficient C_D (dimensionless).
    pub cd: DragCoefficient,
    /// Effective area-to-mass ratio (m²/kg).
    pub area_to_mass: AreaToMass,
    /// Atmosphere density provider.
    pub atmosphere: D,
}

impl DragForce<ExponentialAtmosphere> {
    /// Build a drag model using the [`ExponentialAtmosphere::LEO_500KM`] profile.
    pub fn leo_500km(cd: f64, area_to_mass_m2_kg: f64) -> Self {
        Self {
            cd: DragCoefficient::new(cd),
            area_to_mass: AreaToMass::new(area_to_mass_m2_kg),
            atmosphere: ExponentialAtmosphere::LEO_500KM,
        }
    }
}

impl<D: DensityProvider + Send + Sync> ForceModel for DragForce<D> {
    #[inline]
    fn acceleration(
        &self,
        s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        let r = s.position.distance().value();
        let h = r - R_EARTH.value();
        if h < 0.0 {
            return Err(DynamicsError::AltitudeBelowSurface { altitude_km: h });
        }
        let rho = self.atmosphere.density(Kilometers::new(h)).value();

        let rx = s.position.x().value();
        let ry = s.position.y().value();
        let vx = s.velocity.x().value();
        let vy = s.velocity.y().value();
        let vz = s.velocity.z().value();

        let omega_cross_r = [-OMEGA_EARTH_RAD_S * ry, OMEGA_EARTH_RAD_S * rx, 0.0_f64];
        let v_rel = [
            vx - omega_cross_r[0],
            vy - omega_cross_r[1],
            vz - omega_cross_r[2],
        ];

        let v_mag_m_s =
            (v_rel[0].powi(2) + v_rel[1].powi(2) + v_rel[2].powi(2)).sqrt() * 1_000.0;
        let pre = -0.5 * self.cd.value() * self.area_to_mass.value() * rho * v_mag_m_s;
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            pre * v_rel[0],
            pre * v_rel[1],
            pre * v_rel[2],
        ))
    }
}

/// Type alias: drag model with the built-in exponential atmosphere.
pub type ExponentialDrag = DragForce<ExponentialAtmosphere>;

// =============================================================================
// Third-body perturbations (Sun + Moon)
// =============================================================================

// ---- Private helpers --------------------------------------------------------

/// Compute the geocentric Sun displacement (GCRS, km) from an ephemeris.
fn sun_geocentric(
    eph: &Arc<dyn DynEphemeris + Send + Sync>,
    jd: JulianDate,
) -> Result<Displacement<GCRS, Kilometer>, DynamicsError> {
    let sun_b = eph.try_sun_barycentric(jd).map_err(|e| {
        DynamicsError::EphemerisUnavailable {
            body: "Sun",
            source: Some(Box::new(e)),
        }
    })?;
    let earth_b = eph.try_earth_barycentric(jd).map_err(|e| {
        DynamicsError::EphemerisUnavailable {
            body: "Earth",
            source: Some(Box::new(e)),
        }
    })?;
    let d_ecl: Displacement<EclipticMeanJ2000, AstronomicalUnit> = sun_b - earth_b;
    let rot = ecliptic_of_date_to_mean_equatorial_matrix(JulianDate::J2000);
    let d_gcrs: Displacement<GCRS, AstronomicalUnit> = rot.apply_vec(d_ecl);
    Ok(d_gcrs.to_unit::<Kilometer>())
}

/// Compute the geocentric Moon displacement (GCRS, km) from an ephemeris.
fn moon_geocentric(
    eph: &Arc<dyn DynEphemeris + Send + Sync>,
    jd: JulianDate,
) -> Result<Displacement<GCRS, Kilometer>, DynamicsError> {
    let m = eph.try_moon_geocentric(jd).map_err(|e| {
        DynamicsError::EphemerisUnavailable {
            body: "Moon",
            source: Some(Box::new(e)),
        }
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
        let bodies: [(f64, Result<Displacement<GCRS, Kilometer>, DynamicsError>); 2] = [
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

// =============================================================================
// Solar Radiation Pressure
// =============================================================================

/// Cannonball solar radiation pressure (SRP).
///
/// ```text
/// a_srp = Cr · P0 · (AU / |r_sun_sat|)² · (A/m) · r̂_sun_sat
/// ```
///
/// The resulting acceleration pushes the spacecraft away from the Sun.
/// Eclipse modelling is not yet included — the Sun is always treated as visible.
/// The Sun geocentric position is fetched from the [`DynamicsContext`] ephemeris.
#[derive(Debug, Clone, Copy)]
pub struct CannonballSrp {
    /// Radiation-pressure coefficient (dimensionless).
    pub cr: SrpCoefficient,
    /// Area-to-mass ratio (m²/kg).
    pub area_to_mass: AreaToMass,
}

impl CannonballSrp {
    /// Build a cannonball SRP model.
    pub fn new(cr: f64, area_to_mass_m2_kg: f64) -> Self {
        Self {
            cr: SrpCoefficient::new(cr),
            area_to_mass: AreaToMass::new(area_to_mass_m2_kg),
        }
    }
}

impl ForceModel for CannonballSrp {
    #[inline]
    fn acceleration(
        &self,
        s: &OrbitState,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        let eph = ctx.require_ephemeris()?;
        let sun = sun_geocentric(eph, s.epoch_jd())?;
        let r_sun_sat = Displacement::<GCRS, Kilometer>::new(
            s.position.x().value() - sun.x().value(),
            s.position.y().value() - sun.y().value(),
            s.position.z().value() - sun.z().value(),
        );
        let r = r_sun_sat.magnitude().value();
        if r == 0.0 {
            return Ok(Acceleration::<GCRS, AccelerationUnit>::new(0.0, 0.0, 0.0));
        }
        let r2 = r * r;
        // N/m² · m²/kg = m/s²; convert to km/s² by dividing by 1000.
        let mag_km_s2 =
            self.cr.value() * P0.value() * (AU_IN_KM * AU_IN_KM / r2) * self.area_to_mass.value()
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
    use super::*;
    use crate::astro::dynamics::context::DynamicsContextBuilder;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::qtty::{KilogramsPerCubicMeter, Second};
    use crate::time::JulianDate;

    fn leo() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    fn leo_at(epoch: JulianDate) -> OrbitState {
        OrbitState::new_at_jd(
            epoch,
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn two_body_points_radially_inward() {
        let ctx = DynamicsContext::empty();
        let a = TwoBody::earth().acceleration(&leo(), &ctx).unwrap();
        assert!(a.x().value() < 0.0);
        assert!(a.y().value().abs() < 1e-12);
        assert!(a.z().value().abs() < 1e-12);
    }

    #[test]
    fn two_body_partials_nonzero() {
        let ctx = DynamicsContext::empty();
        let p = TwoBody::earth().partials(&leo(), &ctx).unwrap();
        // A_r diagonal is non-zero; A_v is zero
        let arr = p.d_acc_d_pos.as_array();
        assert!(arr[0][0].abs() > 1e-10 || arr[1][1].abs() > 1e-10);
    }

    #[test]
    fn composite_sums_components() {
        let ctx = DynamicsContext::empty();
        let f = CompositeForce::empty()
            .push(Box::new(TwoBody::earth()))
            .push(Box::new(J2::earth()));
        let a = f.acceleration(&leo(), &ctx).unwrap();
        assert!(a.x().value() < 0.0);
        assert!(a.z().value().abs() < 1e-12);
    }

    #[test]
    fn drag_density_decreases_with_altitude() {
        use crate::qtty::Kilometers;
        let d = DragForce::leo_500km(2.2, 0.02);
        assert!(
            d.atmosphere.density(Kilometers::new(500.0))
                > d.atmosphere.density(Kilometers::new(600.0))
        );
        assert!(
            d.atmosphere.density(Kilometers::new(400.0))
                > d.atmosphere.density(Kilometers::new(500.0))
        );
    }

    #[test]
    fn drag_decays_orbit_altitude() {
        use crate::astro::dynamics::atmosphere::ExponentialAtmosphere;
        use crate::astro::dynamics::integrators::rk4_propagate;
        let mu: f64 = 398_600.441_8;
        let r0 = R_EARTH.value() + 350.0;
        let v0 = (mu / r0).sqrt();
        let s0 = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v0, 0.0),
        );
        let force = CompositeForce::empty()
            .push(Box::new(TwoBody::earth()))
            .push(Box::new(DragForce {
                cd: DragCoefficient::new(2.2),
                area_to_mass: AreaToMass::new(5.0),
                atmosphere: ExponentialAtmosphere {
                    rho0: KilogramsPerCubicMeter::new(1.0e-11),
                    h0: Kilometers::new(350.0),
                    scale_height: Kilometers::new(50.0),
                },
            }));
        let ctx = DynamicsContext::empty();
        // 360 steps × 30 s = 3 h — enough to see drag-driven decay without reaching the surface.
        let s_end = rk4_propagate(&force, s0, Second::new(30.0), 360, &ctx).unwrap();
        let r_end = s_end.position.distance().value();
        assert!(
            r_end < r0,
            "expected drag-driven decay; r0={r0:.3}, r_end={r_end:.3}"
        );
    }

    #[test]
    fn third_body_nonzero_acceleration() {
        use crate::calculus::ephemeris::Vsop87Ephemeris;
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

    #[test]
    fn srp_order_of_magnitude_at_leo() {
        use crate::calculus::ephemeris::Vsop87Ephemeris;
        // Cr=1.5, A/m=0.02 m²/kg ⇒ |a_srp| ≈ 1.5 · 4.56e-6 · 0.02 / 1000 ≈ 1.4e-10 km/s²
        let ctx = DynamicsContextBuilder::new()
            .with_ephemeris(Arc::new(Vsop87Ephemeris))
            .build();
        let srp = CannonballSrp::new(1.5, 0.02);
        let s = leo_at(JulianDate::new(2_451_545.0));
        let a = srp.acceleration(&s, &ctx).unwrap();
        let mag = a.magnitude().value();
        assert!(
            (5e-11..5e-10).contains(&mag),
            "SRP magnitude out of expected band: {mag} km/s²"
        );
    }

    #[test]
    fn srp_zero_when_area_is_zero() {
        use crate::calculus::ephemeris::Vsop87Ephemeris;
        let ctx = DynamicsContextBuilder::new()
            .with_ephemeris(Arc::new(Vsop87Ephemeris))
            .build();
        let srp = CannonballSrp::new(1.5, 0.0);
        let a = srp.acceleration(&leo(), &ctx).unwrap();
        assert!(a.x().value() == 0.0 && a.y().value() == 0.0 && a.z().value() == 0.0);
    }
}
