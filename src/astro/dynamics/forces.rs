// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Force models for spacecraft dynamics.
//!
//! A [`ForceModel`] returns a typed inertial [`Acceleration`] given an
//! [`OrbitState`]. Two textbook models are provided here — central
//! gravity ([`TwoBody`]) and the J2 (Earth oblateness) perturbation
//! ([`J2`]). More advanced perturbations (third body, drag, SRP) live in
//! `siderust-pod-dynamics` until they are matured for upstreaming.
//!
//! ## Example
//!
//! ```rust
//! use siderust::astro::dynamics::forces::{ForceModel, TwoBody};
//! use siderust::astro::dynamics::{OrbitState, Position, Velocity};
//! use siderust::coordinates::frames::GCRS;
//! use siderust::time::JulianDate;
//!
//! let s = OrbitState::new(
//!     JulianDate::new(2_451_545.0),
//!     Position::<GCRS>::new(7000.0, 0.0, 0.0),
//!     Velocity::<GCRS>::new(0.0, 7.5, 0.0),
//! );
//! let a = TwoBody::earth().acceleration(&s);
//! assert!(a.x().value() < 0.0);
//! ```

use crate::astro::dynamics::atmosphere::{DensityProvider, ExponentialAtmosphere};
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::frames::GCRS;
use crate::qtty::Kilometers;

// =============================================================================
// Earth constants (used by multiple force models)
// =============================================================================

/// Earth mean equatorial radius (GRS-80 / WGS-84), km.
pub const R_EARTH_KM: f64 = 6_378.137;

/// Earth mean rotation rate (sidereal), rad/s.
pub const OMEGA_EARTH_RAD_S: f64 = 7.292_115e-5;

// =============================================================================
// Trait
// =============================================================================

/// A force model evaluated on an inertial [`OrbitState`].
///
/// The returned value is the inertial acceleration in the parent frame
/// (currently [`GCRS`]) with canonical units (km/s²).  Internal numerics
/// remain `f64`; the trait surface is typed.
pub trait ForceModel: Send + Sync {
    /// Inertial acceleration acting on `state`.
    fn acceleration(&self, state: &OrbitState) -> Acceleration<GCRS, AccelerationUnit>;
}

// =============================================================================
// CompositeForce
// =============================================================================

/// Sum of N component force models.
pub struct CompositeForce {
    components: Vec<Box<dyn ForceModel>>,
}

impl Default for CompositeForce {
    fn default() -> Self {
        Self::empty()
    }
}

impl CompositeForce {
    /// Empty force (returns zero).
    pub fn empty() -> Self {
        Self { components: Vec::new() }
    }

    /// Append a force component.
    pub fn push(mut self, f: Box<dyn ForceModel>) -> Self {
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

impl ForceModel for CompositeForce {
    fn acceleration(&self, state: &OrbitState) -> Acceleration<GCRS, AccelerationUnit> {
        let mut acc = Acceleration::<GCRS, AccelerationUnit>::new(0.0, 0.0, 0.0);
        for f in &self.components {
            acc = acc + f.acceleration(state);
        }
        acc
    }
}

// =============================================================================
// Two-body central gravity
// =============================================================================

/// Newtonian central-gravity acceleration `−μ r / |r|³`.
#[derive(Debug, Clone, Copy)]
pub struct TwoBody {
    /// Gravitational parameter `GM` in km³/s².
    pub gm: f64,
}

impl TwoBody {
    /// Earth two-body field with EGM2008 GM.
    pub fn earth() -> Self {
        Self { gm: 398_600.441_8 }
    }
}

impl ForceModel for TwoBody {
    #[inline]
    fn acceleration(&self, s: &OrbitState) -> Acceleration<GCRS, AccelerationUnit> {
        let r = s.position.distance().value();
        let r2 = r * r;
        let k = -self.gm / (r2 * r);
        Acceleration::<GCRS, AccelerationUnit>::new(
            k * s.position.x().value(),
            k * s.position.y().value(),
            k * s.position.z().value(),
        )
    }
}

// =============================================================================
// J2 perturbation
// =============================================================================

/// J2 (Earth oblateness) perturbation acceleration in the inertial frame.
///
/// `J2` is the unnormalised coefficient (Earth: `1.082_626_68e-3`).
#[derive(Debug, Clone, Copy)]
pub struct J2 {
    /// `GM` in km³/s².
    pub gm: f64,
    /// Equatorial radius in km.
    pub req_km: f64,
    /// Unnormalised `J2` coefficient.
    pub j2: f64,
}

impl J2 {
    /// Standard Earth values.
    pub fn earth() -> Self {
        Self {
            gm: 398_600.441_8,
            req_km: R_EARTH_KM,
            j2: 1.082_626_68e-3,
        }
    }
}

impl ForceModel for J2 {
    #[inline]
    fn acceleration(&self, s: &OrbitState) -> Acceleration<GCRS, AccelerationUnit> {
        let r = s.position.distance().value();
        let r2 = r * r;
        let rx = s.position.x().value();
        let ry = s.position.y().value();
        let rz = s.position.z().value();
        let z2_over_r2 = (rz * rz) / r2;
        let factor = 1.5 * self.j2 * self.gm * self.req_km * self.req_km / (r2 * r2 * r);
        let cx = 5.0 * z2_over_r2 - 1.0;
        let cz = 5.0 * z2_over_r2 - 3.0;
        Acceleration::<GCRS, AccelerationUnit>::new(
            factor * rx * cx,
            factor * ry * cx,
            factor * rz * cz,
        )
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
/// * `Cd` is the drag coefficient (typical LEO value ≈ 2.2);
/// * `A/m` is the area-to-mass ratio in m² / kg;
/// * `ρ(h)` is the atmospheric mass density provided by a [`DensityProvider`];
/// * `h` is the geocentric altitude `|r| − R_⊕` (km);
/// * `v_rel = v − ω_⊕ × r` is the inertial velocity corrected for
///   the co-rotating atmosphere (Earth angular velocity
///   ω_⊕ = (0, 0, 7.292 115 × 10⁻⁵) rad/s).
///
/// This is a *minimum-viable* atmosphere: density only depends on the
/// magnitude of `r`, no diurnal/latitudinal/solar-activity effects.
#[derive(Debug, Clone)]
pub struct DragForce<D: DensityProvider> {
    /// Drag coefficient (dimensionless).
    pub cd: f64,
    /// Area-to-mass ratio, m² / kg.
    pub area_to_mass_m2_kg: f64,
    /// Atmosphere density provider.
    pub atmosphere: D,
}

impl DragForce<ExponentialAtmosphere> {
    /// Build a drag model using the [`ExponentialAtmosphere::LEO_500KM`] profile.
    pub fn leo_500km(cd: f64, area_to_mass_m2_kg: f64) -> Self {
        Self { cd, area_to_mass_m2_kg, atmosphere: ExponentialAtmosphere::LEO_500KM }
    }
}

impl<D: DensityProvider> ForceModel for DragForce<D> {
    #[inline]
    fn acceleration(&self, s: &OrbitState) -> Acceleration<GCRS, AccelerationUnit> {
        let r = s.position.distance().value();
        let h = r - R_EARTH_KM;
        if h < 0.0 {
            return Acceleration::<GCRS, AccelerationUnit>::new(0.0, 0.0, 0.0);
        }
        let rho = self.atmosphere.density_kg_m3(Kilometers::new(h));

        let rx = s.position.x().value();
        let ry = s.position.y().value();
        let vx = s.velocity.x().value();
        let vy = s.velocity.y().value();
        let vz = s.velocity.z().value();

        let omega_cross_r = [-OMEGA_EARTH_RAD_S * ry, OMEGA_EARTH_RAD_S * rx, 0.0_f64];
        let v_rel = [vx - omega_cross_r[0], vy - omega_cross_r[1], vz - omega_cross_r[2]];

        let v_mag_m_s = (v_rel[0].powi(2) + v_rel[1].powi(2) + v_rel[2].powi(2)).sqrt() * 1_000.0;
        let pre = -0.5 * self.cd * self.area_to_mass_m2_kg * rho * v_mag_m_s;
        Acceleration::<GCRS, AccelerationUnit>::new(pre * v_rel[0], pre * v_rel[1], pre * v_rel[2])
    }
}

/// Type alias: drag model with the built-in exponential atmosphere.
pub type ExponentialDrag = DragForce<ExponentialAtmosphere>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::time::JulianDate;

    fn leo() -> OrbitState {
        OrbitState::new(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn two_body_points_radially_inward() {
        let a = TwoBody::earth().acceleration(&leo());
        assert!(a.x().value() < 0.0);
        assert!(a.y().value().abs() < 1e-12);
        assert!(a.z().value().abs() < 1e-12);
    }

    #[test]
    fn composite_sums_components() {
        let f = CompositeForce::empty()
            .push(Box::new(TwoBody::earth()))
            .push(Box::new(J2::earth()));
        let a = f.acceleration(&leo());
        assert!(a.x().value() < 0.0);
        assert!(a.z().value().abs() < 1e-12);
    }

    #[test]
    fn drag_density_decreases_with_altitude() {
        use crate::astro::dynamics::atmosphere::ExponentialAtmosphere;
        use crate::qtty::Kilometers;
        let d = DragForce::leo_500km(2.2, 0.02);
        assert!(
            d.atmosphere.density_kg_m3(Kilometers::new(500.0))
                > d.atmosphere.density_kg_m3(Kilometers::new(600.0))
        );
        assert!(
            d.atmosphere.density_kg_m3(Kilometers::new(400.0))
                > d.atmosphere.density_kg_m3(Kilometers::new(500.0))
        );
    }

    #[test]
    fn drag_decays_orbit_altitude() {
        use crate::astro::dynamics::atmosphere::ExponentialAtmosphere;
        use crate::astro::dynamics::integrators::rk4_propagate;
        let mu: f64 = 398_600.441_8;
        let r0 = R_EARTH_KM + 350.0;
        let v0 = (mu / r0).sqrt();
        let s0 = OrbitState::new(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v0, 0.0),
        );
        let force = CompositeForce::empty()
            .push(Box::new(TwoBody::earth()))
            .push(Box::new(DragForce {
                cd: 2.2,
                area_to_mass_m2_kg: 5.0,
                atmosphere: ExponentialAtmosphere {
                    rho0_kg_m3: 1.0e-11,
                    h0_km: 350.0,
                    scale_height_km: 50.0,
                },
            }));
        let s_end = rk4_propagate(&force, s0, 30.0, 1440);
        let r_end = s_end.position.distance().value();
        assert!(r_end < r0, "expected drag-driven decay; r0={r0:.3}, r_end={r_end:.3}");
    }
}
