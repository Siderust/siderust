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

use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::frames::GCRS;

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
            req_km: 6_378.137,
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
        // Two-body alone is purely along -x for an equatorial state on +x.
        // J2 adds a small additional radial component but no out-of-plane
        // contribution at z=0.
        assert!(a.x().value() < 0.0);
        assert!(a.z().value().abs() < 1e-12);
    }
}
