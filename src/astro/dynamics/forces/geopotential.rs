// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Full spherical-harmonic geopotential force model.
//!
//! ## Scope
//!
//! Provides [`Geopotential`], a force model that delegates to the
//! [`geopotential_acceleration`] kernel to compute spherical-harmonic
//! gravitational acceleration up to a configurable degree/order.
//!
//! ## Equations
//!
//! Computes `a = ∇U` via [`geopotential_acceleration`] where
//!
//! ```text
//! U = (GM/r) Σ_{n=0}^{N} (R/r)^n  Σ_{m=0}^{n}
//!       [ C̄_{nm} P̄_{nm}(sin φ) cos(mλ)
//!       + S̄_{nm} P̄_{nm}(sin φ) sin(mλ) ]
//! ```
//!
//! The computation is constrained to `degree` and `order` fields, but
//! silently clamped to the provider's limits.
//!
//! ## Units & frames
//!
//! Position km (geocentric).  Acceleration km/s² (same frame as input).
//!
//! ## Failure modes
//!
//! | Error | Condition |
//! |-------|-----------|
//! | [`DynamicsError::GravityFieldUnavailable`] | Context has no gravity provider |
//! | [`DynamicsError::GeopotentialDegreeOutOfRange`] | Requested degree > provider max |
//!
//! ## Validity limits
//!
//! The accuracy is limited by:
//! - The provider's maximum degree/order.
//! - Truncation at the configured `degree` and `order` parameters.
//! - Floating-point precision at very high degrees (n > 180).
//!
//! High-degree models (e.g. EGM2008 full degree 2160) require a full coefficient set
//! and are not included in the built-in providers.
//!
//! ## References
//!
//! * Montenbruck & Gill, *Satellite Orbits*, §3.2.
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, §8.6.

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::gravity::geopotential_acceleration;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::frames::GCRS;

use super::traits::{ForceModel, ForcePartials};

/// Normalised spherical-harmonic geopotential perturbation force model.
///
/// The force model computes the full `(n, m)` harmonic acceleration by
/// querying the [`GravityFieldProvider`][crate::astro::dynamics::gravity::GravityFieldProvider]
/// in the current [`DynamicsContext`].
///
/// The `degree` and `order` fields cap the summation.  Both are silently
/// clamped to the provider's declared limits, but requesting a degree
/// **above** the provider limit returns
/// [`GeopotentialDegreeOutOfRange`][DynamicsError::GeopotentialDegreeOutOfRange].
#[derive(Debug, Clone, Copy)]
pub struct Geopotential {
    /// Maximum degree `N` of the harmonic summation.
    pub degree: usize,
    /// Maximum order `M` of the harmonic summation (≤ degree).
    pub order: usize,
}

impl Geopotential {
    /// Construct a new `Geopotential` with the given truncation limits.
    ///
    /// `order` is clamped to `degree` internally during evaluation.
    pub fn new(degree: usize, order: usize) -> Self {
        Self { degree, order }
    }

    /// Convenience constructor: full order `(n, n)` summation.
    pub fn full(degree: usize) -> Self {
        Self { degree, order: degree }
    }
}

impl ForceModel for Geopotential {
    fn acceleration(
        &self,
        s: &OrbitState,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        let provider = ctx.require_gravity_field()?;

        let pos = [
            s.position.x().value(),
            s.position.y().value(),
            s.position.z().value(),
        ];

        let a = geopotential_acceleration(provider.as_ref(), pos, self.degree, self.order)?;

        Ok(Acceleration::<GCRS, AccelerationUnit>::new(a[0], a[1], a[2]))
    }

    // Analytic partials are not yet implemented; the default impl returns an error.
    fn partials(
        &self,
        _s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<ForcePartials<GCRS>, DynamicsError> {
        Err(DynamicsError::Provider(Box::new(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "Geopotential does not yet provide analytic partials",
        ))))
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::errors::DynamicsError;
    use crate::astro::dynamics::forces::j2::J2;
    use crate::astro::dynamics::forces::traits::ForceModel;
    use crate::astro::dynamics::gravity::LowDegreeEarth;
    use crate::astro::dynamics::state::OrbitState;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::time::JulianDate;
    use std::sync::Arc;

    fn make_state(x: f64, y: f64, z: f64) -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(x, y, z),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    fn ctx_with_gravity() -> DynamicsContext {
        DynamicsContext::builder()
            .with_gravity(Arc::new(LowDegreeEarth))
            .build()
    }

    /// Geopotential(2, 0) should match J2::earth() at several positions.
    #[test]
    fn geopotential_n2m0_matches_j2() {
        let ctx = ctx_with_gravity();
        let geopot = Geopotential::new(2, 0);

        let positions: &[(f64, f64, f64)] = &[
            (7_000.0, 0.0, 0.0),
            (0.0, 6_800.0, 1_500.0),
            (4_500.0, 3_000.0, 3_500.0),
        ];

        for &(x, y, z) in positions {
            let s = make_state(x, y, z);

            let a_geo = geopot.acceleration(&s, &ctx).unwrap();
            let a_j2 = J2::earth().acceleration(&s, &DynamicsContext::empty()).unwrap();

            // Note: Geopotential at (2, 0) includes n=0 (two-body) +  n=2,m=0 (J2 only).
            // J2::earth() gives only the J2 perturbation.  We compare the J2
            // portion by subtracting the two-body portion from Geopotential.
            // Easier: compare just the J2 delta vs Geopotential delta over two-body.
            //
            // Alternative: compare Geopotential(2,0) total against TwoBody + J2 sum.
            use crate::astro::dynamics::forces::two_body::TwoBody;
            let a_tb = TwoBody::earth().acceleration(&s, &DynamicsContext::empty()).unwrap();
            let sum_x = a_tb.x().value() + a_j2.x().value();
            let sum_y = a_tb.y().value() + a_j2.y().value();
            let sum_z = a_tb.z().value() + a_j2.z().value();

            let geo_x = a_geo.x().value();
            let geo_y = a_geo.y().value();
            let geo_z = a_geo.z().value();

            for (got, want, label) in [
                (geo_x, sum_x, "x"),
                (geo_y, sum_y, "y"),
                (geo_z, sum_z, "z"),
            ] {
                let denom = got.abs().max(want.abs()).max(1e-30);
                let rel = (got - want).abs() / denom;
                assert!(
                    rel < 5e-9,
                    "pos=({x},{y},{z}) {label}: geo={got:.9e}, tb+j2={want:.9e}, rel={rel:.2e}"
                );
            }
        }
    }

    /// Without a gravity provider in context, returns GravityFieldUnavailable.
    #[test]
    fn geopotential_no_provider_returns_error() {
        let ctx = DynamicsContext::empty();
        let geopot = Geopotential::new(2, 2);
        let s = make_state(7_000.0, 0.0, 0.0);
        let result = geopot.acceleration(&s, &ctx);
        assert!(
            matches!(result, Err(DynamicsError::GravityFieldUnavailable)),
            "expected GravityFieldUnavailable, got {result:?}"
        );
    }

    /// Requesting degree > provider max returns GeopotentialDegreeOutOfRange.
    #[test]
    fn geopotential_degree_out_of_range() {
        let ctx = ctx_with_gravity();
        let geopot = Geopotential::new(10, 10); // LowDegreeEarth max = 4
        let s = make_state(7_000.0, 0.0, 0.0);
        let result = geopot.acceleration(&s, &ctx);
        assert!(
            matches!(
                result,
                Err(DynamicsError::GeopotentialDegreeOutOfRange { requested: 10, max: 4 })
            ),
            "expected GeopotentialDegreeOutOfRange(10, 4), got {result:?}"
        );
    }
}
