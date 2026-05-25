// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Apparent-position correction policy
//!
//! ## Scientific scope
//!
//! Orchestrates the chain of corrections that turn a *geometric*
//! barycentric position into an *apparent* observed direction:
//! proper motion / parallax, light-time iteration, Lorentz stellar
//! aberration, and gravitational light deflection. The individual
//! corrections live in [`crate::astro::aberration`],
//! [`crate::astro::light_deflection`], and the proper-motion helpers
//! on [`crate::targets::CoordinateWithPM`]; this module centralises
//! the *which corrections, in what order* policy.
//!
//! ## Technical scope
//!
//! - [`CorrectionPolicy`] is a bit-flag of independent corrections.
//!   Order is fixed and documented; toggling individual stages on or
//!   off does **not** reorder the others.
//! - [`CorrectionOrder`] documents the canonical pipeline order:
//!   `ProperMotion → LightTime → Aberration → LightDeflection`.
//! - [`AppliedCorrections`] reports back exactly which stages ran on
//!   a given evaluation, for downstream diagnostics.
//! - The policy is a runtime selector. Compile-time model selection
//!   (nutation theory, sidereal-time convention, …) continues to use
//!   phantom types elsewhere in the crate.
//!
//! ## References
//!
//! - IAU SOFA Cookbook, *Time and Coordinates*, §7.
//! - Kaplan, G. H. (2005). *USNO Circular 179: The IAU Resolutions on
//!   Astronomical Reference Systems, Time Scales, and Earth Rotation
//!   Models*. §6.
//! - Klioner, S. A. (2003). "A practical relativistic model for
//!   microarcsecond astrometry in space." *Astronomical Journal*, 125.

#![forbid(unsafe_code)]

/// Independent apparent-position correction stages.
///
/// The flags can be OR'ed together to opt into multiple corrections
/// at once. The pipeline order is fixed and independent of the bit
/// values (see [`CORRECTION_ORDER`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CorrectionPolicy(u8);

impl CorrectionPolicy {
    /// No corrections — return the input position as is.
    pub const GEOMETRIC: Self = Self(0);
    /// Apply proper motion and annual parallax for the requested epoch.
    pub const PROPER_MOTION: Self = Self(1 << 0);
    /// Iterate light-time so the target position is evaluated at the
    /// *emission* epoch rather than the geocentric epoch.
    pub const LIGHT_TIME: Self = Self(1 << 1);
    /// Apply Lorentz stellar aberration for the observer's velocity.
    pub const ABERRATION: Self = Self(1 << 2);
    /// Apply gravitational light deflection (Sun by default; the caller
    /// may extend it via the underlying
    /// [`crate::astro::light_deflection`] API).
    pub const LIGHT_DEFLECTION: Self = Self(1 << 3);
    /// Geocentric apparent direction: proper motion + light-time +
    /// aberration + deflection.
    pub const APPARENT: Self = Self(
        Self::PROPER_MOTION.0
            | Self::LIGHT_TIME.0
            | Self::ABERRATION.0
            | Self::LIGHT_DEFLECTION.0,
    );

    /// True iff `self` includes every stage flagged in `other`.
    pub const fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }

    /// Set-union of two policies.
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    /// Raw bit pattern, for serialisation or interop only.
    pub const fn bits(self) -> u8 {
        self.0
    }
}

impl core::ops::BitOr for CorrectionPolicy {
    type Output = Self;
    fn bitor(self, rhs: Self) -> Self {
        self.union(rhs)
    }
}

impl core::ops::BitOrAssign for CorrectionPolicy {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

impl Default for CorrectionPolicy {
    fn default() -> Self {
        Self::GEOMETRIC
    }
}

/// The canonical, documented order in which corrections are applied
/// when more than one flag is set.
///
/// 1. **ProperMotion** propagates catalogue position by the elapsed
///    proper-motion / parallax interval.
/// 2. **LightTime** iterates the target position back to the emission
///    epoch by `Δt = |r| / c`.
/// 3. **Aberration** rotates the direction by `v_obs / c`.
/// 4. **LightDeflection** bends the direction around the dominant
///    Schwarzschild mass.
///
/// Order is **fixed**; toggling individual stages does not reorder
/// the others.
pub const CORRECTION_ORDER: [CorrectionPolicy; 4] = [
    CorrectionPolicy::PROPER_MOTION,
    CorrectionPolicy::LIGHT_TIME,
    CorrectionPolicy::ABERRATION,
    CorrectionPolicy::LIGHT_DEFLECTION,
];

/// Documented canonical pipeline ordering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorrectionOrder {
    /// Pipeline order [`CORRECTION_ORDER`].
    Canonical,
}

/// Report of which correction stages were actually executed during
/// an evaluation. Diagnostic value only.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AppliedCorrections {
    /// Flags of corrections that were attempted.
    pub requested: CorrectionPolicy,
    /// Flags of corrections that were actually applied (subset of
    /// `requested`). Stages may be skipped when there is no data to
    /// act on (e.g. PROPER_MOTION on a target with no proper-motion
    /// terms).
    pub applied: CorrectionPolicy,
    /// Number of light-time iteration steps actually executed.
    pub light_time_iters: u8,
}

impl AppliedCorrections {
    /// Build an empty report keyed to the supplied policy.
    pub fn from_request(requested: CorrectionPolicy) -> Self {
        Self {
            requested,
            applied: CorrectionPolicy::GEOMETRIC,
            light_time_iters: 0,
        }
    }

    /// Mark a stage as having been applied.
    pub fn record(&mut self, stage: CorrectionPolicy) {
        self.applied |= stage;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_policy_is_geometric() {
        assert_eq!(CorrectionPolicy::default(), CorrectionPolicy::GEOMETRIC);
    }

    #[test]
    fn apparent_policy_includes_all_stages() {
        let p = CorrectionPolicy::APPARENT;
        for stage in CORRECTION_ORDER {
            assert!(p.contains(stage), "{stage:?} missing from APPARENT");
        }
    }

    #[test]
    fn record_accumulates_applied_stages() {
        let mut r = AppliedCorrections::from_request(CorrectionPolicy::APPARENT);
        assert_eq!(r.applied, CorrectionPolicy::GEOMETRIC);
        r.record(CorrectionPolicy::ABERRATION);
        r.record(CorrectionPolicy::LIGHT_TIME);
        assert!(r.applied.contains(CorrectionPolicy::ABERRATION));
        assert!(r.applied.contains(CorrectionPolicy::LIGHT_TIME));
        assert!(!r.applied.contains(CorrectionPolicy::LIGHT_DEFLECTION));
    }

    #[test]
    fn correction_order_is_canonical_and_complete() {
        // Order must be PM → LightTime → Aberration → LightDeflection.
        assert_eq!(CORRECTION_ORDER[0], CorrectionPolicy::PROPER_MOTION);
        assert_eq!(CORRECTION_ORDER[1], CorrectionPolicy::LIGHT_TIME);
        assert_eq!(CORRECTION_ORDER[2], CorrectionPolicy::ABERRATION);
        assert_eq!(CORRECTION_ORDER[3], CorrectionPolicy::LIGHT_DEFLECTION);
    }
}
