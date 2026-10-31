// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! IAU 2006 frame-bias rotation between GCRS/ICRS and EME2000.
//!
//! Delegates to the same SOFA `iauBp06` J2000 matrix used by all coordinate
//! transform providers. This is the IAU 2006 dynamical mean equator/equinox,
//! not the empirical FK5 catalogue orientation returned by `iauFk5hip`.
//! The fixed bias is about 23 milliarcseconds.
//!
//! # Extension traits
//!
//! The bias rotation is accessible through three extension traits:
//!
//! - [`GcrsFrameBias`] — `GCRS::frame_bias_to_eme2000()`
//! - [`IcrsFrameBias`] — `ICRS::frame_bias_to_eme2000()`
//!   (identical to the GCRS version because GCRS and ICRS share axes for directions)
//! - [`Eme2000FrameBias`] — `EME2000::frame_bias_to_gcrs()` and `EME2000::frame_bias_to_icrs()`
//!
//! Import the relevant trait to use these methods:
//!
//! ```rust
//! use affn::frames::{EME2000, GCRS, ICRS};
//! use siderust::astro::frame_bias::{Eme2000FrameBias, GcrsFrameBias, IcrsFrameBias};
//!
//! let b = GCRS::frame_bias_to_eme2000();
//! let inv = EME2000::frame_bias_to_gcrs();
//! // Round-trip
//! let _ = b * inv;
//! ```
//!
//! # References
//!
//! * IERS Conventions (2010), §5.4.4 and Table 5.2b.
//! * Hilton et al., Cel. Mech. Dyn. Astr. 94 (2006) 351 — IAU 2006 (P03).
//! * SOFA `iauBp06`.

use affn::frames::{EME2000, GCRS, ICRS};
use affn::ops::Rotation3;
#[inline]
fn frame_bias_gcrs_to_eme2000() -> Rotation3 {
    crate::coordinates::transform::frames::bias::frame_bias_icrs_to_j2000()
}

// =============================================================================
// Extension traits
// =============================================================================

/// IAU 2006 frame-bias rotation from [`GCRS`] to [`EME2000`].
pub trait GcrsFrameBias {
    /// Constant IAU 2006 frame-bias rotation `B` from [`GCRS`] to [`EME2000`].
    ///
    /// Built from the IERS Conventions (2010) Table 5.2b angles
    ///
    /// ```text
    /// ξ₀  = −0.0166170″
    /// η₀  = −0.0068192″
    /// dα₀ = −0.0146″
    /// ```
    ///
    /// The rotation angle is approximately 23 milli-arcseconds
    /// (≈ 1.1 × 10⁻⁷ rad), epoch-independent.  It agrees with SOFA
    /// `iauBp06` `rb` at J2000.0 to ≲ 1 × 10⁻¹⁵.
    #[must_use]
    fn frame_bias_to_eme2000() -> Rotation3;
}

impl GcrsFrameBias for GCRS {
    #[inline]
    fn frame_bias_to_eme2000() -> Rotation3 {
        frame_bias_gcrs_to_eme2000()
    }
}

/// IAU 2006 frame-bias rotation from [`ICRS`] to [`EME2000`].
///
/// For *direction* purposes `ICRS` and `GCRS` are bit-identical (per IAU 2000
/// Resolution B1.3), so this is the same matrix as [`GcrsFrameBias::frame_bias_to_eme2000`].
pub trait IcrsFrameBias {
    /// Constant IAU 2006 frame-bias rotation `B` from [`ICRS`] to [`EME2000`].
    ///
    /// Identical to `GCRS::frame_bias_to_eme2000()`.  Magnitude ≈ 23 mas.
    #[must_use]
    fn frame_bias_to_eme2000() -> Rotation3;
}

impl IcrsFrameBias for ICRS {
    #[inline]
    fn frame_bias_to_eme2000() -> Rotation3 {
        frame_bias_gcrs_to_eme2000()
    }
}

/// Inverse IAU 2006 frame-bias rotations from [`EME2000`] to GCRS/ICRS.
pub trait Eme2000FrameBias {
    /// Inverse IAU 2006 frame-bias rotation `Bᵀ` from [`EME2000`] to [`GCRS`].
    ///
    /// Exact algebraic inverse of `GCRS::frame_bias_to_eme2000()`.
    #[must_use]
    fn frame_bias_to_gcrs() -> Rotation3;

    /// Inverse IAU 2006 frame-bias rotation from [`EME2000`] to [`ICRS`].
    ///
    /// Identical to `EME2000::frame_bias_to_gcrs()` for direction purposes,
    /// because GCRS and ICRS share axes.
    #[must_use]
    fn frame_bias_to_icrs() -> Rotation3;
}

impl Eme2000FrameBias for EME2000 {
    #[inline]
    fn frame_bias_to_gcrs() -> Rotation3 {
        frame_bias_gcrs_to_eme2000().inverse()
    }

    #[inline]
    fn frame_bias_to_icrs() -> Rotation3 {
        frame_bias_gcrs_to_eme2000().inverse()
    }
}
