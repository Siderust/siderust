// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! IAU 2006 frame-bias rotation between GCRS/ICRS and EME2000.
//!
//! The three **constant** angles that define the rotation from the
//! GCRS (≡ ICRS for directions) to the mean equator and equinox of
//! J2000.0 (EME2000) are tabulated in IERS Conventions (2010),
//! §5.4.4, Table 5.2b:
//!
//! ```text
//! ξ₀  = −0.0166170″   (offset of CIP at J2000 from mean pole, x)
//! η₀  = −0.0068192″   (offset of CIP at J2000 from mean pole, y)
//! dα₀ = −0.0146″      (offset of ICRS RA origin from mean equinox)
//! ```
//!
//! Using the parametrisation of Eq. (5.32):
//!
//! ```text
//! B = R1(−η₀) · R2(ξ₀) · R3(dα₀)
//! ```
//!
//! where R1, R2, R3 are **passive** elementary rotations.  In `affn`,
//! `Rotation3::r{x,y,z}` are *active*, so `Rᵢ(θ) = Rotation3::r{x,y,z}(−θ)`,
//! yielding:
//!
//! ```text
//! B = Rotation3::rx(η₀) · Rotation3::ry(−ξ₀) · Rotation3::rz(−dα₀)
//! ```
//!
//! The resulting matrix agrees with SOFA `iauBp06` `rb` at J2000.0 to ≲ 1×10⁻¹⁵.
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
use qtty::angular::Radians;

// =============================================================================
// Private constants and helper
// =============================================================================

const FRAME_BIAS_DALPHA0_ARCSEC: f64 = -0.0146;
const FRAME_BIAS_XI0_ARCSEC: f64 = -0.0166170;
const FRAME_BIAS_ETA0_ARCSEC: f64 = -0.0068192;

const ARCSEC_TO_RAD: f64 = std::f64::consts::PI / 648_000.0;

#[inline]
fn frame_bias_gcrs_to_eme2000() -> Rotation3 {
    let xi0 = Radians::new(FRAME_BIAS_XI0_ARCSEC * ARCSEC_TO_RAD);
    let eta0 = Radians::new(FRAME_BIAS_ETA0_ARCSEC * ARCSEC_TO_RAD);
    let da0 = Radians::new(FRAME_BIAS_DALPHA0_ARCSEC * ARCSEC_TO_RAD);
    Rotation3::rx(eta0) * Rotation3::ry(-xi0) * Rotation3::rz(-da0)
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
