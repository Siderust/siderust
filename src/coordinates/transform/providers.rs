// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Transformation Provider Traits
//!
//! This module defines the provider traits for computing time-dependent
//! coordinate transformations. Providers are implemented for specific
//! frame/center pairs and return `affn` operators.
//!
//! ## Architecture
//!
//! The transformation system uses a "hub-and-spoke" model to avoid
//! combinatorial explosion of implementations:
//!
//! - **Frame Hub**: ICRS is the canonical inertial frame. All frame rotations
//!   are computed via ICRS: `F1 → ICRS → F2`.
//!
//! - **Center Hub**: Barycentric is the canonical origin. All center shifts
//!   are computed via Barycentric: `C1 → Barycentric → C2`.
//!
//! ## Provider Traits
//!
//! - [`FrameRotationProvider`]: Computes rotation matrices between frames.
//! - [`CenterShiftProvider`]: Computes translation vectors between centers.
//!
//! ## Sign Conventions
//!
//! ### Frame Rotations
//! - `rotation(F1 → F2)` transforms a vector FROM F1 TO F2.
//! - Applied as: `v_F2 = R * v_F1`
//!
//! ### Center Shifts
//! - `shift(C1 → C2)` is the translation to apply when changing origin.
//! - The shift vector represents the position of C1 as seen from C2.
//! - Applied as: `p_C2 = p_C1 + shift` (shifts the point away from C1 towards C2).
//! - Equivalently: `shift = pos(C1, C2_frame) = -pos(C2, C1_frame)`
//!
//! ## Example
//!
//! ```rust,ignore
//! use siderust::coordinates::transform::providers::*;
//! use siderust::coordinates::frames::{EclipticMeanJ2000, ICRS};
//! use siderust::time::JulianDate;
//! use affn::Rotation3;
//!
//! // Get the rotation from ICRS to EclipticMeanJ2000 at J2000
//! let rot: Rotation3 = FrameRotationProvider::<ICRS, EclipticMeanJ2000>::rotation(
//!     JulianDate::J2000,
//!     &AstroContext::default(),
//! );
//! ```

use crate::astro::eop::EopProvider;
use crate::calculus::ephemeris::Ephemeris;
use crate::coordinates::transform::context::AstroContext;
use crate::time::JulianDate;
use affn::Rotation3;

// =============================================================================
// Frame Rotation Provider
// =============================================================================

/// Trait for computing rotation matrices between reference frames.
///
/// Implementations provide the time-dependent rotation from frame `F1` to
/// frame `F2`. The rotation is computed at a given Julian Date using the
/// provided astronomical context.
///
/// # Type Parameters
///
/// - `F1`: Source reference frame.
/// - `F2`: Target reference frame.
///
/// # Sign Convention
///
/// The returned rotation transforms vectors FROM `F1` TO `F2`:
/// ```text
/// v_F2 = rotation(F1 → F2) * v_F1
/// ```
pub trait FrameRotationProvider<F1, F2> {
    /// Computes the rotation matrix from frame `F1` to frame `F2`.
    ///
    /// # Arguments
    ///
    /// - `jd`: The Julian Date at which to compute the rotation (TT for precession/nutation).
    /// - `ctx`: The astronomical context with model configuration.
    ///
    /// # Returns
    ///
    /// A `Rotation3` that transforms vectors from `F1` to `F2`.
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3;
}

// =============================================================================
// Center Shift Provider
// =============================================================================

/// Trait for computing translation vectors between reference centers.
///
/// Implementations provide the time-dependent translation from center `C1`
/// to center `C2`, expressed in frame `F`. The translation is computed at
/// a given Julian Date using the provided astronomical context.
///
/// # Type Parameters
///
/// - `C1`: Source reference center.
/// - `C2`: Target reference center.
/// - `F`: The frame in which the translation is expressed.
///
/// # Sign Convention
///
/// The returned translation vector represents how to transform a position
/// from center `C1` to center `C2`:
/// ```text
/// p_C2 = p_C1 + shift(C1 → C2, F)
/// ```
///
/// Where `shift(C1 → C2)` is the position of origin C1 as seen from C2,
/// expressed in frame F. This follows the convention:
/// - If C1 is at distance `d` from C2 in direction `u`, then `shift = d * u`.
/// - A point at the C1 origin (p_C1 = 0) becomes `p_C2 = shift`.
///
/// ## Example
///
/// To convert from Heliocentric to Barycentric:
/// - The Sun is at some position relative to the barycenter.
/// - `shift(Helio → Bary) = -pos(Sun in Barycentric) = pos(Bary in Heliocentric)`
/// - A point at the Sun (p_Helio = 0) should map to the Sun's position in Barycentric.
///
/// Wait, let's be precise:
/// - Let `S` = position of Sun in Barycentric frame.
/// - A point P with heliocentric position `p_H` has barycentric position `p_B = p_H + S`.
/// - So `shift(Helio → Bary) = S = pos(Sun in Barycentric)`.
pub trait CenterShiftProvider<C1, C2, F> {
    /// Computes the translation vector from center `C1` to center `C2`.
    ///
    /// # Arguments
    ///
    /// - `jd`: The Julian Date at which to compute the translation.
    /// - `ctx`: The astronomical context with ephemeris configuration.
    ///
    /// # Returns
    ///
    /// A 3-element array `[x, y, z]` representing the shift in frame `F`,
    /// in astronomical units (AU). The caller should convert to the
    /// appropriate unit if needed.
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3];
}

// =============================================================================
// Identity Implementations
// =============================================================================

use crate::astro::nutation;
use crate::astro::precession;
use crate::coordinates::centers::{Barycentric, Geocentric, Heliocentric};
use crate::coordinates::frames::{
    EclipticMeanJ2000, EquatorialMeanJ2000, EquatorialMeanOfDate, EquatorialTrueOfDate, ICRF, ICRS,
};

/// Identity rotation: same frame to same frame.
impl<F> FrameRotationProvider<F, F> for ()
where
    F: affn::ReferenceFrame,
{
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::IDENTITY
    }
}

/// Identity shift: same center to same center.
impl<C, F> CenterShiftProvider<C, C, F> for ()
where
    C: affn::ReferenceCenter,
    F: affn::ReferenceFrame,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        [0.0, 0.0, 0.0]
    }
}

// =============================================================================
// Frame Rotation Implementations (Hub: ICRS)
// =============================================================================

/// Mean obliquity ε₀ at J2000.0 (IAU 2006): 84381.406″, as a typed Radians quantity.
///
/// This is a constant used for the J2000 ecliptic frame, which is fixed
/// and does not vary with epoch.
#[inline]
fn j2000_obliquity() -> qtty::Radians {
    // 84381.406 arcseconds → radians = 84381.406 * π / 648000
    qtty::Radians::new(precession::J2000_MEAN_OBLIQUITY_ARCSEC * std::f64::consts::PI / 648000.0)
}

/// Frame bias rotation from ICRS to mean equator/equinox of J2000.0.
///
/// Values from IERS Conventions (2003), expressed as a small rotation matrix.
const FRAME_BIAS_ICRS_TO_J2000: Rotation3 = Rotation3::from_matrix([
    [
        0.999_999_999_999_994_2,
        0.000_000_070_782_794_8,
        -0.000_000_080_562_171_5,
    ],
    [
        -0.000_000_070_782_797_4,
        0.999_999_999_999_996_9,
        -0.000_000_033_060_408_8,
    ],
    [
        0.000_000_080_562_169_6,
        0.000_000_033_060_414_5,
        0.999_999_999_999_993_2,
    ],
]);

/// ICRS → EclipticMeanJ2000 rotation (J2000 mean ecliptic).
///
/// This composes the ICRS → J2000 mean equator bias with the J2000 obliquity.
/// The `EclipticMeanJ2000` frame represents the ecliptic plane of J2000 (fixed),
/// consistent with VSOP87 ephemeris data.
impl FrameRotationProvider<ICRS, EclipticMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        let bias = FRAME_BIAS_ICRS_TO_J2000;
        let mean_eq_to_ecl = Rotation3::rx(-j2000_obliquity());
        mean_eq_to_ecl * bias
    }
}

/// EclipticMeanJ2000 → ICRS rotation.
impl FrameRotationProvider<EclipticMeanJ2000, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EclipticMeanJ2000>>::rotation(jd, ctx).inverse()
    }
}

/// ICRS → EquatorialMeanJ2000 rotation (frame bias).
impl FrameRotationProvider<ICRS, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        FRAME_BIAS_ICRS_TO_J2000
    }
}

/// EquatorialMeanJ2000 → ICRS rotation.
impl FrameRotationProvider<EquatorialMeanJ2000, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        FRAME_BIAS_ICRS_TO_J2000.inverse()
    }
}

/// EquatorialMeanJ2000 → EclipticMeanJ2000 rotation (J2000 obliquity only).
impl FrameRotationProvider<EquatorialMeanJ2000, EclipticMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::rx(-j2000_obliquity())
    }
}

/// EclipticMeanJ2000 → EquatorialMeanJ2000 rotation.
impl FrameRotationProvider<EclipticMeanJ2000, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::rx(j2000_obliquity())
    }
}

/// EquatorialMeanJ2000 → EquatorialMeanOfDate (IAU 2006 precession).
impl FrameRotationProvider<EquatorialMeanJ2000, EquatorialMeanOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        precession::precession_matrix_iau2006(jd)
    }
}

/// EquatorialMeanOfDate → EquatorialMeanJ2000 (inverse IAU 2006 precession).
impl FrameRotationProvider<EquatorialMeanOfDate, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        precession::precession_matrix_iau2006(jd).inverse()
    }
}

/// EquatorialMeanOfDate → EquatorialTrueOfDate (nutation, IAU 2000B).
impl FrameRotationProvider<EquatorialMeanOfDate, EquatorialTrueOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        // IAU 2000B nutation matrix (ERFA eraNumat convention):
        //   N = Rx(ε+Δε) · Rz(Δψ) · Rx(−ε)
        // Fused 3-rotation constructor: ~17% faster than sequential composition
        let nut = nutation::nutation_iau2000b(jd);
        let eps = nut.mean_obliquity;
        let dpsi = nut.dpsi;
        let deps = nut.deps;
        Rotation3::fused_rx_rz_rx(eps + deps, dpsi, -eps)
    }
}

/// EquatorialTrueOfDate → EquatorialMeanOfDate (inverse nutation, IAU 2000B).
impl FrameRotationProvider<EquatorialTrueOfDate, EquatorialMeanOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EquatorialMeanOfDate, EquatorialTrueOfDate>>::rotation(
            jd, _ctx,
        )
        .inverse()
    }
}

/// EquatorialMeanJ2000 → EquatorialTrueOfDate (IAU 2006 precession + IAU 2000B nutation).
impl FrameRotationProvider<EquatorialMeanJ2000, EquatorialTrueOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        let prec = precession::precession_matrix_iau2006(jd);
        let nut =
            <() as FrameRotationProvider<EquatorialMeanOfDate, EquatorialTrueOfDate>>::rotation(
                jd, _ctx,
            );
        nut * prec
    }
}

/// EquatorialTrueOfDate → EquatorialMeanJ2000 (inverse of precession + nutation).
impl FrameRotationProvider<EquatorialTrueOfDate, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EquatorialMeanJ2000, EquatorialTrueOfDate>>::rotation(jd, _ctx)
            .inverse()
    }
}

/// ICRS → EquatorialMeanOfDate (frame bias + precession).
impl FrameRotationProvider<ICRS, EquatorialMeanOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        let r1: Rotation3 =
            <() as FrameRotationProvider<ICRS, EquatorialMeanJ2000>>::rotation(jd, ctx);
        let r2: Rotation3 = <() as FrameRotationProvider<
            EquatorialMeanJ2000,
            EquatorialMeanOfDate,
        >>::rotation(jd, ctx);
        r2 * r1
    }
}

/// EquatorialMeanOfDate → ICRS.
impl FrameRotationProvider<EquatorialMeanOfDate, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EquatorialMeanOfDate>>::rotation(jd, ctx).inverse()
    }
}

/// ICRS → EquatorialTrueOfDate (IAU 2006/2000B NPB matrix).
///
/// Uses the combined Fukushima-Williams precession-nutation matrix
/// which includes frame bias, IAU 2006 precession, and IAU 2000B nutation.
/// Equivalent to ERFA's `eraPnm06a` (but with IAU 2000B instead of 2000A).
///
/// When the context's EOP provider supplies non-zero celestial pole
/// offsets (`dX`, `dY`), they are folded into the nutation corrections
/// before building the NPB matrix, improving CIP accuracy to ~0.3 mas.
impl FrameRotationProvider<ICRS, EquatorialTrueOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        let nut = nutation::nutation_iau2000b(jd);
        let mut dpsi = nut.dpsi;
        let mut deps = nut.deps;

        // Apply IERS EOP celestial pole offsets (dX, dY → dψ, dε).
        //   dψ_eop = dX / sin(ε_A)
        //   dε_eop = dY
        // Reference: IERS Conventions (2010), §5.5.4
        let eop = ctx.eop_at(jd);
        let dx_rad = qtty::Radians::from(eop.dx);
        let dy_rad = qtty::Radians::from(eop.dy);
        if dx_rad.value() != 0.0 || dy_rad.value() != 0.0 {
            let sin_eps = nut.mean_obliquity.sin();
            if sin_eps.abs() > 1e-15 {
                dpsi += qtty::Radians::new(dx_rad.value() / sin_eps);
            }
            deps += dy_rad;
        }

        precession::precession_nutation_matrix(jd, dpsi, deps)
    }
}

/// EquatorialTrueOfDate → ICRS.
impl FrameRotationProvider<EquatorialTrueOfDate, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EquatorialTrueOfDate>>::rotation(jd, ctx).inverse()
    }
}

/// ICRF → ICRS rotation.
///
/// ICRF is the physical realization of ICRS and is treated as coincident
/// in this provider layer.
impl FrameRotationProvider<ICRF, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::IDENTITY
    }
}

/// ICRS → ICRF rotation.
impl FrameRotationProvider<ICRS, ICRF> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRF, ICRS>>::rotation(jd, ctx).inverse()
    }
}

/// ICRF → EquatorialMeanJ2000 rotation.
impl FrameRotationProvider<ICRF, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EquatorialMeanJ2000>>::rotation(jd, ctx)
    }
}

/// EquatorialMeanJ2000 → ICRF rotation.
impl FrameRotationProvider<EquatorialMeanJ2000, ICRF> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRF, EquatorialMeanJ2000>>::rotation(jd, ctx).inverse()
    }
}

/// ICRF → EclipticMeanJ2000 rotation.
impl FrameRotationProvider<ICRF, EclipticMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EclipticMeanJ2000>>::rotation(jd, ctx)
    }
}

/// EclipticMeanJ2000 → ICRF rotation.
impl FrameRotationProvider<EclipticMeanJ2000, ICRF> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRF, EclipticMeanJ2000>>::rotation(jd, ctx).inverse()
    }
}

/// ICRF → EquatorialMeanOfDate rotation.
impl FrameRotationProvider<ICRF, EquatorialMeanOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EquatorialMeanOfDate>>::rotation(jd, ctx)
    }
}

/// EquatorialMeanOfDate → ICRF rotation.
impl FrameRotationProvider<EquatorialMeanOfDate, ICRF> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRF, EquatorialMeanOfDate>>::rotation(jd, ctx).inverse()
    }
}

/// ICRF → EquatorialTrueOfDate rotation.
impl FrameRotationProvider<ICRF, EquatorialTrueOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EquatorialTrueOfDate>>::rotation(jd, ctx)
    }
}

/// EquatorialTrueOfDate → ICRF rotation.
impl FrameRotationProvider<EquatorialTrueOfDate, ICRF> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRF, EquatorialTrueOfDate>>::rotation(jd, ctx).inverse()
    }
}

// =============================================================================
// Center Shift Implementations (Hub: Barycentric)
// =============================================================================

/// Heliocentric → Barycentric shift.
///
/// Returns the position of the Sun in barycentric coordinates.
/// Uses the active ephemeris backend for the Sun's position.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Heliocentric, Barycentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        // Get Sun's position in barycentric ecliptic coordinates
        let sun_bary = Eph::sun_barycentric(jd);

        // The shift is the Sun's position (in AU)
        [
            sun_bary.x().value(),
            sun_bary.y().value(),
            sun_bary.z().value(),
        ]
    }
}

/// Barycentric → Heliocentric shift.
///
/// Returns the negation of the Sun's barycentric position.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Barycentric, Heliocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let [x, y, z] = <() as CenterShiftProvider<Heliocentric, Barycentric, F>>::shift(jd, ctx);
        [-x, -y, -z]
    }
}

/// Geocentric → Barycentric shift.
///
/// Returns the position of the Earth in barycentric coordinates.
/// Uses the active ephemeris backend for the Earth's position.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Geocentric, Barycentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        // Get Earth's position in barycentric ecliptic coordinates
        let earth_bary = Eph::earth_barycentric(jd);

        // The shift is the Earth's position (in AU)
        [
            earth_bary.x().value(),
            earth_bary.y().value(),
            earth_bary.z().value(),
        ]
    }
}

/// Barycentric → Geocentric shift.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Barycentric, Geocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let [x, y, z] = <() as CenterShiftProvider<Geocentric, Barycentric, F>>::shift(jd, ctx);
        [-x, -y, -z]
    }
}

/// Heliocentric → Geocentric shift.
///
/// Composed via Barycentric hub.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Heliocentric, Geocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        // Helio → Bary → Geo
        let [x1, y1, z1] =
            <() as CenterShiftProvider<Heliocentric, Barycentric, F>>::shift(jd, ctx);
        let [x2, y2, z2] = <() as CenterShiftProvider<Barycentric, Geocentric, F>>::shift(jd, ctx);
        [x1 + x2, y1 + y2, z1 + z2]
    }
}

/// Geocentric → Heliocentric shift.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Geocentric, Heliocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let [x, y, z] = <() as CenterShiftProvider<Heliocentric, Geocentric, F>>::shift(jd, ctx);
        [-x, -y, -z]
    }
}

// =============================================================================
// GCRS ↔ ICRS Frame Rotation
// =============================================================================

use crate::coordinates::frames::GCRS as GCRSFrame;

/// GCRS → ICRS rotation.
///
/// The GCRS is kinematically non-rotating with respect to ICRS. At the
/// precision level of this library, the axes coincide (< 1 mas).
/// This is the identity rotation.
impl FrameRotationProvider<GCRSFrame, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::IDENTITY
    }
}

/// ICRS → GCRS rotation.
impl FrameRotationProvider<ICRS, GCRSFrame> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::IDENTITY
    }
}

/// GCRS → EquatorialMeanJ2000 rotation (via ICRS).
impl FrameRotationProvider<GCRSFrame, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EquatorialMeanJ2000>>::rotation(jd, ctx)
    }
}

/// EquatorialMeanJ2000 → GCRS rotation.
impl FrameRotationProvider<EquatorialMeanJ2000, GCRSFrame> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EquatorialMeanJ2000, ICRS>>::rotation(jd, ctx)
    }
}

/// GCRS → EquatorialTrueOfDate rotation.
impl FrameRotationProvider<GCRSFrame, EquatorialTrueOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EquatorialTrueOfDate>>::rotation(jd, ctx)
    }
}

/// EquatorialTrueOfDate → GCRS rotation.
impl FrameRotationProvider<EquatorialTrueOfDate, GCRSFrame> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EquatorialTrueOfDate, ICRS>>::rotation(jd, ctx)
    }
}

/// GCRS → EclipticMeanJ2000 rotation.
impl FrameRotationProvider<GCRSFrame, EclipticMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EclipticMeanJ2000>>::rotation(jd, ctx)
    }
}

/// EclipticMeanJ2000 → GCRS rotation.
impl FrameRotationProvider<EclipticMeanJ2000, GCRSFrame> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EclipticMeanJ2000, ICRS>>::rotation(jd, ctx)
    }
}

// =============================================================================
// Galactic ↔ ICRS Frame Rotation (IAU 1958)
// =============================================================================

use crate::coordinates::frames::Galactic;

/// Galactic → ICRS rotation matrix (constant).
///
/// This rotation matrix converts vectors from the Galactic frame (IAU 1958)
/// to ICRS. It is based on the equatorial coordinates of the Galactic poles
/// and the Galactic center:
///
/// - North Galactic Pole (ICRS): α = 192.85948°, δ = 27.12825°
/// - Galactic center direction: l = 0°
/// - Position angle of the Galactic equator on the ICRS equator: 122.93192°
///
/// The matrix is R = Rz(αG) · Ry(90° − δG) · Rz(ΘG) where
/// αG, δG are NGP coordinates and ΘG = 122.93192° is the position angle.
///
/// # References
///
/// * Murray, C.A. (1989). "The transformation of coordinates between the
///   systems of B1950 and J2000", Astronomy & Astrophysics, 218, 325-329.
/// * Liu, J.C. et al. (2011). "Reconsidering the Galactic coordinate system",
///   Astronomy & Astrophysics, 536, A102.
// Pre-computed rotation matrix: Galactic → ICRS
// This is the transpose of the equatorial-to-galactic matrix AG from
// Hipparcos Volume 1, Section 1.5.3 (ESA 1997).
const GALACTIC_TO_ICRS: Rotation3 = Rotation3::from_matrix([
    [
        -0.054_875_560_416_215_4,
        -0.873_437_090_234_885_1,
        -0.483_835_015_548_713_2,
    ],
    [
        0.494_109_427_875_583_7,
        -0.444_829_629_960_011_2,
        0.746_982_244_580_286_6,
    ],
    [
        -0.867_666_149_019_004_7,
        -0.198_076_373_431_201_5,
        0.455_983_776_175_066_9,
    ],
]);

/// Galactic → ICRS rotation.
impl FrameRotationProvider<Galactic, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        GALACTIC_TO_ICRS
    }
}

/// ICRS → Galactic rotation.
impl FrameRotationProvider<ICRS, Galactic> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        GALACTIC_TO_ICRS.inverse()
    }
}

/// Galactic → EquatorialMeanJ2000 rotation (via ICRS).
impl FrameRotationProvider<Galactic, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        let r1 = <() as FrameRotationProvider<Galactic, ICRS>>::rotation(jd, ctx);
        let r2 = <() as FrameRotationProvider<ICRS, EquatorialMeanJ2000>>::rotation(jd, ctx);
        r2 * r1
    }
}

/// EquatorialMeanJ2000 → Galactic rotation.
impl FrameRotationProvider<EquatorialMeanJ2000, Galactic> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<Galactic, EquatorialMeanJ2000>>::rotation(jd, ctx).inverse()
    }
}

/// Galactic → EclipticMeanJ2000 rotation (via ICRS).
impl FrameRotationProvider<Galactic, EclipticMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        let r1 = <() as FrameRotationProvider<Galactic, ICRS>>::rotation(jd, ctx);
        let r2 = <() as FrameRotationProvider<ICRS, EclipticMeanJ2000>>::rotation(jd, ctx);
        r2 * r1
    }
}

/// EclipticMeanJ2000 → Galactic rotation.
impl FrameRotationProvider<EclipticMeanJ2000, Galactic> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<Galactic, EclipticMeanJ2000>>::rotation(jd, ctx).inverse()
    }
}

// =============================================================================
// FK4 B1950 ↔ ICRS Frame Rotation
// =============================================================================

use crate::coordinates::frames::FK4B1950;

/// FK4 B1950 → ICRS rotation matrix (constant).
///
/// This is the combined rotation from FK4/B1950 to ICRS, accounting for:
/// 1. The FK4 → FK5/J2000 rotation (Standish 1982, Aoki et al. 1983)
/// 2. The FK5/J2000 → ICRS frame bias
///
/// E-terms of aberration are NOT removed by this matrix; they should be
/// handled separately when converting catalog positions.
///
/// The matrix elements are derived from the standard Standish (1982)
/// rotation angles:
///   ε₀ = −0.525″ (η), z₀ = −0.552″ (ξ), θ₀ = −1.145″ (da)
/// combined with the frame bias to ICRS.
///
/// # References
///
/// * Standish, E.M. (1982). A&A 115, 20–22.
/// * Aoki et al. (1983). A&A 128, 263–267.
/// * Kaplan et al. (1989). AJ 97, 1197.
// Pre-computed FK4 B1950 → FK5 J2000 rotation matrix from ERFA eraFk524.
// This is the Standish (1982) rotation without E-term removal.
// Note: FK5/J2000 ≈ ICRS to within 20 mas; the ICRS frame bias is composed
// at runtime in the FK4→ICRS provider below.
const FK4_TO_FK5: Rotation3 = Rotation3::from_matrix([
    [
        0.999_925_679_495_687_7,
        -0.011_181_483_220_466_2,
        -0.004_859_003_815_359_2,
    ],
    [
        0.011_181_483_239_171_7,
        0.999_937_484_893_313_5,
        -0.000_027_162_594_714_2,
    ],
    [
        0.004_859_003_772_314_3,
        -0.000_027_170_293_744_0,
        0.999_988_194_602_374_2,
    ],
]);

/// FK4B1950 → ICRS rotation (FK4→FK5 then FK5→ICRS via frame bias inverse).
impl FrameRotationProvider<FK4B1950, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        // FK4→FK5(J2000) then J2000→ICRS (bias inverse)
        FRAME_BIAS_ICRS_TO_J2000.inverse() * FK4_TO_FK5
    }
}

/// ICRS → FK4B1950 rotation.
impl FrameRotationProvider<ICRS, FK4B1950> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        FK4_TO_FK5.inverse() * FRAME_BIAS_ICRS_TO_J2000
    }
}

/// FK4B1950 → EquatorialMeanJ2000 rotation (direct Standish matrix).
impl FrameRotationProvider<FK4B1950, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        FK4_TO_FK5
    }
}

/// EquatorialMeanJ2000 → FK4B1950 rotation.
impl FrameRotationProvider<EquatorialMeanJ2000, FK4B1950> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<FK4B1950, EquatorialMeanJ2000>>::rotation(jd, ctx).inverse()
    }
}

// =============================================================================
// TEME ↔ EquatorialTrueOfDate (TOD) Frame Rotation
// =============================================================================

use crate::coordinates::frames::TEME;

/// TEME → EquatorialTrueOfDate (TOD) rotation.
///
/// TEME and TOD share the same pole (true celestial pole / CIP), but
/// differ in the origin of right ascension:
/// - TOD uses the **true equinox** (precession + nutation in longitude).
/// - TEME uses the **mean equinox** (no nutation in longitude applied).
///
/// The rotation from TEME to TOD is therefore:
/// ```text
/// R(TEME → TOD) = Rz(Δψ · cos(ε_A))
/// ```
/// where Δψ is the nutation in longitude and ε_A is the mean obliquity.
/// The quantity Δψ·cos(ε_A) is the **equation of the equinoxes**.
///
/// # References
///
/// * Vallado et al. (2006), AIAA 2006-6753, §3.
/// * IERS Conventions (2010), §5.5.
impl FrameRotationProvider<TEME, EquatorialTrueOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        let nut = nutation::nutation_iau2000b(jd);
        let eq_eq = nut.dpsi * nut.mean_obliquity.cos();
        Rotation3::rz(eq_eq)
    }
}

/// EquatorialTrueOfDate (TOD) → TEME rotation.
impl FrameRotationProvider<EquatorialTrueOfDate, TEME> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<TEME, EquatorialTrueOfDate>>::rotation(jd, ctx).inverse()
    }
}

/// TEME → ICRS rotation (via TOD → EquatorialMeanJ2000 → ICRS).
impl FrameRotationProvider<TEME, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        let teme_to_tod =
            <() as FrameRotationProvider<TEME, EquatorialTrueOfDate>>::rotation(jd, ctx);
        let tod_to_icrs =
            <() as FrameRotationProvider<EquatorialTrueOfDate, ICRS>>::rotation(jd, ctx);
        tod_to_icrs * teme_to_tod
    }
}

/// ICRS → TEME rotation.
impl FrameRotationProvider<ICRS, TEME> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<TEME, ICRS>>::rotation(jd, ctx).inverse()
    }
}

// =============================================================================
// Planetary Body-Fixed ↔ ICRS Frame Rotations
// =============================================================================

use crate::astro::{HasIauRotation, IauRotationParams};
use crate::coordinates::frames::planetary::{
    JupiterSystemIII, MarsFixed, MercuryFixed, MoonPrincipalAxes, NeptuneFixed, PlutoFixed,
    SaturnFixed, UranusFixed, VenusFixed,
};

/// Computes the body-fixed → ICRS rotation matrix from IAU rotation parameters.
///
/// The transformation from a body-fixed frame to ICRS is:
/// ```text
/// R = Rz(−(α₀ + 90°)) · Rx(−(90° − δ₀)) · Rz(−W)
/// ```
///
/// This implements the standard IAU convention where:
/// - α₀, δ₀ define the north pole direction in ICRS
/// - W is the prime meridian angle
fn iau_body_fixed_to_icrs(params: &IauRotationParams, jd: JulianDate) -> Rotation3 {
    let alpha0 = params.alpha0(jd).to::<qtty::Radian>();
    let delta0 = params.delta0(jd).to::<qtty::Radian>();
    let w = params.w(jd).to::<qtty::Radian>();
    let right_angle = qtty::Degrees::new(90.0).to::<qtty::Radian>();

    // R = Rz(-(alpha0 + 90°)) · Rx(-(90° - delta0)) · Rz(-W)
    // Mapping to affn's fused ZXZ constructor:
    // from_euler_zxz(z1, x, z2) = Rz(z2) · Rx(x) · Rz(z1)
    let z1 = -w;
    let x = -(right_angle - delta0);
    let z2 = -(alpha0 + right_angle);
    Rotation3::from_euler_zxz(z1, x, z2)
}

/// Generates FrameRotationProvider impls for a body-fixed frame ↔ ICRS.
macro_rules! impl_body_fixed_rotation {
    ($frame:ty, $body_ty:ty) => {
        impl FrameRotationProvider<$frame, ICRS> for () {
            #[inline]
            fn rotation<Eph, Eop: EopProvider, Nut>(
                jd: JulianDate,
                _ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> Rotation3 {
                iau_body_fixed_to_icrs(&<$body_ty as HasIauRotation>::ROTATION, jd)
            }
        }

        impl FrameRotationProvider<ICRS, $frame> for () {
            #[inline]
            fn rotation<Eph, Eop: EopProvider, Nut>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> Rotation3 {
                <() as FrameRotationProvider<$frame, ICRS>>::rotation(jd, ctx).inverse()
            }
        }
    };
}

impl_body_fixed_rotation!(MercuryFixed, crate::bodies::solar_system::Mercury);
impl_body_fixed_rotation!(VenusFixed, crate::bodies::solar_system::Venus);
impl_body_fixed_rotation!(MarsFixed, crate::bodies::solar_system::Mars);
impl_body_fixed_rotation!(MoonPrincipalAxes, crate::bodies::solar_system::Moon);
impl_body_fixed_rotation!(JupiterSystemIII, crate::bodies::solar_system::Jupiter);
impl_body_fixed_rotation!(SaturnFixed, crate::bodies::solar_system::Saturn);
impl_body_fixed_rotation!(UranusFixed, crate::bodies::solar_system::Uranus);
impl_body_fixed_rotation!(NeptuneFixed, crate::bodies::solar_system::Neptune);
impl_body_fixed_rotation!(PlutoFixed, crate::bodies::solar_system::Pluto);

// =============================================================================
// Planetocentric Center Shift Implementations
// =============================================================================

use crate::coordinates::centers::{
    Jovicentric, Marscentric, Mercurycentric, Neptunocentric, Plutocentric, Saturnocentric,
    Selenocentric, Uranocentric, Venuscentric,
};

/// Generates center shift providers for a planetocentric center via the
/// Barycentric hub. Planet positions are computed directly from VSOP87E
/// barycentric solutions exposed on body marker types.
macro_rules! impl_planet_center_shift_vsop {
    ($center:ty, $planet_ty:path) => {
        /// Planetocentric → Barycentric shift.
        impl<F: affn::ReferenceFrame> CenterShiftProvider<$center, Barycentric, F> for () {
            fn shift<Eph: Ephemeris, Eop, Nut>(
                jd: JulianDate,
                _ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> [f64; 3] {
                let bary_pos = <$planet_ty>::vsop87e(jd);
                [
                    bary_pos.x().value(),
                    bary_pos.y().value(),
                    bary_pos.z().value(),
                ]
            }
        }

        /// Barycentric → Planetocentric shift.
        impl<F: affn::ReferenceFrame> CenterShiftProvider<Barycentric, $center, F> for () {
            fn shift<Eph: Ephemeris, Eop, Nut>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> [f64; 3] {
                let [x, y, z] =
                    <() as CenterShiftProvider<$center, Barycentric, F>>::shift(jd, ctx);
                [-x, -y, -z]
            }
        }

        /// Heliocentric → Planetocentric shift (via Barycentric hub).
        impl<F: affn::ReferenceFrame> CenterShiftProvider<Heliocentric, $center, F> for () {
            fn shift<Eph: Ephemeris, Eop, Nut>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> [f64; 3] {
                let [x1, y1, z1] =
                    <() as CenterShiftProvider<Heliocentric, Barycentric, F>>::shift(jd, ctx);
                let [x2, y2, z2] =
                    <() as CenterShiftProvider<Barycentric, $center, F>>::shift(jd, ctx);
                [x1 + x2, y1 + y2, z1 + z2]
            }
        }

        /// Planetocentric → Heliocentric shift.
        impl<F: affn::ReferenceFrame> CenterShiftProvider<$center, Heliocentric, F> for () {
            fn shift<Eph: Ephemeris, Eop, Nut>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> [f64; 3] {
                let [x, y, z] =
                    <() as CenterShiftProvider<Heliocentric, $center, F>>::shift(jd, ctx);
                [-x, -y, -z]
            }
        }

        /// Geocentric → Planetocentric shift (via Barycentric hub).
        impl<F: affn::ReferenceFrame> CenterShiftProvider<Geocentric, $center, F> for () {
            fn shift<Eph: Ephemeris, Eop, Nut>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> [f64; 3] {
                let [x1, y1, z1] =
                    <() as CenterShiftProvider<Geocentric, Barycentric, F>>::shift(jd, ctx);
                let [x2, y2, z2] =
                    <() as CenterShiftProvider<Barycentric, $center, F>>::shift(jd, ctx);
                [x1 + x2, y1 + y2, z1 + z2]
            }
        }

        /// Planetocentric → Geocentric shift.
        impl<F: affn::ReferenceFrame> CenterShiftProvider<$center, Geocentric, F> for () {
            fn shift<Eph: Ephemeris, Eop, Nut>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> [f64; 3] {
                let [x, y, z] = <() as CenterShiftProvider<Geocentric, $center, F>>::shift(jd, ctx);
                [-x, -y, -z]
            }
        }
    };
}

impl_planet_center_shift_vsop!(Mercurycentric, crate::bodies::solar_system::Mercury);
impl_planet_center_shift_vsop!(Venuscentric, crate::bodies::solar_system::Venus);
impl_planet_center_shift_vsop!(Marscentric, crate::bodies::solar_system::Mars);
impl_planet_center_shift_vsop!(Jovicentric, crate::bodies::solar_system::Jupiter);
impl_planet_center_shift_vsop!(Saturnocentric, crate::bodies::solar_system::Saturn);
impl_planet_center_shift_vsop!(Uranocentric, crate::bodies::solar_system::Uranus);
impl_planet_center_shift_vsop!(Neptunocentric, crate::bodies::solar_system::Neptune);

/// Plutocentric → Barycentric shift.
///
/// Pluto currently uses the orbital elements in `solar_system::PLUTO` to
/// produce heliocentric coordinates, then shifts by the Sun's barycentric
/// offset. This preserves existing behavior until a dedicated
/// `solar_system::Pluto` barycentric ephemeris method is added.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Plutocentric, Barycentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        use crate::bodies::solar_system;

        let helio_pos = solar_system::PLUTO.orbit.kepler_position(jd);
        let [sx, sy, sz] =
            <() as CenterShiftProvider<Heliocentric, Barycentric, F>>::shift(jd, ctx);

        [
            helio_pos.x().value() + sx,
            helio_pos.y().value() + sy,
            helio_pos.z().value() + sz,
        ]
    }
}

/// Barycentric → Plutocentric shift.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Barycentric, Plutocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let [x, y, z] = <() as CenterShiftProvider<Plutocentric, Barycentric, F>>::shift(jd, ctx);
        [-x, -y, -z]
    }
}

/// Heliocentric → Plutocentric shift (via Barycentric hub).
impl<F: affn::ReferenceFrame> CenterShiftProvider<Heliocentric, Plutocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let [x1, y1, z1] =
            <() as CenterShiftProvider<Heliocentric, Barycentric, F>>::shift(jd, ctx);
        let [x2, y2, z2] =
            <() as CenterShiftProvider<Barycentric, Plutocentric, F>>::shift(jd, ctx);
        [x1 + x2, y1 + y2, z1 + z2]
    }
}

/// Plutocentric → Heliocentric shift.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Plutocentric, Heliocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let [x, y, z] = <() as CenterShiftProvider<Heliocentric, Plutocentric, F>>::shift(jd, ctx);
        [-x, -y, -z]
    }
}

/// Geocentric → Plutocentric shift (via Barycentric hub).
impl<F: affn::ReferenceFrame> CenterShiftProvider<Geocentric, Plutocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let [x1, y1, z1] = <() as CenterShiftProvider<Geocentric, Barycentric, F>>::shift(jd, ctx);
        let [x2, y2, z2] =
            <() as CenterShiftProvider<Barycentric, Plutocentric, F>>::shift(jd, ctx);
        [x1 + x2, y1 + y2, z1 + z2]
    }
}

/// Plutocentric → Geocentric shift.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Plutocentric, Geocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let [x, y, z] = <() as CenterShiftProvider<Geocentric, Plutocentric, F>>::shift(jd, ctx);
        [-x, -y, -z]
    }
}

/// Selenocentric (Moon) → Barycentric shift.
///
/// The Moon's position is obtained from the ephemeris (geocentric ecliptic),
/// then combined with Earth's barycentric position to get the Moon's
/// barycentric position.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Selenocentric, Barycentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        // Moon's geocentric position (in km, ecliptic)
        let moon_geo = Eph::moon_geocentric(jd);
        // Convert km to AU
        let km_per_au = 149_597_870.7;
        let moon_x = moon_geo.x().value() / km_per_au;
        let moon_y = moon_geo.y().value() / km_per_au;
        let moon_z = moon_geo.z().value() / km_per_au;

        // Earth's barycentric position
        let [ex, ey, ez] = <() as CenterShiftProvider<Geocentric, Barycentric, F>>::shift(jd, ctx);

        // Moon barycentric = Moon geocentric + Earth barycentric
        [moon_x + ex, moon_y + ey, moon_z + ez]
    }
}

/// Barycentric → Selenocentric shift.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Barycentric, Selenocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let [x, y, z] = <() as CenterShiftProvider<Selenocentric, Barycentric, F>>::shift(jd, ctx);
        [-x, -y, -z]
    }
}

/// Heliocentric → Selenocentric shift (via Barycentric hub).
impl<F: affn::ReferenceFrame> CenterShiftProvider<Heliocentric, Selenocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let [x1, y1, z1] =
            <() as CenterShiftProvider<Heliocentric, Barycentric, F>>::shift(jd, ctx);
        let [x2, y2, z2] =
            <() as CenterShiftProvider<Barycentric, Selenocentric, F>>::shift(jd, ctx);
        [x1 + x2, y1 + y2, z1 + z2]
    }
}

/// Selenocentric → Heliocentric shift.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Selenocentric, Heliocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let [x, y, z] = <() as CenterShiftProvider<Heliocentric, Selenocentric, F>>::shift(jd, ctx);
        [-x, -y, -z]
    }
}

/// Geocentric → Selenocentric shift.
///
/// Direct computation from the ephemeris: the Moon's geocentric position
/// gives us the Earth→Moon vector. The shift from Geo→Seleno is the
/// position of Earth's centre in Selenocentric (Moon-centred) coordinates,
/// which is the negative of the Moon's geocentric position.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Geocentric, Selenocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let moon_geo = Eph::moon_geocentric(jd);
        let km_per_au = 149_597_870.7;
        // Shift = position of Geocentric origin (Earth) in Selenocentric (Moon) frame
        // = -(Moon position relative to Earth)
        [
            -moon_geo.x().value() / km_per_au,
            -moon_geo.y().value() / km_per_au,
            -moon_geo.z().value() / km_per_au,
        ]
    }
}

/// Selenocentric → Geocentric shift.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Selenocentric, Geocentric, F> for () {
    fn shift<Eph: Ephemeris, Eop, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let [x, y, z] = <() as CenterShiftProvider<Geocentric, Selenocentric, F>>::shift(jd, ctx);
        [-x, -y, -z]
    }
}

// =============================================================================
// Convenience Functions
// =============================================================================

/// Computes the rotation matrix from frame `F1` to frame `F2`.
///
/// This is a convenience function that dispatches to the appropriate
/// [`FrameRotationProvider`] implementation.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::transform::providers::frame_rotation;
/// use siderust::coordinates::transform::context::AstroContext;
/// use siderust::coordinates::frames::{ICRS, EclipticMeanJ2000};
/// use siderust::time::JulianDate;
///
/// let rot = frame_rotation::<ICRS, EclipticMeanJ2000>(JulianDate::J2000, &AstroContext::default());
/// ```
#[inline]
pub fn frame_rotation<F1, F2>(jd: JulianDate, ctx: &AstroContext) -> Rotation3
where
    (): FrameRotationProvider<F1, F2>,
{
    <() as FrameRotationProvider<F1, F2>>::rotation(jd, ctx)
}

/// Computes the center shift from `C1` to `C2` in frame `F`.
///
/// Returns the shift as a 3-element array in AU.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::transform::providers::center_shift;
/// use siderust::coordinates::transform::context::AstroContext;
/// use siderust::coordinates::centers::{Heliocentric, Geocentric};
/// use siderust::coordinates::frames::EclipticMeanJ2000;
/// use siderust::time::JulianDate;
///
/// let shift = center_shift::<Heliocentric, Geocentric, EclipticMeanJ2000>(
///     JulianDate::J2000,
///     &AstroContext::default(),
/// );
/// ```
#[inline]
pub fn center_shift<C1, C2, F>(jd: JulianDate, ctx: &AstroContext) -> [f64; 3]
where
    (): CenterShiftProvider<C1, C2, F>,
{
    <() as CenterShiftProvider<C1, C2, F>>::shift(jd, ctx)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-10;

    #[test]
    fn test_identity_frame_rotation() {
        let rot = frame_rotation::<ICRS, ICRS>(JulianDate::J2000, &AstroContext::default());
        let v = [1.0, 2.0, 3.0];
        let result = rot.apply_array(v);
        assert!((result[0] - v[0]).abs() < EPSILON);
        assert!((result[1] - v[1]).abs() < EPSILON);
        assert!((result[2] - v[2]).abs() < EPSILON);
    }

    #[test]
    fn test_icrs_to_ecliptic_rotation() {
        let rot =
            frame_rotation::<ICRS, EclipticMeanJ2000>(JulianDate::J2000, &AstroContext::default());

        // Includes a small frame-bias (ICRS↔J2000), so don't assume a pure X-axis rotation.
        // Instead: verify it behaves like a proper rotation (finite + length-preserving).
        let v = [1.0, 2.0, 3.0];
        let w = rot.apply_array(v);

        assert!(w[0].is_finite() && w[1].is_finite() && w[2].is_finite());

        let nv = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let nw = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
        assert!((nv - nw).abs() < 1e-12);
    }

    #[test]
    fn test_ecliptic_icrs_roundtrip() {
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;

        let r1 = frame_rotation::<ICRS, EclipticMeanJ2000>(jd, &ctx);
        let r2 = frame_rotation::<EclipticMeanJ2000, ICRS>(jd, &ctx);

        let v = [1.0, 2.0, 3.0];
        let roundtrip = r2.apply_array(r1.apply_array(v));

        assert!((roundtrip[0] - v[0]).abs() < EPSILON);
        assert!((roundtrip[1] - v[1]).abs() < EPSILON);
        assert!((roundtrip[2] - v[2]).abs() < EPSILON);
    }

    #[test]
    fn test_identity_center_shift() {
        let shift = center_shift::<Barycentric, Barycentric, EclipticMeanJ2000>(
            JulianDate::J2000,
            &AstroContext::default(),
        );
        assert!((shift[0]).abs() < EPSILON);
        assert!((shift[1]).abs() < EPSILON);
        assert!((shift[2]).abs() < EPSILON);
    }

    #[test]
    fn test_helio_bary_geo_composition() {
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;

        // Helio → Geo should equal Helio → Bary + Bary → Geo
        let helio_geo = center_shift::<Heliocentric, Geocentric, EclipticMeanJ2000>(jd, &ctx);
        let helio_bary = center_shift::<Heliocentric, Barycentric, EclipticMeanJ2000>(jd, &ctx);
        let bary_geo = center_shift::<Barycentric, Geocentric, EclipticMeanJ2000>(jd, &ctx);

        let composed = [
            helio_bary[0] + bary_geo[0],
            helio_bary[1] + bary_geo[1],
            helio_bary[2] + bary_geo[2],
        ];

        assert!((helio_geo[0] - composed[0]).abs() < EPSILON);
        assert!((helio_geo[1] - composed[1]).abs() < EPSILON);
        assert!((helio_geo[2] - composed[2]).abs() < EPSILON);
    }

    #[test]
    fn test_center_shift_antisymmetry() {
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;

        let forward = center_shift::<Heliocentric, Geocentric, EclipticMeanJ2000>(jd, &ctx);
        let backward = center_shift::<Geocentric, Heliocentric, EclipticMeanJ2000>(jd, &ctx);

        assert!((forward[0] + backward[0]).abs() < EPSILON);
        assert!((forward[1] + backward[1]).abs() < EPSILON);
        assert!((forward[2] + backward[2]).abs() < EPSILON);
    }

    #[test]
    fn test_frame_bias_is_non_identity() {
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;

        let rot = frame_rotation::<ICRS, EquatorialMeanJ2000>(jd, &ctx);
        let v = [0.0, 1.0, 0.0];
        let out = rot.apply_array(v);

        let delta = (out[0] - v[0]).abs() + (out[1] - v[1]).abs() + (out[2] - v[2]).abs();
        assert!(delta > 1e-12, "frame bias should not be identity");
    }

    #[test]
    fn test_icrs_ecliptic_roundtrip_is_identity() {
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;

        let r = frame_rotation::<ICRS, EclipticMeanJ2000>(jd, &ctx);
        let rinv = frame_rotation::<EclipticMeanJ2000, ICRS>(jd, &ctx);

        let v = [1.0, 0.0, 0.0];
        let round = rinv.apply_array(r.apply_array(v));
        let err = (round[0] - v[0]).abs() + (round[1] - v[1]).abs() + (round[2] - v[2]).abs();
        assert!(
            err < 1e-12,
            "ICRS↔EclipticMeanJ2000 roundtrip should be identity, got {:?}",
            round
        );
    }

    #[test]
    fn test_precession_identity_at_j2000() {
        // The IAU 2006 precession matrix includes frame bias (ICRS→J2000,
        // ~23 mas ≈ 1.1e-7 rad), so at J2000 it's not exactly identity.
        // Tolerance must accommodate the frame bias offset.
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;

        let rot = frame_rotation::<EquatorialMeanJ2000, EquatorialMeanOfDate>(jd, &ctx);
        let v = [1.0, 0.0, 0.0];
        let out = rot.apply_array(v);

        assert!((out[0] - v[0]).abs() < 1e-6);
        assert!((out[1] - v[1]).abs() < 1e-6);
        assert!((out[2] - v[2]).abs() < 1e-6);
    }

    #[test]
    fn test_nutation_rotation_roundtrip() {
        let ctx = AstroContext::default();
        let jd = JulianDate::new(2_460_000.5);

        let rot = frame_rotation::<EquatorialMeanOfDate, EquatorialTrueOfDate>(jd, &ctx);
        let inv = frame_rotation::<EquatorialTrueOfDate, EquatorialMeanOfDate>(jd, &ctx);

        let v = [0.3, 0.4, 0.5];
        let roundtrip = inv.apply_array(rot.apply_array(v));

        assert!((roundtrip[0] - v[0]).abs() < 1e-12);
        assert!((roundtrip[1] - v[1]).abs() < 1e-12);
        assert!((roundtrip[2] - v[2]).abs() < 1e-12);
    }

    #[test]
    fn test_icrf_icrs_identity_rotation() {
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;
        let rot = frame_rotation::<ICRF, ICRS>(jd, &ctx);
        let v = [0.1, -0.2, 0.3];
        let out = rot.apply_array(v);
        assert!((out[0] - v[0]).abs() < 1e-15);
        assert!((out[1] - v[1]).abs() < 1e-15);
        assert!((out[2] - v[2]).abs() < 1e-15);
    }

    #[test]
    fn test_icrf_ecliptic_roundtrip() {
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;
        let r = frame_rotation::<ICRF, EclipticMeanJ2000>(jd, &ctx);
        let rinv = frame_rotation::<EclipticMeanJ2000, ICRF>(jd, &ctx);

        let v = [1.0, 2.0, 3.0];
        let round = rinv.apply_array(r.apply_array(v));
        assert!((round[0] - v[0]).abs() < 1e-12);
        assert!((round[1] - v[1]).abs() < 1e-12);
        assert!((round[2] - v[2]).abs() < 1e-12);
    }

    #[test]
    fn test_icrf_to_ecliptic_matches_icrs_to_ecliptic() {
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;
        let via_icrf = frame_rotation::<ICRF, EclipticMeanJ2000>(jd, &ctx);
        let via_icrs = frame_rotation::<ICRS, EclipticMeanJ2000>(jd, &ctx);

        let v = [0.3, -0.1, 0.8];
        let a = via_icrf.apply_array(v);
        let b = via_icrs.apply_array(v);
        assert!((a[0] - b[0]).abs() < 1e-15);
        assert!((a[1] - b[1]).abs() < 1e-15);
        assert!((a[2] - b[2]).abs() < 1e-15);
    }
}
