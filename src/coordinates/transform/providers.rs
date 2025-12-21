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
//! use siderust::coordinates::frames::{Ecliptic, ICRS};
//! use siderust::astro::JulianDate;
//! use affn::Rotation3;
//!
//! // Get the rotation from ICRS to Ecliptic at J2000
//! let rot: Rotation3 = FrameRotationProvider::<ICRS, Ecliptic>::rotation(
//!     JulianDate::J2000,
//!     &AstroContext::default(),
//! );
//! ```

use crate::astro::JulianDate;
use crate::coordinates::transform::context::AstroContext;
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
    /// - `jd`: The Julian Date at which to compute the rotation.
    /// - `ctx`: The astronomical context with model configuration.
    ///
    /// # Returns
    ///
    /// A `Rotation3` that transforms vectors from `F1` to `F2`.
    fn rotation<Eph, Eop, Nut>(jd: JulianDate, ctx: &AstroContext<Eph, Eop, Nut>) -> Rotation3;
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
    fn shift<Eph, Eop, Nut>(jd: JulianDate, ctx: &AstroContext<Eph, Eop, Nut>) -> [f64; 3];
}

// =============================================================================
// Identity Implementations
// =============================================================================

use crate::coordinates::frames::{Ecliptic, Equatorial, ICRS};
use crate::coordinates::centers::{Barycentric, Geocentric, Heliocentric};

/// Identity rotation: same frame to same frame.
impl<F> FrameRotationProvider<F, F> for ()
where
    F: affn::ReferenceFrame,
{
    #[inline]
    fn rotation<Eph, Eop, Nut>(_jd: JulianDate, _ctx: &AstroContext<Eph, Eop, Nut>) -> Rotation3 {
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
    fn shift<Eph, Eop, Nut>(_jd: JulianDate, _ctx: &AstroContext<Eph, Eop, Nut>) -> [f64; 3] {
        [0.0, 0.0, 0.0]
    }
}

// =============================================================================
// Frame Rotation Implementations (Hub: ICRS)
// =============================================================================

/// Mean obliquity of the ecliptic at J2000.0 (radians).
///
/// This is the angle between the ecliptic and equatorial planes.
/// Value: 23.439291111° = 84381.448" (IAU 1976)
const OBLIQUITY_J2000: f64 = 0.409092804222329; // 23.4392911° in radians

/// ICRS → Ecliptic rotation (J2000 mean ecliptic).
///
/// ICRS is closely aligned with the J2000 equatorial frame, so we use
/// a simple rotation around the X-axis by the obliquity.
impl FrameRotationProvider<ICRS, Ecliptic> for () {
    #[inline]
    fn rotation<Eph, Eop, Nut>(_jd: JulianDate, _ctx: &AstroContext<Eph, Eop, Nut>) -> Rotation3 {
        // Rotate from equatorial (ICRS) to ecliptic around X-axis
        Rotation3::from_x_rotation(OBLIQUITY_J2000)
    }
}

/// Ecliptic → ICRS rotation.
impl FrameRotationProvider<Ecliptic, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop, Nut>(_jd: JulianDate, _ctx: &AstroContext<Eph, Eop, Nut>) -> Rotation3 {
        // Inverse of ICRS → Ecliptic
        Rotation3::from_x_rotation(-OBLIQUITY_J2000)
    }
}

/// ICRS → Equatorial rotation.
///
/// ICRS is very close to J2000 Equatorial (differences are at the
/// milliarcsecond level), so we use identity for now.
/// For high-precision work, a frame bias matrix would be applied.
impl FrameRotationProvider<ICRS, Equatorial> for () {
    #[inline]
    fn rotation<Eph, Eop, Nut>(_jd: JulianDate, _ctx: &AstroContext<Eph, Eop, Nut>) -> Rotation3 {
        // ICRS ≈ J2000 Equatorial (within ~20 mas)
        Rotation3::IDENTITY
    }
}

/// Equatorial → ICRS rotation.
impl FrameRotationProvider<Equatorial, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop, Nut>(_jd: JulianDate, _ctx: &AstroContext<Eph, Eop, Nut>) -> Rotation3 {
        Rotation3::IDENTITY
    }
}

/// Equatorial → Ecliptic rotation.
impl FrameRotationProvider<Equatorial, Ecliptic> for () {
    #[inline]
    fn rotation<Eph, Eop, Nut>(jd: JulianDate, ctx: &AstroContext<Eph, Eop, Nut>) -> Rotation3 {
        // Equatorial → ICRS → Ecliptic
        let r1: Rotation3 = <() as FrameRotationProvider<Equatorial, ICRS>>::rotation(jd, ctx);
        let r2: Rotation3 = <() as FrameRotationProvider<ICRS, Ecliptic>>::rotation(jd, ctx);
        r2 * r1
    }
}

/// Ecliptic → Equatorial rotation.
impl FrameRotationProvider<Ecliptic, Equatorial> for () {
    #[inline]
    fn rotation<Eph, Eop, Nut>(jd: JulianDate, ctx: &AstroContext<Eph, Eop, Nut>) -> Rotation3 {
        // Ecliptic → ICRS → Equatorial
        let r1: Rotation3 = <() as FrameRotationProvider<Ecliptic, ICRS>>::rotation(jd, ctx);
        let r2: Rotation3 = <() as FrameRotationProvider<ICRS, Equatorial>>::rotation(jd, ctx);
        r2 * r1
    }
}

// =============================================================================
// Center Shift Implementations (Hub: Barycentric)
// =============================================================================

/// Heliocentric → Barycentric shift.
///
/// Returns the position of the Sun in barycentric coordinates.
/// Uses VSOP87E for the Sun's position.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Heliocentric, Barycentric, F> for () {
    fn shift<Eph, Eop, Nut>(jd: JulianDate, _ctx: &AstroContext<Eph, Eop, Nut>) -> [f64; 3] {
        use crate::bodies::solar_system::Sun;
        
        // Get Sun's position in barycentric ecliptic coordinates
        let sun_bary = Sun::vsop87e(jd);
        let pos = sun_bary.get_position();
        
        // The shift is the Sun's position (in AU)
        [pos.x().value(), pos.y().value(), pos.z().value()]
    }
}

/// Barycentric → Heliocentric shift.
///
/// Returns the negation of the Sun's barycentric position.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Barycentric, Heliocentric, F> for () {
    fn shift<Eph, Eop, Nut>(jd: JulianDate, ctx: &AstroContext<Eph, Eop, Nut>) -> [f64; 3] {
        let [x, y, z] = <() as CenterShiftProvider<Heliocentric, Barycentric, F>>::shift(jd, ctx);
        [-x, -y, -z]
    }
}

/// Geocentric → Barycentric shift.
///
/// Returns the position of the Earth in barycentric coordinates.
/// Uses VSOP87E for the Earth's position.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Geocentric, Barycentric, F> for () {
    fn shift<Eph, Eop, Nut>(jd: JulianDate, _ctx: &AstroContext<Eph, Eop, Nut>) -> [f64; 3] {
        use crate::bodies::solar_system::Earth;
        
        // Get Earth's position in barycentric ecliptic coordinates
        let earth_bary = Earth::vsop87e(jd);
        let pos = earth_bary.get_position();
        
        // The shift is the Earth's position (in AU)
        [pos.x().value(), pos.y().value(), pos.z().value()]
    }
}

/// Barycentric → Geocentric shift.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Barycentric, Geocentric, F> for () {
    fn shift<Eph, Eop, Nut>(jd: JulianDate, ctx: &AstroContext<Eph, Eop, Nut>) -> [f64; 3] {
        let [x, y, z] = <() as CenterShiftProvider<Geocentric, Barycentric, F>>::shift(jd, ctx);
        [-x, -y, -z]
    }
}

/// Heliocentric → Geocentric shift.
///
/// Composed via Barycentric hub.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Heliocentric, Geocentric, F> for () {
    fn shift<Eph, Eop, Nut>(jd: JulianDate, ctx: &AstroContext<Eph, Eop, Nut>) -> [f64; 3] {
        // Helio → Bary → Geo
        let [x1, y1, z1] = <() as CenterShiftProvider<Heliocentric, Barycentric, F>>::shift(jd, ctx);
        let [x2, y2, z2] = <() as CenterShiftProvider<Barycentric, Geocentric, F>>::shift(jd, ctx);
        [x1 + x2, y1 + y2, z1 + z2]
    }
}

/// Geocentric → Heliocentric shift.
impl<F: affn::ReferenceFrame> CenterShiftProvider<Geocentric, Heliocentric, F> for () {
    fn shift<Eph, Eop, Nut>(jd: JulianDate, ctx: &AstroContext<Eph, Eop, Nut>) -> [f64; 3] {
        let [x, y, z] = <() as CenterShiftProvider<Heliocentric, Geocentric, F>>::shift(jd, ctx);
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
/// use siderust::coordinates::frames::{ICRS, Ecliptic};
/// use siderust::astro::JulianDate;
///
/// let rot = frame_rotation::<ICRS, Ecliptic>(JulianDate::J2000, &AstroContext::default());
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
/// use siderust::coordinates::frames::Ecliptic;
/// use siderust::astro::JulianDate;
///
/// let shift = center_shift::<Heliocentric, Geocentric, Ecliptic>(
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
        let rot = frame_rotation::<ICRS, Ecliptic>(JulianDate::J2000, &AstroContext::default());
        
        // X-axis should be unchanged (rotation is around X)
        let x = rot.apply_array([1.0, 0.0, 0.0]);
        assert!((x[0] - 1.0).abs() < EPSILON);
        assert!(x[1].abs() < EPSILON);
        assert!(x[2].abs() < EPSILON);
        
        // Y and Z should be rotated
        let y = rot.apply_array([0.0, 1.0, 0.0]);
        assert!(y[1].abs() < 1.0); // Y component reduced
        assert!(y[2].abs() > 0.0); // Z component introduced
    }

    #[test]
    fn test_ecliptic_icrs_roundtrip() {
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;
        
        let r1 = frame_rotation::<ICRS, Ecliptic>(jd, &ctx);
        let r2 = frame_rotation::<Ecliptic, ICRS>(jd, &ctx);
        
        let v = [1.0, 2.0, 3.0];
        let roundtrip = r2.apply_array(r1.apply_array(v));
        
        assert!((roundtrip[0] - v[0]).abs() < EPSILON);
        assert!((roundtrip[1] - v[1]).abs() < EPSILON);
        assert!((roundtrip[2] - v[2]).abs() < EPSILON);
    }

    #[test]
    fn test_identity_center_shift() {
        let shift = center_shift::<Barycentric, Barycentric, Ecliptic>(
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
        let helio_geo = center_shift::<Heliocentric, Geocentric, Ecliptic>(jd, &ctx);
        let helio_bary = center_shift::<Heliocentric, Barycentric, Ecliptic>(jd, &ctx);
        let bary_geo = center_shift::<Barycentric, Geocentric, Ecliptic>(jd, &ctx);
        
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
        
        let forward = center_shift::<Heliocentric, Geocentric, Ecliptic>(jd, &ctx);
        let backward = center_shift::<Geocentric, Heliocentric, Ecliptic>(jd, &ctx);
        
        assert!((forward[0] + backward[0]).abs() < EPSILON);
        assert!((forward[1] + backward[1]).abs() < EPSILON);
        assert!((forward[2] + backward[2]).abs() < EPSILON);
    }
}
