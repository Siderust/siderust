// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Polar Motion, IERS Conventions
//!
//! Builds the polar-motion matrix **W** that links the Terrestrial
//! Intermediate Reference System (TIRS) to the International Terrestrial
//! Reference System (ITRS).
//!
//! ## Scientific scope
//!
//! Polar motion describes the deviation of Earth's rotation axis from its
//! crust-fixed reference position. It is characterised by the pole
//! coordinates `(xₚ, yₚ)` published by the IERS, which combine the Chandler
//! wobble, an annual term, and irregular components driven by mass
//! redistribution in the atmosphere, oceans and core. Even a few hundred
//! milliarcseconds of polar motion produce metre-level shifts at Earth's
//! surface, so polar motion must be modelled for any precise terrestrial-to-
//! celestial transformation.
//!
//! ## Technical scope
//!
//! The polar-motion matrix is assembled as
//!
//! ```text
//! W = R₃(−s') · R₂(xₚ) · R₁(yₚ)
//! ```
//!
//! where `s'` is the TIO locator, here implemented as the linear-in-time
//! IERS form `s' ≈ −47 μas × t` (Julian centuries from J2000 on TT).
//! Convenience helpers convert IERS [`Arcseconds`] inputs to [`Radians`]
//! and assemble the rotation via `affn::Rotation3` factors.
//!
//! ## References
//!
//! * IERS Conventions (2010), §5.4.2, eq. 5.13
//! * SOFA routines `iauSp00`, `iauPom00`

use crate::qtty::*;
use crate::time::JulianDate;
use affn::Rotation3;

/// TIO locator s'.
///
/// ```text
/// s' ≈ −47 μas × t
/// ```
///
/// where t is Julian centuries from J2000 on the TT scale.
///
/// ## References
/// * IERS Conventions (2010), eq. 5.13
/// * SOFA routine `iauSp00`
#[inline]
pub fn tio_locator_sp(jd_tt: JulianDate) -> Radians {
    let t = jd_tt.julian_centuries();
    MicroArcseconds::new(-47.0 * t).to::<Radian>()
}

/// Polar motion matrix **W**.
///
/// Transforms from TIRS to ITRS (or equivalently, from the intermediate
/// frame to the Earth-fixed frame):
///
/// ```text
/// W = R₃(−s') · R₂(xₚ) · R₁(yₚ)
/// ```
///
/// ## References
/// * IERS Conventions (2010), §5.4.2
/// * SOFA routine `iauPom00`
pub fn polar_motion_matrix(xp: Radians, yp: Radians, sp: Radians) -> Rotation3 {
    Rotation3::rz(-sp) * Rotation3::ry(xp) * Rotation3::rx(yp)
}

/// Convenience: compute W from pole coordinates as [`Arcseconds`] and Julian Date.
///
/// This handles the unit conversion and TIO locator computation.
#[inline]
pub fn polar_motion_matrix_from_eop(
    xp: Arcseconds,
    yp: Arcseconds,
    jd_tt: JulianDate,
) -> Rotation3 {
    let sp = tio_locator_sp(jd_tt);
    polar_motion_matrix(xp.to::<Radian>(), yp.to::<Radian>(), sp)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tio_locator_at_j2000_is_zero() {
        let sp = tio_locator_sp(crate::J2000);
        assert!(
            sp.value().abs() < 1e-15,
            "s' at J2000 should be ~0, got {}",
            sp
        );
    }

    #[test]
    fn tio_locator_is_small() {
        // After 25 years: s' ≈ −47e-6 × 0.25 arcsec ≈ −11.75 μas
        let jd = crate::time::JulianDate::new(2_460_000.5); // ~2023
        let sp = tio_locator_sp(jd);
        let sp_uas = sp.to::<MicroArcsecond>().value();
        assert!(
            sp_uas.abs() < 50.0,
            "s' should be ~12 μas, got {} μas",
            sp_uas
        );
    }

    #[test]
    fn polar_motion_identity_with_zero_poles() {
        let zero = Radians::new(0.0);
        let w = polar_motion_matrix(zero, zero, zero);
        let m = w.as_matrix();
        for (i, row) in m.iter().enumerate().take(3) {
            assert!(
                (row[i] - 1.0).abs() < 1e-15,
                "W[{}][{}] = {}, expected 1",
                i,
                i,
                row[i]
            );
        }
    }

    #[test]
    fn polar_motion_is_proper_rotation() {
        let xp = Arcseconds::new(0.1).to::<Radian>();
        let yp = Arcseconds::new(0.3).to::<Radian>();
        let sp = Radians::new(-1e-9);
        let w = polar_motion_matrix(xp, yp, sp);
        let m = w.as_matrix();

        // Check det ≈ 1
        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        assert!((det - 1.0).abs() < 1e-14, "det(W) = {}, expected ≈ 1", det);
    }

    #[test]
    fn polar_motion_small_angles() {
        // For small pole coordinates (~0.3″), the rotation should be tiny
        let w = polar_motion_matrix_from_eop(
            Arcseconds::new(0.3),
            Arcseconds::new(0.2),
            crate::time::JulianDate::new(2_460_000.5),
        );
        let m = w.as_matrix();

        // Diagonal elements should be very close to 1
        for (i, row) in m.iter().enumerate().take(3) {
            assert!(
                (row[i] - 1.0).abs() < 1e-10,
                "W[{}][{}] = {}, should be ≈ 1 for small poles",
                i,
                i,
                row[i]
            );
        }
    }
}
