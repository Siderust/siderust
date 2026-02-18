// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # ITRS → Equatorial rotation provider
//!
//! This module centralises the composite rotation matrix that transforms
//! Earth-fixed (ITRS/ECEF) Cartesian vectors into the mean equatorial frame
//! of J2000.0.
//!
//! ## Rotation chain
//!
//! ```text
//! ITRS  ──W⁻¹──▶  TIRS  ──ERA──▶  CIRS  ──Q──▶  GCRS/ICRS  ──P──▶  EquatorialMeanJ2000
//! ```
//!
//! Where:
//! - **W** = polar-motion matrix (xₚ, yₚ, s′)  — from EOP
//! - **ERA** = Earth Rotation Angle (from UT1)  — from EOP `dUT1`
//! - **Q** = CIO/CIP matrix (X, Y, s)  — from nutation IAU 2000B ± EOP `dX,dY`
//! - **P** = precession from ICRS to EquatorialMeanJ2000  — frame rotation provider
//!
//! ## Time-scale contract
//!
//! The `jd` argument passed to [`itrs_to_equatorial_mean_j2000_rotation`]
//! must be a **TT Julian Date** (Terrestrial Time).  Internally the function
//! converts to UT1 using EOP `dUT1`.  Do not pass UTC or UT1 directly.

use crate::astro::cio;
use crate::astro::earth_rotation::jd_ut1_from_tt_eop;
use crate::astro::eop::{EopProvider, EopValues};
use crate::astro::era::earth_rotation_angle;
use crate::astro::nutation::nutation_iau2000b;
use crate::astro::polar_motion::polar_motion_matrix_from_eop;
use crate::coordinates::frames::{EquatorialMeanJ2000, ICRS};
use crate::coordinates::transform::context::AstroContext;
use crate::coordinates::transform::providers::FrameRotationProvider;
use crate::time::JulianDate;

/// Apply IERS celestial pole offsets dX, dY to raw IAU 2000B nutation values.
///
/// Returns corrected `(dpsi, deps)` in radians.
#[inline]
pub(crate) fn nutation_with_celestial_pole_offsets(
    jd: JulianDate,
    eop: EopValues,
) -> (qtty::Radians, qtty::Radians) {
    let nut = nutation_iau2000b(jd);
    let mut dpsi = nut.dpsi;
    let mut deps = nut.deps;

    // Apply IERS celestial pole offsets dX,dY as first-order corrections.
    // dψ_eop = dX/sin(εA), dε_eop = dY.
    let dx_rad = qtty::Radians::from(eop.dx);
    let dy_rad = qtty::Radians::from(eop.dy);
    if dx_rad.value() != 0.0 || dy_rad.value() != 0.0 {
        let sin_eps = nut.mean_obliquity.sin();
        if sin_eps.abs() > 1e-15 {
            dpsi += qtty::Radians::new(dx_rad.value() / sin_eps);
        }
        deps += dy_rad;
    }

    (dpsi, deps)
}

/// Composite rotation from ITRF/ECEF to EquatorialMeanJ2000 at epoch `jd`.
///
/// This is the single choke-point for Earth-rotation in siderust:
/// every code path that converts an Earth-fixed vector to an equatorial
/// Cartesian frame calls this function.
///
/// # Arguments
///
/// * `jd`  — observation epoch in **Terrestrial Time** (TT Julian Date).
/// * `ctx` — [`AstroContext`] providing the EOP and ephemeris providers.
///
/// # Time-scale contract
///
/// `jd` **must** be Terrestrial Time (TT). The function internally calls
/// `ctx.eop_at(jd)` which expects a UTC Julian Date — because TT and UTC
/// differ by only ≈ 69 s (accumulated leap seconds + 32.184 s), an epoch
/// passed as TT will introduce at most ~0.5 ms EOP interpolation error,
/// which is below the precision of the daily IERS table.  The subsequent
/// UT1 conversion (`jd_ut1_from_tt_eop`) is then exact.
///
/// # Returns
///
/// A 3×3 rotation matrix `R` such that `v_eq = R * v_itrs`.
pub(crate) fn itrs_to_equatorial_mean_j2000_rotation<Eph, Eop: EopProvider, Nut>(
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop, Nut>,
) -> affn::Rotation3 {
    let eop = ctx.eop_at(jd);
    let jd_ut1 = jd_ut1_from_tt_eop(jd, &eop);

    let (dpsi, deps) = nutation_with_celestial_pole_offsets(jd, eop);
    let cip = cio::cip_cio(jd, dpsi, deps);
    let q = cio::gcrs_to_cirs_matrix(cip.x, cip.y, cip.s);
    let w = polar_motion_matrix_from_eop(eop.xp, eop.yp, jd);
    let era = earth_rotation_angle(jd_ut1);

    // NOTE: `cio::gcrs_to_cirs_matrix` and `Rotation3::rz` conventions are
    // aligned such that this composition matches the ERFA/SOFA chain for
    // ITRS -> GCRS in this codebase.
    let itrs_to_gcrs = q * affn::Rotation3::rz(-era) * w.inverse();
    let icrs_to_j2000 = <() as FrameRotationProvider<ICRS, EquatorialMeanJ2000>>::rotation(jd, ctx);

    icrs_to_j2000 * itrs_to_gcrs
}
