// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # ITRS → Equatorial rotation provider
//!
//! Centralises the composite rotation matrix that transforms Earth-fixed
//! (ITRS/ECEF) Cartesian vectors into the mean equatorial frame of J2000.0,
//! using EOP data and the configured nutation model.
//!
//! ## Scientific scope
//!
//! Earth's instantaneous orientation in inertial space is the product of
//! polar motion, diurnal rotation, precession and nutation, all of which
//! must be assembled to project a terrestrial vector onto a celestial frame.
//! This module owns the canonical assembly of those matrices for siderust,
//! so that every topocentric, horizontal or stellar transform path uses the
//! same IAU-compliant chain and therefore yields self-consistent positions.
//!
//! ## Technical scope
//!
//! The chain implemented is
//!
//! ```text
//! ITRS  ──W⁻¹──▶  TIRS  ──ERA──▶  CIRS  ──Q──▶  GCRS/ICRS  ──P──▶  EquatorialMeanJ2000
//! ```
//!
//! with **W** the polar-motion matrix `(xₚ, yₚ, s′)` from EOP, **ERA** the
//! Earth Rotation Angle from UT1 (derived from EOP `dUT1`), **Q** the
//! CIO/CIP matrix `(X, Y, s)` from the active nutation model with optional
//! IERS celestial-pole offsets `dX, dY`, and **P** the ICRS → mean J2000.0
//! frame rotation. Inputs are **TT Julian Dates**; the conversion to UT1 is
//! performed internally using the EOP `dUT1`.
//!
//! ## References
//!
//! * IERS Conventions (2010), Chapter 5
//! * SOFA Earth-rotation cookbook (`iauC2t06a` and friends)

use crate::astro::cio;
use crate::astro::earth_rotation::{jd_ut1_from_tt_eop, try_jd_utc_from_tt};
use crate::astro::eop::{EopProvider, EopValues};
use crate::astro::era::earth_rotation_angle;
use crate::astro::nutation::NutationModel;
use crate::astro::polar_motion::polar_motion_matrix_from_eop;
use crate::coordinates::frames::{EquatorialMeanJ2000, ICRS};
use crate::coordinates::transform::context::AstroContext;
use crate::coordinates::transform::providers::FrameRotationProvider;
use crate::time::JulianDate;

/// Apply IERS celestial pole offsets dX, dY to raw IAU 2000B nutation values.
///
/// Returns corrected `(dpsi, deps)` in radians.
#[inline]
pub(crate) fn nutation_with_celestial_pole_offsets<Nut: NutationModel>(
    jd: JulianDate,
    eop: EopValues,
) -> (crate::qtty::Radians, crate::qtty::Radians) {
    let nut = Nut::nutation(jd);
    let mut dpsi = nut.dpsi;
    let mut deps = nut.deps;

    // Apply IERS celestial pole offsets dX,dY as first-order corrections.
    // dψ_eop = dX/sin(εA), dε_eop = dY.
    let dx_rad = crate::qtty::Radians::from(eop.dx);
    let dy_rad = crate::qtty::Radians::from(eop.dy);
    let zero = crate::qtty::Radians::new(0.0);
    if dx_rad != zero || dy_rad != zero {
        let sin_eps = nut.mean_obliquity.sin();
        if sin_eps.abs() > 1e-15 {
            dpsi += dx_rad / sin_eps;
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
/// * `jd` , observation epoch in **Terrestrial Time** (TT Julian Date).
/// * `ctx`, [`AstroContext`] providing the EOP and ephemeris providers.
///
/// # Time-scale contract
///
/// `jd` **must** be Terrestrial Time (TT). The function converts it to UTC
/// before consulting the context's EOP provider, because IERS EOP tables are
/// indexed by UTC civil date.
///
/// # Returns
///
/// A 3×3 rotation matrix `R` such that `v_eq = R * v_itrs`.
pub(crate) fn itrs_to_equatorial_mean_j2000_rotation<Eph, Eop: EopProvider, Nut: NutationModel>(
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop>,
) -> affn::Rotation3 {
    let eop = match try_jd_utc_from_tt(jd) {
        Ok(jd_utc) => ctx.eop_at(jd_utc),
        Err(_) => EopValues::default(),
    };
    let jd_ut1 = jd_ut1_from_tt_eop(jd, &eop);

    let (dpsi, deps) = nutation_with_celestial_pole_offsets::<Nut>(jd, eop);
    let cip = cio::cip_cio(jd, dpsi, deps);
    let q = cio::gcrs_to_cirs_matrix(cip.x, cip.y, cip.s);
    let w = polar_motion_matrix_from_eop(eop.xp, eop.yp, jd);
    let era = earth_rotation_angle(jd_ut1);

    // NOTE: `cio::gcrs_to_cirs_matrix` and `Rotation3::rz` conventions are
    // aligned such that this composition matches the SOFA chain for
    // ITRS -> GCRS in this codebase.
    let itrs_to_gcrs = q * affn::Rotation3::rz(-era) * w.inverse();
    let icrs_to_j2000 =
        <() as FrameRotationProvider<ICRS, EquatorialMeanJ2000>>::rotation::<Eph, Eop, Nut>(
            jd, ctx,
        );

    icrs_to_j2000 * itrs_to_gcrs
}
