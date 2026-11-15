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
//! The observer state and center translations obtain that assembly from the
//! same frame provider used by the context-aware horizon projection.
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

use crate::astro::eop::EopProvider;
use crate::astro::nutation::NutationModel;
use crate::coordinates::frames::{EquatorialMeanJ2000, ECEF};
use crate::coordinates::transform::context::AstroContext;
use crate::coordinates::transform::providers::FrameRotationProvider;
use crate::time::JulianDate;

/// Composite rotation from ITRF/ECEF to EquatorialMeanJ2000 at epoch `jd`.
///
/// Delegates to the context-aware ECEF frame provider.
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
    <() as FrameRotationProvider<ECEF, EquatorialMeanJ2000>>::rotation::<Eph, Eop, Nut>(jd, ctx)
}
