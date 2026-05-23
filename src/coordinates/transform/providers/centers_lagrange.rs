// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Center-shift providers for dynamic Sun-Earth Lagrange centers.
//!
//! Each center is represented by its barycentric ecliptic J2000 position from
//! the embedded Chebyshev archive when available. The committed placeholder
//! archive is empty, so transform providers fall back to the N-body solver with
//! the caller-selected ephemeris backend.

use super::*;
use crate::ephemeris::lagrange::{fallback_or_solve, SunEarthLagrangePoint};

macro_rules! impl_lagrange_center_shift {
    ($center:ty, $point:expr) => {
        impl<F> CenterShiftProvider<$center, Barycentric, F> for ()
        where
            F: affn::ReferenceFrame,
            (): FrameRotationProvider<EclipticMeanJ2000, F>,
        {
            #[inline]
            fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop>,
            ) -> AuShift {
                let bary_pos = fallback_or_solve::<Eph>($point, jd);
                rotate_shift_from_ecliptic::<_, F, Eph, Eop, Nut>(bary_pos, jd, ctx)
            }
        }

        impl_reverse_center_shifts!($center);
    };
}

impl_lagrange_center_shift!(SunEarthL1, SunEarthLagrangePoint::L1);
impl_lagrange_center_shift!(SunEarthL2, SunEarthLagrangePoint::L2);
impl_lagrange_center_shift!(SunEarthL3, SunEarthLagrangePoint::L3);
impl_lagrange_center_shift!(SunEarthL4, SunEarthLagrangePoint::L4);
impl_lagrange_center_shift!(SunEarthL5, SunEarthLagrangePoint::L5);
