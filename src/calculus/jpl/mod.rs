// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Shared JPL DE4xx implementation.
//!
//! This module provides the common infrastructure shared by all JPL DE4xx
//! ephemeris backends (DE440, DE441, etc.):
//!
//! - [`eval`] — Chebyshev polynomial evaluation and [`SegmentDescriptor`](eval::SegmentDescriptor).
//! - [`bodies`] — Generic body-chain resolution (Sun, Earth, Moon positions/velocities).
//! - [`DeData`] — Trait abstracting over the per-version coefficient data.
//! - [`DeEphemeris`] — Generic zero-sized ephemeris backend for any `DeData` impl.

pub mod bodies;
pub mod eval;

#[cfg(feature = "de440")]
pub mod de440;
#[cfg(feature = "de441")]
pub mod de441;

use eval::SegmentDescriptor;

use crate::calculus::ephemeris::{AuPerDay, Ephemeris};
use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::EclipticMeanJ2000,
};
use crate::targets::Target;
use crate::time::JulianDate;
use qtty::{AstronomicalUnit, Kilometer};

use core::marker::PhantomData;

/// Trait abstracting over DE4xx coefficient data sources.
///
/// Each JPL DE version (DE440, DE441, …) implements this trait with
/// its own embedded segment descriptors. The descriptors are associated
/// constants so the compiler can fully monomorphize all evaluation paths.
pub trait DeData: 'static {
    /// Segment descriptor for the Sun (NAIF 10 → SSB).
    const SUN: SegmentDescriptor;
    /// Segment descriptor for the Earth-Moon Barycenter (NAIF 3 → SSB).
    const EMB: SegmentDescriptor;
    /// Segment descriptor for the Moon (NAIF 301 → EMB).
    const MOON: SegmentDescriptor;
}

/// Generic zero-sized ephemeris backend for any JPL DE4xx dataset.
///
/// `D` selects the coefficient data at compile time. All dispatch is
/// monomorphized away — there is no runtime cost vs. a hand-written backend.
///
/// # Example
/// ```rust,ignore
/// use siderust::calculus::jpl::DeEphemeris;
/// use siderust::calculus::de440::De440Data;
///
/// type De440Ephemeris = DeEphemeris<De440Data>;
/// ```
#[derive(Debug, Clone, Copy, Default)]
pub struct DeEphemeris<D: DeData>(PhantomData<D>);

impl<D: DeData> Ephemeris for DeEphemeris<D> {
    #[inline]
    fn sun_barycentric(
        jd: JulianDate,
    ) -> Target<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>> {
        bodies::sun_barycentric(jd, &D::SUN)
    }

    #[inline]
    fn earth_barycentric(
        jd: JulianDate,
    ) -> Target<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>> {
        bodies::earth_barycentric(jd, &D::EMB, &D::MOON)
    }

    #[inline]
    fn earth_heliocentric(
        jd: JulianDate,
    ) -> Target<Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>> {
        bodies::earth_heliocentric(jd, &D::SUN, &D::EMB, &D::MOON)
    }

    #[inline]
    fn earth_barycentric_velocity(jd: JulianDate) -> Velocity<EclipticMeanJ2000, AuPerDay> {
        bodies::earth_barycentric_velocity(jd, &D::EMB, &D::MOON)
    }

    #[inline]
    fn moon_geocentric(jd: JulianDate) -> Position<Geocentric, EclipticMeanJ2000, Kilometer> {
        bodies::moon_geocentric(jd, &D::MOON)
    }
}
