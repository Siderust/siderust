// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Parameterized surface-observer reference centre.

use super::body::{BodySurfaceModel, LatitudeConvention};
use super::reference_surface::ReferenceSurface;
use super::site::SurfaceSite;
use crate::coordinates::centers::ReferenceCenter;
use crate::qtty::{Degrees, Meters};
use affn::AffineCenter;
use core::marker::PhantomData;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Runtime site parameters stored on [`SurfaceObserver`] coordinates.
///
/// The celestial body is encoded in the centre type parameter `B` of
/// [`SurfaceObserver<B>`]; these fields carry the site geometry only.
/// Parameter equality is checked at runtime (same pattern as
/// [`crate::coordinates::centers::Topocentric`]).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SurfaceObserverParams {
    /// Longitude, positive eastward (degrees).
    pub lon: Degrees,
    /// Latitude per [`Self::latitude_convention`] (degrees).
    pub lat: Degrees,
    /// Height above the reference surface (metres).
    pub height: Meters,
    /// Reference surface used for Cartesian / ENU conversion.
    pub surface: ReferenceSurface,
    /// Latitude convention used by `lat`.
    pub latitude_convention: LatitudeConvention,
}

impl SurfaceObserverParams {
    /// Build params from a typed [`SurfaceSite`].
    pub fn from_site<B: BodySurfaceModel>(site: &SurfaceSite<B>) -> Self {
        Self {
            lon: site.lon(),
            lat: site.lat(),
            height: site.height(),
            surface: *site.surface(),
            latitude_convention: site.latitude_convention(),
        }
    }
}

impl Default for SurfaceObserverParams {
    /// Origin on the Moon IAU sphere (placeholder for `ReferenceCenter::Params`).
    ///
    /// Prefer constructing params from a concrete [`SurfaceSite`]. This default
    /// exists only to satisfy the `affn::ReferenceCenter::Params: Default` bound.
    fn default() -> Self {
        use super::body::{LatitudeConvention, MOON_IAU_SPHERE};
        Self {
            lon: Degrees::new(0.0),
            lat: Degrees::new(0.0),
            height: Meters::new(0.0),
            surface: MOON_IAU_SPHERE,
            latitude_convention: LatitudeConvention::Planetocentric,
        }
    }
}

/// Observer on the surface of celestial body `B`.
///
/// Analogous to [`crate::coordinates::centers::Topocentric`] for Earth, but
/// generic over any [`BodySurfaceModel`]. Horizontal positions that use this
/// centre carry [`SurfaceObserverParams`] so the observation site remains
/// attached to the coordinate value.
///
/// Earth-based APIs that already use [`crate::coordinates::centers::Topocentric`]
/// remain unchanged; prefer `Topocentric` when the site is a WGS84
/// [`crate::coordinates::centers::Geodetic`]`::<ECEF>` and
/// [`SurfaceObserver`]`<Earth>` / [`SurfaceSite`]`<Earth>` when sharing the
/// generic body-surface pipeline with Moon/Mars.
#[derive(Debug, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SurfaceObserver<B: BodySurfaceModel>(PhantomData<B>);

impl<B: BodySurfaceModel> ReferenceCenter for SurfaceObserver<B> {
    type Params = SurfaceObserverParams;

    fn center_name() -> &'static str {
        // Associated function cannot return a dynamic body name without
        // allocation; the body is recoverable from `B::body_name()` at the
        // call site. Keep a stable kind name for diagnostics.
        "SurfaceObserver"
    }
}

impl<B: BodySurfaceModel> AffineCenter for SurfaceObserver<B> {}

impl<B: BodySurfaceModel> SurfaceObserver<B> {
    /// Body name associated with this centre kind.
    #[inline]
    pub fn body_name() -> &'static str {
        B::body_name()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::surface::MoonSurface;

    #[test]
    fn center_name_is_stable() {
        assert_eq!(
            <SurfaceObserver<MoonSurface> as ReferenceCenter>::center_name(),
            "SurfaceObserver"
        );
        assert_eq!(SurfaceObserver::<MoonSurface>::body_name(), "Moon");
    }
}
