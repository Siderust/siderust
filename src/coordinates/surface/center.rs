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
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize))]
#[cfg_attr(feature = "serde", serde(bound(serialize = "")))]
pub struct SurfaceObserverParams<B: BodySurfaceModel> {
    /// Longitude, positive eastward (degrees).
    lon: Degrees,
    /// Latitude per [`Self::latitude_convention`] (degrees).
    lat: Degrees,
    /// Height above the reference surface (metres).
    height: Meters,
    /// Reference surface used for Cartesian / ENU conversion.
    surface: ReferenceSurface,
    /// Latitude convention used by `lat`.
    latitude_convention: LatitudeConvention,
    #[cfg_attr(feature = "serde", serde(skip))]
    _body: PhantomData<B>,
}

#[cfg(feature = "serde")]
impl<'de, B: BodySurfaceModel> Deserialize<'de> for SurfaceObserverParams<B> {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        struct Fields {
            lon: Degrees,
            lat: Degrees,
            height: Meters,
            surface: ReferenceSurface,
            latitude_convention: LatitudeConvention,
        }
        let site = Fields::deserialize(deserializer)?;
        let site = SurfaceSite::<B>::with_surface(
            site.lon,
            site.lat,
            site.height,
            site.surface,
            site.latitude_convention,
        )
        .map_err(serde::de::Error::custom)?;
        Ok(Self::from_site(&site))
    }
}

impl<B: BodySurfaceModel> PartialEq for SurfaceObserverParams<B> {
    fn eq(&self, other: &Self) -> bool {
        self.lon == other.lon
            && self.lat == other.lat
            && self.height == other.height
            && self.surface == other.surface
            && self.latitude_convention == other.latitude_convention
    }
}

impl<B: BodySurfaceModel> SurfaceObserverParams<B> {
    /// Observer longitude, positive east.
    pub fn lon(&self) -> Degrees {
        self.lon
    }

    /// Observer latitude, in the selected convention.
    pub fn lat(&self) -> Degrees {
        self.lat
    }

    /// Height above the reference surface (negative below).
    pub fn height(&self) -> Meters {
        self.height
    }

    /// The validated reference surface.
    pub fn surface(&self) -> &ReferenceSurface {
        &self.surface
    }

    /// Latitude convention for this observer.
    pub fn latitude_convention(&self) -> LatitudeConvention {
        self.latitude_convention
    }

    /// Build params from a typed [`SurfaceSite`].
    pub fn from_site(site: &SurfaceSite<B>) -> Self {
        Self {
            lon: site.lon(),
            lat: site.lat(),
            height: site.height(),
            surface: *site.surface(),
            latitude_convention: site.latitude_convention(),
            _body: PhantomData,
        }
    }
}

impl<B: BodySurfaceModel> Default for SurfaceObserverParams<B> {
    /// Default body-specific site parameters at lon=lat=height=0.
    ///
    /// The body, reference surface, and latitude convention always match `B`.
    /// Prefer building params from a validated [`SurfaceSite`].
    fn default() -> Self {
        let site = SurfaceSite::<B>::new(Degrees::new(0.0), Degrees::new(0.0), Meters::new(0.0))
            .expect("BodySurfaceModel default site must be geometrically valid");
        Self::from_site(&site)
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
    type Params = SurfaceObserverParams<B>;

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

    #[cfg(feature = "serde")]
    #[test]
    fn observer_params_deserialization_enforces_site_invariants() {
        use crate::qtty::{DEG, M};
        let site = SurfaceSite::<MoonSurface>::new(10.0 * DEG, 20.0 * DEG, 0.0 * M).unwrap();
        let mut payload = serde_json::to_value(site.to_params()).unwrap();
        payload["lat"] = serde_json::to_value(95.0 * DEG).unwrap();
        assert!(serde_json::from_value::<SurfaceObserverParams<MoonSurface>>(payload).is_err());
    }

    #[test]
    fn center_name_is_stable() {
        assert_eq!(
            <SurfaceObserver<MoonSurface> as ReferenceCenter>::center_name(),
            "SurfaceObserver"
        );
        assert_eq!(SurfaceObserver::<MoonSurface>::body_name(), "Moon");
    }

    #[test]
    fn observer_params_getters_preserve_site_and_cannot_mutate_invariants() {
        use crate::qtty::{DEG, M};
        let site = SurfaceSite::<MoonSurface>::new(30.0 * DEG, 20.0 * DEG, -10.0 * M)
            .unwrap();
        let params = SurfaceObserverParams::from_site(&site);
        assert_eq!(params.lon(), site.lon());
        assert_eq!(params.lat(), site.lat());
        assert_eq!(params.height(), site.height());
        assert_eq!(params.surface(), site.surface());
        assert_eq!(params.latitude_convention(), site.latitude_convention());
        assert_eq!(params, site.to_params());
    }
}
