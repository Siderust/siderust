// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed surface observer sites.

use super::body::{BodySurfaceModel, LatitudeConvention};
use super::center::SurfaceObserverParams;
use super::enu::{body_fixed_position, local_enu_basis, EnuBasis};
use super::error::SurfaceError;
use super::reference_surface::ReferenceSurface;
use crate::coordinates::cartesian::Position;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::coordinates::surface::EarthSurface;
#[cfg(not(feature = "std"))]
use crate::qtty::Scalar;
use crate::qtty::{Degrees, LengthUnit, Meter, Meters, Quantity};
use core::marker::PhantomData;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Observer site on the reference surface of celestial body `B`.
///
/// Longitude is positive **east**. Latitude follows
/// [`LatitudeConvention`] (body default unless overridden). Height is
/// measured above the configured [`ReferenceSurface`] along the local
/// vertical of that convention.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::surface::{MoonSurface, SurfaceSite};
/// use siderust::qtty::*;
///
/// let site = SurfaceSite::<MoonSurface>::new(
///     20.0 * DEG,   // lon east
///     -3.0 * DEG,   // lat
///     0.0 * M,      // height
/// ).unwrap();
/// assert_eq!(site.body_name(), "Moon");
/// ```
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SurfaceSite<B: BodySurfaceModel> {
    lon: Degrees,
    lat: Degrees,
    height: Meters,
    surface: ReferenceSurface,
    latitude_convention: LatitudeConvention,
    #[cfg_attr(feature = "serde", serde(skip))]
    _body: PhantomData<B>,
}

impl<B: BodySurfaceModel> PartialEq for SurfaceSite<B> {
    fn eq(&self, other: &Self) -> bool {
        self.lon == other.lon
            && self.lat == other.lat
            && self.height == other.height
            && self.surface == other.surface
            && self.latitude_convention == other.latitude_convention
    }
}

impl<B: BodySurfaceModel> SurfaceSite<B> {
    /// Construct a site using the body's default surface and latitude convention.
    pub fn new(lon: Degrees, lat: Degrees, height: Meters) -> Result<Self, SurfaceError> {
        B::ensure_supported()?;
        Self::with_surface(
            lon,
            lat,
            height,
            B::default_surface(),
            B::default_latitude_convention(),
        )
    }

    /// Construct a site with an explicit reference surface and latitude convention.
    pub fn with_surface(
        lon: Degrees,
        lat: Degrees,
        height: Meters,
        surface: ReferenceSurface,
        latitude_convention: LatitudeConvention,
    ) -> Result<Self, SurfaceError> {
        B::ensure_supported()?;
        // Re-validate surface parameters (const defaults are trusted; custom ones may not be).
        match surface {
            ReferenceSurface::Sphere { radius } => {
                let _ = ReferenceSurface::sphere(radius)?;
            }
            ReferenceSurface::Ellipsoid { a, f } => {
                let _ = ReferenceSurface::ellipsoid(a, f)?;
            }
        }
        if !lon.value().is_finite() {
            return Err(SurfaceError::NonFiniteInput { field: "longitude" });
        }
        if !lat.value().is_finite() {
            return Err(SurfaceError::NonFiniteInput { field: "latitude" });
        }
        if !height.value().is_finite() {
            return Err(SurfaceError::NonFiniteInput { field: "height" });
        }
        let lat_n = normalize_latitude(lat)?;
        let lon_n = Degrees::new(lon.value().rem_euclid(360.0));
        // Map to (−180, 180] for stable display while preserving east-positive.
        let lon_n = if lon_n.value() > 180.0 {
            Degrees::new(lon_n.value() - 360.0)
        } else {
            lon_n
        };

        Ok(Self {
            lon: lon_n,
            lat: lat_n,
            height,
            surface,
            latitude_convention,
            _body: PhantomData,
        })
    }

    /// Body name.
    #[inline]
    pub fn body_name(&self) -> &'static str {
        B::body_name()
    }

    /// East-positive longitude.
    #[inline]
    pub fn lon(&self) -> Degrees {
        self.lon
    }

    /// Latitude in this site's convention.
    #[inline]
    pub fn lat(&self) -> Degrees {
        self.lat
    }

    /// Height above the reference surface.
    #[inline]
    pub fn height(&self) -> Meters {
        self.height
    }

    /// Reference surface.
    #[inline]
    pub fn surface(&self) -> &ReferenceSurface {
        &self.surface
    }

    /// Latitude convention.
    #[inline]
    pub fn latitude_convention(&self) -> LatitudeConvention {
        self.latitude_convention
    }

    /// Body-fixed Cartesian position of the observer (metres).
    pub fn body_fixed_cartesian(&self) -> Result<[Meters; 3], SurfaceError> {
        body_fixed_position(
            &self.surface,
            self.lon,
            self.lat,
            self.height,
            self.latitude_convention,
        )
    }

    /// Local ENU basis in the body-fixed frame.
    pub fn enu_basis(&self) -> Result<EnuBasis, SurfaceError> {
        local_enu_basis(&self.surface, self.lon, self.lat, self.latitude_convention)
    }

    /// Observer position as a typed body-centred body-fixed Cartesian vector.
    pub fn observer_position<U>(
        &self,
    ) -> Result<Position<B::Center, B::FixedFrame, U>, SurfaceError>
    where
        U: LengthUnit,
        Quantity<U>: From<Quantity<Meter>>,
    {
        let [x, y, z] = self.body_fixed_cartesian()?;
        Ok(Position::<B::Center, B::FixedFrame, U>::new(
            Quantity::<U>::from(x),
            Quantity::<U>::from(y),
            Quantity::<U>::from(z),
        ))
    }

    /// Convert to centre parameters for [`super::SurfaceObserver`].
    #[inline]
    pub fn to_params(&self) -> SurfaceObserverParams {
        SurfaceObserverParams::from_site(self)
    }
}

impl SurfaceSite<EarthSurface> {
    /// Build from an existing WGS84 [`Geodetic`]`<ECEF>` site.
    ///
    /// Uses the WGS84 ellipsoid and planetographic (geodetic) latitude.
    pub fn from_geodetic(site: Geodetic<ECEF>) -> Result<Self, SurfaceError> {
        Self::with_surface(
            site.lon,
            site.lat,
            site.height.to::<Meter>(),
            EarthSurface::default_surface(),
            LatitudeConvention::Planetographic,
        )
    }

    /// View this Earth site as [`Geodetic`]`<ECEF>` (WGS84).
    ///
    /// Only valid when this site uses the default WGS84 ellipsoid and
    /// planetographic latitude; otherwise returns [`SurfaceError::DegenerateGeometry`].
    pub fn to_geodetic(&self) -> Result<Geodetic<ECEF>, SurfaceError> {
        if self.latitude_convention != LatitudeConvention::Planetographic {
            return Err(SurfaceError::DegenerateGeometry {
                reason: "Earth Geodetic requires planetographic latitude",
            });
        }
        if self.surface != EarthSurface::default_surface() {
            return Err(SurfaceError::DegenerateGeometry {
                reason: "Earth Geodetic requires the WGS84 reference surface",
            });
        }
        Ok(Geodetic::<ECEF>::new(self.lon, self.lat, self.height))
    }
}

fn normalize_latitude(lat: Degrees) -> Result<Degrees, SurfaceError> {
    let mut v = lat.value();
    // Fold into [-90, 90] with longitude reflection omitted (caller owns lon).
    if !(-90.0..=90.0).contains(&v) {
        // Allow mild overshoot from floating conversions only.
        if (-90.000_000_1..=90.000_000_1).contains(&v) {
            v = v.clamp(-90.0, 90.0);
        } else {
            return Err(SurfaceError::InvalidLatitude { degrees: v });
        }
    }
    Ok(Degrees::new(v))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::surface::{MarsSurface, MoonSurface};
    use crate::qtty::{DEG, M};

    #[test]
    fn moon_site_defaults() {
        let site = SurfaceSite::<MoonSurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
        assert_eq!(site.body_name(), "Moon");
        let p = site.body_fixed_cartesian().unwrap();
        assert!((p[0].value() - 1_737_400.0).abs() < 1e-6);
    }

    #[test]
    fn rejects_bad_latitude() {
        let err = SurfaceSite::<MarsSurface>::new(0.0 * DEG, 100.0 * DEG, 0.0 * M).unwrap_err();
        assert!(matches!(err, SurfaceError::InvalidLatitude { .. }));
    }

    #[test]
    fn earth_geodetic_roundtrip() {
        let g = Geodetic::<ECEF>::new((-17.89) * DEG, 28.75 * DEG, 2200.0 * M);
        let site = SurfaceSite::<EarthSurface>::from_geodetic(g).unwrap();
        let back = site.to_geodetic().unwrap();
        assert!((back.lat.value() - g.lat.value()).abs() < 1e-12);
        assert!((back.lon.value() - g.lon.value()).abs() < 1e-12);
    }
}
