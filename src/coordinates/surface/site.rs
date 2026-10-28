// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed surface observer sites.

use super::body::{BodySurfaceModel, LatitudeConvention};
use super::center::SurfaceObserverParams;
use super::enu::{body_fixed_position, local_enu_basis, EnuBasis};
use super::error::SurfaceError;
use super::reference_surface::ReferenceSurface;
use crate::bodies::solar_system::Earth;
use crate::coordinates::cartesian::Position;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
#[cfg(not(feature = "std"))]
use crate::qtty::Scalar;
use crate::qtty::{Degrees, LengthUnit, Meter, Meters, Quantity};
use core::marker::PhantomData;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Observer site on the reference surface of celestial body `B`.
///
/// Longitude is positive **east** and is never silently remapped from an
/// external planetographic west-positive value. Convert first with
/// [`BodySurfaceModel::east_positive_longitude`] when ingesting NAIF SPICE
/// planetographic coordinates. Latitude follows [`LatitudeConvention`]
/// (body default unless overridden). Height is measured above the configured
/// [`ReferenceSurface`] along its outward normal (radial on a sphere). On an
/// ellipsoid, the local horizon is independent of whether the same site is
/// given planetocentric or planetographic latitude.
///
/// # Example
///
/// ```rust
/// use siderust::bodies::solar_system::Moon;
/// use siderust::coordinates::surface::SurfaceSite;
/// use siderust::qtty::*;
///
/// let site = SurfaceSite::<Moon>::new(
///     20.0 * DEG,   // lon east
///     -3.0 * DEG,   // lat
///     0.0 * M,      // height
/// ).unwrap();
/// assert_eq!(site.body_name(), "Moon");
/// ```
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize))]
#[cfg_attr(feature = "serde", serde(bound(serialize = "")))]
pub struct SurfaceSite<B: BodySurfaceModel> {
    lon: Degrees,
    lat: Degrees,
    height: Meters,
    surface: ReferenceSurface,
    latitude_convention: LatitudeConvention,
    #[cfg_attr(feature = "serde", serde(skip))]
    _body: PhantomData<B>,
}

#[cfg(feature = "serde")]
impl<'de, B: BodySurfaceModel> Deserialize<'de> for SurfaceSite<B> {
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
        Self::with_surface(
            site.lon,
            site.lat,
            site.height,
            site.surface,
            site.latitude_convention,
        )
        .map_err(serde::de::Error::custom)
    }
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
    ///
    /// Longitude must already be east-positive.
    pub fn new(lon: Degrees, lat: Degrees, height: Meters) -> Result<Self, SurfaceError> {
        Self::with_surface(
            lon,
            lat,
            height,
            B::default_surface(),
            B::default_latitude_convention(),
        )
    }

    /// Construct a site from an external planetographic longitude.
    ///
    /// Applies [`BodySurfaceModel::east_positive_longitude`] explicitly, then
    /// uses the body's default surface and latitude convention. Prefer
    /// [`Self::new`] when the longitude is already east-positive.
    pub fn from_planetographic_longitude(
        lon_planetographic: Degrees,
        lat: Degrees,
        height: Meters,
    ) -> Result<Self, SurfaceError> {
        Self::new(B::east_positive_longitude(lon_planetographic), lat, height)
    }

    /// Construct a site with an explicit reference surface and latitude convention.
    ///
    /// Negative heights are valid for a sphere when `h > -radius`. For an
    /// ellipsoid, require `h > -min(M, s)`, where `M` is the meridional
    /// curvature radius at the selected surface point and `s` is its support
    /// distance along the normal. This excludes focal singularities and
    /// positions whose local vertical points towards the body centre. All
    /// Cartesian components must remain finite. Deserialization uses the
    /// same validation.
    pub fn with_surface(
        lon: Degrees,
        lat: Degrees,
        height: Meters,
        surface: ReferenceSurface,
        latitude_convention: LatitudeConvention,
    ) -> Result<Self, SurfaceError> {
        // Public enum variants can bypass constructors; always enforce invariants.
        surface.validate()?;
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

        // Validate geometry immediately (including negative-height domain),
        // not just when the observer is first converted to Cartesian.
        // Deserialization uses this constructor and therefore shares the
        // same invariants.
        body_fixed_position(&surface, lon_n, lat_n, height, latitude_convention)?;

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

    /// This site's longitude expressed in the body's external planetographic
    /// sense (for SPICE / cartographic export). Storage remains east-positive.
    #[inline]
    pub fn planetographic_longitude(&self) -> Degrees {
        B::planetographic_longitude(self.lon)
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
    ///
    /// Returns an error if the requested unit overflows or rounds the entire
    /// nonzero position to zero.
    pub fn observer_position<U>(
        &self,
    ) -> Result<Position<B::Center, B::FixedFrame, U>, SurfaceError>
    where
        U: LengthUnit,
        Quantity<U>: From<Quantity<Meter>>,
    {
        let [x, y, z] = self.body_fixed_cartesian()?;
        let xyz = [
            Quantity::<U>::from(x),
            Quantity::<U>::from(y),
            Quantity::<U>::from(z),
        ];
        if !xyz.iter().all(|component| component.value().is_finite()) {
            return Err(SurfaceError::NonFiniteInput {
                field: "observer position in output unit",
            });
        }
        if xyz.iter().all(|component| component.value() == 0.0) {
            return Err(SurfaceError::DegenerateGeometry {
                reason: "observer position is not representable in output unit",
            });
        }
        Ok(Position::<B::Center, B::FixedFrame, U>::new(
            xyz[0], xyz[1], xyz[2],
        ))
    }

    /// Convert to centre parameters for [`super::SurfaceObserver`].
    #[inline]
    pub fn to_params(&self) -> SurfaceObserverParams<B> {
        SurfaceObserverParams::from_site(self)
    }
}

impl SurfaceSite<Earth> {
    /// Build from an existing WGS84 [`Geodetic`]`<ECEF>` site.
    ///
    /// Uses the WGS84 ellipsoid and planetographic (geodetic) latitude.
    /// Longitude remains east-positive (Earth NAIF exception).
    pub fn from_geodetic(site: Geodetic<ECEF>) -> Result<Self, SurfaceError> {
        Self::with_surface(
            site.lon,
            site.lat,
            site.height.to::<Meter>(),
            Earth::default_surface(),
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
        if self.surface != Earth::default_surface() {
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
    use crate::bodies::solar_system::{Mars, Moon};
    use crate::qtty::{DEG, M};

    #[test]
    fn moon_site_defaults() {
        let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
        assert_eq!(site.body_name(), "Moon");
        let p = site.body_fixed_cartesian().unwrap();
        assert!((p[0].value() - 1_737_400.0).abs() < 1e-6);
    }

    #[test]
    fn rejects_bad_latitude() {
        let err = SurfaceSite::<Mars>::new(0.0 * DEG, 100.0 * DEG, 0.0 * M).unwrap_err();
        assert!(matches!(err, SurfaceError::InvalidLatitude { .. }));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn deserialization_enforces_site_invariants() {
        let site = SurfaceSite::<Moon>::new(20.0 * DEG, 30.0 * DEG, 0.0 * M).unwrap();
        let mut payload = serde_json::to_value(site).unwrap();
        payload["lat"] = serde_json::to_value(120.0 * DEG).unwrap();
        assert!(serde_json::from_value::<SurfaceSite<Moon>>(payload).is_err());

        let mut payload = serde_json::to_value(site).unwrap();
        payload["surface"] = serde_json::json!({"Sphere": {
            "radius": serde_json::to_value(-1.0 * M).unwrap()
        }});
        assert!(serde_json::from_value::<SurfaceSite<Moon>>(payload).is_err());
    }

    #[test]
    fn earth_geodetic_roundtrip() {
        let g = Geodetic::<ECEF>::new((-17.89) * DEG, 28.75 * DEG, 2200.0 * M);
        let site = SurfaceSite::<Earth>::from_geodetic(g).unwrap();
        let back = site.to_geodetic().unwrap();
        assert!((back.lat.value() - g.lat.value()).abs() < 1e-12);
        assert!((back.lon.value() - g.lon.value()).abs() < 1e-12);
    }

    #[test]
    fn longitude_normalization_and_checked_geodetic_conversion() {
        use crate::coordinates::surface::{LatitudeConvention, ReferenceSurface};
        let wrapped = SurfaceSite::<Moon>::new(370.0 * DEG, 0.0 * DEG, -100.0 * M).unwrap();
        let negative = SurfaceSite::<Moon>::new(-190.0 * DEG, 0.0 * DEG, -100.0 * M).unwrap();
        assert_eq!(wrapped.lon().value(), 10.0);
        assert_eq!(negative.lon().value(), 170.0);
        assert_eq!(wrapped.height().value(), -100.0);
        assert_ne!(wrapped, negative);

        let custom = SurfaceSite::<Earth>::with_surface(
            0.0 * DEG,
            10.0 * DEG,
            0.0 * M,
            ReferenceSurface::sphere(6_371_000.0 * M).unwrap(),
            LatitudeConvention::Planetographic,
        )
        .unwrap();
        assert!(custom.to_geodetic().is_err());

        let centric = SurfaceSite::<Earth>::with_surface(
            0.0 * DEG,
            10.0 * DEG,
            0.0 * M,
            Earth::default_surface(),
            LatitudeConvention::Planetocentric,
        )
        .unwrap();
        assert!(centric.to_geodetic().is_err());
    }

    #[test]
    fn mars_from_planetographic_west_longitude() {
        let site =
            SurfaceSite::<Mars>::from_planetographic_longitude(137.4 * DEG, (-4.6) * DEG, 0.0 * M)
                .unwrap();
        assert!((site.lon().value() + 137.4).abs() < 1e-12);
        assert!((site.planetographic_longitude().value() - 137.4).abs() < 1e-12);
    }

    #[test]
    fn site_nonfinite_input_is_rejected_at_construction() {
        let longitude = SurfaceSite::<Moon>::new(f64::INFINITY * DEG, 0.0 * DEG, 0.0 * M);
        let latitude = SurfaceSite::<Moon>::new(0.0 * DEG, f64::NAN * DEG, 0.0 * M);
        let height = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, f64::NAN * M);
        assert!(matches!(
            longitude,
            Err(SurfaceError::NonFiniteInput { .. })
        ));
        assert!(matches!(latitude, Err(SurfaceError::NonFiniteInput { .. })));
        assert!(matches!(height, Err(SurfaceError::NonFiniteInput { .. })));
    }
}
