// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Spherical and ellipsoidal reference surfaces for body-surface observers.

use super::error::SurfaceError;
use crate::qtty::{Meter, Meters, Quantity};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Explicit reference surface used to convert surface coordinates to
/// body-fixed Cartesian positions and local verticals.
///
/// Surfaces are **runtime** parameters so callers can override defaults
/// (e.g. a mission-specific lunar radius) without proliferating marker types.
/// Body defaults are documented on [`super::BodySurfaceModel`].
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize))]
pub enum ReferenceSurface {
    /// Sphere of mean/reference radius `radius` (metres).
    ///
    /// Planetocentric and planetographic latitudes coincide on a sphere.
    Sphere {
        /// Reference radius in metres (must be finite and positive).
        radius: Meters,
    },
    /// Oblate reference ellipsoid with equatorial radius `a` and flattening `f`.
    ///
    /// Latitude interpretation follows [`super::LatitudeConvention`]:
    /// planetographic specifies the surface-normal direction; planetocentric
    /// specifies the centre-to-surface direction. Height and ENU up are
    /// always along the ellipsoid normal.
    Ellipsoid {
        /// Semi-major (equatorial) axis in metres.
        a: Meters,
        /// Flattening `f = (a − b) / a`, in `(0, 1)`.
        f: f64,
    },
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for ReferenceSurface {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        // Preserve the externally-tagged enum representation but forbid invalid
        // radii and flattening values from entering the typed geometry model.
        #[derive(Deserialize)]
        enum RawSurface {
            Sphere { radius: Meters },
            Ellipsoid { a: Meters, f: f64 },
        }
        match RawSurface::deserialize(deserializer)? {
            RawSurface::Sphere { radius } => Self::sphere(radius),
            RawSurface::Ellipsoid { a, f } => Self::ellipsoid(a, f),
        }
        .map_err(serde::de::Error::custom)
    }
}

impl ReferenceSurface {
    /// Validate a reference surface, including values constructed via public variants.
    pub fn validate(&self) -> Result<(), SurfaceError> {
        match *self {
            Self::Sphere { radius } => Self::sphere(radius).map(|_| ()),
            Self::Ellipsoid { a, f } => Self::ellipsoid(a, f).map(|_| ()),
        }
    }

    /// Construct a validated spherical surface.
    pub fn sphere(radius: Meters) -> Result<Self, SurfaceError> {
        let r = radius.value();
        if !r.is_finite() || r <= 0.0 {
            return Err(SurfaceError::InvalidDimension {
                parameter: "radius",
                value: r,
            });
        }
        Ok(Self::Sphere { radius })
    }

    /// Construct a validated ellipsoidal surface.
    pub fn ellipsoid(a: Meters, f: f64) -> Result<Self, SurfaceError> {
        let a_v = a.value();
        if !a_v.is_finite() || a_v <= 0.0 {
            return Err(SurfaceError::InvalidDimension {
                parameter: "semi_major_axis",
                value: a_v,
            });
        }
        if !f.is_finite() || f <= 0.0 || f >= 1.0 {
            return Err(SurfaceError::InvalidFlattening { value: f });
        }
        Ok(Self::Ellipsoid { a, f })
    }

    /// Equatorial / reference radius in metres.
    #[inline]
    pub fn equatorial_radius(&self) -> Meters {
        match *self {
            Self::Sphere { radius } => radius,
            Self::Ellipsoid { a, .. } => a,
        }
    }

    /// Flattening (`0` for a sphere).
    #[inline]
    pub fn flattening(&self) -> f64 {
        match *self {
            Self::Sphere { .. } => 0.0,
            Self::Ellipsoid { f, .. } => f,
        }
    }

    /// First eccentricity squared (`0` for a sphere).
    #[inline]
    pub fn e2(&self) -> f64 {
        let f = self.flattening();
        2.0 * f - f * f
    }

    /// Semi-minor axis in metres.
    #[inline]
    pub fn polar_radius(&self) -> Meters {
        match *self {
            Self::Sphere { radius } => radius,
            Self::Ellipsoid { a, f } => Quantity::<Meter>::new(a.value() * (1.0 - f)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::M;

    #[test]
    fn rejects_non_positive_sphere() {
        assert!(ReferenceSurface::sphere(0.0 * M).is_err());
        assert!(ReferenceSurface::sphere((-1.0) * M).is_err());
    }

    #[test]
    fn rejects_bad_flattening() {
        assert!(ReferenceSurface::ellipsoid(1.0 * M, 0.0).is_err());
        assert!(ReferenceSurface::ellipsoid(1.0 * M, 1.0).is_err());
        assert!(ReferenceSurface::ellipsoid(1.0 * M, f64::NAN).is_err());
    }

    #[cfg(feature = "serde")]
    #[test]
    fn deserialization_rejects_invalid_surface() {
        let surface = ReferenceSurface::sphere(1_000.0 * M).unwrap();
        let mut payload = serde_json::to_value(surface).unwrap();
        payload["Sphere"]["radius"] = serde_json::to_value(-1.0 * M).unwrap();
        assert!(serde_json::from_value::<ReferenceSurface>(payload).is_err());

        let surface = ReferenceSurface::ellipsoid(1_000.0 * M, 0.1).unwrap();
        let mut payload = serde_json::to_value(surface).unwrap();
        payload["Ellipsoid"]["f"] = serde_json::json!(1.0);
        assert!(serde_json::from_value::<ReferenceSurface>(payload).is_err());
    }

    #[test]
    fn sphere_e2_is_zero() {
        let s = ReferenceSurface::sphere(1_737_400.0 * M).unwrap();
        assert_eq!(s.e2(), 0.0);
        assert_eq!(s.polar_radius(), s.equatorial_radius());
    }
}
