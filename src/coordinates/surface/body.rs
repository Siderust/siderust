// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Body → surface / frame / center associations.

use super::error::SurfaceError;
use super::reference_surface::ReferenceSurface;
use crate::bodies::solar_system::{Earth, Mars, Mercury, Moon, Pluto, Venus};
use crate::coordinates::centers::{
    Geocentric, Marscentric, Mercurycentric, Plutocentric, ReferenceCenter, Selenocentric,
    Venuscentric,
};
use crate::coordinates::frames::{
    MarsFixed, MercuryFixed, MoonPrincipalAxes, PlutoFixed, ReferenceFrame, VenusFixed, ECEF,
};
use crate::qtty::Meters;
use affn::AffineCenter;
use core::fmt::Debug;

/// Latitude convention for a surface site.
///
/// Planetocentric and planetographic latitudes coincide on a sphere and
/// differ on an ellipsoid. They must **not** be treated as interchangeable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum LatitudeConvention {
    /// Angle between the equatorial plane and the radius vector from the
    /// body centre to the point (also called planetocentric latitude).
    Planetocentric,
    /// Angle between the equatorial plane and the surface normal
    /// (geodetic / planetographic latitude).
    Planetographic,
}

impl LatitudeConvention {
    /// Human-readable label.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Planetocentric => "planetocentric",
            Self::Planetographic => "planetographic",
        }
    }
}

/// Compile-time association of a celestial body with its surface geometry,
/// body-centred origin, and body-fixed frame.
///
/// Implementors must document the reference-surface provenance and the
/// mathematical meaning of the body-fixed frame (e.g. IAU WGCCRE analytical
/// orientation versus a JPL DE principal-axes realisation).
pub trait BodySurfaceModel: Copy + Clone + Debug + Send + Sync + 'static {
    /// Body-centred reference centre (`Params = ()`).
    type Center: ReferenceCenter<Params = ()> + AffineCenter;
    /// Body-fixed reference frame used for surface Cartesian coordinates.
    type FixedFrame: ReferenceFrame;

    /// Canonical body name for diagnostics and documentation.
    fn body_name() -> &'static str;

    /// Default reference surface for this body.
    fn default_surface() -> ReferenceSurface;

    /// Default latitude convention for [`Self::default_surface`].
    fn default_latitude_convention() -> LatitudeConvention;

    /// Whether a solid or cartographic reference surface is supported.
    ///
    /// Gas giants return `false` unless an explicit 1-bar (or similar) model
    /// is opted into by a specialised implementor.
    fn has_reference_surface() -> bool {
        true
    }

    /// Validate that this body exposes a usable reference surface.
    fn ensure_supported() -> Result<(), SurfaceError> {
        if Self::has_reference_surface() {
            Ok(())
        } else {
            Err(SurfaceError::UnsupportedSurface {
                body: Self::body_name(),
                reason: "no documented solid/reference surface; provide an explicit model",
            })
        }
    }
}

// =============================================================================
// Default surfaces (IAU / WGS84 provenance)
// =============================================================================

/// WGS 84 ellipsoid used for terrestrial [`ECEF`] geodesy in siderust.
///
/// `a = 6 378 137 m`, `f = 1/298.257 223 563` (NIMA TR8350.2).
pub const EARTH_WGS84_SURFACE: ReferenceSurface = ReferenceSurface::Ellipsoid {
    a: Meters::new(6_378_137.0),
    f: 1.0 / 298.257_223_563,
};

/// IAU mean lunar radius as a sphere: `R = 1 737.4 km`.
///
/// Source: Archinal et al. (2018) / IAU WGCCRE. This is a **spherical**
/// cartographic convenience model, not a lunar topographic DEM. Selenodetic
/// heights are measured above this sphere.
///
/// The associated body-fixed frame is [`MoonPrincipalAxes`] as realised by
/// siderust's IAU 2015 analytical `(α₀, δ₀, W)` polynomials — **not** the
/// JPL DE-series principal-axes frame and **not** the Mean Earth/Polar Axis
/// (ME) cartographic frame. Differences between these realisations can reach
/// arcminutes and must be accounted for when comparing to SPICE `MOON_PA` or
/// `MOON_ME` kernels.
pub const MOON_IAU_SPHERE: ReferenceSurface = ReferenceSurface::Sphere {
    radius: Meters::new(1_737_400.0),
};

/// IAU 2015 Mars reference ellipsoid: `a = 3 396.19 km`, `f ≈ 1/169.894`.
///
/// Derived from IAU equatorial radius 3396.19 km and polar radius 3376.20 km
/// (`f = (a−b)/a`). Latitude default is **planetographic** (areographic).
pub const MARS_IAU_ELLIPSOID: ReferenceSurface = ReferenceSurface::Ellipsoid {
    a: Meters::new(3_396_190.0),
    f: (3_396_190.0 - 3_376_200.0) / 3_396_190.0,
};

/// IAU mean Mercury radius (sphere), 2439.7 km.
pub const MERCURY_IAU_SPHERE: ReferenceSurface = ReferenceSurface::Sphere {
    radius: Meters::new(2_439_700.0),
};

/// IAU mean Venus radius (sphere), 6051.8 km.
pub const VENUS_IAU_SPHERE: ReferenceSurface = ReferenceSurface::Sphere {
    radius: Meters::new(6_051_800.0),
};

/// IAU mean Pluto radius (sphere), 1188.3 km.
pub const PLUTO_IAU_SPHERE: ReferenceSurface = ReferenceSurface::Sphere {
    radius: Meters::new(1_188_300.0),
};

// =============================================================================
// Marker wrappers (keep solar-system markers free of orphan-rule conflicts)
// =============================================================================

/// Earth surface-observer body tag (WGS84 / ECEF).
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Default)]
pub struct EarthSurface;

/// Moon surface-observer body tag (IAU sphere / MoonPrincipalAxes).
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Default)]
pub struct MoonSurface;

/// Mars surface-observer body tag (IAU ellipsoid / MarsFixed).
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Default)]
pub struct MarsSurface;

/// Mercury surface-observer body tag.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Default)]
pub struct MercurySurface;

/// Venus surface-observer body tag.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Default)]
pub struct VenusSurface;

/// Pluto surface-observer body tag.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Default)]
pub struct PlutoSurface;

impl BodySurfaceModel for EarthSurface {
    type Center = Geocentric;
    type FixedFrame = ECEF;

    fn body_name() -> &'static str {
        "Earth"
    }
    fn default_surface() -> ReferenceSurface {
        EARTH_WGS84_SURFACE
    }
    fn default_latitude_convention() -> LatitudeConvention {
        LatitudeConvention::Planetographic
    }
}

impl BodySurfaceModel for MoonSurface {
    type Center = Selenocentric;
    type FixedFrame = MoonPrincipalAxes;

    fn body_name() -> &'static str {
        "Moon"
    }
    fn default_surface() -> ReferenceSurface {
        MOON_IAU_SPHERE
    }
    fn default_latitude_convention() -> LatitudeConvention {
        // Sphere: planetocentric ≡ planetographic; document as planetocentric.
        LatitudeConvention::Planetocentric
    }
}

impl BodySurfaceModel for MarsSurface {
    type Center = Marscentric;
    type FixedFrame = MarsFixed;

    fn body_name() -> &'static str {
        "Mars"
    }
    fn default_surface() -> ReferenceSurface {
        MARS_IAU_ELLIPSOID
    }
    fn default_latitude_convention() -> LatitudeConvention {
        LatitudeConvention::Planetographic
    }
}

impl BodySurfaceModel for MercurySurface {
    type Center = Mercurycentric;
    type FixedFrame = MercuryFixed;

    fn body_name() -> &'static str {
        "Mercury"
    }
    fn default_surface() -> ReferenceSurface {
        MERCURY_IAU_SPHERE
    }
    fn default_latitude_convention() -> LatitudeConvention {
        LatitudeConvention::Planetocentric
    }
}

impl BodySurfaceModel for VenusSurface {
    type Center = Venuscentric;
    type FixedFrame = VenusFixed;

    fn body_name() -> &'static str {
        "Venus"
    }
    fn default_surface() -> ReferenceSurface {
        VENUS_IAU_SPHERE
    }
    fn default_latitude_convention() -> LatitudeConvention {
        LatitudeConvention::Planetocentric
    }
}

impl BodySurfaceModel for PlutoSurface {
    type Center = Plutocentric;
    type FixedFrame = PlutoFixed;

    fn body_name() -> &'static str {
        "Pluto"
    }
    fn default_surface() -> ReferenceSurface {
        PLUTO_IAU_SPHERE
    }
    fn default_latitude_convention() -> LatitudeConvention {
        LatitudeConvention::Planetocentric
    }
}

/// Convenience: treat the solar-system [`Earth`] marker as [`EarthSurface`].
impl BodySurfaceModel for Earth {
    type Center = Geocentric;
    type FixedFrame = ECEF;
    fn body_name() -> &'static str {
        EarthSurface::body_name()
    }
    fn default_surface() -> ReferenceSurface {
        EarthSurface::default_surface()
    }
    fn default_latitude_convention() -> LatitudeConvention {
        EarthSurface::default_latitude_convention()
    }
}

/// Convenience: treat the solar-system [`Moon`] marker as [`MoonSurface`].
impl BodySurfaceModel for Moon {
    type Center = Selenocentric;
    type FixedFrame = MoonPrincipalAxes;
    fn body_name() -> &'static str {
        MoonSurface::body_name()
    }
    fn default_surface() -> ReferenceSurface {
        MoonSurface::default_surface()
    }
    fn default_latitude_convention() -> LatitudeConvention {
        MoonSurface::default_latitude_convention()
    }
}

/// Convenience: treat the solar-system [`Mars`] marker as [`MarsSurface`].
impl BodySurfaceModel for Mars {
    type Center = Marscentric;
    type FixedFrame = MarsFixed;
    fn body_name() -> &'static str {
        MarsSurface::body_name()
    }
    fn default_surface() -> ReferenceSurface {
        MarsSurface::default_surface()
    }
    fn default_latitude_convention() -> LatitudeConvention {
        MarsSurface::default_latitude_convention()
    }
}

/// Convenience: [`Mercury`] → [`MercurySurface`].
impl BodySurfaceModel for Mercury {
    type Center = Mercurycentric;
    type FixedFrame = MercuryFixed;
    fn body_name() -> &'static str {
        MercurySurface::body_name()
    }
    fn default_surface() -> ReferenceSurface {
        MercurySurface::default_surface()
    }
    fn default_latitude_convention() -> LatitudeConvention {
        MercurySurface::default_latitude_convention()
    }
}

/// Convenience: [`Venus`] → [`VenusSurface`].
impl BodySurfaceModel for Venus {
    type Center = Venuscentric;
    type FixedFrame = VenusFixed;
    fn body_name() -> &'static str {
        VenusSurface::body_name()
    }
    fn default_surface() -> ReferenceSurface {
        VenusSurface::default_surface()
    }
    fn default_latitude_convention() -> LatitudeConvention {
        VenusSurface::default_latitude_convention()
    }
}

/// Convenience: [`Pluto`] → [`PlutoSurface`].
impl BodySurfaceModel for Pluto {
    type Center = Plutocentric;
    type FixedFrame = PlutoFixed;
    fn body_name() -> &'static str {
        PlutoSurface::body_name()
    }
    fn default_surface() -> ReferenceSurface {
        PlutoSurface::default_surface()
    }
    fn default_latitude_convention() -> LatitudeConvention {
        PlutoSurface::default_latitude_convention()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn support_matrix_rocky_bodies() {
        assert!(EarthSurface::has_reference_surface());
        assert!(MoonSurface::has_reference_surface());
        assert!(MarsSurface::has_reference_surface());
        assert!(EarthSurface::ensure_supported().is_ok());
    }

    #[test]
    fn mars_flattening_matches_iau() {
        let s = MarsSurface::default_surface();
        let f = s.flattening();
        assert!((f - 0.005_886_007_9).abs() < 1e-9, "f={f}");
    }
}
