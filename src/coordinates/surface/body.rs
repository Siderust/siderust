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
    MarsFixed, MercuryFixed, MoonIauMeanEarth, PlutoFixed, ReferenceFrame, VenusFixed, ECEF,
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
    /// A specialised gas-giant implementor can opt in after defining a
    /// physical reference surface (e.g. a 1-bar level) and fixed frame.
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
                reason: "no documented reference surface; implement a body model with an explicit surface",
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
/// The associated body-fixed frame uses the complete NAIF text-PCK
/// `IAU_MOON` analytical orientation (IAU 2009): Mean-Earth/Polar-Axis
/// approximation, including 13 periodic terms and quadratic W. The
/// `MoonPrincipalAxes` type is retained for compatibility and can also be
/// named `MoonIauMeanEarth`; neither is the JPL DE `MOON_PA` frame nor
/// the DE-specific binary-PCK `MOON_ME` realization. Model evaluation
/// uses numerical J2000 TDB days; a TT epoch is an approximation.
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
// Canonical solar-system body markers (legacy *Surface names are type aliases)
// =============================================================================

/// Earth observer marker, identical to [`Earth`].
pub type EarthSurface = Earth;
/// Moon observer marker, identical to [`Moon`].
pub type MoonSurface = Moon;
/// Mars observer marker, identical to [`Mars`].
pub type MarsSurface = Mars;
/// Mercury observer marker, identical to [`Mercury`].
pub type MercurySurface = Mercury;
/// Venus observer marker, identical to [`Venus`].
pub type VenusSurface = Venus;
/// Pluto observer marker, identical to [`Pluto`].
pub type PlutoSurface = Pluto;

impl BodySurfaceModel for Earth {
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

impl BodySurfaceModel for Moon {
    type Center = Selenocentric;
    type FixedFrame = MoonIauMeanEarth;

    fn body_name() -> &'static str {
        "Moon"
    }
    fn default_surface() -> ReferenceSurface {
        MOON_IAU_SPHERE
    }
    fn default_latitude_convention() -> LatitudeConvention {
        LatitudeConvention::Planetocentric
    }
}

impl BodySurfaceModel for Mars {
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

impl BodySurfaceModel for Mercury {
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

impl BodySurfaceModel for Venus {
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

impl BodySurfaceModel for Pluto {
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

    #[test]
    fn default_surface_and_latitude_conventions_for_supported_rocky_bodies() {
        use crate::coordinates::surface::SurfaceSite;
        use crate::qtty::{DEG, M};

        assert_eq!(EarthSurface::body_name(), "Earth");
        assert_eq!(MoonSurface::body_name(), "Moon");
        assert_eq!(MarsSurface::body_name(), "Mars");
        assert_eq!(MercurySurface::body_name(), "Mercury");
        assert_eq!(VenusSurface::body_name(), "Venus");
        assert_eq!(PlutoSurface::body_name(), "Pluto");

        assert_eq!(
            EarthSurface::default_latitude_convention(),
            LatitudeConvention::Planetographic
        );
        assert_eq!(
            MarsSurface::default_latitude_convention(),
            LatitudeConvention::Planetographic
        );
        for latitude in [
            MoonSurface::default_latitude_convention(),
            MercurySurface::default_latitude_convention(),
            VenusSurface::default_latitude_convention(),
            PlutoSurface::default_latitude_convention(),
        ] {
            assert_eq!(latitude, LatitudeConvention::Planetocentric);
        }

        let mercury = SurfaceSite::<MercurySurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M)
            .unwrap();
        let venus = SurfaceSite::<VenusSurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M)
            .unwrap();
        let pluto = SurfaceSite::<PlutoSurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M)
            .unwrap();
        assert_eq!(mercury.body_fixed_cartesian().unwrap()[0].value(), 2_439_700.0);
        assert_eq!(venus.body_fixed_cartesian().unwrap()[0].value(), 6_051_800.0);
        assert_eq!(pluto.body_fixed_cartesian().unwrap()[0].value(), 1_188_300.0);
    }
}
