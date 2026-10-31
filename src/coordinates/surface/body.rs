// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Body → surface / frame / center associations and longitude conventions.

use super::reference_surface::ReferenceSurface;
use crate::bodies::solar_system::{Earth, Mars, Mercury, Moon, Pluto, Venus};
use crate::coordinates::centers::{
    Geocentric, Marscentric, Mercurycentric, Plutocentric, ReferenceCenter, Selenocentric,
    Venuscentric,
};
use crate::coordinates::frames::{
    MarsFixed, MercuryFixed, MoonIauMeanEarth, PlutoFixed, ReferenceFrame, VenusFixed, ECEF,
};
use crate::qtty::{Degrees, Meters};
use affn::AffineCenter;
use core::fmt::Debug;

/// Latitude convention for a surface site.
///
/// Planetocentric and planetographic latitudes coincide on a sphere and
/// differ on an ellipsoid. They must **not** be treated as interchangeable.
///
/// This enum covers **latitude only**. Longitude sense is independent: see
/// [`PlanetographicLongitudeSense`] and [`BodySurfaceModel::planetographic_longitude_sense`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum LatitudeConvention {
    /// Angle between the equatorial plane and the radius vector from the
    /// body centre to the reference-surface point (planetocentric latitude).
    /// Height then offsets that point along the surface normal; at nonzero
    /// height this is not the latitude of the displaced Cartesian position.
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

/// Sense of **planetographic** longitude in external cartographic systems
/// (notably NAIF SPICE `pgrrec` / `recpgr`).
///
/// Siderust surface sites always store longitude **east-positive**, matching
/// planetocentric / body-fixed Cartesian conventions (`x = ρ cos φ cos λ`,
/// `y = ρ cos φ sin λ` with λ east). Never apply a hidden sign flip inside
/// [`super::SurfaceSite`] constructors.
///
/// ## Default NAIF / WGCCRE rules (what Siderust encodes)
///
/// Distinct from latitude. Under the usual SPICE defaults:
/// - **West-positive** for most prograde rotators (e.g. Mars, Mercury).
/// - **East-positive** for retrograde rotators (e.g. Venus; negative `W` rate).
/// - **East-positive** exceptions for Earth, Moon, and Sun.
///
/// Pluto is treated here as west-positive because its published `W` rate is
/// positive under NAIF’s spin-sense rule. Cartographic practice for dwarf
/// planets can vary; always confirm the convention of your external source.
///
/// Planetocentric longitude is always east-positive in SPICE (`latrec` /
/// `reclat`) and needs no conversion for Siderust.
///
/// ## SPICE kernel-pool overrides are not automatic
///
/// SPICE may redefine the planetographic longitude sense via kernel-pool
/// assignments such as `BODY<body ID>_PGR_POSITIVE_LON = 'EAST'|'WEST'`
/// (see NAIF `pgrrec` / `plnsns` documentation). Siderust’s
/// [`BodySurfaceModel::planetographic_longitude_sense`] reflects the **default**
/// body rules above; it does **not** read SPICE kernel pools or track such
/// overrides. If an external dataset was produced under a non-default
/// `PGR_POSITIVE_LON` setting, convert with an explicit
/// [`PlanetographicLongitudeSense`] rather than assuming the body default.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum PlanetographicLongitudeSense {
    /// Longitude increases eastward (Earth, Moon, Venus, …).
    EastPositive,
    /// Longitude increases westward (Mars, Mercury, Pluto, …).
    WestPositive,
}

impl PlanetographicLongitudeSense {
    /// Human-readable label.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::EastPositive => "east-positive",
            Self::WestPositive => "west-positive",
        }
    }
}

/// Convert an external planetographic longitude into Siderust's east-positive
/// internal longitude.
///
/// For [`PlanetographicLongitudeSense::EastPositive`] this is the identity.
/// For [`PlanetographicLongitudeSense::WestPositive`], `λ_east = −λ_west`
/// (then normalised by [`super::SurfaceSite`] construction).
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::surface::{
///     east_positive_from_planetographic, PlanetographicLongitudeSense,
/// };
/// use siderust::qtty::DEG;
///
/// // SPICE Mars planetographic west longitude 137.4° → Siderust −137.4° east.
/// let east = east_positive_from_planetographic(
///     137.4 * DEG,
///     PlanetographicLongitudeSense::WestPositive,
/// );
/// assert!((east.value() + 137.4).abs() < 1e-12);
/// ```
#[inline]
pub fn east_positive_from_planetographic(
    lon_planetographic: Degrees,
    sense: PlanetographicLongitudeSense,
) -> Degrees {
    match sense {
        PlanetographicLongitudeSense::EastPositive => lon_planetographic,
        PlanetographicLongitudeSense::WestPositive => Degrees::new(-lon_planetographic.value()),
    }
}

/// Convert a Siderust east-positive longitude into an external planetographic
/// longitude with the given sense.
#[inline]
pub fn planetographic_from_east_positive(
    lon_east_positive: Degrees,
    sense: PlanetographicLongitudeSense,
) -> Degrees {
    east_positive_from_planetographic(lon_east_positive, sense)
}

/// Compile-time association of a celestial body with its surface geometry,
/// body-centred origin, and body-fixed frame.
///
/// Implementors must document the reference-surface provenance and the
/// mathematical meaning of the body-fixed frame (e.g. IAU WGCCRE analytical
/// orientation versus a JPL DE principal-axes realisation).
///
/// # Trait contracts
///
/// - [`Self::default_surface`] must return a surface that passes
///   [`ReferenceSurface::validate`] and admits a geometrically valid site at
///   longitude = latitude = height = 0. This underpins
///   [`super::SurfaceObserverParams`]`::default`, required by affn's
///   `ReferenceCenter::Params: Default` bound.
/// - Bodies without a documented solid/cartographic reference surface (e.g.
///   gas giants without an explicit 1-bar model) simply **do not** implement
///   this trait. There is no runtime “unsupported” flag.
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

    /// NAIF/WGCCRE planetographic **longitude** sense for this body.
    ///
    /// Used only when converting external planetographic longitudes into
    /// Siderust's east-positive internal convention. It does not alter
    /// [`super::SurfaceSite`] storage.
    fn planetographic_longitude_sense() -> PlanetographicLongitudeSense;

    /// Convert an external planetographic longitude into east-positive form
    /// for this body.
    #[inline]
    fn east_positive_longitude(lon_planetographic: Degrees) -> Degrees {
        east_positive_from_planetographic(
            lon_planetographic,
            Self::planetographic_longitude_sense(),
        )
    }

    /// Convert a Siderust east-positive longitude into this body's external
    /// planetographic longitude sense.
    #[inline]
    fn planetographic_longitude(lon_east_positive: Degrees) -> Degrees {
        planetographic_from_east_positive(lon_east_positive, Self::planetographic_longitude_sense())
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
/// The associated body-fixed frame is the NAIF text-PCK `IAU_MOON` analytical
/// Mean-Earth/Polar-Axis orientation (IAU 2009), including 13 periodic terms
/// and quadratic W. The [`MoonIauMeanEarth`] marker is distinct from
/// [`crate::coordinates::frames::MoonPrincipalAxes`]. DE-specific PA and ME
/// realizations require binary PCK and frame kernels.
/// Model evaluation converts the public TT epoch to TDB using tempoch.
pub const MOON_IAU_SPHERE: ReferenceSurface = ReferenceSurface::Sphere {
    radius: Meters::new(1_737_400.0),
};

/// IAU 2015 Mars reference ellipsoid: `a = 3 396.19 km`, `f ≈ 1/169.894`.
///
/// Derived from IAU equatorial radius 3396.19 km and polar radius 3376.20 km
/// (`f = (a−b)/a`). Latitude default is **planetographic** (areographic).
/// Mars planetographic **longitude** in SPICE is west-positive; Siderust still
/// stores east-positive longitudes (see [`PlanetographicLongitudeSense`]).
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
    fn planetographic_longitude_sense() -> PlanetographicLongitudeSense {
        // NAIF exception: Earth planetographic longitude is east-positive.
        PlanetographicLongitudeSense::EastPositive
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
    fn planetographic_longitude_sense() -> PlanetographicLongitudeSense {
        // NAIF exception: Moon planetographic longitude is east-positive.
        PlanetographicLongitudeSense::EastPositive
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
    fn planetographic_longitude_sense() -> PlanetographicLongitudeSense {
        // Prograde rotator: SPICE planetographic longitude is west-positive.
        PlanetographicLongitudeSense::WestPositive
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
    fn planetographic_longitude_sense() -> PlanetographicLongitudeSense {
        PlanetographicLongitudeSense::WestPositive
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
    fn planetographic_longitude_sense() -> PlanetographicLongitudeSense {
        // Retrograde: SPICE planetographic longitude is east-positive.
        PlanetographicLongitudeSense::EastPositive
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
    fn planetographic_longitude_sense() -> PlanetographicLongitudeSense {
        // Positive W-dot ⇒ prograde under NAIF's spin-sense rule ⇒ west-positive.
        PlanetographicLongitudeSense::WestPositive
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::surface::SurfaceSite;
    use crate::qtty::{DEG, M};

    #[test]
    fn default_surfaces_admit_origin_sites() {
        // Trait contract underpinning SurfaceObserverParams::default.
        assert!(SurfaceSite::<Earth>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).is_ok());
        assert!(SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).is_ok());
        assert!(SurfaceSite::<Mars>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).is_ok());
        assert!(SurfaceSite::<Mercury>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).is_ok());
        assert!(SurfaceSite::<Venus>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).is_ok());
        assert!(SurfaceSite::<Pluto>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).is_ok());
    }

    #[test]
    fn mars_flattening_matches_iau() {
        let s = Mars::default_surface();
        let f = s.flattening();
        assert!((f - 0.005_886_007_9).abs() < 1e-9, "f={f}");
    }

    #[test]
    fn longitude_senses_match_naif_cartographic_rules() {
        assert_eq!(
            Earth::planetographic_longitude_sense(),
            PlanetographicLongitudeSense::EastPositive
        );
        assert_eq!(
            Moon::planetographic_longitude_sense(),
            PlanetographicLongitudeSense::EastPositive
        );
        assert_eq!(
            Venus::planetographic_longitude_sense(),
            PlanetographicLongitudeSense::EastPositive
        );
        assert_eq!(
            Mars::planetographic_longitude_sense(),
            PlanetographicLongitudeSense::WestPositive
        );
        assert_eq!(
            Mercury::planetographic_longitude_sense(),
            PlanetographicLongitudeSense::WestPositive
        );
        assert_eq!(
            Pluto::planetographic_longitude_sense(),
            PlanetographicLongitudeSense::WestPositive
        );
    }

    #[test]
    fn mars_spice_west_longitude_converts_to_east_positive_without_hidden_flips() {
        // Independent CSPICE pgrrec check (pck00011.tpc): Mars west-positive
        // planetographic (137.4°W, −4.6°) equals east-positive lon = −137.4°.
        let lon_west = 137.4 * DEG;
        let lon_east = Mars::east_positive_longitude(lon_west);
        assert!((lon_east.value() + 137.4).abs() < 1e-12);
        assert!((Mars::planetographic_longitude(lon_east).value() - 137.4).abs() < 1e-12);

        let site = SurfaceSite::<Mars>::new(lon_east, (-4.6) * DEG, 0.0 * M).unwrap();
        let p = site.body_fixed_cartesian().unwrap();
        // Metres; CSPICE pgrrec with a=3396.19 km, b=3376.20 km:
        // [-2491.967067275659, -2291.481184471932, -269.1840291055111] km
        assert_near(p[0].value() / 1000.0, -2491.967067275659, 1e-9);
        assert_near(p[1].value() / 1000.0, -2291.481184471932, 1e-9);
        assert_near(p[2].value() / 1000.0, -269.1840291055111, 1e-9);

        // Using the west value as if it were east places the site in the
        // opposite hemisphere — the conversion must stay explicit.
        let wrong = SurfaceSite::<Mars>::new(lon_west, (-4.6) * DEG, 0.0 * M).unwrap();
        let w = wrong.body_fixed_cartesian().unwrap();
        assert!((w[1].value() - p[1].value()).abs() > 1e6);
    }

    #[test]
    fn earth_and_venus_planetographic_longitude_is_identity() {
        let lon = 28.75 * DEG;
        assert_eq!(Earth::east_positive_longitude(lon), lon);
        assert_eq!(Venus::east_positive_longitude(lon), lon);
        assert_eq!(Moon::east_positive_longitude(lon), lon);
    }

    #[test]
    fn default_surface_and_latitude_conventions_for_rocky_bodies() {
        assert_eq!(Earth::body_name(), "Earth");
        assert_eq!(Moon::body_name(), "Moon");
        assert_eq!(Mars::body_name(), "Mars");
        assert_eq!(Mercury::body_name(), "Mercury");
        assert_eq!(Venus::body_name(), "Venus");
        assert_eq!(Pluto::body_name(), "Pluto");

        assert_eq!(
            Earth::default_latitude_convention(),
            LatitudeConvention::Planetographic
        );
        assert_eq!(
            Mars::default_latitude_convention(),
            LatitudeConvention::Planetographic
        );
        for latitude in [
            Moon::default_latitude_convention(),
            Mercury::default_latitude_convention(),
            Venus::default_latitude_convention(),
            Pluto::default_latitude_convention(),
        ] {
            assert_eq!(latitude, LatitudeConvention::Planetocentric);
        }

        let mercury = SurfaceSite::<Mercury>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
        let venus = SurfaceSite::<Venus>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
        let pluto = SurfaceSite::<Pluto>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
        assert_eq!(
            mercury.body_fixed_cartesian().unwrap()[0].value(),
            2_439_700.0
        );
        assert_eq!(
            venus.body_fixed_cartesian().unwrap()[0].value(),
            6_051_800.0
        );
        assert_eq!(
            pluto.body_fixed_cartesian().unwrap()[0].value(),
            1_188_300.0
        );
    }

    fn assert_near(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() < tol,
            "got {actual}, expected {expected} (tol {tol})"
        );
    }

    #[test]
    fn unsupported_surface_error_remains_available_for_diagnostics() {
        // Gas giants do not implement BodySurfaceModel; the error variant is
        // retained for custom models that reject construction explicitly.
        use super::super::error::SurfaceError;
        let err = SurfaceError::UnsupportedSurface {
            body: "Jupiter",
            reason: "no default 1-bar reference surface",
        };
        assert!(alloc::string::ToString::to_string(&err).contains("Jupiter"));
    }
}
