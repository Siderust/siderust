// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! C-compatible struct definitions for siderust-ffi.
//!
//! These `#[repr(C)]` types are the data bridge between C/C++ callers and
//! the rich generic Rust types in siderust.
//!
//! # Shared Coordinate Model
//!
//! The coordinate types in this module implement the `affn` coordinate semantics:
//!
//! - **Frames**: [`SiderustFrame`] enumerates supported reference frames
//! - **Centers**: [`SiderustCenter`] enumerates supported reference centers
//! - **Units**: [`SiderustLengthUnit`] enumerates supported length units
//!
//! Adapters (C++, Python, JavaScript) should use these enum values rather than
//! runtime strings when possible. For string-based adapter APIs:
//!
//! | Frame String          | `SiderustFrame` Variant  |
//! |-----------------------|--------------------------|
//! | `"ICRS"`              | `ICRS = 1`               |
//! | `"EclipticMeanJ2000"` | `EclipticMeanJ2000 = 2`  |
//! | `"EquatorialMeanJ2000"`| `EquatorialMeanJ2000 = 3`|
//! | `"EquatorialMeanOfDate"`| `EquatorialMeanOfDate = 4`|
//! | `"EquatorialTrueOfDate"`| `EquatorialTrueOfDate = 5`|
//! | `"Horizontal"`        | `Horizontal = 6`         |
//! | `"ECEF"`              | `ECEF = 7`               |
//! | `"Galactic"`          | `Galactic = 8`           |
//! | `"GCRS"`              | `GCRS = 9`               |
//! | `"ICRF"`              | `ICRF = 15`              |
//!
//! | Center String      | `SiderustCenter` Variant |
//! |--------------------|--------------------------|
//! | `"Barycentric"`    | `Barycentric = 1`        |
//! | `"Heliocentric"`   | `Heliocentric = 2`       |
//! | `"Geocentric"`     | `Geocentric = 3`         |
//! | `"Topocentric"`    | `Topocentric = 4`        |
//! | `"Bodycentric"`    | `Bodycentric = 5`        |
//!
//! # Coordinate Arithmetic
//!
//! Following `affn` semantics:
//! - `Position - Position` yields a displacement/vector (not another position)
//! - `Position + Position` is **not valid**
//! - `Position + Displacement` yields a position
//! - `Displacement + Displacement` yields a displacement
//!
//! Adapters should enforce these rules at their language boundary.

use crate::ffi_utils::FfiFrom;
use qtty::*;
use siderust::calculus::azimuth::{
    AzimuthCrossingDirection, AzimuthCrossingEvent, AzimuthExtremum, AzimuthExtremumKind,
};
use siderust::coordinates::centers::{
    BodycentricParams as RustBodycentricParams, Geodetic,
    OrbitReferenceCenter as RustOrbitRefCenter,
};
use siderust::coordinates::frames::ECEF;
use tempoch::{Interval, JulianDate, ModifiedJulianDate, Period, MJD};

// Re-export tempoch-ffi types so the generated header can reference them.
// The extern crate declaration is needed because the dep name maps through a hyphen.
pub use ::tempoch_ffi::TempochPeriodMjd;

// ═══════════════════════════════════════════════════════════════════════════
// Enumerations
// ═══════════════════════════════════════════════════════════════════════════

/// Reference frame identifier for C interop.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustFrame {
    /// International Celestial Reference System.
    ICRS = 1,
    /// Mean ecliptic of J2000.0.
    EclipticMeanJ2000 = 2,
    /// Mean equatorial of J2000.0.
    EquatorialMeanJ2000 = 3,
    /// Mean equatorial of date.
    EquatorialMeanOfDate = 4,
    /// True equatorial of date (includes nutation).
    EquatorialTrueOfDate = 5,
    /// Local horizontal (azimuth/altitude).
    Horizontal = 6,
    /// Earth-Centred Earth-Fixed.
    ECEF = 7,
    /// Galactic coordinate system.
    Galactic = 8,
    /// Geocentric Celestial Reference System.
    GCRS = 9,
    /// Ecliptic of date.
    EclipticOfDate = 10,
    /// True ecliptic of date.
    EclipticTrueOfDate = 11,
    /// Celestial Intermediate Reference System.
    CIRS = 12,
    /// Terrestrial Intermediate Reference System.
    TIRS = 13,
    /// International Terrestrial Reference Frame.
    ITRF = 14,
    /// International Celestial Reference Frame.
    ICRF = 15,
}

impl SiderustFrame {
    /// Parse a frame name string (case-insensitive).
    ///
    /// Accepts the canonical names documented in the module header:
    /// `"icrs"`, `"ecliptic_mean_j2000"`, `"equatorial_mean_j2000"`, etc.
    /// Also accepts common aliases like `"horizontal"`, `"galactic"`, `"gcrs"`.
    ///
    /// Returns `None` if the string does not match any known frame.
    pub fn from_str(s: &str) -> Option<Self> {
        let lower = s.to_ascii_lowercase();
        match lower.as_str() {
            "icrs" => Some(Self::ICRS),
            "ecliptic_mean_j2000" | "eclipticmeanj2000" | "mean_ecliptic_j2000" => {
                Some(Self::EclipticMeanJ2000)
            }
            "equatorial_mean_j2000" | "equatorialmeanj2000" | "mean_equatorial_j2000" => {
                Some(Self::EquatorialMeanJ2000)
            }
            "equatorial_mean_of_date" | "equatorialmeanofdate" | "mean_of_date" => {
                Some(Self::EquatorialMeanOfDate)
            }
            "equatorial_true_of_date" | "equatorialtrueofdate" | "true_of_date" => {
                Some(Self::EquatorialTrueOfDate)
            }
            "horizontal" | "altaz" | "local_horizontal" => Some(Self::Horizontal),
            "ecef" | "earth_fixed" => Some(Self::ECEF),
            "galactic" | "gal" => Some(Self::Galactic),
            "gcrs" | "geocentric_celestial" => Some(Self::GCRS),
            "ecliptic_of_date" | "eclipticofdate" => Some(Self::EclipticOfDate),
            "ecliptic_true_of_date" | "ecliptictrueofdate" => Some(Self::EclipticTrueOfDate),
            "cirs" | "celestial_intermediate" => Some(Self::CIRS),
            "tirs" | "terrestrial_intermediate" => Some(Self::TIRS),
            "itrf" | "terrestrial_reference" => Some(Self::ITRF),
            "icrf" | "celestial_reference" => Some(Self::ICRF),
            _ => None,
        }
    }

    /// Return the canonical string name of this frame.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::ICRS => "icrs",
            Self::EclipticMeanJ2000 => "ecliptic_mean_j2000",
            Self::EquatorialMeanJ2000 => "equatorial_mean_j2000",
            Self::EquatorialMeanOfDate => "equatorial_mean_of_date",
            Self::EquatorialTrueOfDate => "equatorial_true_of_date",
            Self::Horizontal => "horizontal",
            Self::ECEF => "ecef",
            Self::Galactic => "galactic",
            Self::GCRS => "gcrs",
            Self::EclipticOfDate => "ecliptic_of_date",
            Self::EclipticTrueOfDate => "ecliptic_true_of_date",
            Self::CIRS => "cirs",
            Self::TIRS => "tirs",
            Self::ITRF => "itrf",
            Self::ICRF => "icrf",
        }
    }
}

/// Reference center identifier for C interop.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustCenter {
    /// Solar-system barycentre.
    Barycentric = 1,
    /// Sun centre.
    Heliocentric = 2,
    /// Earth centre.
    Geocentric = 3,
    /// Observer site on the Earth's surface.
    Topocentric = 4,
    /// Centre of a specific body.
    Bodycentric = 5,
}

impl SiderustCenter {
    /// Parse a center name string (case-insensitive).
    ///
    /// Accepts the canonical names documented in the module header:
    /// `"barycentric"`, `"heliocentric"`, `"geocentric"`, `"topocentric"`, `"bodycentric"`.
    /// Also accepts common aliases like `"ssb"`, `"solar"`, `"earth"`, `"observer"`.
    ///
    /// Returns `None` if the string does not match any known center.
    pub fn from_str(s: &str) -> Option<Self> {
        let lower = s.to_ascii_lowercase();
        match lower.as_str() {
            "barycentric" | "ssb" | "solar_system_barycenter" => Some(Self::Barycentric),
            "heliocentric" | "solar" | "sun" => Some(Self::Heliocentric),
            "geocentric" | "earth" => Some(Self::Geocentric),
            "topocentric" | "observer" | "site" => Some(Self::Topocentric),
            "bodycentric" | "body" => Some(Self::Bodycentric),
            _ => None,
        }
    }

    /// Return the canonical string name of this center.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Barycentric => "barycentric",
            Self::Heliocentric => "heliocentric",
            Self::Geocentric => "geocentric",
            Self::Topocentric => "topocentric",
            Self::Bodycentric => "bodycentric",
        }
    }
}

/// Crossing event direction.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustCrossingDirection {
    /// The body is crossing upward through the threshold.
    Rising = 0,
    /// The body is crossing downward through the threshold.
    Setting = 1,
}

/// Culmination event kind.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustCulminationKind {
    /// Upper culmination (maximum altitude).
    Max = 0,
    /// Lower culmination (minimum altitude).
    Min = 1,
}

/// Asteroid taxonomic class.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustAsteroidClass {
    /// Main-belt asteroid.
    MainBelt = 0,
    /// Near-Earth asteroid.
    NearEarth = 1,
    /// Trojan asteroid (co-orbital with a planet).
    Trojan = 2,
    /// Centaur (orbiting between Jupiter and Neptune).
    Centaur = 3,
    /// Trans-Neptunian object.
    TransNeptunian = 4,
    /// Dwarf planet.
    DwarfPlanet = 5,
}

/// Orbit reference frame for comets.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustOrbitFrame {
    /// Orbit defined relative to the Sun.
    Heliocentric = 0,
    /// Orbit defined relative to the solar-system barycentre.
    Barycentric = 1,
}

/// Length unit for coordinate positions.
///
/// Specifies the unit of measure for coordinate distances (x, y, z components
/// in Cartesian positions, or distance in spherical positions).
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SiderustLengthUnit {
    /// Astronomical Units (≈149.6 million km).
    #[default]
    AU = 0,
    /// Kilometres.
    Km = 1,
    /// Light-years.
    LightYear = 2,
    /// Parsecs.
    Parsec = 3,
    /// Metres.
    Meter = 4,
}

/// Solar-system body identifier for generic altitude/azimuth dispatch.
///
/// Each variant maps to a concrete unit type in `siderust::bodies::solar_system`.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustBody {
    /// The Sun.
    Sun = 0,
    /// Earth's Moon.
    Moon = 1,
    /// Mercury.
    Mercury = 2,
    /// Venus.
    Venus = 3,
    /// Mars.
    Mars = 4,
    /// Jupiter.
    Jupiter = 5,
    /// Saturn.
    Saturn = 6,
    /// Uranus.
    Uranus = 7,
    /// Neptune.
    Neptune = 8,
}

/// Subject kind discriminant for [`SiderustSubject`].
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustSubjectKind {
    /// Solar-system body (the `body` field is valid).
    Body = 0,
    /// Star opaque handle (the `star_handle` field is valid).
    Star = 1,
    /// Fixed ICRS direction (the `icrs_dir` field is valid).
    Icrs = 2,
    /// Target opaque handle (the `target_handle` field is valid).
    Target = 3,
    /// Generic target opaque handle (the `generic_target_handle` field is valid).
    GenericTarget = 4,
}

/// Unified subject for altitude / azimuth / tracking computations.
///
/// A tagged struct that can represent any entity on which altitude and
/// azimuth queries can be performed: solar-system bodies, catalog stars,
/// fixed ICRS directions, or opaque Target handles.
///
/// Only the field corresponding to `kind` is valid:
///
/// | `kind`           | Valid field(s)           |
/// |------------------|--------------------------|
/// | `Body`           | `body`                   |
/// | `Star`           | `star_handle`            |
/// | `Icrs`           | `icrs_dir`               |
/// | `Target`         | `target_handle`          |
/// | `GenericTarget`  | `generic_target_handle`  |
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustSubject {
    /// Discriminant selecting which field is active.
    pub kind: SiderustSubjectKind,
    /// Solar-system body discriminant.  Valid when `kind == Body`.
    pub body: SiderustBody,
    /// Opaque star handle (non-null).  Valid when `kind == Star`.
    pub star_handle: *const crate::bodies::SiderustStar,
    /// ICRS spherical direction.  Valid when `kind == Icrs`.
    pub icrs_dir: SiderustSphericalDir,
    /// Opaque target handle (non-null).  Valid when `kind == Target`.
    pub target_handle: *const crate::target::SiderustTarget,
    /// Opaque generic target handle.  Valid when `kind == GenericTarget`.
    pub generic_target_handle: *const crate::target::SiderustGenericTarget,
}

// SAFETY: `SiderustSubject` contains raw pointers (`star_handle`,
// `target_handle`) which prevent the auto-derived `Send`/`Sync` impls.
// These impls are sound because:
//
//  1. The pointees (`SiderustStar`, `SiderustTarget`) are heap-allocated
//     immutable objects created via `Box::into_raw`.  They are never
//     mutated after creation, so sharing across threads is safe.
//
//  2. The pointers are only dereferenced inside `dispatch_subject!` which
//     runs within an `ffi_guard!` (= `catch_unwind`) block.  A null check
//     precedes every dereference.
//
//  3. Ownership of the pointed-to object is held by the C caller.  The
//     caller must ensure the handle remains valid for the lifetime of the
//     `SiderustSubject` that borrows it.  This is an invariant of the
//     extern "C" API contract, documented on each function.
unsafe impl Send for SiderustSubject {}
unsafe impl Sync for SiderustSubject {}

impl SiderustSubject {
    /// Construct a `Body` subject.
    pub(crate) fn body(body: SiderustBody) -> Self {
        Self {
            kind: SiderustSubjectKind::Body,
            body,
            star_handle: std::ptr::null(),
            icrs_dir: SiderustSphericalDir::zeroed(),
            target_handle: std::ptr::null(),
            generic_target_handle: std::ptr::null(),
        }
    }

    /// Construct a `Star` subject (borrows the opaque handle).
    pub(crate) fn star(handle: *const crate::bodies::SiderustStar) -> Self {
        Self {
            kind: SiderustSubjectKind::Star,
            body: SiderustBody::Sun,
            star_handle: handle,
            icrs_dir: SiderustSphericalDir::zeroed(),
            target_handle: std::ptr::null(),
            generic_target_handle: std::ptr::null(),
        }
    }

    /// Construct an `Icrs` subject from a spherical direction.
    pub(crate) fn icrs(dir: SiderustSphericalDir) -> Self {
        Self {
            kind: SiderustSubjectKind::Icrs,
            body: SiderustBody::Sun,
            star_handle: std::ptr::null(),
            icrs_dir: dir,
            target_handle: std::ptr::null(),
            generic_target_handle: std::ptr::null(),
        }
    }

    /// Construct a `Target` subject (borrows the opaque handle).
    pub(crate) fn target(handle: *const crate::target::SiderustTarget) -> Self {
        Self {
            kind: SiderustSubjectKind::Target,
            body: SiderustBody::Sun,
            star_handle: std::ptr::null(),
            icrs_dir: SiderustSphericalDir::zeroed(),
            target_handle: handle,
            generic_target_handle: std::ptr::null(),
        }
    }

    /// Construct a `GenericTarget` subject (borrows the opaque handle).
    #[allow(dead_code)]
    pub(crate) fn generic_target(handle: *const crate::target::SiderustGenericTarget) -> Self {
        Self {
            kind: SiderustSubjectKind::GenericTarget,
            body: SiderustBody::Sun,
            star_handle: std::ptr::null(),
            icrs_dir: SiderustSphericalDir::zeroed(),
            target_handle: std::ptr::null(),
            generic_target_handle: handle,
        }
    }
}

/// Proper motion RA convention.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustRaConvention {
    /// True RA rate µα (deg/yr).
    MuAlpha = 0,
    /// Catalog rate µα★ = µα cos(δ) (deg/yr).
    MuAlphaStar = 1,
}

/// Coordinate kind discriminant for [`SiderustTargetCoord`].
///
/// Specifies which coordinate representation is stored in the target.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustTargetCoordKind {
    /// Spherical direction (RA/Dec) without distance.
    SphericalDir = 0,
    /// Spherical position with distance.
    SphericalPos = 1,
    /// Cartesian position (x, y, z).
    CartesianPos = 2,
}

/// Union of coordinate types for [`SiderustGenericTargetData`].
///
/// Callers must check `kind` to determine which field is valid.
/// Only the field corresponding to `kind` contains valid data:
///
/// | `kind`        | Valid field    |
/// |---------------|----------------|
/// | `SphericalDir`| `spherical_dir`|
/// | `SphericalPos`| `spherical_pos`|
/// | `CartesianPos`| `cartesian_pos`|
#[repr(C)]
#[derive(Clone, Copy)]
pub union SiderustTargetCoordUnion {
    /// Spherical direction (valid when kind == SphericalDir).
    pub spherical_dir: SiderustSphericalDir,
    /// Spherical position (valid when kind == SphericalPos).
    pub spherical_pos: SiderustSphericalPos,
    /// Cartesian position (valid when kind == CartesianPos).
    pub cartesian_pos: SiderustCartesianPos,
}

/// Generic target data, stored as a flat struct for FFI.
///
/// Represents a point in celestial coordinates at a specific epoch,
/// optionally with proper motion. This struct mirrors the Rust
/// `CoordinateWithPM<T>` type but in a C-compatible form.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct SiderustGenericTargetData {
    /// Which coordinate representation is stored.
    pub kind: SiderustTargetCoordKind,
    /// Padding for alignment.
    pub _pad1: [u8; 4],
    /// The coordinate data (union, check `kind` to determine which field).
    pub coord: SiderustTargetCoordUnion,
    /// Epoch as a Julian Date.
    pub epoch_jd: f64,
    /// Whether proper motion is present.
    pub has_proper_motion: bool,
    /// Padding for alignment.
    pub _pad2: [u8; 7],
    /// Proper motion (valid only when `has_proper_motion` is true).
    pub proper_motion: SiderustProperMotion,
}

// ═══════════════════════════════════════════════════════════════════════════
// Structs
// ═══════════════════════════════════════════════════════════════════════════

/// Geodetic position (WGS84) for C interop.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustGeodetict {
    /// Longitude in degrees (east positive).
    pub lon_deg: f64,
    /// Latitude in degrees (north positive).
    pub lat_deg: f64,
    /// Height above ellipsoid in metres.
    pub height_m: f64,
}

impl SiderustGeodetict {
    /// Convert to the Rust domain type.
    pub fn to_rust(&self) -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(
            Degrees::new(self.lon_deg),
            Degrees::new(self.lat_deg),
            Meters::new(self.height_m),
        )
    }

    /// Create from the Rust domain type.
    pub fn from_rust(g: &Geodetic<ECEF>) -> Self {
        Self {
            lon_deg: g.lon.value(),
            lat_deg: g.lat.value(),
            height_m: g.height.value(),
        }
    }
}

/// Keplerian orbital elements.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustOrbit {
    /// Semi-major axis in astronomical units.
    pub semi_major_axis_au: f64,
    /// Orbital eccentricity.
    pub eccentricity: f64,
    /// Orbital inclination in degrees.
    pub inclination_deg: f64,
    /// Longitude of the ascending node in degrees.
    pub lon_ascending_node_deg: f64,
    /// Argument of perihelion in degrees.
    pub arg_perihelion_deg: f64,
    /// Mean anomaly at epoch in degrees.
    pub mean_anomaly_deg: f64,
    /// Epoch as a Julian Date.
    pub epoch_jd: f64,
}

impl SiderustOrbit {
    /// Convert to the Rust domain type.
    pub fn to_rust(&self) -> siderust::astro::orbit::Orbit {
        siderust::astro::orbit::Orbit::new(
            AstronomicalUnits::new(self.semi_major_axis_au),
            self.eccentricity,
            Degrees::new(self.inclination_deg),
            Degrees::new(self.lon_ascending_node_deg),
            Degrees::new(self.arg_perihelion_deg),
            Degrees::new(self.mean_anomaly_deg),
            JulianDate::new(self.epoch_jd),
        )
    }

    /// Create from the Rust domain type.
    pub fn from_rust(o: &siderust::astro::orbit::Orbit) -> Self {
        Self {
            semi_major_axis_au: o.semi_major_axis.value(),
            eccentricity: o.eccentricity,
            inclination_deg: o.inclination.value(),
            lon_ascending_node_deg: o.longitude_of_ascending_node.value(),
            arg_perihelion_deg: o.argument_of_perihelion.value(),
            mean_anomaly_deg: o.mean_anomaly_at_epoch.value(),
            epoch_jd: o.epoch.value(),
        }
    }
}

/// Bodycentric reference center: which standard center the orbit is relative to.
/// Must match `OrbitReferenceCenter` in siderust: Barycentric=0, Heliocentric=1, Geocentric=2.
pub type SiderustOrbitRefCenter = u8;

/// Parameters for a body-centric coordinate system.
///
/// Pairs Keplerian orbital elements with the reference center of those elements.
/// Corresponds to `siderust::coordinates::centers::BodycentricParams`.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustBodycentricParams {
    /// Keplerian orbital elements of the body.
    pub orbit: SiderustOrbit,
    /// Reference center: 0=Barycentric, 1=Heliocentric, 2=Geocentric.
    pub orbit_center: SiderustOrbitRefCenter,
    /// Padding bytes for alignment; must be zeroed.
    pub _pad: [u8; 7],
}

impl SiderustBodycentricParams {
    /// Convert to the Rust domain type.
    pub fn to_rust(&self) -> RustBodycentricParams {
        let orbit = self.orbit.to_rust();
        let orbit_center = match self.orbit_center {
            0 => RustOrbitRefCenter::Barycentric,
            1 => RustOrbitRefCenter::Heliocentric,
            2 => RustOrbitRefCenter::Geocentric,
            _ => RustOrbitRefCenter::Heliocentric, // safe default
        };
        RustBodycentricParams::new(orbit, orbit_center)
    }

    /// Create from the Rust domain type.
    pub fn from_rust(p: &RustBodycentricParams) -> Self {
        let orbit_center = match p.orbit_center {
            RustOrbitRefCenter::Barycentric => 0u8,
            RustOrbitRefCenter::Heliocentric => 1u8,
            RustOrbitRefCenter::Geocentric => 2u8,
        };
        Self {
            orbit: SiderustOrbit::from_rust(&p.orbit),
            orbit_center,
            _pad: [0u8; 7],
        }
    }
}

/// Search options for altitude computations.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustSearchOpts {
    /// Time tolerance in days (default: ~1 µs = 1e-9 days).
    pub time_tolerance_days: f64,
    /// Scan step in days. Set to 0 or negative to use the body's default.
    pub scan_step_days: f64,
    /// Whether `scan_step_days` is valid (non-zero).
    pub has_scan_step: bool,
}

impl SiderustSearchOpts {
    /// Convert to the Rust domain type.
    pub fn to_rust(&self) -> siderust::SearchOpts {
        let mut opts = siderust::SearchOpts::default();
        if self.time_tolerance_days > 0.0 {
            opts.time_tolerance = Days::new(self.time_tolerance_days);
        }
        if self.has_scan_step && self.scan_step_days > 0.0 {
            opts.scan_step_days = Some(Days::new(self.scan_step_days));
        }
        opts
    }
}

impl Default for SiderustSearchOpts {
    fn default() -> Self {
        Self {
            time_tolerance_days: 1e-9,
            scan_step_days: 0.0,
            has_scan_step: false,
        }
    }
}

/// A threshold-crossing event.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustCrossingEvent {
    /// Time of the crossing (Modified Julian Date).
    pub mjd: f64,
    /// Direction of the crossing.
    pub direction: SiderustCrossingDirection,
}

impl SiderustCrossingEvent {
    /// Create from the Rust domain type.
    pub fn from_rust(e: &siderust::CrossingEvent) -> Self {
        Self {
            mjd: e.mjd.value(),
            direction: match e.direction {
                siderust::CrossingDirection::Rising => SiderustCrossingDirection::Rising,
                siderust::CrossingDirection::Setting => SiderustCrossingDirection::Setting,
            },
        }
    }
}

/// A culmination (local extremum in altitude) event.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustCulminationEvent {
    /// Time of the culmination (Modified Julian Date).
    pub mjd: f64,
    /// Altitude at the extremum in degrees.
    pub altitude_deg: f64,
    /// Kind of extremum.
    pub kind: SiderustCulminationKind,
}

impl SiderustCulminationEvent {
    /// Create from the Rust domain type.
    pub fn from_rust(e: &siderust::CulminationEvent) -> Self {
        Self {
            mjd: e.mjd.value(),
            altitude_deg: e.altitude.value(),
            kind: match e.kind {
                siderust::CulminationKind::Max => SiderustCulminationKind::Max,
                siderust::CulminationKind::Min => SiderustCulminationKind::Min,
            },
        }
    }
}

/// Altitude computation query parameters.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustAltitudeQuery {
    /// Observer location.
    pub observer: SiderustGeodetict,
    /// Start of the search window (Modified Julian Date).
    pub start_mjd: f64,
    /// End of the search window (Modified Julian Date).
    pub end_mjd: f64,
    /// Minimum altitude threshold in degrees.
    pub min_altitude_deg: f64,
    /// Maximum altitude threshold in degrees.
    pub max_altitude_deg: f64,
}

impl SiderustAltitudeQuery {
    /// Convert to the Rust domain type.
    pub fn to_rust(&self) -> siderust::AltitudeQuery {
        siderust::AltitudeQuery {
            observer: self.observer.to_rust(),
            window: Interval::new(
                ModifiedJulianDate::new(self.start_mjd),
                ModifiedJulianDate::new(self.end_mjd),
            ),
            min_altitude: Degrees::new(self.min_altitude_deg),
            max_altitude: Degrees::new(self.max_altitude_deg),
        }
    }
}

/// Proper motion of a star (equatorial).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustProperMotion {
    /// RA proper motion in degrees per Julian year.
    pub pm_ra_deg_yr: f64,
    /// Dec proper motion in degrees per Julian year.
    pub pm_dec_deg_yr: f64,
    /// Convention for `pm_ra`.
    pub ra_convention: SiderustRaConvention,
}

/// Planet data (value type, copyable).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustPlanet {
    /// Planet mass in kilograms.
    pub mass_kg: f64,
    /// Mean equatorial radius in kilometres.
    pub radius_km: f64,
    /// Keplerian orbital elements.
    pub orbit: SiderustOrbit,
}

impl SiderustPlanet {
    /// Create from the Rust domain type.
    pub fn from_rust(p: &siderust::bodies::Planet) -> Self {
        Self {
            mass_kg: p.mass.value(),
            radius_km: p.radius.value(),
            orbit: SiderustOrbit::from_rust(&p.orbit),
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Coordinate types
// ═══════════════════════════════════════════════════════════════════════════

/// Spherical direction (lon/lat or RA/Dec) in degrees.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustSphericalDir {
    /// Polar angle / latitude in degrees.
    pub polar_deg: f64,
    /// Azimuth / longitude / right ascension in degrees.
    pub azimuth_deg: f64,
    /// Reference frame.
    pub frame: SiderustFrame,
}

impl SiderustSphericalDir {
    /// Return a zeroed direction (used as placeholder for unused fields).
    pub(crate) fn zeroed() -> Self {
        Self {
            polar_deg: 0.0,
            azimuth_deg: 0.0,
            frame: SiderustFrame::ICRS,
        }
    }
}

/// Cartesian unit-direction vector.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustCartesianDir {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
    /// Reference frame.
    pub frame: SiderustFrame,
}

/// Spherical position (direction + distance).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustSphericalPos {
    /// Longitude in degrees.
    pub lon_deg: f64,
    /// Latitude in degrees.
    pub lat_deg: f64,
    /// Distance in the unit appropriate for the context (AU, km, ly, etc.).
    pub distance: f64,
    /// Reference frame.
    pub frame: SiderustFrame,
    /// Reference centre.
    pub center: SiderustCenter,
}

/// Cartesian position (x, y, z + metadata).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustCartesianPos {
    /// X coordinate in the specified length unit.
    pub x: f64,
    /// Y coordinate in the specified length unit.
    pub y: f64,
    /// Z coordinate in the specified length unit.
    pub z: f64,
    /// Reference frame.
    pub frame: SiderustFrame,
    /// Reference centre.
    pub center: SiderustCenter,
    /// Length unit for x, y, z coordinates.
    pub length_unit: SiderustLengthUnit,
}

/// Cartesian velocity (vx, vy, vz + frame metadata).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustCartesianVel {
    /// X velocity component (AU/day).
    pub vx: f64,
    /// Y velocity component (AU/day).
    pub vy: f64,
    /// Z velocity component (AU/day).
    pub vz: f64,
    /// Reference frame.
    pub frame: SiderustFrame,
}

// ═══════════════════════════════════════════════════════════════════════════
// Azimuth types
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth extremum kind (maximum or minimum bearing).
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustAzimuthExtremumKind {
    /// Maximum azimuth (easternmost bearing).
    Max = 0,
    /// Minimum azimuth (westernmost bearing).
    Min = 1,
}

/// An azimuth bearing-crossing event.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustAzimuthCrossingEvent {
    /// Time of the event (Modified Julian Date).
    pub mjd: f64,
    /// Crossing direction.
    pub direction: SiderustCrossingDirection,
    /// Padding bytes for alignment; must be zeroed.
    pub _pad: [u8; 4],
}

/// An azimuth extremum (minimum or maximum bearing).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustAzimuthExtremum {
    /// Time of the extremum (Modified Julian Date).
    pub mjd: f64,
    /// Azimuth at the extremum, in degrees.
    pub azimuth_deg: f64,
    /// Kind of extremum.
    pub kind: SiderustAzimuthExtremumKind,
    /// Padding bytes for alignment; must be zeroed.
    pub _pad: [u8; 4],
}

// ═══════════════════════════════════════════════════════════════════════════
// Lunar phase types
// ═══════════════════════════════════════════════════════════════════════════

/// Principal lunar phase kind.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustPhaseKind {
    /// New Moon.
    NewMoon = 0,
    /// First Quarter.
    FirstQuarter = 1,
    /// Full Moon.
    FullMoon = 2,
    /// Last Quarter.
    LastQuarter = 3,
}

/// Named phase label (includes intermediate phases).
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustMoonPhaseLabel {
    /// New Moon.
    NewMoon = 0,
    /// Waxing Crescent.
    WaxingCrescent = 1,
    /// First Quarter.
    FirstQuarter = 2,
    /// Waxing Gibbous.
    WaxingGibbous = 3,
    /// Full Moon.
    FullMoon = 4,
    /// Waning Gibbous.
    WaningGibbous = 5,
    /// Last Quarter.
    LastQuarter = 6,
    /// Waning Crescent.
    WaningCrescent = 7,
}

/// Moon phase geometry at a given instant.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustMoonPhaseGeometry {
    /// Phase angle in radians.
    pub phase_angle_rad: f64,
    /// Fraction of the Moon's disk that is illuminated (0–1).
    pub illuminated_fraction: f64,
    /// Elongation in radians.
    pub elongation_rad: f64,
    /// Non-zero if waxing, zero if waning.
    pub waxing: u8,
    /// Padding bytes for alignment; must be zeroed.
    pub _pad: [u8; 7],
}

/// A principal lunar phase event.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustPhaseEvent {
    /// Time of the event (Modified Julian Date).
    pub mjd: f64,
    /// Phase kind.
    pub kind: SiderustPhaseKind,
    /// Padding bytes for alignment; must be zeroed.
    pub _pad: [u8; 4],
}

// ═══════════════════════════════════════════════════════════════════════════
// FfiFrom implementations
// ═══════════════════════════════════════════════════════════════════════════

impl FfiFrom<Period<MJD>> for TempochPeriodMjd {
    fn ffi_from(p: &Period<MJD>) -> Self {
        TempochPeriodMjd {
            start_mjd: p.start.value(),
            end_mjd: p.end.value(),
        }
    }
}

impl FfiFrom<siderust::CrossingEvent> for SiderustCrossingEvent {
    fn ffi_from(e: &siderust::CrossingEvent) -> Self {
        SiderustCrossingEvent::from_rust(e)
    }
}

impl FfiFrom<siderust::CulminationEvent> for SiderustCulminationEvent {
    fn ffi_from(e: &siderust::CulminationEvent) -> Self {
        SiderustCulminationEvent::from_rust(e)
    }
}

impl FfiFrom<AzimuthCrossingEvent> for SiderustAzimuthCrossingEvent {
    fn ffi_from(e: &AzimuthCrossingEvent) -> Self {
        SiderustAzimuthCrossingEvent {
            mjd: e.mjd.value(),
            direction: match e.direction {
                AzimuthCrossingDirection::Rising => SiderustCrossingDirection::Rising,
                AzimuthCrossingDirection::Setting => SiderustCrossingDirection::Setting,
            },
            _pad: [0; 4],
        }
    }
}

impl FfiFrom<AzimuthExtremum> for SiderustAzimuthExtremum {
    fn ffi_from(e: &AzimuthExtremum) -> Self {
        SiderustAzimuthExtremum {
            mjd: e.mjd.value(),
            azimuth_deg: e.azimuth.value(),
            kind: match e.kind {
                AzimuthExtremumKind::Max => SiderustAzimuthExtremumKind::Max,
                AzimuthExtremumKind::Min => SiderustAzimuthExtremumKind::Min,
            },
            _pad: [0; 4],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── SiderustGeodetict ────────────────────────────────────────────────

    #[test]
    fn geodetic_roundtrip() {
        let orig = SiderustGeodetict {
            lon_deg: 10.0,
            lat_deg: -20.0,
            height_m: 500.0,
        };
        let rust = orig.to_rust();
        let back = SiderustGeodetict::from_rust(&rust);
        assert!((back.lon_deg - 10.0).abs() < 1e-10);
        assert!((back.lat_deg - (-20.0)).abs() < 1e-10);
        assert!((back.height_m - 500.0).abs() < 1e-10);
    }

    // ── SiderustOrbit ────────────────────────────────────────────────────

    #[test]
    fn orbit_roundtrip() {
        let orig = SiderustOrbit {
            semi_major_axis_au: 1.0,
            eccentricity: 0.017,
            inclination_deg: 7.25,
            lon_ascending_node_deg: 48.3,
            arg_perihelion_deg: 102.9,
            mean_anomaly_deg: 100.0,
            epoch_jd: 2_451_545.0,
        };
        let rust = orig.to_rust();
        let back = SiderustOrbit::from_rust(&rust);
        assert!((back.semi_major_axis_au - 1.0).abs() < 1e-10);
        assert!((back.eccentricity - 0.017).abs() < 1e-10);
        assert!((back.inclination_deg - 7.25).abs() < 1e-10);
        assert!((back.epoch_jd - 2_451_545.0).abs() < 1e-10);
    }

    // ── SiderustSearchOpts ───────────────────────────────────────────────

    #[test]
    fn search_opts_default() {
        let d = SiderustSearchOpts::default();
        assert!((d.time_tolerance_days - 1e-9).abs() < 1e-15);
        assert!(!d.has_scan_step);
    }

    #[test]
    fn search_opts_to_rust_without_scan_step() {
        let opts = SiderustSearchOpts {
            time_tolerance_days: 1e-6,
            scan_step_days: 0.0,
            has_scan_step: false,
        };
        let rust = opts.to_rust();
        assert!((rust.time_tolerance.value() - 1e-6).abs() < 1e-12);
        assert!(rust.scan_step_days.is_none());
    }

    #[test]
    fn search_opts_to_rust_with_scan_step() {
        let opts = SiderustSearchOpts {
            time_tolerance_days: 1e-9,
            scan_step_days: 0.1,
            has_scan_step: true,
        };
        let rust = opts.to_rust();
        assert!(rust.scan_step_days.is_some());
        assert!((rust.scan_step_days.unwrap().value() - 0.1).abs() < 1e-12);
    }

    #[test]
    fn search_opts_zero_tolerance_uses_default() {
        let opts = SiderustSearchOpts {
            time_tolerance_days: 0.0, // ≤ 0 → keep default
            scan_step_days: 0.0,
            has_scan_step: false,
        };
        let rust = opts.to_rust();
        // With 0 tolerance the default is kept
        let default_rust = SiderustSearchOpts::default().to_rust();
        assert_eq!(
            rust.time_tolerance.value(),
            default_rust.time_tolerance.value()
        );
    }

    // ── SiderustAltitudeQuery ────────────────────────────────────────────

    #[test]
    fn altitude_query_to_rust() {
        let q = SiderustAltitudeQuery {
            observer: SiderustGeodetict {
                lon_deg: 0.0,
                lat_deg: 51.5,
                height_m: 10.0,
            },
            start_mjd: 60000.0,
            end_mjd: 60001.0,
            min_altitude_deg: 10.0,
            max_altitude_deg: 90.0,
        };
        let rust = q.to_rust();
        assert!((rust.min_altitude.value() - 10.0).abs() < 1e-10);
        assert!((rust.max_altitude.value() - 90.0).abs() < 1e-10);
    }

    // ── SiderustPlanet ───────────────────────────────────────────────────

    #[test]
    fn planet_from_rust_earth() {
        let p = SiderustPlanet::from_rust(&siderust::bodies::EARTH);
        // Earth mass ~5.97e24 kg
        assert!(
            p.mass_kg > 5e24 && p.mass_kg < 7e24,
            "Earth mass out of range: {}",
            p.mass_kg
        );
        // Earth radius ~6371 km
        assert!(
            p.radius_km > 6000.0 && p.radius_km < 7000.0,
            "Earth radius out of range: {}",
            p.radius_km
        );
    }

    // ── Enum PartialEq/Clone ─────────────────────────────────────────────

    #[test]
    fn enum_clone_and_eq() {
        assert_eq!(SiderustFrame::ICRS, SiderustFrame::ICRS.clone());
        assert_ne!(SiderustFrame::ICRS, SiderustFrame::GCRS);
        assert_eq!(
            SiderustCenter::Geocentric,
            SiderustCenter::Geocentric.clone()
        );
        assert_eq!(
            SiderustCrossingDirection::Rising,
            SiderustCrossingDirection::Rising
        );
        assert_eq!(SiderustCulminationKind::Max, SiderustCulminationKind::Max);
        assert_eq!(
            SiderustAzimuthExtremumKind::Max,
            SiderustAzimuthExtremumKind::Max
        );
        assert_eq!(SiderustPhaseKind::FullMoon, SiderustPhaseKind::FullMoon);
        assert_eq!(
            SiderustMoonPhaseLabel::WaxingCrescent,
            SiderustMoonPhaseLabel::WaxingCrescent
        );
    }

    #[test]
    fn enum_debug() {
        let s = format!("{:?}", SiderustFrame::Horizontal);
        assert!(s.contains("Horizontal"));
        let s2 = format!("{:?}", SiderustAsteroidClass::MainBelt);
        assert!(s2.contains("MainBelt"));
        let s3 = format!("{:?}", SiderustOrbitFrame::Heliocentric);
        assert!(s3.contains("Heliocentric"));
        let s4 = format!("{:?}", SiderustRaConvention::MuAlphaStar);
        assert!(s4.contains("MuAlphaStar"));
    }

    // ── FfiFrom for Period<MJD> ──────────────────────────────────────────

    #[test]
    fn ffi_from_period_mjd() {
        use tempoch::{Interval, ModifiedJulianDate, MJD};
        let p: tempoch::Period<MJD> = Interval::new(
            ModifiedJulianDate::new(60000.0),
            ModifiedJulianDate::new(60001.0),
        );
        let ffi = TempochPeriodMjd::ffi_from(&p);
        assert!((ffi.start_mjd - 60000.0).abs() < 1e-10);
        assert!((ffi.end_mjd - 60001.0).abs() < 1e-10);
    }

    // ── SiderustCrossingEvent / SiderustCulminationEvent ─────────────────

    #[test]
    fn crossing_event_from_rust_rising() {
        let e = siderust::CrossingEvent {
            mjd: tempoch::ModifiedJulianDate::new(60000.5),
            direction: siderust::CrossingDirection::Rising,
        };
        let ffi = SiderustCrossingEvent::from_rust(&e);
        assert!((ffi.mjd - 60000.5).abs() < 1e-10);
        assert_eq!(ffi.direction, SiderustCrossingDirection::Rising);
    }

    #[test]
    fn crossing_event_from_rust_setting() {
        let e = siderust::CrossingEvent {
            mjd: tempoch::ModifiedJulianDate::new(60000.5),
            direction: siderust::CrossingDirection::Setting,
        };
        let ffi = SiderustCrossingEvent::ffi_from(&e);
        assert_eq!(ffi.direction, SiderustCrossingDirection::Setting);
    }

    #[test]
    fn culmination_event_from_rust_max() {
        let e = siderust::CulminationEvent {
            mjd: tempoch::ModifiedJulianDate::new(60000.0),
            altitude: qtty::Degrees::new(45.0),
            kind: siderust::CulminationKind::Max,
        };
        let ffi = SiderustCulminationEvent::from_rust(&e);
        assert!((ffi.altitude_deg - 45.0).abs() < 1e-10);
        assert_eq!(ffi.kind, SiderustCulminationKind::Max);
    }

    #[test]
    fn culmination_event_from_rust_min() {
        let e = siderust::CulminationEvent {
            mjd: tempoch::ModifiedJulianDate::new(60000.0),
            altitude: qtty::Degrees::new(-10.0),
            kind: siderust::CulminationKind::Min,
        };
        let ffi = SiderustCulminationEvent::ffi_from(&e);
        assert_eq!(ffi.kind, SiderustCulminationKind::Min);
    }

    // ── AzimuthExtremum FfiFrom ──────────────────────────────────────────

    #[test]
    fn azimuth_extremum_ffi_from_max() {
        use siderust::calculus::azimuth::{AzimuthExtremum, AzimuthExtremumKind};
        let e = AzimuthExtremum {
            mjd: tempoch::ModifiedJulianDate::new(60000.0),
            azimuth: qtty::Degrees::new(180.0),
            kind: AzimuthExtremumKind::Max,
        };
        let ffi = SiderustAzimuthExtremum::ffi_from(&e);
        assert!((ffi.azimuth_deg - 180.0).abs() < 1e-10);
        assert_eq!(ffi.kind, SiderustAzimuthExtremumKind::Max);
    }

    #[test]
    fn azimuth_extremum_ffi_from_min() {
        use siderust::calculus::azimuth::{AzimuthExtremum, AzimuthExtremumKind};
        let e = AzimuthExtremum {
            mjd: tempoch::ModifiedJulianDate::new(60000.0),
            azimuth: qtty::Degrees::new(0.0),
            kind: AzimuthExtremumKind::Min,
        };
        let ffi = SiderustAzimuthExtremum::ffi_from(&e);
        assert_eq!(ffi.kind, SiderustAzimuthExtremumKind::Min);
    }

    // ── AzimuthCrossingEvent FfiFrom ─────────────────────────────────────

    #[test]
    fn azimuth_crossing_event_ffi_from_rising() {
        use siderust::calculus::azimuth::{AzimuthCrossingDirection, AzimuthCrossingEvent};
        let e = AzimuthCrossingEvent {
            mjd: tempoch::ModifiedJulianDate::new(60000.0),
            direction: AzimuthCrossingDirection::Rising,
        };
        let ffi = SiderustAzimuthCrossingEvent::ffi_from(&e);
        assert_eq!(ffi.direction, SiderustCrossingDirection::Rising);
    }

    #[test]
    fn azimuth_crossing_event_ffi_from_setting() {
        use siderust::calculus::azimuth::{AzimuthCrossingDirection, AzimuthCrossingEvent};
        let e = AzimuthCrossingEvent {
            mjd: tempoch::ModifiedJulianDate::new(60000.0),
            direction: AzimuthCrossingDirection::Setting,
        };
        let ffi = SiderustAzimuthCrossingEvent::ffi_from(&e);
        assert_eq!(ffi.direction, SiderustCrossingDirection::Setting);
    }

    // ── ABI layout assertions ────────────────────────────────────────────
    // These tests lock down the size and alignment of every repr(C) type
    // so accidental field additions/reorderings break CI before they
    // reach consumers.

    macro_rules! assert_layout {
        ($ty:ty, size = $size:expr, align = $align:expr) => {
            assert_eq!(
                std::mem::size_of::<$ty>(),
                $size,
                concat!("size_of::<", stringify!($ty), ">() mismatch")
            );
            assert_eq!(
                std::mem::align_of::<$ty>(),
                $align,
                concat!("align_of::<", stringify!($ty), ">() mismatch")
            );
        };
    }

    #[test]
    fn layout_status_enum() {
        assert_layout!(crate::error::SiderustStatus, size = 4, align = 4);
    }

    #[test]
    fn layout_frame_enum() {
        assert_layout!(SiderustFrame, size = 4, align = 4);
    }

    #[test]
    fn layout_center_enum() {
        assert_layout!(SiderustCenter, size = 4, align = 4);
    }

    #[test]
    fn layout_body_enum() {
        assert_layout!(SiderustBody, size = 4, align = 4);
    }

    #[test]
    fn layout_geodetic() {
        assert_layout!(SiderustGeodetict, size = 24, align = 8);
    }

    #[test]
    fn layout_orbit() {
        assert_layout!(SiderustOrbit, size = 56, align = 8);
    }

    #[test]
    fn layout_search_opts() {
        // 2 × f64 + bool + padding = 24
        assert_eq!(std::mem::size_of::<SiderustSearchOpts>(), 24);
    }

    #[test]
    fn layout_crossing_event() {
        // f64 + i32 + padding = 16
        assert_eq!(std::mem::size_of::<SiderustCrossingEvent>(), 16);
        assert_eq!(std::mem::align_of::<SiderustCrossingEvent>(), 8);
    }

    #[test]
    fn layout_culmination_event() {
        // 2 × f64 + i32 + padding = 24
        assert_eq!(std::mem::size_of::<SiderustCulminationEvent>(), 24);
        assert_eq!(std::mem::align_of::<SiderustCulminationEvent>(), 8);
    }

    #[test]
    fn layout_spherical_dir() {
        // 2 × f64 + i32 + padding = 24
        assert_eq!(std::mem::size_of::<SiderustSphericalDir>(), 24);
        assert_eq!(std::mem::align_of::<SiderustSphericalDir>(), 8);
    }

    #[test]
    fn layout_cartesian_dir() {
        // 3 × f64 + i32 + padding = 32
        assert_eq!(std::mem::size_of::<SiderustCartesianDir>(), 32);
        assert_eq!(std::mem::align_of::<SiderustCartesianDir>(), 8);
    }

    #[test]
    fn layout_cartesian_pos() {
        // 3 × f64 + i32 (frame) + i32 (center) + i32 (length_unit) + padding = 40
        assert_eq!(std::mem::size_of::<SiderustCartesianPos>(), 40);
        assert_eq!(std::mem::align_of::<SiderustCartesianPos>(), 8);
    }

    #[test]
    fn layout_cartesian_vel() {
        // 3 × f64 + i32 + padding = 32
        assert_eq!(std::mem::size_of::<SiderustCartesianVel>(), 32);
        assert_eq!(std::mem::align_of::<SiderustCartesianVel>(), 8);
    }

    #[test]
    fn layout_azimuth_crossing_event() {
        // f64 + i32 + 4 pad = 16
        assert_eq!(std::mem::size_of::<SiderustAzimuthCrossingEvent>(), 16);
        assert_eq!(std::mem::align_of::<SiderustAzimuthCrossingEvent>(), 8);
    }

    #[test]
    fn layout_azimuth_extremum() {
        // 2 × f64 + i32 + 4 pad = 24
        assert_eq!(std::mem::size_of::<SiderustAzimuthExtremum>(), 24);
        assert_eq!(std::mem::align_of::<SiderustAzimuthExtremum>(), 8);
    }

    #[test]
    fn layout_moon_phase_geometry() {
        // 3 × f64 + u8 + 7 pad = 32
        assert_eq!(std::mem::size_of::<SiderustMoonPhaseGeometry>(), 32);
        assert_eq!(std::mem::align_of::<SiderustMoonPhaseGeometry>(), 8);
    }

    #[test]
    fn layout_phase_event() {
        // f64 + i32 + 4 pad = 16
        assert_eq!(std::mem::size_of::<SiderustPhaseEvent>(), 16);
        assert_eq!(std::mem::align_of::<SiderustPhaseEvent>(), 8);
    }

    #[test]
    fn layout_planet() {
        // 2 × f64 + SiderustOrbit(56) = 72
        assert_eq!(std::mem::size_of::<SiderustPlanet>(), 72);
        assert_eq!(std::mem::align_of::<SiderustPlanet>(), 8);
    }

    // ───────────────────────────────────────────────────────────
    // Frame/Center parsing tests
    // ───────────────────────────────────────────────────────────

    #[test]
    fn frame_from_str_canonical() {
        assert_eq!(SiderustFrame::from_str("icrs"), Some(SiderustFrame::ICRS));
        assert_eq!(
            SiderustFrame::from_str("horizontal"),
            Some(SiderustFrame::Horizontal)
        );
        assert_eq!(
            SiderustFrame::from_str("galactic"),
            Some(SiderustFrame::Galactic)
        );
        assert_eq!(SiderustFrame::from_str("gcrs"), Some(SiderustFrame::GCRS));
    }

    #[test]
    fn frame_from_str_case_insensitive() {
        assert_eq!(SiderustFrame::from_str("ICRS"), Some(SiderustFrame::ICRS));
        assert_eq!(SiderustFrame::from_str("Icrs"), Some(SiderustFrame::ICRS));
        assert_eq!(
            SiderustFrame::from_str("HORIZONTAL"),
            Some(SiderustFrame::Horizontal)
        );
    }

    #[test]
    fn frame_from_str_aliases() {
        assert_eq!(
            SiderustFrame::from_str("altaz"),
            Some(SiderustFrame::Horizontal)
        );
        assert_eq!(
            SiderustFrame::from_str("gal"),
            Some(SiderustFrame::Galactic)
        );
        assert_eq!(
            SiderustFrame::from_str("mean_of_date"),
            Some(SiderustFrame::EquatorialMeanOfDate)
        );
    }

    #[test]
    fn frame_from_str_unknown() {
        assert_eq!(SiderustFrame::from_str("unknown"), None);
        assert_eq!(SiderustFrame::from_str(""), None);
    }

    #[test]
    fn frame_roundtrip() {
        for frame in [
            SiderustFrame::ICRS,
            SiderustFrame::Horizontal,
            SiderustFrame::Galactic,
            SiderustFrame::GCRS,
            SiderustFrame::ECEF,
        ] {
            let name = frame.as_str();
            assert_eq!(SiderustFrame::from_str(name), Some(frame));
        }
    }

    #[test]
    fn center_from_str_canonical() {
        assert_eq!(
            SiderustCenter::from_str("barycentric"),
            Some(SiderustCenter::Barycentric)
        );
        assert_eq!(
            SiderustCenter::from_str("heliocentric"),
            Some(SiderustCenter::Heliocentric)
        );
        assert_eq!(
            SiderustCenter::from_str("geocentric"),
            Some(SiderustCenter::Geocentric)
        );
        assert_eq!(
            SiderustCenter::from_str("topocentric"),
            Some(SiderustCenter::Topocentric)
        );
    }

    #[test]
    fn center_from_str_aliases() {
        assert_eq!(
            SiderustCenter::from_str("ssb"),
            Some(SiderustCenter::Barycentric)
        );
        assert_eq!(
            SiderustCenter::from_str("sun"),
            Some(SiderustCenter::Heliocentric)
        );
        assert_eq!(
            SiderustCenter::from_str("earth"),
            Some(SiderustCenter::Geocentric)
        );
        assert_eq!(
            SiderustCenter::from_str("observer"),
            Some(SiderustCenter::Topocentric)
        );
    }

    #[test]
    fn center_roundtrip() {
        for center in [
            SiderustCenter::Barycentric,
            SiderustCenter::Heliocentric,
            SiderustCenter::Geocentric,
            SiderustCenter::Topocentric,
            SiderustCenter::Bodycentric,
        ] {
            let name = center.as_str();
            assert_eq!(SiderustCenter::from_str(name), Some(center));
        }
    }
}
