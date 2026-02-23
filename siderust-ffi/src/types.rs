// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! C-compatible struct definitions for siderust-ffi.
//!
//! These `#[repr(C)]` types are the data bridge between C/C++ callers and
//! the rich generic Rust types in siderust.

use qtty::*;
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use tempoch::{Interval, JulianDate, ModifiedJulianDate};

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
    ICRS = 1,
    EclipticMeanJ2000 = 2,
    EquatorialMeanJ2000 = 3,
    EquatorialMeanOfDate = 4,
    EquatorialTrueOfDate = 5,
    Horizontal = 6,
    ECEF = 7,
    Galactic = 8,
    GCRS = 9,
    EclipticOfDate = 10,
    EclipticTrueOfDate = 11,
    CIRS = 12,
    TIRS = 13,
    ITRF = 14,
    ICRF = 15,
}

/// Reference center identifier for C interop.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustCenter {
    Barycentric = 1,
    Heliocentric = 2,
    Geocentric = 3,
    Topocentric = 4,
    Bodycentric = 5,
}

/// Crossing event direction.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustCrossingDirection {
    Rising = 0,
    Setting = 1,
}

/// Culmination event kind.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustCulminationKind {
    Max = 0,
    Min = 1,
}

/// Asteroid taxonomic class.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustAsteroidClass {
    MainBelt = 0,
    NearEarth = 1,
    Trojan = 2,
    Centaur = 3,
    TransNeptunian = 4,
    DwarfPlanet = 5,
}

/// Orbit reference frame for comets.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustOrbitFrame {
    Heliocentric = 0,
    Barycentric = 1,
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
    pub fn to_rust(&self) -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(
            Degrees::new(self.lon_deg),
            Degrees::new(self.lat_deg),
            Meters::new(self.height_m),
        )
    }

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
    pub semi_major_axis_au: f64,
    pub eccentricity: f64,
    pub inclination_deg: f64,
    pub lon_ascending_node_deg: f64,
    pub arg_perihelion_deg: f64,
    pub mean_anomaly_deg: f64,
    pub epoch_jd: f64,
}

impl SiderustOrbit {
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
    pub observer: SiderustGeodetict,
    pub start_mjd: f64,
    pub end_mjd: f64,
    pub min_altitude_deg: f64,
    pub max_altitude_deg: f64,
}

impl SiderustAltitudeQuery {
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
    pub mass_kg: f64,
    pub radius_km: f64,
    pub orbit: SiderustOrbit,
}

impl SiderustPlanet {
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
    /// Longitude / Right Ascension in degrees.
    pub lon_deg: f64,
    /// Latitude / Declination in degrees.
    pub lat_deg: f64,
    /// Reference frame.
    pub frame: SiderustFrame,
}

/// Cartesian unit-direction vector.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustCartesianDir {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub frame: SiderustFrame,
}

/// Spherical position (direction + distance).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustSphericalPos {
    pub lon_deg: f64,
    pub lat_deg: f64,
    /// Distance in the unit appropriate for the context (AU, km, ly, etc.).
    pub distance: f64,
    pub frame: SiderustFrame,
    pub center: SiderustCenter,
}

/// Cartesian position (x, y, z + metadata).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustCartesianPos {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub frame: SiderustFrame,
    pub center: SiderustCenter,
}

// ═══════════════════════════════════════════════════════════════════════════
// Azimuth types
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth extremum kind (maximum or minimum bearing).
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustAzimuthExtremumKind {
    Max = 0,
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
    pub _pad: [u8; 4],
}

// ═══════════════════════════════════════════════════════════════════════════
// Lunar phase types
// ═══════════════════════════════════════════════════════════════════════════

/// Principal lunar phase kind.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustPhaseKind {
    NewMoon = 0,
    FirstQuarter = 1,
    FullMoon = 2,
    LastQuarter = 3,
}

/// Named phase label (includes intermediate phases).
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustMoonPhaseLabel {
    NewMoon = 0,
    WaxingCrescent = 1,
    FirstQuarter = 2,
    WaxingGibbous = 3,
    FullMoon = 4,
    WaningGibbous = 5,
    LastQuarter = 6,
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
    pub _pad: [u8; 4],
}
