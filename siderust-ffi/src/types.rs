// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! C-compatible struct definitions for siderust-ffi.
//!
//! These `#[repr(C)]` types are the data bridge between C/C++ callers and
//! the rich generic Rust types in siderust.

use crate::ffi_utils::FfiFrom;
use qtty::*;
use siderust::calculus::azimuth::{
    AzimuthCrossingDirection, AzimuthCrossingEvent, AzimuthExtremum, AzimuthExtremumKind,
};
use siderust::coordinates::centers::Geodetic;
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

/// Solar-system body identifier for generic altitude/azimuth dispatch.
///
/// Each variant maps to a concrete unit type in `siderust::bodies::solar_system`.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustBody {
    Sun = 0,
    Moon = 1,
    Mercury = 2,
    Venus = 3,
    Mars = 4,
    Jupiter = 5,
    Saturn = 6,
    Uranus = 7,
    Neptune = 8,
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
    /// Polar angle / latitude in degrees.
    pub polar_deg: f64,
    /// Azimuth / longitude / right ascension in degrees.
    pub azimuth_deg: f64,
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
}
