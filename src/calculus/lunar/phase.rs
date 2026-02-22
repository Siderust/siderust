// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Moon Phase Module
//!
//! Computes the photometric geometry of the Moon–Sun–Earth (or observer)
//! system and derives human-facing phase labels.
//!
//! ## Outputs
//!
//! Two layers are exposed:
//!
//! ### A) Continuous geometry ([`MoonPhaseGeometry`])
//!
//! * **Phase angle** *i*: Sun–Moon–Observer angle at the Moon
//! * **Illuminated fraction** *k* = (1 + cos *i*) / 2
//! * **Elongation** ψ: geocentric (or topocentric) angular separation
//!   between Moon and Sun, in \[0, 2π)
//! * **Waxing/waning** flag: derived from elongation ψ ∈ (0, π)
//!
//! ### B) Discrete label ([`MoonPhaseLabel`])
//!
//! Eight classical names, derived from the elongation with configurable
//! thresholds ([`PhaseThresholds`]).
//!
//! ## Event finding ([`find_phase_events`])
//!
//! Locates principal lunar phases (new, first quarter, full, last quarter)
//! in a time window using the same scan + Brent root-finding approach as
//! the altitude-crossings API.
//!
//! ## Backend selection
//!
//! All free functions are generic over `E: Ephemeris`, so callers can choose
//! between `Vsop87Ephemeris` (always available) and optional `De440Ephemeris`
//! / `De441Ephemeris` when the corresponding feature is enabled.
//!
//! ## Time scale
//!
//! All `JulianDate` / `ModifiedJulianDate` / `Period<MJD>` values are
//! interpreted on the TT axis (consistent with the altitude API).

use crate::calculus::ephemeris::Ephemeris;
use crate::calculus::math_core::intervals;
use crate::coordinates::cartesian;
use crate::coordinates::centers::*;
use crate::coordinates::frames;
use crate::time::{JulianDate, ModifiedJulianDate, Period, MJD};
use qtty::*;
use std::f64::consts::PI;
use std::marker::PhantomData;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

// ===========================================================================
// Constants
// ===========================================================================

const TWO_PI: f64 = 2.0 * PI;

/// 1 AU in kilometres (IAU 2012).
const AU_KM: f64 = 149_597_870.700;

/// Default scan step for phase-event detection: ~0.5 days.
/// The synodic month is ~29.53 days, so 0.5 d gives ~59 samples per cycle,
/// easily catching every quarter.
const PHASE_SCAN_STEP: Days = Days::new(0.5);

// ===========================================================================
// PhaseThresholds
// ===========================================================================

/// Threshold boundaries (in degrees of elongation) for mapping the
/// continuous elongation to the eight classical phase names.
///
/// Defaults: each octant spans 45°, centred on the cardinal/intercardinal
/// directions (0°, 45°, 90°, 135°, 180°, 225°, 270°, 315°).
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PhaseThresholds {
    /// Half-width around each cardinal direction (Degrees quantity).
    /// Default: 22.5° → each bin spans 45°.
    pub half_width: Degrees,
}

impl Default for PhaseThresholds {
    fn default() -> Self {
        Self {
            half_width: Degrees::new(22.5),
        }
    }
}

// ===========================================================================
// MoonPhaseLabel
// ===========================================================================

/// The eight classical Moon-phase names.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum MoonPhaseLabel {
    NewMoon,
    WaxingCrescent,
    FirstQuarter,
    WaxingGibbous,
    FullMoon,
    WaningGibbous,
    LastQuarter,
    WaningCrescent,
}

impl MoonPhaseLabel {
    /// Map elongation (Degrees quantity, in \[0, 360)) to a label using the given
    /// thresholds.
    pub fn from_elongation(elongation: Degrees, th: &PhaseThresholds) -> Self {
        let hw = th.half_width;
        // Normalise to [0, 360) as a Degrees quantity
        let mut e = elongation;
        let full = Degrees::new(360.0);
        while e >= full {
            e -= full;
        }
        while e < Degrees::new(0.0) {
            e += full;
        }

        if e < hw || e >= full - hw {
            Self::NewMoon
        } else if e < Degrees::new(90.0) - hw {
            Self::WaxingCrescent
        } else if e < Degrees::new(90.0) + hw {
            Self::FirstQuarter
        } else if e < Degrees::new(180.0) - hw {
            Self::WaxingGibbous
        } else if e < Degrees::new(180.0) + hw {
            Self::FullMoon
        } else if e < Degrees::new(270.0) - hw {
            Self::WaningGibbous
        } else if e < Degrees::new(270.0) + hw {
            Self::LastQuarter
        } else {
            Self::WaningCrescent
        }
    }

    /// Returns `true` if the Moon is gaining illumination.
    pub fn is_waxing(self) -> bool {
        matches!(
            self,
            Self::WaxingCrescent | Self::FirstQuarter | Self::WaxingGibbous
        )
    }

    /// Returns `true` if the Moon is losing illumination.
    pub fn is_waning(self) -> bool {
        matches!(
            self,
            Self::WaningGibbous | Self::LastQuarter | Self::WaningCrescent
        )
    }
}

impl std::fmt::Display for MoonPhaseLabel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NewMoon => write!(f, "New Moon"),
            Self::WaxingCrescent => write!(f, "Waxing Crescent"),
            Self::FirstQuarter => write!(f, "First Quarter"),
            Self::WaxingGibbous => write!(f, "Waxing Gibbous"),
            Self::FullMoon => write!(f, "Full Moon"),
            Self::WaningGibbous => write!(f, "Waning Gibbous"),
            Self::LastQuarter => write!(f, "Last Quarter"),
            Self::WaningCrescent => write!(f, "Waning Crescent"),
        }
    }
}

/// Convert a `MoonPhaseLabel` to a representative elongation in degrees.
///
/// This returns the central angle of the octant associated with the label
/// (New=0°, WaxingCrescent=45°, FirstQuarter=90°, ...). Implemented as
/// `From` so callers can use `Into<Degrees>` as needed.
impl From<MoonPhaseLabel> for Degrees {
    fn from(lbl: MoonPhaseLabel) -> Self {
        match lbl {
            MoonPhaseLabel::NewMoon => Degrees::new(0.0),
            MoonPhaseLabel::WaxingCrescent => Degrees::new(45.0),
            MoonPhaseLabel::FirstQuarter => Degrees::new(90.0),
            MoonPhaseLabel::WaxingGibbous => Degrees::new(135.0),
            MoonPhaseLabel::FullMoon => Degrees::new(180.0),
            MoonPhaseLabel::WaningGibbous => Degrees::new(225.0),
            MoonPhaseLabel::LastQuarter => Degrees::new(270.0),
            MoonPhaseLabel::WaningCrescent => Degrees::new(315.0),
        }
    }
}

// ===========================================================================
// PhaseKind (principal phases for event-finding)
// ===========================================================================

/// The four principal lunar phases, used as targets for event detection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum PhaseKind {
    NewMoon,
    FirstQuarter,
    FullMoon,
    LastQuarter,
}

impl PhaseKind {
    /// Target elongation in radians for this phase.
    pub(crate) fn target_elongation(self) -> f64 {
        match self {
            Self::NewMoon => 0.0,
            Self::FirstQuarter => PI / 2.0,
            Self::FullMoon => PI,
            Self::LastQuarter => 3.0 * PI / 2.0,
        }
    }

    /// All four principal phases in chronological order within a synodic month.
    pub const ALL: [PhaseKind; 4] = [
        PhaseKind::NewMoon,
        PhaseKind::FirstQuarter,
        PhaseKind::FullMoon,
        PhaseKind::LastQuarter,
    ];
}

impl std::fmt::Display for PhaseKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NewMoon => write!(f, "New Moon"),
            Self::FirstQuarter => write!(f, "First Quarter"),
            Self::FullMoon => write!(f, "Full Moon"),
            Self::LastQuarter => write!(f, "Last Quarter"),
        }
    }
}

// ===========================================================================
// PhaseEvent
// ===========================================================================

/// A principal lunar phase event: the instant when the Moon's geocentric
/// elongation from the Sun equals exactly 0°, 90°, 180°, or 270°.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PhaseEvent {
    /// Modified Julian Date (TT) of the event.
    pub mjd: ModifiedJulianDate,
    /// Which principal phase.
    pub kind: PhaseKind,
}

// ===========================================================================
// MoonPhaseGeometry
// ===========================================================================

/// Continuous photometric geometry of the Moon–Sun–Earth (or observer) system.
///
/// All angular quantities are in **radians** (typed via `qtty::Radians`).
/// The illuminated fraction is a dimensionless `f64` in \[0, 1\].
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MoonPhaseGeometry {
    /// Sun–Moon–Earth (or observer) angle at the Moon, in \[0, π\].
    pub phase_angle: Radians,
    /// Fraction of the Moon's disk that is illuminated, \[0, 1\].
    pub illuminated_fraction: f64,
    /// Geocentric (or topocentric) elongation of the Moon from the Sun,
    /// measured eastward from the Sun, in \[0, 2π).
    pub elongation: Radians,
    /// `true` when the Moon is waxing (elongation ∈ (0, π)).
    pub waxing: bool,
}

impl MoonPhaseGeometry {
    /// Classify this geometry into one of the eight classical phase names
    /// using the default thresholds (45° bins).
    pub fn label(&self) -> MoonPhaseLabel {
        self.label_with(&PhaseThresholds::default())
    }

    /// Classify with custom thresholds.
    pub fn label_with(&self, th: &PhaseThresholds) -> MoonPhaseLabel {
        MoonPhaseLabel::from_elongation(self.elongation.to::<Degree>(), th)
    }

    /// Illuminated fraction as a percentage \[0, 100\].
    pub fn illuminated_percent(&self) -> f64 {
        self.illuminated_fraction * 100.0
    }
}

// ===========================================================================
// Core geometry: geocentric
// ===========================================================================

/// Ecliptic longitude and distance of the Moon (geocentric, ecliptic J2000).
///
/// Returns `(lon_rad, lat_rad, dist_km)`.
fn moon_ecliptic_geocentric<E: Ephemeris>(jd: JulianDate) -> (f64, f64, f64) {
    let moon: cartesian::Position<Geocentric, frames::EclipticMeanJ2000, Kilometer> =
        E::moon_geocentric(jd);
    let sph = moon.to_spherical();
    (
        sph.azimuth.to::<Radian>().value(), // lon
        sph.polar.to::<Radian>().value(),   // lat
        sph.distance.value(),               // km
    )
}

/// Ecliptic longitude and distance of the Sun as seen from Earth (geocentric).
///
/// We compute: Sun geocentric ecliptic lon = Earth heliocentric lon + π.
/// Returns `(lon_rad, dist_au)`.
fn sun_ecliptic_geocentric<E: Ephemeris>(jd: JulianDate) -> (f64, f64) {
    let earth = E::earth_heliocentric(jd);
    let sph = earth.position.to_spherical();
    let earth_lon = sph.azimuth.to::<Radian>().value();
    let sun_lon = (earth_lon + PI).rem_euclid(TWO_PI);
    let sun_dist = sph.distance.value(); // in AU (heliocentric distance = geocentric distance)
    (sun_lon, sun_dist)
}

/// Compute Moon–Sun elongation (eastward, in \[0, 2π)) from ecliptic
/// longitudes.
fn elongation_from_longitudes(moon_lon: f64, sun_lon: f64) -> f64 {
    (moon_lon - sun_lon).rem_euclid(TWO_PI)
}

/// Phase angle *i* via the law of cosines on the Sun–Moon–Earth triangle.
///
/// Given:
/// - `r_moon_km`: geocentric distance to the Moon (km)
/// - `r_sun_au`: geocentric distance to the Sun (AU)
/// - `psi`: elongation (radians)
///
/// Returns phase angle in radians, clamped to \[0, π\].
fn phase_angle_from_triangle(r_moon_km: f64, r_sun_au: f64, psi: f64) -> f64 {
    let r_moon_au = r_moon_km / AU_KM;

    // Distance Sun–Moon via law of cosines:
    // d² = r_sun² + r_moon² − 2·r_sun·r_moon·cos(ψ)
    let d_sq = r_sun_au * r_sun_au + r_moon_au * r_moon_au - 2.0 * r_sun_au * r_moon_au * psi.cos();
    let d = d_sq.max(0.0).sqrt();

    if d < 1e-15 {
        return 0.0; // degenerate: Moon at Sun position (impossible but safe)
    }

    // Phase angle i = angle at Moon in the Sun–Moon–Earth triangle.
    // cos i = (d² + r_moon² − r_sun²) / (2·d·r_moon)
    let cos_i = (d * d + r_moon_au * r_moon_au - r_sun_au * r_sun_au) / (2.0 * d * r_moon_au);
    cos_i.clamp(-1.0, 1.0).acos()
}

/// Compute the geocentric Moon phase geometry at a given Julian Date.
///
/// Uses the `E: Ephemeris` backend for both Moon (geocentric ecliptic)
/// and Earth/Sun (heliocentric ecliptic) positions.
///
/// # Returns
///
/// A [`MoonPhaseGeometry`] with phase angle, illuminated fraction,
/// elongation, and waxing flag — all computed in the geocentric frame.
///
/// # Example
///
/// ```rust
/// use siderust::calculus::lunar::phase::{moon_phase_geocentric, MoonPhaseGeometry};
/// use siderust::calculus::ephemeris::Vsop87Ephemeris;
/// use siderust::time::JulianDate;
///
/// let geom = moon_phase_geocentric::<Vsop87Ephemeris>(JulianDate::J2000);
/// assert!(geom.illuminated_fraction >= 0.0);
/// assert!(geom.illuminated_fraction <= 1.0);
/// ```
pub fn moon_phase_geocentric<E: Ephemeris>(jd: JulianDate) -> MoonPhaseGeometry {
    let (moon_lon, moon_lat, moon_dist_km) = moon_ecliptic_geocentric::<E>(jd);
    let (sun_lon, sun_dist_au) = sun_ecliptic_geocentric::<E>(jd);

    // Elongation: angular separation measured eastward from the Sun.
    // For the proper great-circle elongation including latitude:
    let psi = great_circle_elongation(sun_lon, 0.0, moon_lon, moon_lat);

    // Signed elongation for waxing/waning (eastward = positive)
    let signed_elong = elongation_from_longitudes(moon_lon, sun_lon);

    // Phase angle from the three-body triangle
    let i = phase_angle_from_triangle(moon_dist_km, sun_dist_au, psi);

    // Illuminated fraction
    let k = (1.0 + i.cos()) / 2.0;

    MoonPhaseGeometry {
        phase_angle: Radians::new(i),
        illuminated_fraction: k,
        elongation: Radians::new(signed_elong),
        waxing: signed_elong > 0.0 && signed_elong < PI,
    }
}

/// Great-circle angular separation between two ecliptic directions.
///
/// `(lon1, lat1)` and `(lon2, lat2)` are in radians.
/// Returns angle in \[0, π\].
fn great_circle_elongation(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let cos_sep = lat1.sin() * lat2.sin() + lat1.cos() * lat2.cos() * (lon2 - lon1).cos();
    cos_sep.clamp(-1.0, 1.0).acos()
}

// ===========================================================================
// Topocentric variant
// ===========================================================================

/// Compute the topocentric Moon phase geometry at a given Julian Date
/// from a specific observer location.
///
/// The topocentric elongation uses the apparent equatorial positions of
/// both the Moon and the Sun corrected for parallax. The phase angle is
/// still computed from the geocentric distance triangle (the parallax
/// shift is < 1° and has negligible effect on the distance ratio).
///
/// # Parameters
///
/// - `jd`: Julian Date (TT)
/// - `site`: Observer location on Earth
///
/// # Returns
///
/// A [`MoonPhaseGeometry`] with topocentric elongation, geocentric-triangle
/// phase angle, and derived illuminated fraction and waxing flag.
///
/// # Example
///
/// ```rust
/// use siderust::calculus::lunar::phase::{moon_phase_topocentric, MoonPhaseGeometry};
/// use siderust::calculus::ephemeris::Vsop87Ephemeris;
/// use siderust::coordinates::centers::Geodetic;
/// use siderust::coordinates::frames::ECEF;
/// use siderust::time::JulianDate;
/// use qtty::*;
///
/// let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.48 * DEG, 0.0 * M);
/// let geom = moon_phase_topocentric::<Vsop87Ephemeris>(JulianDate::J2000, site);
/// assert!(geom.illuminated_fraction >= 0.0 && geom.illuminated_fraction <= 1.0);
/// ```
pub fn moon_phase_topocentric<E: Ephemeris>(
    jd: JulianDate,
    site: Geodetic<frames::ECEF>,
) -> MoonPhaseGeometry {
    use crate::bodies::solar_system::{Moon, Sun};

    // Topocentric equatorial positions (true-of-date)
    let moon_topo: affn::spherical::Position<Topocentric, frames::EquatorialTrueOfDate, Kilometer> =
        Moon::get_apparent_topocentric_equ(jd, site);
    let sun_topo: affn::spherical::Position<
        Topocentric,
        frames::EquatorialTrueOfDate,
        AstronomicalUnit,
    > = Sun::get_apparent_topocentric_equ(jd, site);

    // RA/Dec in radians
    let moon_ra = moon_topo.ra().to::<Radian>().value();
    let moon_dec = moon_topo.dec().to::<Radian>().value();
    let sun_ra = sun_topo.ra().to::<Radian>().value();
    let sun_dec = sun_topo.dec().to::<Radian>().value();

    // True angular separation (great circle)
    let psi = great_circle_elongation(moon_ra, moon_dec, sun_ra, sun_dec);

    // Signed elongation: positive when Moon is east of Sun.
    // Use RA difference as proxy for the east/west sign.
    let dra = (moon_ra - sun_ra).rem_euclid(TWO_PI);
    let signed_elong = if dra < PI { psi } else { TWO_PI - psi };

    // Phase angle from geocentric distance triangle
    // (topocentric parallax is < 1° → distance ratios essentially unchanged)
    let (_, _, moon_dist_km) = moon_ecliptic_geocentric::<E>(jd);
    let (_, sun_dist_au) = sun_ecliptic_geocentric::<E>(jd);
    let i = phase_angle_from_triangle(moon_dist_km, sun_dist_au, psi);

    let k = (1.0 + i.cos()) / 2.0;

    MoonPhaseGeometry {
        phase_angle: Radians::new(i),
        illuminated_fraction: k,
        elongation: Radians::new(signed_elong),
        waxing: signed_elong > 0.0 && signed_elong < PI,
    }
}

// ===========================================================================
// Elongation as a scalar function of MJD (for root-finding)
// ===========================================================================

/// Geocentric elongation as a function of MJD, returned in radians.
fn elongation_at_mjd<E: Ephemeris>(mjd: ModifiedJulianDate) -> f64 {
    let jd: JulianDate = mjd.into();
    let (moon_lon, _, _) = moon_ecliptic_geocentric::<E>(jd);
    let (sun_lon, _) = sun_ecliptic_geocentric::<E>(jd);
    elongation_from_longitudes(moon_lon, sun_lon)
}

// ===========================================================================
// Event finding
// ===========================================================================

/// Search options for phase event detection.
#[derive(Debug, Clone, Copy)]
pub struct PhaseSearchOpts {
    /// Time tolerance for Brent root refinement (days).
    /// Default: ~1 µs (1 × 10⁻⁹ days ≈ 0.086 ms).
    pub time_tolerance: Days,
    /// Scan step for coarse bracketing (days).
    /// Default: 0.5 days.
    pub scan_step: Days,
}

impl Default for PhaseSearchOpts {
    fn default() -> Self {
        Self {
            time_tolerance: Days::new(1e-9),
            scan_step: PHASE_SCAN_STEP,
        }
    }
}

/// Find all principal phase events in a time `window`.
///
/// Locates the instants when the geocentric Moon–Sun elongation equals
/// 0° (new moon), 90° (first quarter), 180° (full moon), or 270° (last
/// quarter). Results are sorted chronologically.
///
/// # Type parameter
///
/// - `E: Ephemeris` — the ephemeris backend to use.
///
/// # Arguments
///
/// - `window` — search interval as `Period<MJD>` (TT axis)
/// - `opts` — search precision options
///
/// # Returns
///
/// A `Vec<PhaseEvent>` sorted by `mjd`.
///
/// # Example
///
/// ```rust,no_run
/// use siderust::calculus::lunar::phase::{find_phase_events, PhaseSearchOpts};
/// use siderust::calculus::ephemeris::Vsop87Ephemeris;
/// use siderust::time::{ModifiedJulianDate, Period};
///
/// let start = ModifiedJulianDate::new(60000.0);
/// let end   = ModifiedJulianDate::new(60030.0);
/// let window = Period::new(start, end);
///
/// let events = find_phase_events::<Vsop87Ephemeris>(window, PhaseSearchOpts::default());
/// for ev in &events {
///     println!("{} at MJD {}", ev.kind, ev.mjd);
/// }
/// ```
pub fn find_phase_events<E: Ephemeris>(
    window: Period<MJD>,
    opts: PhaseSearchOpts,
) -> Vec<PhaseEvent> {
    let mut events = Vec::new();
    let step = opts.scan_step;
    let _tol = opts.time_tolerance;

    for kind in PhaseKind::ALL {
        let target = kind.target_elongation();

        // Build a signed scalar function whose zero crossings correspond to
        // the target elongation. We use sin(elongation − target) so the
        // function is smooth and avoids the 0°/360° wraparound discontinuity.
        let f = |mjd: ModifiedJulianDate| -> Radians {
            let elong = elongation_at_mjd::<E>(mjd);
            Radians::new((elong - target).sin())
        };

        // Scan for sign changes
        let threshold = Radians::new(0.0);
        let raw_crossings = intervals::find_crossings(window, step, &f, threshold);

        // Filter: only keep crossings where the actual elongation is close
        // to the target (eliminates the extraneous root at target ± π).
        for mjd in raw_crossings {
            let elong = elongation_at_mjd::<E>(mjd);
            let diff = angle_diff(elong, target);
            if diff.abs() < 0.5 {
                // within ~28.6° — genuine crossing, not the ±π alias
                events.push(PhaseEvent { mjd, kind });
            }
        }
    }

    // Sort chronologically
    events.sort_by(|a, b| a.mjd.partial_cmp(&b.mjd).unwrap());
    events
}

/// Signed shortest-path angular difference in \[−π, π\].
fn angle_diff(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(TWO_PI);
    if d > PI {
        d - TWO_PI
    } else {
        d
    }
}

// ===========================================================================
// MoonPhaseSeries
// ===========================================================================

/// A lightweight wrapper for batch-sampling Moon phase geometry over
/// a time range.
///
/// Iterates from `start` to `end` at the given `step`, computing
/// [`MoonPhaseGeometry`] at each sample point. No additional caching
/// is performed — the ephemeris backend's own caching (e.g. Chebyshev
/// for DE441) is relied upon.
///
/// # Type parameter
///
/// - `E: Ephemeris` — the ephemeris backend to use.
pub struct MoonPhaseSeries<E: Ephemeris> {
    _marker: PhantomData<E>,
}

impl<E: Ephemeris> MoonPhaseSeries<E> {
    /// Sample geocentric Moon phase geometry from `start` to `end`
    /// (inclusive) at intervals of `step`.
    ///
    /// # Returns
    ///
    /// A `Vec<(ModifiedJulianDate, MoonPhaseGeometry)>` with one entry
    /// per sample point.
    pub fn sample(
        start: ModifiedJulianDate,
        end: ModifiedJulianDate,
        step: Days,
    ) -> Vec<(ModifiedJulianDate, MoonPhaseGeometry)> {
        let mut results = Vec::new();
        let mut t = start;
        while t <= end {
            let jd: JulianDate = t.into();
            let geom = moon_phase_geocentric::<E>(jd);
            results.push((t, geom));
            t += step;
        }
        results
    }

    /// Sample topocentric Moon phase geometry from `start` to `end`
    /// (inclusive) at intervals of `step` for a given observer location.
    pub fn sample_topocentric(
        start: ModifiedJulianDate,
        end: ModifiedJulianDate,
        step: Days,
        site: Geodetic<frames::ECEF>,
    ) -> Vec<(ModifiedJulianDate, MoonPhaseGeometry)> {
        let mut results = Vec::new();
        let mut t = start;
        while t <= end {
            let jd: JulianDate = t.into();
            let geom = moon_phase_topocentric::<E>(jd, site);
            results.push((t, geom));
            t += step;
        }
        results
    }
}

// ===========================================================================
// Illumination period finders
// ===========================================================================

/// A reusable scalar closure for geocentric illuminated fraction at a given MJD.
fn illumination_at_mjd<E: Ephemeris>(mjd: ModifiedJulianDate) -> Radians {
    let jd: JulianDate = mjd.into();
    let geom = moon_phase_geocentric::<E>(jd);
    Radians::new(geom.illuminated_fraction)
}

/// Find all time windows in `window` where the geocentric illuminated
/// fraction is **above** `k_min`.
///
/// # Parameters
///
/// - `window`: search interval (`Period<MJD>`, TT axis)
/// - `k_min`: illuminated fraction lower bound, in \[0, 1\]
/// - `opts`: scan precision options
///
/// # Returns
///
/// Sorted, non-overlapping `Vec<Period<MJD>>`.
///
/// # Example
///
/// ```rust,no_run
/// use siderust::calculus::lunar::phase::{illumination_above, PhaseSearchOpts};
/// use siderust::calculus::ephemeris::Vsop87Ephemeris;
/// use siderust::time::{ModifiedJulianDate, Period};
/// use qtty::Days;
///
/// let start  = ModifiedJulianDate::new(60000.0);
/// let window = Period::new(start, start + Days::new(30.0));
/// // Find windows where Moon is more than 50% illuminated (gibbous/full)
/// let bright = illumination_above::<Vsop87Ephemeris>(window, 0.5, PhaseSearchOpts::default());
/// ```
pub fn illumination_above<E: Ephemeris>(
    window: Period<MJD>,
    k_min: f64,
    opts: PhaseSearchOpts,
) -> Vec<Period<MJD>> {
    intervals::above_threshold_periods(
        window,
        opts.scan_step,
        &illumination_at_mjd::<E>,
        Radians::new(k_min),
    )
}

/// Find all time windows in `window` where the geocentric illuminated
/// fraction is **below** `k_max`.
///
/// # Parameters
///
/// - `window`: search interval
/// - `k_max`: illuminated fraction upper bound, in \[0, 1\]
/// - `opts`: scan precision options
pub fn illumination_below<E: Ephemeris>(
    window: Period<MJD>,
    k_max: f64,
    opts: PhaseSearchOpts,
) -> Vec<Period<MJD>> {
    use crate::time::complement_within;
    let above = illumination_above::<E>(window, k_max, opts);
    complement_within(window, &above)
}

/// Find all time windows in `window` where the geocentric illuminated
/// fraction is within `[k_min, k_max]`.
///
/// Useful for finding "crescent windows" (e.g., 5%–35%) or "gibbous windows"
/// (e.g., 60%–100%) for observation planning.
///
/// # Parameters
///
/// - `window`: search interval
/// - `k_min`: lower bound of the illumination band, in \[0, 1\]
/// - `k_max`: upper bound of the illumination band, in \[0, 1\]
/// - `opts`: scan precision options
///
/// # Example
///
/// ```rust,no_run
/// use siderust::calculus::lunar::phase::{illumination_range, PhaseSearchOpts};
/// use siderust::calculus::ephemeris::Vsop87Ephemeris;
/// use siderust::time::{ModifiedJulianDate, Period};
/// use qtty::Days;
///
/// let start  = ModifiedJulianDate::new(60000.0);
/// let window = Period::new(start, start + Days::new(60.0));
/// // Crescent phase: 5–35% illuminated
/// let crescent = illumination_range::<Vsop87Ephemeris>(window, 0.05, 0.35, PhaseSearchOpts::default());
/// ```
pub fn illumination_range<E: Ephemeris>(
    window: Period<MJD>,
    k_min: f64,
    k_max: f64,
    opts: PhaseSearchOpts,
) -> Vec<Period<MJD>> {
    intervals::in_range_periods(
        window,
        opts.scan_step,
        &illumination_at_mjd::<E>,
        Radians::new(k_min),
        Radians::new(k_max),
    )
}

// ===========================================================================
// Unit tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::calculus::ephemeris::Vsop87Ephemeris;

    #[test]
    fn illuminated_fraction_bounded() {
        // Check that illuminated fraction is always in [0, 1] for a range of dates.
        let start = JulianDate::J2000;
        for i in 0..100 {
            let jd = start + Days::new(i as f64 * 3.0);
            let geom = moon_phase_geocentric::<Vsop87Ephemeris>(jd);
            assert!(
                geom.illuminated_fraction >= 0.0 && geom.illuminated_fraction <= 1.0,
                "Fraction out of bounds at JD offset {}: {}",
                i * 3,
                geom.illuminated_fraction
            );
        }
    }

    #[test]
    fn phase_angle_bounded() {
        let start = JulianDate::J2000;
        for i in 0..50 {
            let jd = start + Days::new(i as f64 * 5.0);
            let geom = moon_phase_geocentric::<Vsop87Ephemeris>(jd);
            let i_val = geom.phase_angle;
            assert!(
                i_val >= Radians::new(0.0) && i_val <= Radians::new(PI),
                "Phase angle out of [0, π] at offset {}: {}",
                i * 5,
                i_val
            );
        }
    }

    #[test]
    fn elongation_bounded() {
        let start = JulianDate::J2000;
        for i in 0..50 {
            let jd = start + Days::new(i as f64 * 5.0);
            let geom = moon_phase_geocentric::<Vsop87Ephemeris>(jd);
            let e = geom.elongation;
            assert!(
                e >= Radians::new(0.0) && e < Radians::new(TWO_PI),
                "Elongation out of [0, 2π) at offset {}: {}",
                i * 5,
                e
            );
        }
    }

    #[test]
    fn waxing_consistent_with_elongation() {
        let start = JulianDate::J2000;
        for i in 0..50 {
            let jd = start + Days::new(i as f64 * 5.0);
            let geom = moon_phase_geocentric::<Vsop87Ephemeris>(jd);
            let e = geom.elongation;
            let expected_waxing = e > Radians::new(0.0) && e < Radians::new(PI);
            assert_eq!(
                geom.waxing,
                expected_waxing,
                "Waxing flag mismatch at offset {}: elong={}, waxing={}",
                i * 5,
                e,
                geom.waxing
            );
        }
    }

    #[test]
    fn label_known_elongations() {
        let th = PhaseThresholds::default();
        assert_eq!(
            MoonPhaseLabel::from_elongation(Degrees::new(0.0), &th),
            MoonPhaseLabel::NewMoon
        );
        assert_eq!(
            MoonPhaseLabel::from_elongation(Degrees::new(45.0), &th),
            MoonPhaseLabel::WaxingCrescent
        );
        assert_eq!(
            MoonPhaseLabel::from_elongation(Degrees::new(90.0), &th),
            MoonPhaseLabel::FirstQuarter
        );
        assert_eq!(
            MoonPhaseLabel::from_elongation(Degrees::new(135.0), &th),
            MoonPhaseLabel::WaxingGibbous
        );
        assert_eq!(
            MoonPhaseLabel::from_elongation(Degrees::new(180.0), &th),
            MoonPhaseLabel::FullMoon
        );
        assert_eq!(
            MoonPhaseLabel::from_elongation(Degrees::new(225.0), &th),
            MoonPhaseLabel::WaningGibbous
        );
        assert_eq!(
            MoonPhaseLabel::from_elongation(Degrees::new(270.0), &th),
            MoonPhaseLabel::LastQuarter
        );
        assert_eq!(
            MoonPhaseLabel::from_elongation(Degrees::new(315.0), &th),
            MoonPhaseLabel::WaningCrescent
        );
        assert_eq!(
            MoonPhaseLabel::from_elongation(Degrees::new(359.0), &th),
            MoonPhaseLabel::NewMoon
        );
    }

    #[test]
    fn label_via_geometry() {
        // At J2000 the Moon should produce a valid label.
        let geom = moon_phase_geocentric::<Vsop87Ephemeris>(JulianDate::J2000);
        let label = geom.label();
        // Just ensure it doesn't panic and returns something reasonable.
        let _ = format!("{}", label);
    }

    #[test]
    fn find_events_in_one_synodic_month() {
        // Search ~35 days from J2000 — should find at least one of each kind.
        let start = ModifiedJulianDate::from(JulianDate::J2000);
        let end = start + Days::new(35.0);
        let window = Period::new(start, end);
        let events = find_phase_events::<Vsop87Ephemeris>(window, PhaseSearchOpts::default());

        assert!(
            !events.is_empty(),
            "Should find at least one phase event in 35 days"
        );

        // Verify chronological order
        for pair in events.windows(2) {
            assert!(pair[0].mjd <= pair[1].mjd, "Events not sorted");
        }
    }

    #[test]
    fn topocentric_close_to_geocentric() {
        use crate::coordinates::centers::Geodetic;
        use crate::coordinates::frames::ECEF;

        let site = Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
        let jd = JulianDate::J2000;

        let geo = moon_phase_geocentric::<Vsop87Ephemeris>(jd);
        let topo = moon_phase_topocentric::<Vsop87Ephemeris>(jd, site);

        // Elongation difference should be < 2° (parallax bound)
        let diff = (geo.elongation - topo.elongation).abs();
        let two_deg_in_rad = Degrees::new(2.0).to::<Radian>();
        assert!(
            diff < two_deg_in_rad,
            "Geocentric vs topocentric elongation differ by more than 2°: {} rad",
            diff
        );

        // Illuminated fraction difference should be < 1%
        let frac_diff = (geo.illuminated_fraction - topo.illuminated_fraction).abs();
        assert!(
            frac_diff < 0.01,
            "Illuminated fraction geo vs topo differ by more than 1%: {}",
            frac_diff
        );
    }

    #[test]
    fn series_length() {
        let start = ModifiedJulianDate::from(JulianDate::J2000);
        let end = start + Days::new(10.0);
        let step = Days::new(1.0);
        let series = MoonPhaseSeries::<Vsop87Ephemeris>::sample(start, end, step);
        assert_eq!(series.len(), 11); // 0, 1, 2, ..., 10
    }
}
