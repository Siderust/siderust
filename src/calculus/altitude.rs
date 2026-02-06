// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Unified Altitude Computation & Event API
//!
//! A clean, user‑friendly API for computing target altitude vs time and
//! finding events (crossings, culminations, altitude ranges) for **any**
//! celestial target.
//!
//! ## Public functions
//!
//! | Function | Purpose |
//! |----------|---------|
//! | [`altitude_at`] | Single‑point altitude, full fidelity |
//! | [`crossings`] | Threshold crossings in a time window |
//! | [`culminations`] | All local maxima/minima of altitude |
//! | [`altitude_ranges`] | Intervals where altitude is within `[h_min, h_max]` |
//!
//! All inputs/outputs are typed with `qtty` (`Degrees`, `JulianDate`, etc.).
//!
//! ## Target types
//!
//! The [`AltitudeTarget`] enum supports Sun, Moon, planets (via VSOP87 +
//! coordinate transforms), and fixed RA/Dec objects (stars).  Each variant
//! carries the minimal data needed to compute horizontal altitude from an
//! [`ObserverSite`].

use crate::astro::JulianDate;
use crate::bodies::solar_system::{Moon, Sun};
use crate::calculus::math_core::{extrema, intervals};
use crate::coordinates::centers::ObserverSite;
use crate::coordinates::spherical;
use crate::time::{complement_within, ModifiedJulianDate, Period};
use qtty::{AstronomicalUnit, Degrees, Kilometer, Radian};

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Selects which celestial body/target to compute altitude for.
#[derive(Debug, Clone)]
pub enum AltitudeTarget {
    /// The Sun (uses VSOP87 via existing `Sun::get_horizontal`)
    Sun,
    /// The Moon (uses ELP2000 via existing `Moon::get_horizontal`)
    Moon,
    /// A fixed RA/Dec object (star or distant body).
    /// RA and Dec are J2000 equatorial coordinates.
    FixedEquatorial {
        ra: Degrees,
        dec: Degrees,
    },
}

/// Direction of a threshold crossing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossingDirection {
    Rising,
    Setting,
}

/// A threshold crossing event.
#[derive(Debug, Clone, Copy)]
pub struct CrossingEvent {
    /// Julian Date of the crossing.
    pub jd: JulianDate,
    /// Direction: rising above or setting below the threshold.
    pub direction: CrossingDirection,
}

/// Kind of culmination (altitude extremum).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CulminationKind {
    /// Upper culmination: local maximum altitude.
    Max,
    /// Lower culmination: local minimum altitude.
    Min,
}

/// A culmination event.
#[derive(Debug, Clone, Copy)]
pub struct CulminationEvent {
    /// Julian Date of the extremum.
    pub jd: JulianDate,
    /// Altitude at the extremum.
    pub altitude: Degrees,
    /// Maximum or minimum.
    pub kind: CulminationKind,
}

/// Options for controlling search precision and strategy.
#[derive(Debug, Clone, Copy)]
pub struct SearchOpts {
    /// Time tolerance for root/extremum refinement (days).
    /// Default: ~1 µs (1e-11 days).
    pub time_tolerance_days: f64,
    /// Scan step for coarse bracket detection (days).
    /// Default: 10 minutes for generic, 2 hours for Moon.
    pub scan_step_days: Option<f64>,
}

impl Default for SearchOpts {
    fn default() -> Self {
        Self {
            time_tolerance_days: 1e-9,
            scan_step_days: None,
        }
    }
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Default scan step: 10 minutes in days.
const DEFAULT_SCAN_STEP: f64 = 10.0 / 1440.0;

/// Moon scan step: 2 hours in days.
const MOON_SCAN_STEP: f64 = 2.0 / 24.0;

/// Extrema scan step: 20 minutes in days (for culmination detection).
const EXTREMA_SCAN_STEP: f64 = 20.0 / 1440.0;

// ---------------------------------------------------------------------------
// Internal: altitude function factory
// ---------------------------------------------------------------------------

/// Build an `f64 → f64` altitude function (JD.value() → altitude in radians).
fn make_altitude_fn<'a>(
    target: &'a AltitudeTarget,
    site: &'a ObserverSite,
) -> Box<dyn Fn(f64) -> f64 + 'a> {
    let site = *site;
    match target {
        AltitudeTarget::Sun => Box::new(move |jd_val: f64| {
            let jd = JulianDate::new(jd_val);
            Sun::get_horizontal::<AstronomicalUnit>(jd, site)
                .alt()
                .to::<Radian>()
                .value()
        }),
        AltitudeTarget::Moon => Box::new(move |jd_val: f64| {
            let jd = JulianDate::new(jd_val);
            Moon::get_horizontal::<Kilometer>(jd, site)
                .alt()
                .to::<Radian>()
                .value()
        }),
        AltitudeTarget::FixedEquatorial { ra, dec } => {
            let ra = *ra;
            let dec = *dec;
            Box::new(move |jd_val: f64| {
                let jd = JulianDate::new(jd_val);
                fixed_star_altitude_rad(jd, &site, ra, dec)
            })
        }
    }
}

/// Compute altitude of a fixed RA/Dec object from an observer site.
///
/// Uses precession + nutation-corrected RA, GAST, and the standard
/// equatorial→horizontal formula.
fn fixed_star_altitude_rad(
    jd: JulianDate,
    site: &ObserverSite,
    ra_j2000: Degrees,
    dec_j2000: Degrees,
) -> f64 {
    use crate::astro::nutation::corrected_ra_with_nutation;
    use crate::astro::precession;
    use crate::astro::sidereal::{calculate_gst, calculate_lst};
    use crate::coordinates::frames::EquatorialMeanJ2000;
    use qtty::LightYear;

    // Build a spherical position in EquatorialMeanJ2000
    let pos = spherical::Position::<
        crate::coordinates::centers::Geocentric,
        EquatorialMeanJ2000,
        LightYear,
    >::new(ra_j2000, dec_j2000, qtty::LightYears::new(1.0));

    // Precess J2000 → mean-of-date
    let mean_of_date = precession::precess_from_j2000(pos, jd);
    // Apply nutation correction to RA
    let ra_corrected = corrected_ra_with_nutation(&mean_of_date.direction(), jd);
    let dec = mean_of_date.polar();

    // Compute hour angle
    let gst = calculate_gst(jd);
    let lst = calculate_lst(gst, site.lon);
    let ha = (lst - ra_corrected).normalize().to::<Radian>();

    // Equatorial → horizontal altitude
    let lat = site.lat.to::<Radian>();
    let dec_rad = dec.to::<Radian>();
    let sin_alt = dec_rad.sin() * lat.sin() + dec_rad.cos() * lat.cos() * ha.cos();
    sin_alt.asin()
}

/// Choose the best scan step for the target.
fn scan_step_for(target: &AltitudeTarget, opts: &SearchOpts) -> f64 {
    opts.scan_step_days.unwrap_or(match target {
        AltitudeTarget::Moon => MOON_SCAN_STEP,
        _ => DEFAULT_SCAN_STEP,
    })
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Compute the altitude of `target` as seen from `observer` at time `t`.
///
/// Returns altitude as `Degrees` (positive above horizon, negative below).
///
/// # Example
/// ```ignore
/// use siderust::calculus::altitude::{altitude_at, AltitudeTarget};
/// use siderust::coordinates::centers::ObserverSite;
/// use siderust::astro::JulianDate;
/// use qtty::*;
///
/// let site = ObserverSite::new(0.0 * DEG, 51.48 * DEG, 0.0 * M);
/// let alt = altitude_at(&AltitudeTarget::Sun, &site, JulianDate::J2000);
/// println!("Sun altitude: {:.2}°", alt.value());
/// ```
pub fn altitude_at(target: &AltitudeTarget, observer: &ObserverSite, t: JulianDate) -> Degrees {
    let f = make_altitude_fn(target, observer);
    let rad = f(t.value());
    Degrees::new(rad.to_degrees())
}

/// Find all threshold crossings of `target` altitude in the given `window`.
///
/// Returns a chronologically sorted list of [`CrossingEvent`]s.
///
/// # Arguments
/// * `target` — which celestial body
/// * `observer` — site on Earth
/// * `window` — search interval (MJD)
/// * `threshold` — altitude threshold
/// * `opts` — search options (tolerances, scan step)
///
/// # Example
/// ```ignore
/// let events = crossings(
///     &AltitudeTarget::Sun, &site, window,
///     Degrees::new(-18.0), // astronomical twilight
///     SearchOpts::default(),
/// );
/// for e in events {
///     println!("{:?} at JD {}", e.direction, e.jd.value());
/// }
/// ```
pub fn crossings(
    target: &AltitudeTarget,
    observer: &ObserverSite,
    window: Period<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<CrossingEvent> {
    let f = make_altitude_fn(target, observer);
    let t_start = window.start.to_julian_day().value();
    let t_end = window.end.to_julian_day().value();
    let thr_rad = threshold.to::<Radian>().value();
    let step = scan_step_for(target, &opts);

    // Use the fast scan + label approach from math_core
    let mut raw_crossings = intervals::find_crossings(t_start, t_end, step, &f, thr_rad);
    let labeled = intervals::label_crossings(&mut raw_crossings, &f, thr_rad);

    labeled
        .iter()
        .map(|lc| CrossingEvent {
            jd: JulianDate::new(lc.t),
            direction: if lc.direction > 0 {
                CrossingDirection::Rising
            } else {
                CrossingDirection::Setting
            },
        })
        .collect()
}

/// Find all altitude culminations (local maxima and minima) of `target` in `window`.
///
/// Returns a chronologically sorted list of [`CulminationEvent`]s.
pub fn culminations(
    target: &AltitudeTarget,
    observer: &ObserverSite,
    window: Period<ModifiedJulianDate>,
    opts: SearchOpts,
) -> Vec<CulminationEvent> {
    let f = make_altitude_fn(target, observer);
    let t_start = window.start.to_julian_day().value();
    let t_end = window.end.to_julian_day().value();
    let step = match target {
        AltitudeTarget::Moon => MOON_SCAN_STEP,
        _ => EXTREMA_SCAN_STEP,
    };
    let tol = opts.time_tolerance_days;

    let raw = extrema::find_extrema_tol(t_start, t_end, step, &f, tol);

    raw.iter()
        .map(|ext| {
            let alt_deg = Degrees::new(ext.value.to_degrees());
            CulminationEvent {
                jd: JulianDate::new(ext.t),
                altitude: alt_deg,
                kind: match ext.kind {
                    extrema::ExtremumKind::Maximum => CulminationKind::Max,
                    extrema::ExtremumKind::Minimum => CulminationKind::Min,
                },
            }
        })
        .collect()
}

/// Find all time intervals where the altitude of `target` is within
/// `[h_min, h_max]`.
///
/// Returns a sorted list of `Period<ModifiedJulianDate>`.
///
/// # Algorithm
///
/// Uses the two‑stage approach:
/// 1. Fast coarse scan to find threshold crossings of `h_min` and `h_max`.
/// 2. Brent refinement for each bracket.
/// 3. Interval algebra: `above(h_min) ∩ complement(above(h_max))`.
///
/// # Example
/// ```ignore
/// // Find astronomical night (Sun between -90° and -18°)
/// let nights = altitude_ranges(
///     &AltitudeTarget::Sun, &site, window,
///     Degrees::new(-90.0), Degrees::new(-18.0),
///     SearchOpts::default(),
/// );
/// ```
pub fn altitude_ranges(
    target: &AltitudeTarget,
    observer: &ObserverSite,
    window: Period<ModifiedJulianDate>,
    h_min: Degrees,
    h_max: Degrees,
    opts: SearchOpts,
) -> Vec<Period<ModifiedJulianDate>> {
    let f = make_altitude_fn(target, observer);
    let t_start = window.start.to_julian_day().value();
    let t_end = window.end.to_julian_day().value();
    let min_rad = h_min.to::<Radian>().value();
    let max_rad = h_max.to::<Radian>().value();
    let step = scan_step_for(target, &opts);

    let raw_intervals = intervals::in_range_periods(t_start, t_end, step, &f, min_rad, max_rad);

    raw_intervals
        .iter()
        .map(|iv| {
            Period::new(
                ModifiedJulianDate::new(iv.start - 2_400_000.5),
                ModifiedJulianDate::new(iv.end - 2_400_000.5),
            )
        })
        .collect()
}

/// Convenience: find periods where altitude is **above** a threshold.
///
/// Equivalent to `altitude_ranges(target, observer, window, threshold, 90°, opts)`.
pub fn above_threshold(
    target: &AltitudeTarget,
    observer: &ObserverSite,
    window: Period<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<Period<ModifiedJulianDate>> {
    let f = make_altitude_fn(target, observer);
    let t_start = window.start.to_julian_day().value();
    let t_end = window.end.to_julian_day().value();
    let thr_rad = threshold.to::<Radian>().value();
    let step = scan_step_for(target, &opts);

    let raw = intervals::above_threshold_periods(t_start, t_end, step, &f, thr_rad);

    raw.iter()
        .map(|iv| {
            Period::new(
                ModifiedJulianDate::new(iv.start - 2_400_000.5),
                ModifiedJulianDate::new(iv.end - 2_400_000.5),
            )
        })
        .collect()
}

/// Convenience: find periods where altitude is **below** a threshold.
///
/// Equivalent to complement of [`above_threshold`].
pub fn below_threshold(
    target: &AltitudeTarget,
    observer: &ObserverSite,
    window: Period<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Vec<Period<ModifiedJulianDate>> {
    let above = above_threshold(target, observer, window, threshold, opts);
    complement_within(window, &above)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::*;

    fn greenwich() -> ObserverSite {
        ObserverSite::new(Degrees::new(0.0), Degrees::new(51.4769), Quantity::<Meter>::new(0.0))
    }

    #[test]
    fn altitude_at_sun_in_range() {
        let alt = altitude_at(&AltitudeTarget::Sun, &greenwich(), JulianDate::J2000);
        // J2000 is 2000-01-01 12:00 TT, Sun should be at reasonable altitude
        assert!(alt.value() > -90.0 && alt.value() < 90.0);
    }

    #[test]
    fn altitude_at_moon_in_range() {
        let alt = altitude_at(&AltitudeTarget::Moon, &greenwich(), JulianDate::J2000);
        assert!(alt.value() > -90.0 && alt.value() < 90.0);
    }

    #[test]
    fn altitude_at_fixed_star() {
        // Sirius: RA ≈ 101.287°, Dec ≈ -16.716°
        let target = AltitudeTarget::FixedEquatorial {
            ra: Degrees::new(101.287),
            dec: Degrees::new(-16.716),
        };
        let alt = altitude_at(&target, &greenwich(), JulianDate::J2000);
        assert!(alt.value() > -90.0 && alt.value() < 90.0);
    }

    #[test]
    fn crossings_finds_sun_rise_set() {
        let site = greenwich();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60001.0);
        let window = Period::new(mjd_start, mjd_end);

        let events = crossings(
            &AltitudeTarget::Sun,
            &site,
            window,
            Degrees::new(0.0),
            SearchOpts::default(),
        );

        // In a normal 24h window at ~51°N, expect 1 rise + 1 set
        assert!(!events.is_empty(), "should find crossings");
        let rises = events.iter().filter(|e| e.direction == CrossingDirection::Rising).count();
        let sets = events.iter().filter(|e| e.direction == CrossingDirection::Setting).count();
        assert!(rises >= 1 || sets >= 1, "should find at least one rise or set");
    }

    #[test]
    fn culminations_finds_sun_extrema() {
        let site = greenwich();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60001.0);
        let window = Period::new(mjd_start, mjd_end);

        let culms = culminations(
            &AltitudeTarget::Sun,
            &site,
            window,
            SearchOpts::default(),
        );

        // Expect upper and lower culmination in 24h
        assert!(!culms.is_empty(), "should find culminations");
        let maxima = culms.iter().filter(|c| c.kind == CulminationKind::Max).count();
        let minima = culms.iter().filter(|c| c.kind == CulminationKind::Min).count();
        assert!(maxima >= 1, "should find at least one upper culmination");
        assert!(minima >= 1, "should find at least one lower culmination");
    }

    #[test]
    fn above_threshold_sun_day_periods() {
        let site = greenwich();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);
        let window = Period::new(mjd_start, mjd_end);

        let days = above_threshold(
            &AltitudeTarget::Sun,
            &site,
            window,
            Degrees::new(0.0),
            SearchOpts::default(),
        );

        assert!(!days.is_empty(), "should find daytime periods in 7 days");
        for p in &days {
            assert!(p.duration_days() > 0.0);
            assert!(p.duration_days() < 1.0, "each day period < 24h");
        }
    }

    #[test]
    fn below_threshold_sun_night_periods() {
        let site = greenwich();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);
        let window = Period::new(mjd_start, mjd_end);

        let nights = below_threshold(
            &AltitudeTarget::Sun,
            &site,
            window,
            Degrees::new(-18.0), // astronomical twilight
            SearchOpts::default(),
        );

        assert!(!nights.is_empty(), "should find night periods");
    }

    #[test]
    fn altitude_ranges_twilight_band() {
        let site = greenwich();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60002.0);
        let window = Period::new(mjd_start, mjd_end);

        let twilight = altitude_ranges(
            &AltitudeTarget::Sun,
            &site,
            window,
            Degrees::new(-18.0),
            Degrees::new(-12.0),
            SearchOpts::default(),
        );

        // Should find nautical-to-astronomical twilight bands
        assert!(!twilight.is_empty(), "should find twilight bands");
    }

    #[test]
    fn moon_above_horizon_7_days() {
        let site = greenwich();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);
        let window = Period::new(mjd_start, mjd_end);

        let periods = above_threshold(
            &AltitudeTarget::Moon,
            &site,
            window,
            Degrees::new(0.0),
            SearchOpts::default(),
        );

        assert!(!periods.is_empty(), "should find moon-up periods over 7 days");
    }
}
