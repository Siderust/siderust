//! # Altitude Window Periods
//!
//! This module provides generic tools for finding time intervals where a celestial
//! body's altitude is within a specified range. This is useful for:
//! - Finding astronomical night periods (Sun altitude < -18°)
//! - Finding civil/nautical/astronomical twilight windows
//! - Planning observations when a target is above a minimum altitude
//!
//! ## Algorithm
//!
//! The module uses a bracket-and-refine approach:
//! 1. Coarse scan to find brackets where altitude crosses threshold boundaries
//! 2. Newton-Raphson refinement with bisection fallback for high accuracy
//! 3. Classification of crossings (entering vs exiting the altitude range)
//! 4. Pairing of crossings to form contiguous intervals
//!
//! ## Accuracy
//!
//! Timing precision is ~1 microsecond using Newton-Raphson with a tight
//! convergence criterion (`1e-11` days ≈ 0.86 µs).

use crate::astro::JulianDate;
use crate::coordinates::centers::ObserverSite;
use crate::time::{ModifiedJulianDate, Period};
use qtty::{Day, Days, Degrees, Minutes, Radian};

// =============================================================================
// Constants
// =============================================================================

/// Scan step for coarse bracket detection (10 minutes).
const SCAN_STEP: Minutes = Minutes::new(10.0);

// Root-finding constants and implementations have been moved to
// `crate::calculus::root_finding` for reuse.

// =============================================================================
// Public Types
// =============================================================================

/// Type alias for altitude periods using Modified Julian Date.
///
/// This provides backward compatibility with the previous `AltitudePeriod` struct
/// while using the generic `Period<T>` implementation.
pub type AltitudePeriod = Period<ModifiedJulianDate>;


/// Altitude condition specification for period finding.
///
/// Supports three types of conditions:
/// - `Below`: altitude < threshold (e.g., night periods)
/// - `Above`: altitude > threshold (e.g., daytime periods)
/// - `Between`: min ≤ altitude ≤ max (e.g., twilight zones)
#[derive(Debug, Clone, Copy)]
pub enum AltitudeCondition {
    /// Find periods where altitude is below a threshold (altitude < threshold).
    Below(Degrees),
    /// Find periods where altitude is above a threshold (altitude > threshold).
    Above(Degrees),
    /// Find periods where altitude is within a range (min ≤ altitude ≤ max).
    Between { min: Degrees, max: Degrees },
}

impl AltitudeCondition {
    /// Creates a condition for finding periods where altitude is BELOW a value.
    ///
    /// # Example
    /// ```ignore
    /// // Find periods where Sun is below -18° (astronomical night)
    /// let condition = AltitudeCondition::below(Degrees::new(-18.0));
    /// ```
    pub fn below(degrees: Degrees) -> Self {
        Self::Below(degrees)
    }

    /// Creates a condition for finding periods where altitude is ABOVE a value.
    ///
    /// # Example
    /// ```ignore
    /// // Find periods where a star is above 30°
    /// let condition = AltitudeCondition::above(Degrees::new(30.0));
    /// ```
    pub fn above(degrees: Degrees) -> Self {
        Self::Above(degrees)
    }

    /// Creates a condition for finding periods where altitude is within a range.
    ///
    /// # Arguments
    /// - `min`: Minimum altitude (inclusive)
    /// - `max`: Maximum altitude (inclusive)
    ///
    /// # Example
    /// ```ignore
    /// // Find astronomical night: altitude between -90° and -18°
    /// let condition = AltitudeCondition::between(
    ///     Degrees::new(-90.0),
    ///     Degrees::new(-18.0)
    /// );
    /// ```
    pub fn between(min: Degrees, max: Degrees) -> Self {
        Self::Between { min, max }
    }

    /// Checks if an altitude value (in radians) satisfies this condition.
    fn is_inside(&self, altitude_rad: f64) -> bool {
        match self {
            Self::Below(threshold) => altitude_rad < threshold.to::<Radian>().value(),
            Self::Above(threshold) => altitude_rad > threshold.to::<Radian>().value(),
            Self::Between { min, max } => {
                let alt = altitude_rad;
                let min_rad = min.to::<Radian>().value();
                let max_rad = max.to::<Radian>().value();
                alt >= min_rad && alt <= max_rad
            }
        }
    }
}


// =============================================================================
// Main API: Find Altitude Periods
// =============================================================================

/// Finds all time periods where altitude satisfies the given condition.
///
/// This is a generic function that works with any altitude evaluator closure.
/// For specific use cases like Sun altitude, use the convenience functions.
///
/// # Arguments
/// - `altitude_fn`: Closure that returns altitude in **radians** at a given JD
/// - `mjd_start`: Start of search interval (MJD)
/// - `mjd_end`: End of search interval (MJD)
/// - `condition`: The altitude condition specification (Below, Above, or Between)
///
/// # Returns
/// - `Some(Vec<AltitudePeriod>)`: Periods where condition is satisfied
/// - `None`: If the altitude never satisfies the condition
///
/// # Algorithm
/// 1. Coarse scan the interval at `SCAN_STEP` intervals
/// 2. Detect sign changes crossing the boundary/boundaries
/// 3. Refine each crossing using Newton-Raphson + bisection fallback
/// 4. Classify crossings as "entering" or "exiting" the valid range
/// 5. Pair crossings to form contiguous intervals
pub fn find_altitude_periods<F>(
    altitude_fn: F,
    mjd_start: ModifiedJulianDate,
    mjd_end: ModifiedJulianDate,
    condition: AltitudeCondition,
) -> Option<Vec<AltitudePeriod>>
where
    F: Fn(JulianDate) -> f64,
{
    let jd_start = mjd_start.to_julian_day();
    let jd_end = mjd_end.to_julian_day();

    // Collect all boundary crossings (may be 1 or 2 boundaries depending on condition)
    let boundaries = match condition {
        AltitudeCondition::Below(threshold) | AltitudeCondition::Above(threshold) => {
            vec![threshold.to::<Radian>().value()]
        }
        AltitudeCondition::Between { min, max } => {
            vec![min.to::<Radian>().value(), max.to::<Radian>().value()]
        }
    };

    let mut all_crossings: Vec<JulianDate> = Vec::new();

    // Find crossings for each boundary
    for &boundary_rad in &boundaries {
        let step: Days = SCAN_STEP.to::<Day>();
        let mut jd = jd_start;
        let mut prev_f = altitude_fn(jd) - boundary_rad;

        while jd < jd_end {
            let next_jd = (jd + step).min(jd_end);
            let next_f = altitude_fn(next_jd) - boundary_rad;

            // Sign change indicates a crossing
            if prev_f * next_f < 0.0 {
                if let Some(root) = crate::calculus::root_finding::find_crossing(
                    jd,
                    next_jd,
                    &altitude_fn,
                    boundary_rad,
                ) {
                    // Only include if within our interval
                    if root >= jd_start && root <= jd_end {
                        all_crossings.push(root);
                    }
                }
            }

            jd = next_jd;
            prev_f = next_f;
        }
    }

    // Sort crossings chronologically
    all_crossings.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Deduplicate crossings that are very close (can happen at interval boundaries)
    const DEDUPE_EPS: f64 = 1e-8; // ~1 ms
    all_crossings.dedup_by(|a, b| (a.value() - b.value()).abs() < DEDUPE_EPS);

    // Classify each crossing: +1 = entering valid range, -1 = exiting valid range
    let mut labeled: Vec<(JulianDate, i32)> = Vec::new();
    for &root in &all_crossings {
        let dt = Days::new(10.0 * crate::calculus::root_finding::FD_STEP_DAYS);
        let alt_before = altitude_fn(root - dt);
        let alt_after = altitude_fn(root + dt);

        let inside_before = condition.is_inside(alt_before);
        let inside_after = condition.is_inside(alt_after);

        if !inside_before && inside_after {
            labeled.push((root, 1)); // entering
        } else if inside_before && !inside_after {
            labeled.push((root, -1)); // exiting
        }
        // Ignore crossings that don't change inside/outside status
    }

    // Check if we start inside the valid range
    let start_altitude = altitude_fn(jd_start);
    let start_inside = condition.is_inside(start_altitude);

    // Build intervals by pairing enter/exit crossings
    let mut periods: Vec<AltitudePeriod> = Vec::new();

    if labeled.is_empty() {
        // No crossings found
        if start_inside {
            // Entire interval is valid
            return Some(vec![AltitudePeriod::new(mjd_start, mjd_end)]);
        } else {
            // Never valid
            return None;
        }
    }

    // Handle intervals
    let mut i = 0;

    // If we start inside and first crossing is an exit, add initial interval
    if start_inside && labeled[0].1 == -1 {
        let exit_mjd = ModifiedJulianDate::new(labeled[0].0.value() - 2400000.5);
        // Verify midpoint is actually inside the range
        let mid = JulianDate::new((jd_start.value() + labeled[0].0.value()) * 0.5);
        let mid_inside = condition.is_inside(altitude_fn(mid));
        if mid_inside {
            periods.push(AltitudePeriod::new(mjd_start, exit_mjd));
        }
        i = 1;
    }

    // Process remaining crossings as enter/exit pairs
    while i < labeled.len() {
        if labeled[i].1 == 1 {
            // This is an entry crossing
            let enter_jd = labeled[i].0;
            let enter_mjd = ModifiedJulianDate::new(enter_jd.value() - 2400000.5);

            // Find the corresponding exit
            let exit_mjd = if i + 1 < labeled.len() && labeled[i + 1].1 == -1 {
                let exit_jd = labeled[i + 1].0;
                i += 2;
                ModifiedJulianDate::new(exit_jd.value() - 2400000.5)
            } else {
                // No exit found, extends to end of interval
                i += 1;
                mjd_end
            };

            // Verify midpoint is actually inside the range
            let mid = JulianDate::new((enter_jd.value() + exit_mjd.to_julian_day().value()) * 0.5);
            let mid_inside = condition.is_inside(altitude_fn(mid));
            if mid_inside {
                periods.push(AltitudePeriod::new(enter_mjd, exit_mjd));
            }
        } else {
            // Unexpected exit without entry (should not happen if start_inside handled correctly)
            i += 1;
        }
    }

    if periods.is_empty() {
        None
    } else {
        Some(periods)
    }
}

// =============================================================================
// Convenience: Sun Altitude Periods
// =============================================================================

use crate::astro::sidereal::{calculate_gst, calculate_lst};
use crate::bodies::solar_system::Sun;
use qtty::AstronomicalUnit;

/// Computes the Sun's altitude in radians at a given Julian Date and observer site.
///
/// This function provides a reusable building block for any Sun-altitude-based
/// calculations (twilight, night periods, solar elevation studies, etc.).
///
/// # Arguments
/// - `jd`: Julian Date
/// - `site`: Observer site (latitude, longitude, height)
///
/// # Returns
/// Sun altitude in radians, positive above horizon, negative below.
pub fn sun_altitude_rad(jd: JulianDate, site: &ObserverSite) -> f64 {
    let sun = Sun::get_apparent_topocentric_equ::<AstronomicalUnit>(jd, *site);
    // For Topocentric positions, use azimuth() for RA and polar() for Dec
    let ra = sun.azimuth();
    let dec = sun.polar();
    
    let gst = calculate_gst(jd);
    let lst = calculate_lst(gst, site.lon);
    let ha = (lst - ra).normalize().to::<Radian>().value();
    
    let lat_rad = site.lat.to::<Radian>().value();
    let dec_rad = dec.to::<Radian>().value();
    
    (dec_rad.sin() * lat_rad.sin() + dec_rad.cos() * lat_rad.cos() * ha.cos()).asin()
}

/// Finds periods of astronomical night (Sun altitude below a threshold).
///
/// This is a convenience wrapper around [`find_altitude_periods`] specifically
/// for Sun-based calculations.
///
/// # Arguments
/// - `site`: Observer location
/// - `mjd_start`: Start of search interval (MJD)
/// - `mjd_end`: End of search interval (MJD)
/// - `threshold_deg`: Altitude threshold in degrees (e.g., -18.0 for astronomical twilight)
///
/// # Returns
/// - `Some(Vec<AltitudePeriod>)`: Night periods when Sun is below threshold
/// - `None`: Sun never goes below threshold in this interval (polar day)
///
/// # Example
/// ```ignore
/// use siderust::calculus::events::altitude_periods::*;
/// use siderust::coordinates::centers::ObserverSite;
/// use siderust::astro::ModifiedJulianDate;
/// use qtty::*;
///
/// let site = ObserverSite::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
/// let start = ModifiedJulianDate::new(60000.0);
/// let end = ModifiedJulianDate::new(60007.0);
///
/// if let Some(nights) = find_sun_below_altitude(site, start, end, Degrees::new(-18.0)) {
///     for night in nights {
///         println!("Night: MJD {} to {}", night.start.value(), night.end.value());
///     }
/// }
/// ```
pub fn find_sun_below_altitude(
    site: ObserverSite,
    mjd_start: ModifiedJulianDate,
    mjd_end: ModifiedJulianDate,
    threshold_deg: Degrees,
) -> Option<Vec<AltitudePeriod>> {
    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site);
    find_altitude_periods(altitude_fn, mjd_start, mjd_end, AltitudeCondition::below(threshold_deg))
}

/// Finds periods when Sun altitude is above a threshold.
///
/// # Arguments
/// - `site`: Observer location
/// - `mjd_start`: Start of search interval (MJD)
/// - `mjd_end`: End of search interval (MJD)
/// - `threshold_deg`: Altitude threshold in degrees (e.g., 0.0 for daytime)
///
/// # Returns
/// - `Some(Vec<AltitudePeriod>)`: Daylight periods when Sun is above threshold
/// - `None`: Sun never goes above threshold in this interval (polar night)
pub fn find_sun_above_altitude(
    site: ObserverSite,
    mjd_start: ModifiedJulianDate,
    mjd_end: ModifiedJulianDate,
    threshold_deg: Degrees,
) -> Option<Vec<AltitudePeriod>> {
    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site);
    find_altitude_periods(altitude_fn, mjd_start, mjd_end, AltitudeCondition::above(threshold_deg))
}

/// Finds periods when Sun altitude is within a specified range.
///
/// This is useful for finding specific twilight zones or other altitude windows.
///
/// # Arguments
/// - `site`: Observer location
/// - `mjd_start`: Start of search interval (MJD)
/// - `mjd_end`: End of search interval (MJD)
/// - `min_altitude`: Minimum altitude in degrees (inclusive)
/// - `max_altitude`: Maximum altitude in degrees (inclusive)
///
/// # Returns
/// - `Some(Vec<AltitudePeriod>)`: Periods when Sun is within the altitude range
/// - `None`: Sun never enters the specified range in this interval
///
/// # Example
/// ```ignore
/// // Find astronomical night periods (Sun between -90° and -18°)
/// let nights = find_sun_in_altitude_range(
///     site,
///     start,
///     end,
///     Degrees::new(-90.0),
///     Degrees::new(-18.0),
/// );
/// ```
pub fn find_sun_in_altitude_range(
    site: ObserverSite,
    mjd_start: ModifiedJulianDate,
    mjd_end: ModifiedJulianDate,
    min_altitude: Degrees,
    max_altitude: Degrees,
) -> Option<Vec<AltitudePeriod>> {
    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site);
    find_altitude_periods(
        altitude_fn,
        mjd_start,
        mjd_end,
        AltitudeCondition::between(min_altitude, max_altitude),
    )
}

// =============================================================================
// Twilight Definitions
// =============================================================================

/// Standard twilight threshold definitions.
pub mod twilight {
    use qtty::Degrees;

    /// Civil twilight: Sun center 6° below horizon (-6°)
    pub const CIVIL: Degrees = Degrees::new(-6.0);
    
    /// Nautical twilight: Sun center 12° below horizon (-12°)
    pub const NAUTICAL: Degrees = Degrees::new(-12.0);
    
    /// Astronomical twilight: Sun center 18° below horizon (-18°)
    pub const ASTRONOMICAL: Degrees = Degrees::new(-18.0);
    
    /// Sunrise/sunset: Sun center at geometric horizon (0°)
    /// Note: For apparent sunrise/sunset, use -0.833° to account for refraction
    pub const HORIZON: Degrees = Degrees::new(0.0);
    
    /// Apparent sunrise/sunset accounting for atmospheric refraction (-0.833°)
    pub const APPARENT_HORIZON: Degrees = Degrees::new(-0.833);
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::*;

    fn greenwich_site() -> ObserverSite {
        ObserverSite::new(Degrees::new(0.0), Degrees::new(51.4769), Quantity::<Meter>::new(0.0))
    }

    #[test]
    fn test_sun_altitude_basic() {
        let site = greenwich_site();
        let jd = JulianDate::J2000; // Noon-ish at J2000
        let alt = sun_altitude_rad(jd, &site);
        // At J2000 (Jan 1, 2000, 12:00 TT), Sun should be low in sky for Greenwich
        // Just verify we get a reasonable value
        assert!(alt > -std::f64::consts::FRAC_PI_2 && alt < std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn test_condition_creation() {
        let below = AltitudeCondition::below(Degrees::new(-18.0));
        assert!(matches!(below, AltitudeCondition::Below(_)));
        assert!(below.is_inside(-19.0_f64.to_radians()));
        assert!(!below.is_inside(-17.0_f64.to_radians()));

        let above = AltitudeCondition::above(Degrees::new(30.0));
        assert!(matches!(above, AltitudeCondition::Above(_)));
        assert!(above.is_inside(31.0_f64.to_radians()));
        assert!(!above.is_inside(29.0_f64.to_radians()));

        let between = AltitudeCondition::between(Degrees::new(-90.0), Degrees::new(-18.0));
        assert!(matches!(between, AltitudeCondition::Between { .. }));
        assert!(between.is_inside(-50.0_f64.to_radians()));
        assert!(between.is_inside(-18.0_f64.to_radians())); // max boundary inclusive
        assert!(between.is_inside(-90.0_f64.to_radians())); // min boundary inclusive
        assert!(!between.is_inside(-10.0_f64.to_radians()));
        assert!(!between.is_inside(-100.0_f64.to_radians()));
    }

    #[test]
    fn test_find_night_periods() {
        let site = greenwich_site();
        // Search a week around a reasonable date
        let mjd_start = ModifiedJulianDate::new(60000.0); // Around 2023
        let mjd_end = ModifiedJulianDate::new(60007.0);

        let nights = find_sun_below_altitude(site, mjd_start, mjd_end, twilight::ASTRONOMICAL);
        
        // Greenwich at 51° latitude should have astronomical night in any week
        assert!(nights.is_some(), "Should find night periods at 51° latitude");
        
        let nights = nights.unwrap();
        assert!(!nights.is_empty(), "Should have at least one night period");
        
        // Each night should have positive duration
        for night in &nights {
            assert!(night.duration_days() > 0.0, "Night duration should be positive");
            assert!(night.duration_days() < 1.0, "Night should be less than 24 hours");
        }
    }

    #[test]
    fn test_find_altitude_range_periods() {
        let site = greenwich_site();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);

        // Find astronomical night using range: [-90°, -18°]
        let nights = find_sun_in_altitude_range(
            site,
            mjd_start,
            mjd_end,
            Degrees::new(-90.0),
            Degrees::new(-18.0),
        );

        assert!(nights.is_some(), "Should find night periods using range");
        let nights = nights.unwrap();
        assert!(!nights.is_empty(), "Should have at least one night period");

        // Find nautical twilight zone: [-18°, -12°]
        let nautical = find_sun_in_altitude_range(
            site,
            mjd_start,
            mjd_end,
            Degrees::new(-18.0),
            Degrees::new(-12.0),
        );

        // Should find some nautical twilight periods
        assert!(nautical.is_some(), "Should find nautical twilight periods");
    }
}
