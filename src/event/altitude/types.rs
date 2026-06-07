// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Altitude Type Definitions
//!
//! ## Scientific scope
//!
//! Light data structures used to express the *result* of an altitude
//! computation: the kinematic event types (rising/setting threshold
//! crossings, upper/lower culminations) and the geometric query
//! description (observer, time window, altitude band). No physics is
//! performed here; the types are unit‑typed via `qtty` so callers cannot
//! accidentally mix degrees, radians, or time scales.
//!
//! ## Technical scope
//!
//! Pure data types (no functions). Defines:
//! - [`CrossingDirection`] / [`CrossingEvent`], threshold transit results.
//! - [`CulminationKind`] / [`CulminationEvent`], altitude extrema.
//!
//! ## References
//! None.

use crate::qtty::*;
use crate::time::ModifiedJulianDate;

// ---------------------------------------------------------------------------
// Crossing Types
// ---------------------------------------------------------------------------

/// Direction of a threshold crossing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossingDirection {
    /// Altitude is increasing through the threshold (object rising).
    Rising,
    /// Altitude is decreasing through the threshold (object setting).
    Setting,
}

impl std::fmt::Display for CrossingDirection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Rising => write!(f, "Rising"),
            Self::Setting => write!(f, "Setting"),
        }
    }
}

/// A threshold crossing event.
#[derive(Debug, Clone, Copy)]
pub struct CrossingEvent {
    /// Modified Julian Date of the crossing.
    pub mjd: ModifiedJulianDate,
    /// Direction: rising above or setting below the threshold.
    pub direction: CrossingDirection,
}

impl std::fmt::Display for CrossingEvent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} at {}", self.direction, self.mjd)
    }
}

// ---------------------------------------------------------------------------
// Culmination Types
// ---------------------------------------------------------------------------

/// Kind of culmination (altitude extremum).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CulminationKind {
    /// Upper culmination: local maximum altitude.
    Max,
    /// Lower culmination: local minimum altitude.
    Min,
}

impl std::fmt::Display for CulminationKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Max => write!(f, "Upper Transit"),
            Self::Min => write!(f, "Lower Transit"),
        }
    }
}

/// A culmination event.
#[derive(Debug, Clone, Copy)]
pub struct CulminationEvent {
    /// Modified Julian Date of the extremum.
    pub mjd: ModifiedJulianDate,
    /// Altitude at the extremum.
    pub altitude: Degrees,
    /// Maximum or minimum.
    pub kind: CulminationKind,
}

impl std::fmt::Display for CulminationEvent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} at {} (alt: {})", self.kind, self.mjd, self.altitude)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::Degrees;

    #[test]
    fn crossing_direction_display() {
        assert_eq!(CrossingDirection::Rising.to_string(), "Rising");
        assert_eq!(CrossingDirection::Setting.to_string(), "Setting");
    }

    #[test]
    fn crossing_event_display_includes_direction_and_mjd() {
        let event = CrossingEvent {
            mjd: ModifiedJulianDate::new(60_000.0),
            direction: CrossingDirection::Rising,
        };
        let text = event.to_string();
        assert!(text.contains("Rising"));
        assert!(text.contains("60000"));
    }

    #[test]
    fn culmination_kind_display() {
        assert_eq!(CulminationKind::Max.to_string(), "Upper Transit");
        assert_eq!(CulminationKind::Min.to_string(), "Lower Transit");
    }

    #[test]
    fn culmination_event_display_includes_kind_altitude_and_mjd() {
        let event = CulminationEvent {
            mjd: ModifiedJulianDate::new(60_001.5),
            altitude: Degrees::new(45.0),
            kind: CulminationKind::Max,
        };
        let text = event.to_string();
        assert!(text.contains("Upper Transit"));
        assert!(text.contains("45"));
        assert!(text.contains("60001.5"));
    }
}
