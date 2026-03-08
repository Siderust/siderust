// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Altitude Type Definitions
//!
//! Core types for altitude computation, crossings, and culminations.

use crate::time::{ModifiedJulianDate, Period, MJD};
use qtty::*;

// ---------------------------------------------------------------------------
// Crossing Types
// ---------------------------------------------------------------------------

/// Direction of a threshold crossing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossingDirection {
    Rising,
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
// Query & Period Types (from bodies/altitude.rs)
// ---------------------------------------------------------------------------

/// Describes *what* to search for: the observer, the time window, and the
/// altitude band of interest.
///
/// All fields use strongly‑typed `qtty` quantities.
#[derive(Debug, Clone, Copy)]
pub struct AltitudeQuery {
    /// Observer location on Earth.
    pub observer: crate::coordinates::centers::Geodetic<crate::coordinates::frames::ECEF>,
    /// Time window to search (MJD on the TT axis).
    pub window: Period<MJD>,
    /// Lower bound of the altitude band (inclusive).
    pub min_altitude: Degrees,
    /// Upper bound of the altitude band (inclusive).
    pub max_altitude: Degrees,
}
