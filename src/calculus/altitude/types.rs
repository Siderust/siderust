// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Altitude Type Definitions
//!
//! Core types for altitude computation, crossings, and culminations.

use crate::time::{ModifiedJulianDate, Period};
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

/// A threshold crossing event.
#[derive(Debug, Clone, Copy)]
pub struct CrossingEvent {
    /// Modified Julian Date of the crossing.
    pub mjd: ModifiedJulianDate,
    /// Direction: rising above or setting below the threshold.
    pub direction: CrossingDirection,
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

/// A culmination event.
#[derive(Debug, Clone, Copy)]
pub struct CulminationEvent {
    /// Modified Julian Date of the extremum.
    pub jd: ModifiedJulianDate,
    /// Altitude at the extremum.
    pub altitude: Degrees,
    /// Maximum or minimum.
    pub kind: CulminationKind,
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
    pub observer: crate::coordinates::centers::ObserverSite,
    /// Time window to search (Modified Julian Date).
    pub window: Period<ModifiedJulianDate>,
    /// Lower bound of the altitude band (inclusive).
    pub min_altitude: Degrees,
    /// Upper bound of the altitude band (inclusive).
    pub max_altitude: Degrees,
}
