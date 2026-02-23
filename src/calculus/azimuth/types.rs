// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Azimuth Type Definitions
//!
//! Core types for azimuth computation, bearing crossings, and azimuth extrema.

use crate::calculus::altitude::CrossingDirection;
use crate::time::{ModifiedJulianDate, Period, MJD};
use qtty::*;

// Re-export CrossingDirection so consumers only need to import from this module.
pub use crate::calculus::altitude::CrossingDirection as AzimuthCrossingDirection;

// ---------------------------------------------------------------------------
// Bearing Crossing Types
// ---------------------------------------------------------------------------

/// A bearing-crossing event: the moment a body's azimuth passes through a
/// specific compass bearing.
///
/// `direction` is [`CrossingDirection::Rising`] when the azimuth is increasing
/// (moving clockwise / eastward through the bearing), and
/// [`CrossingDirection::Setting`] when decreasing (westward).
#[derive(Debug, Clone, Copy)]
pub struct AzimuthCrossingEvent {
    /// Modified Julian Date of the crossing (TT axis).
    pub mjd: ModifiedJulianDate,
    /// Whether azimuth was increasing (Rising) or decreasing (Setting).
    pub direction: CrossingDirection,
}

// ---------------------------------------------------------------------------
// Azimuth Extremum Types
// ---------------------------------------------------------------------------

/// Kind of azimuth extremum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AzimuthExtremumKind {
    /// Local maximum azimuth (southernmost / easternmost bearing).
    Max,
    /// Local minimum azimuth (northernmost / westernmost bearing).
    Min,
}

/// An azimuth extremum event.
#[derive(Debug, Clone, Copy)]
pub struct AzimuthExtremum {
    /// Modified Julian Date of the extremum (TT axis).
    pub mjd: ModifiedJulianDate,
    /// Azimuth at the extremum (North-clockwise, `[0°, 360°)`).
    pub azimuth: Degrees,
    /// Maximum or minimum.
    pub kind: AzimuthExtremumKind,
}

// ---------------------------------------------------------------------------
// Query Type
// ---------------------------------------------------------------------------

/// Describes what to search for: observer, time window, and the azimuth
/// range of interest.
///
/// All fields use strongly-typed `qtty` quantities; time is MJD on the TT axis.
///
/// ## Wrap-around ranges
///
/// If `min_azimuth > max_azimuth` the query is interpreted as a
/// **wrap-around** (North-crossing) range: the body is "in range" when
/// `az ≥ min_azimuth  OR  az ≤ max_azimuth`.
///
/// For example `{ min: 350°, max: 10° }` matches azimuths from 350° through
/// North (0°) to 10°.
#[derive(Debug, Clone, Copy)]
pub struct AzimuthQuery {
    /// Observer location on Earth.
    pub observer: crate::coordinates::centers::Geodetic<crate::coordinates::frames::ECEF>,
    /// Time window to search (MJD on the TT axis).
    pub window: Period<MJD>,
    /// Lower (or start-of-wrap) bound of the azimuth band.
    pub min_azimuth: Degrees,
    /// Upper (or end-of-wrap) bound of the azimuth band.
    pub max_azimuth: Degrees,
}
