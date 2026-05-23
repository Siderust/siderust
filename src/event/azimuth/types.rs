// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Azimuth Type Definitions
//!
//! ## Scientific scope
//!
//! Pure data structures expressing the *result* of an azimuth analysis:
//! bearing‑crossing events (instants when *A(t)* sweeps through a fixed
//! compass bearing) and azimuth extrema (turning points of *A(t)*),
//! together with a query descriptor for range searches over the circular
//! `[0°, 360°)` domain. Wrap‑around ranges spanning North are encoded by
//! `min_azimuth > max_azimuth`; this is a convention, not a constraint
//! enforced at the type level.
//!
//! ## Technical scope
//!
//! No functions. Defines:
//! - [`AzimuthCrossingDirection`] (re‑exported from `altitude`),
//! - [`AzimuthCrossingEvent`],
//! - [`AzimuthExtremumKind`] / [`AzimuthExtremum`],
//! - [`AzimuthQuery`].
//!
//! ## References
//! None.

use crate::event::altitude::CrossingDirection;
use crate::qtty::*;
use crate::time::{ModifiedJulianDate, Period};

// Re-export CrossingDirection so consumers only need to import from this module.
pub use crate::event::altitude::CrossingDirection as AzimuthCrossingDirection;

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

impl std::fmt::Display for AzimuthCrossingEvent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Azimuth {} at {}", self.direction, self.mjd)
    }
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

impl std::fmt::Display for AzimuthExtremumKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Max => write!(f, "Max Azimuth"),
            Self::Min => write!(f, "Min Azimuth"),
        }
    }
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

impl std::fmt::Display for AzimuthExtremum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} at {} (az: {})", self.kind, self.mjd, self.azimuth)
    }
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
    pub window: Period<ModifiedJulianDate>,
    /// Lower (or start-of-wrap) bound of the azimuth band.
    pub min_azimuth: Degrees,
    /// Upper (or end-of-wrap) bound of the azimuth band.
    pub max_azimuth: Degrees,
}
