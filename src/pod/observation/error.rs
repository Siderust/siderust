// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Error type for the observation-model crate.

/// Unified error type for all observation-model failures.
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::PodObservationsError;
///
/// let e = PodObservationsError::MissingSatelliteState {
///     prn: "G01".to_string(),
///     epoch_jd: 2_451_545.0,
/// };
/// assert!(e.to_string().contains("G01"));
/// ```
#[derive(Debug, thiserror::Error)]
pub enum PodObservationsError {
    /// Satellite state (position/velocity) is not available for the requested
    /// PRN and epoch.
    #[error("satellite state not available for PRN {prn} at JD {epoch_jd:.6}")]
    MissingSatelliteState {
        /// GNSS PRN or body identifier.
        prn: String,
        /// Julian date of the request.
        epoch_jd: f64,
    },

    /// Station position is not available for the requested station identifier.
    #[error("station position not available for station '{station_id}'")]
    MissingStationPosition {
        /// Station identifier (e.g. DOMES number or 4-char code).
        station_id: String,
    },

    /// Computed elevation is below the atmospheric model cutoff.
    #[error("elevation {el_deg:.1}° is below the atmospheric model cutoff")]
    ElevationTooLow {
        /// Elevation in degrees.
        el_deg: f64,
    },

    /// Light-time iteration failed to converge within the allowed iterations.
    #[error("light-time iteration did not converge")]
    LightTimeNotConverged,
}
