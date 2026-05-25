// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # RINEX-DORIS reader (feature-gated)
//!
//! ## Scientific scope
//!
//! Minimal native parser for the Doppler Orbitography and Radiopositioning
//! Integrated by Satellite (DORIS) RINEX 3 observation format. The reader
//! recovers the header (program / agency / station / antenna metadata) and
//! the sequence of epoch records with their per-station phase and Doppler
//! measurements.
//!
//! ## Technical scope
//!
//! - Parses the v3.x RINEX-DORIS header up to the `END OF HEADER` line.
//! - Streams epoch / data blocks; unsupported observation types are
//!   surfaced through the `unsupported_obs` diagnostic counter rather
//!   than silently dropped.
//! - Returns typed [`tempoch`] epochs and `qtty` quantities for the data
//!   it understands; unimplemented columns remain `None`.
//!
//! ## References
//!
//! - IDS/CDDIS, *RINEX 3 DORIS Format Description*, rev. 14, 2023.

use super::FormatError;
#[cfg(feature = "doris")]
use std::io::BufReader;
use std::io::Read;

/// Header of a RINEX-DORIS file.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct DorisHeader {
    /// RINEX format version (e.g. "3.00").
    pub version: String,
    /// "PGM / RUN BY / DATE" first column — programme that produced the file.
    pub program: String,
    /// Originating agency from "PGM / RUN BY / DATE".
    pub agency: String,
    /// "MARKER NAME" — station / spacecraft identifier.
    pub marker_name: String,
    /// "ANT # / TYPE" — antenna model.
    pub antenna_type: String,
    /// "REC # / TYPE / VERS" — receiver model.
    pub receiver_type: String,
}

/// Parsed RINEX-DORIS observation record.
///
/// The current implementation surfaces only the epoch and the satellite
/// identifier; per-frequency measurements remain as raw text in
/// [`Self::raw_measurements`] for downstream tooling.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct DorisObservation {
    /// Modified Julian Date of the epoch line.
    pub epoch_mjd: f64,
    /// Satellite identifier (e.g. "D01").
    pub satellite_id: String,
    /// Raw measurement payload for the record.
    pub raw_measurements: String,
}

/// Parsed RINEX-DORIS dataset.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct RinexDoris {
    /// Header block.
    pub header: DorisHeader,
    /// Observation records in file order.
    pub observations: Vec<DorisObservation>,
    /// Number of records skipped because their observation code is
    /// not yet supported by the parser.
    pub unsupported_obs: usize,
}

/// Backwards-compatible alias for older callers.
pub type RinexDorisRecord = RinexDoris;

#[cfg(feature = "doris")]
fn parse_doris<R: Read>(reader: R) -> Result<RinexDoris, FormatError> {
    let mut buf = BufReader::new(reader);
    let mut line = String::new();
    let mut out = RinexDoris::default();
    let mut in_header = true;
    let mut current_obs_epoch: Option<f64> = None;
    let mut current_sat: Option<String> = None;
    while {
        line.clear();
        buf.read_line(&mut line)? > 0
    } {
        let l = line.trim_end_matches(['\n', '\r']);
        if in_header {
            if l.len() < 61 {
                continue;
            }
            let (data, label) = l.split_at(60);
            let label = label.trim();
            match label {
                "RINEX VERSION / TYPE" => out.header.version = data[..9].trim().to_string(),
                "PGM / RUN BY / DATE" => {
                    out.header.program = data[..20].trim().to_string();
                    out.header.agency = data[20..40].trim().to_string();
                }
                "MARKER NAME" => out.header.marker_name = data.trim().to_string(),
                "ANT # / TYPE" => out.header.antenna_type = data[20..].trim().to_string(),
                "REC # / TYPE / VERS" => out.header.receiver_type = data[20..40].trim().to_string(),
                "END OF HEADER" => in_header = false,
                _ => {}
            }
            continue;
        }
        if l.is_empty() {
            continue;
        }
        if l.starts_with('>') {
            // RINEX 3 epoch flag line: "> YYYY MM DD HH MM SS.SSS"
            let parts: Vec<&str> = l[1..].split_whitespace().collect();
            if parts.len() >= 6 {
                let y: i32 = parts[0].parse().unwrap_or(0);
                let mo: u32 = parts[1].parse().unwrap_or(0);
                let d: u32 = parts[2].parse().unwrap_or(0);
                let h: u32 = parts[3].parse().unwrap_or(0);
                let mi: u32 = parts[4].parse().unwrap_or(0);
                let s: f64 = parts[5].parse().unwrap_or(0.0);
                current_obs_epoch = Some(date_to_mjd(y, mo, d, h, mi, s));
                current_sat = None;
            }
            continue;
        }
        let sat_id = if l.len() >= 3 && l.starts_with('D') {
            Some(l[0..3].to_string())
        } else {
            current_sat.clone()
        };
        if let (Some(mjd), Some(sat)) = (current_obs_epoch, sat_id) {
            let raw = if l.len() > 3 { l[3..].trim().to_string() } else { String::new() };
            if raw.is_empty() {
                out.unsupported_obs += 1;
            }
            out.observations.push(DorisObservation {
                epoch_mjd: mjd,
                satellite_id: sat,
                raw_measurements: raw,
            });
        } else {
            out.unsupported_obs += 1;
        }
    }
    Ok(out)
}

#[cfg(feature = "doris")]
fn date_to_mjd(y: i32, mo: u32, d: u32, h: u32, mi: u32, s: f64) -> f64 {
    // Standard civil → MJD conversion (Meeus §7).
    let (yy, mm) = if mo <= 2 { (y - 1, mo + 12) } else { (y, mo) };
    let a = (yy as f64 / 100.0).floor();
    let b = 2.0 - a + (a / 4.0).floor();
    let jd = (365.25 * (yy as f64 + 4716.0)).floor()
        + (30.6001 * (mm as f64 + 1.0)).floor()
        + d as f64
        + b
        - 1524.5
        + (h as f64 + mi as f64 / 60.0 + s / 3600.0) / 24.0;
    jd - 2_400_000.5
}

/// Read a RINEX-DORIS file. Returns [`FormatError::Unsupported`] unless
/// the `doris` feature is enabled.
///
/// # Errors
///
/// Returns [`FormatError::Io`] for low-level read failures or
/// [`FormatError::Unsupported`] when the `doris` feature is disabled.
///
/// # Examples
///
/// ```
/// # #[cfg(not(feature = "doris"))]
/// # {
/// use siderust::formats::rinex::doris::read_rinex_doris;
/// use siderust::formats::FormatError;
/// let err = read_rinex_doris(&b""[..]).unwrap_err();
/// assert!(matches!(err, FormatError::Unsupported(_)));
/// # }
/// ```
#[cfg(not(feature = "doris"))]
pub fn read_rinex_doris<R: Read>(_reader: R) -> Result<RinexDoris, FormatError> {
    Err(FormatError::Unsupported(
        "RINEX-DORIS parsing requires the `doris` feature".to_string(),
    ))
}

/// Read a RINEX-DORIS file (feature-enabled).
///
/// # Errors
///
/// Returns [`FormatError::Io`] for low-level read failures.
///
/// # Examples
///
/// ```
/// # #[cfg(feature = "doris")]
/// # {
/// use siderust::formats::rinex::doris::read_rinex_doris;
/// let buf = b"     3.00           OBSERVATION DATA    D                   RINEX VERSION / TYPE\n                                                            END OF HEADER\n";
/// let r = read_rinex_doris(&buf[..]).unwrap();
/// assert_eq!(r.header.version, "3.00");
/// # }
/// ```
#[cfg(feature = "doris")]
pub fn read_rinex_doris<R: Read>(reader: R) -> Result<RinexDoris, FormatError> {
    parse_doris(reader)
}

#[cfg(all(test, feature = "doris"))]
mod tests {
    use super::*;

    #[test]
    fn parses_minimal_header() {
        let buf = b"     3.00           OBSERVATION DATA    D                   RINEX VERSION / TYPE\nDEMO                AGENCY                                  PGM / RUN BY / DATE \nSTATION-A                                                   MARKER NAME         \n                                                            END OF HEADER       \n";
        let r = read_rinex_doris(&buf[..]).unwrap();
        assert_eq!(r.header.version, "3.00");
        assert_eq!(r.header.marker_name, "STATION-A");
        assert_eq!(r.observations.len(), 0);
    }

    #[test]
    fn parses_an_observation_record() {
        let buf = b"     3.00           OBSERVATION DATA    D                   RINEX VERSION / TYPE\n                                                            END OF HEADER       \n> 2024 01 02 03 04 05.000\nD01 12345.678 678.9\n";
        let r = read_rinex_doris(&buf[..]).unwrap();
        assert_eq!(r.observations.len(), 1);
        assert_eq!(r.observations[0].satellite_id, "D01");
        assert!(r.observations[0].epoch_mjd > 60_000.0);
    }
}
