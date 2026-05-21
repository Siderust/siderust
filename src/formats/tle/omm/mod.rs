// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! CCSDS Orbit Mean-elements Message (OMM) parsing/writing.
//!
//! Supports the three encodings emitted by Celestrak and other catalog
//! providers:
//!
//! * [`kvn`] — keyword-value notation (CCSDS 502.0-B-2 §4.1).
//! * [`xml`] — XML Schema for OMM (CCSDS 502.0-B-2 Annex E).
//! * [`json`] — Celestrak's JSON encoding (`OBJECT_NAME`, `MEAN_MOTION`, …).
//!   Requires the `serde` feature.
//!
//! All three round-trip through the same in-memory [`Omm`] record. Use
//! [`Omm::from_tle`] / [`Omm::to_tle`] to convert from/to the classic
//! 2LE/3LE [`crate::formats::tle::Tle`] representation.

#[cfg(feature = "serde")]
pub mod json;
pub mod kvn;
pub mod xml;

use chrono::{Datelike, TimeZone, Timelike, Utc};
use qtty::angular::Degrees;
use qtty::angular::Turn;
use qtty::angular_rate::AngularRate;
use qtty::time::Day;
use tempoch::{Time, UTC};

use crate::formats::tle::TleError;
use crate::formats::tle::{Classification, InternationalDesignator, SatelliteNumber, Tle};

/// Strongly-typed OMM record.
///
/// Mirrors the CCSDS 502.0 OMM/SGP4 schema. Fields not relevant to SGP4
/// propagation (`COMMENT`, `CENTER_NAME`, etc.) are preserved as
/// `metadata` for round-trip round-tripping but are otherwise opaque.
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::{parse_3le, omm::Omm};
/// let tle = parse_3le(
///     "ISS (ZARYA)",
///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
/// ).unwrap();
/// let omm = Omm::from_tle(&tle);
/// assert_eq!(omm.norad_id.0, 25544);
/// ```
#[derive(Clone, Debug)]
pub struct Omm {
    /// `OBJECT_NAME` — satellite name (free-form).
    pub object_name: String,
    /// `OBJECT_ID` — COSPAR international designator (e.g. `"1998-067A"`).
    pub object_id: String,
    /// Originator-supplied epoch (UTC).
    pub epoch: Time<UTC>,
    /// Mean motion (revolutions per day, typed).
    pub mean_motion: AngularRate<Turn, Day>,
    /// Eccentricity (dimensionless, 0 ≤ e < 1).
    pub eccentricity: f64,
    /// Inclination.
    pub inclination: Degrees,
    /// Right ascension of the ascending node.
    pub ra_of_asc_node: Degrees,
    /// Argument of perigee.
    pub arg_of_pericenter: Degrees,
    /// Mean anomaly.
    pub mean_anomaly: Degrees,
    /// `EPHEMERIS_TYPE` — usually `0` for SGP4.
    pub ephemeris_type: u8,
    /// `CLASSIFICATION_TYPE` — `U`/`C`/`S`.
    pub classification: Classification,
    /// `NORAD_CAT_ID` — catalog id (Alpha-5-decoded).
    pub norad_id: SatelliteNumber,
    /// `ELEMENT_SET_NO`.
    pub element_set_no: u16,
    /// `REV_AT_EPOCH`.
    pub rev_at_epoch: u32,
    /// `BSTAR` drag term (1 / Earth radii).
    pub bstar: f64,
    /// `MEAN_MOTION_DOT` — first time derivative (revs/day²).
    pub mean_motion_dot: f64,
    /// `MEAN_MOTION_DDOT` — second time derivative (revs/day³).
    pub mean_motion_ddot: f64,
}

impl Omm {
    /// Build an OMM record from a parsed [`Tle`].
    ///
    /// `OBJECT_NAME` falls back to the catalog id rendering when the TLE
    /// has no name (i.e. parsed via [`crate::formats::tle::parse_tle`] rather than
    /// [`crate::formats::tle::parse_3le`]). `OBJECT_ID` is reconstructed from the
    /// 6-character TLE international designator into the canonical OMM
    /// `YYYY-NNNP` shape (with the 2-digit launch year expanded using the
    /// same 1957/2000 cutover as the TLE epoch).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::{parse_3le, omm::Omm};
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let omm = Omm::from_tle(&tle);
    /// assert_eq!(omm.object_id, "1998-067A");
    /// ```
    pub fn from_tle(tle: &Tle) -> Self {
        let object_name = tle
            .name
            .clone()
            .unwrap_or_else(|| format!("NORAD {}", tle.norad_id.0));
        let object_id = expand_intl_designator(&tle.international_designator.0);
        Self {
            object_name,
            object_id,
            epoch: tle.epoch,
            mean_motion: tle.mean_motion,
            eccentricity: tle.eccentricity,
            inclination: tle.inclination,
            ra_of_asc_node: tle.raan,
            arg_of_pericenter: tle.argument_of_perigee,
            mean_anomaly: tle.mean_anomaly,
            ephemeris_type: 0,
            classification: tle.classification,
            norad_id: tle.norad_id,
            element_set_no: tle.element_set_number,
            rev_at_epoch: tle.revolution_number_at_epoch,
            bstar: tle.bstar,
            mean_motion_dot: tle.mean_motion_dot,
            mean_motion_ddot: tle.mean_motion_ddot,
        }
    }

    /// Convert the OMM into a classic [`Tle`].
    ///
    /// `name` is set from `OBJECT_NAME`, and the international designator
    /// is contracted from `YYYY-NNNP` back into the 6/7-character TLE
    /// field.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::{parse_3le, omm::Omm};
    /// let l1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927";
    /// let l2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
    /// let tle1 = parse_3le("ISS (ZARYA)", l1, l2).unwrap();
    /// let tle2 = Omm::from_tle(&tle1).to_tle();
    /// assert_eq!(tle2.norad_id, tle1.norad_id);
    /// ```
    pub fn to_tle(&self) -> Tle {
        Tle {
            name: Some(self.object_name.clone()),
            norad_id: self.norad_id,
            classification: self.classification,
            international_designator: InternationalDesignator(contract_intl_designator(
                &self.object_id,
            )),
            epoch: self.epoch,
            mean_motion_dot: self.mean_motion_dot,
            mean_motion_ddot: self.mean_motion_ddot,
            bstar: self.bstar,
            element_set_number: self.element_set_no,
            revolution_number_at_epoch: self.rev_at_epoch,
            inclination: self.inclination,
            raan: self.ra_of_asc_node,
            eccentricity: self.eccentricity,
            argument_of_perigee: self.arg_of_pericenter,
            mean_anomaly: self.mean_anomaly,
            mean_motion: self.mean_motion,
        }
    }
}

/// Convert a TLE 6/8-character international designator (e.g. `"98067A"`)
/// into the OMM `YYYY-NNNP` canonical form (e.g. `"1998-067A"`).
pub(crate) fn expand_intl_designator(field: &str) -> String {
    let s = field.trim();
    if s.len() < 5 || !s.is_ascii() {
        return s.to_string();
    }
    let bytes = s.as_bytes();
    if !bytes[0].is_ascii_digit() || !bytes[1].is_ascii_digit() {
        return s.to_string();
    }
    let yy: i32 = s[..2].parse().unwrap_or(0);
    let yyyy = if yy >= 57 { 1900 + yy } else { 2000 + yy };
    let launch = &s[2..5];
    let piece = &s[5..];
    if piece.is_empty() {
        format!("{yyyy:04}-{launch}")
    } else {
        format!("{yyyy:04}-{launch}{piece}")
    }
}

/// Inverse of [`expand_intl_designator`]: contract `YYYY-NNNP` back to
/// `YYNNNP`. Tolerant of input that is already in TLE format.
pub(crate) fn contract_intl_designator(s: &str) -> String {
    let s = s.trim();
    if let Some(dash) = s.find('-') {
        let (yyyy, rest) = s.split_at(dash);
        let rest = &rest[1..];
        if yyyy.len() == 4 && yyyy.bytes().all(|b| b.is_ascii_digit()) {
            let yy = &yyyy[2..];
            return format!("{yy}{rest}");
        }
    }
    s.to_string()
}

/// Format a [`Time<UTC>`] as the ISO-8601 fractional-second OMM epoch.
pub(crate) fn format_epoch(t: Time<UTC>) -> Result<String, TleError> {
    let dt = t
        .try_to_chrono()
        .map_err(|e| TleError::EpochConversion(format!("{e:?}")))?;
    let micros = dt.timestamp_subsec_micros();
    Ok(format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}.{:06}",
        dt.year(),
        dt.month(),
        dt.day(),
        dt.hour(),
        dt.minute(),
        dt.second(),
        micros
    ))
}

/// Parse an ISO-8601 OMM epoch (with optional fractional seconds and
/// optional `Z`/`+00:00` zone) into a UTC [`Time`].
pub(crate) fn parse_epoch(raw: &str) -> Result<Time<UTC>, TleError> {
    let s = raw.trim();
    let trimmed = s.trim_end_matches('Z');
    let trimmed = trimmed.split('+').next().unwrap_or(trimmed);
    let dt = chrono::NaiveDateTime::parse_from_str(trimmed, "%Y-%m-%dT%H:%M:%S%.f")
        .or_else(|_| chrono::NaiveDateTime::parse_from_str(trimmed, "%Y-%m-%dT%H:%M:%S"))
        .map_err(|e| TleError::OmmInvalidEpoch {
            raw: raw.to_string(),
            reason: match e.kind() {
                chrono::format::ParseErrorKind::Invalid => "invalid component",
                chrono::format::ParseErrorKind::OutOfRange => "out of range",
                _ => "malformed ISO-8601",
            },
        })?;
    let utc_dt = Utc.from_utc_datetime(&dt);
    Time::<UTC>::try_from_chrono(utc_dt).map_err(|e| TleError::EpochConversion(format!("{e:?}")))
}
