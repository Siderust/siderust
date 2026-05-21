// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Core typed record produced by every parser in `siderust-tle`.

use qtty::angular::Degrees;
use qtty::angular::Turn;
use qtty::angular_rate::AngularRate;
use qtty::time::Day;
use tempoch::{Time, UTC};

use super::TleError;

/// Strongly typed NORAD catalog number with optional Alpha-5 encoding.
///
/// The wrapped `u32` is the *expanded* catalog id (≥ 100 000 for Alpha-5
/// inputs).
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::SatelliteNumber;
/// assert_eq!(SatelliteNumber::parse("25544").unwrap(), SatelliteNumber(25_544));
/// assert_eq!(SatelliteNumber::parse("T0001").unwrap(), SatelliteNumber(270_001));
/// ```
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct SatelliteNumber(pub u32);

impl SatelliteNumber {
    /// Parse a 5-character TLE catalog field, supporting Alpha-5 for ids ≥ 100 000.
    ///
    /// Alpha-5 encoding (per Celestrak): the leading character is a letter
    /// from `A` (10) … `Z` (33), excluding `I` and `O`, mapped via:
    ///
    /// ```text
    ///   value = digit_index * 10000 + remaining_digits
    /// ```
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::SatelliteNumber;
    /// assert_eq!(SatelliteNumber::parse("A0000").unwrap(), SatelliteNumber(100_000));
    /// assert!(SatelliteNumber::parse("I0001").is_err());
    /// ```
    pub fn parse(field: &str) -> Result<Self, TleError> {
        let s = field.trim_start();
        if s.is_empty() || s.len() > 5 {
            return Err(TleError::InvalidAlpha5 { raw: field.into() });
        }
        let bytes = s.as_bytes();
        let first = bytes[0];
        if first.is_ascii_digit() {
            let n: u32 = s.parse().map_err(|_| TleError::InvalidNumber {
                field: "satellite_number",
                raw: field.into(),
            })?;
            return Ok(SatelliteNumber(n));
        }
        if !first.is_ascii_uppercase() || first == b'I' || first == b'O' {
            return Err(TleError::InvalidAlpha5 { raw: field.into() });
        }
        let value = match first {
            b'A'..=b'H' => first - b'A' + 10,
            b'J'..=b'N' => first - b'J' + 18,
            b'P'..=b'Z' => first - b'P' + 23,
            _ => return Err(TleError::InvalidAlpha5 { raw: field.into() }),
        } as u32;
        let tail = &s[1..];
        if tail.len() != 4 || !tail.bytes().all(|b| b.is_ascii_digit()) {
            return Err(TleError::InvalidAlpha5 { raw: field.into() });
        }
        let tail_n: u32 = tail
            .parse()
            .map_err(|_| TleError::InvalidAlpha5 { raw: field.into() })?;
        Ok(SatelliteNumber(value * 10_000 + tail_n))
    }

    /// Format this catalog number into the 5-character TLE field, applying
    /// Alpha-5 encoding for ids ≥ 100 000.
    ///
    /// Returns an error for catalog numbers that cannot be represented in
    /// Alpha-5 (≥ 340 000).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::SatelliteNumber;
    /// assert_eq!(SatelliteNumber(25_544).format_alpha5().unwrap(), "25544");
    /// assert_eq!(SatelliteNumber(270_001).format_alpha5().unwrap(), "T0001");
    /// ```
    pub fn format_alpha5(self) -> Result<String, TleError> {
        let n = self.0;
        if n < 100_000 {
            return Ok(format!("{n:05}"));
        }
        if n >= 340_000 {
            return Err(TleError::InvalidAlpha5 {
                raw: format!("{n}"),
            });
        }
        let high = (n / 10_000) as u8;
        let tail = n % 10_000;
        let letter = match high {
            10..=17 => b'A' + (high - 10),
            18..=22 => b'J' + (high - 18),
            23..=33 => b'P' + (high - 23),
            _ => unreachable!(),
        };
        Ok(format!("{}{:04}", letter as char, tail))
    }
}

/// TLE classification character.
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::Classification;
/// assert_eq!(Classification::from_char('U').unwrap(), Classification::Unclassified);
/// ```
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Classification {
    /// Unclassified satellite (`U`).
    Unclassified,
    /// Classified satellite (`C`).
    Classified,
    /// Secret satellite (`S`).
    Secret,
}

impl Classification {
    /// Parse a TLE classification character (`U`, `C`, or `S`).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::Classification;
    /// assert!(Classification::from_char('Z').is_err());
    /// ```
    pub fn from_char(c: char) -> Result<Self, TleError> {
        match c {
            'U' => Ok(Self::Unclassified),
            'C' => Ok(Self::Classified),
            'S' => Ok(Self::Secret),
            _ => Err(TleError::InvalidClassification { raw: c }),
        }
    }

    /// Render the classification back to its TLE character.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::Classification;
    /// assert_eq!(Classification::Unclassified.as_char(), 'U');
    /// ```
    pub fn as_char(self) -> char {
        match self {
            Self::Unclassified => 'U',
            Self::Classified => 'C',
            Self::Secret => 'S',
        }
    }
}

/// COSPAR (international) designator, stored verbatim from the TLE field.
///
/// Format: 2-digit launch year, 3-digit launch number, up to 3-char piece
/// (e.g. `"98067A"`).
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::InternationalDesignator;
/// let d = InternationalDesignator("98067A".to_string());
/// assert_eq!(d.0, "98067A");
/// ```
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct InternationalDesignator(pub String);

/// Parsed TLE record with strongly-typed orbital elements.
///
/// Public fields use typed `qtty`/`tempoch` quantities wherever a typed
/// equivalent exists. The drag-related fields ([`mean_motion_dot`],
/// [`mean_motion_ddot`], [`bstar`]) and [`eccentricity`] remain `f64`:
///
// rationale: SGP4's BSTAR is a fitted, dimensionally-mixed pseudo-coefficient
// (units of 1/Earth-radii), not a clean physical drag coefficient — there is
// no canonical typed `qtty` representation. Likewise the *as-printed* time
// derivatives of mean motion and the dimensionless eccentricity have no
// typed wrapper that adds safety; we therefore expose them as plain `f64`.
//
/// [`mean_motion_dot`]: Tle::mean_motion_dot
/// [`mean_motion_ddot`]: Tle::mean_motion_ddot
/// [`bstar`]: Tle::bstar
/// [`eccentricity`]: Tle::eccentricity
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::parse_3le;
/// let tle = parse_3le(
///     "ISS (ZARYA)",
///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
/// ).unwrap();
/// assert_eq!(tle.norad_id.0, 25544);
/// ```
#[derive(Clone, Debug)]
pub struct Tle {
    /// Optional satellite name (from 3LE format, line 0).
    pub name: Option<String>,
    /// NORAD catalog number (satellite ID).
    pub norad_id: SatelliteNumber,
    /// Classification level (`U`, `C`, or `S`).
    pub classification: Classification,
    /// COSPAR (international) designator: 2-digit launch year, 3-digit launch
    /// number, piece.
    pub international_designator: InternationalDesignator,
    /// Epoch of the orbital elements (UTC).
    pub epoch: Time<UTC>,
    /// First time derivative of mean motion (revolutions per day²), as printed.
    pub mean_motion_dot: f64,
    /// Second time derivative of mean motion (revolutions per day³).
    pub mean_motion_ddot: f64,
    /// SGP4 BSTAR drag term (1 / Earth radii).
    pub bstar: f64,
    /// Element set number, incremented by the originating data source.
    pub element_set_number: u16,
    /// Number of orbits completed by the satellite at the epoch.
    pub revolution_number_at_epoch: u32,
    /// Inclination of the orbit.
    pub inclination: Degrees,
    /// Right ascension of the ascending node.
    pub raan: Degrees,
    /// Eccentricity of the orbit (dimensionless, 0 ≤ e < 1).
    pub eccentricity: f64,
    /// Argument of perigee.
    pub argument_of_perigee: Degrees,
    /// Mean anomaly.
    pub mean_anomaly: Degrees,
    /// Mean motion (revolutions per day, typed).
    pub mean_motion: AngularRate<Turn, Day>,
}
