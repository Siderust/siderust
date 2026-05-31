// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Unified error taxonomy for `siderust-tle`.
//!
//! The crate exposes exactly one error enum, [`TleError`], covering every
//! parser branch (classic 2LE/3LE, OMM-KVN, OMM-XML, OMM-JSON) as well as
//! programmatic-construction failures from [`crate::formats::tle::TleBuilder`].

/// Errors produced by every parser, writer and builder in `siderust-tle`.
///
/// All public APIs return `Result<_, TleError>`; the crate never depends on
/// `anyhow`. Variants are deliberately structural so callers can react
/// programmatically (e.g. tolerate `BadChecksum` while rejecting `BadLength`).
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::{parse_tle, TleError};
///
/// let too_short = "1 25544";
/// let other     = "2 25544";
/// match parse_tle(too_short, other) {
///     Err(TleError::BadLength { line, .. }) => assert_eq!(line, 1),
///     other => panic!("unexpected: {other:?}"),
/// }
/// ```
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum TleError {
    /// A TLE line did not have the required 69-character length.
    #[error("line {line} length is {got}, expected 69")]
    BadLength {
        /// TLE line number (1 or 2).
        line: u8,
        /// Actual length in characters.
        got: usize,
    },
    /// A TLE line's first character does not match the expected line number marker.
    #[error("line {line} starts with {found:?}, expected leading marker {expected:?}")]
    BadLeadingChar {
        /// TLE line number (1 or 2).
        line: u8,
        /// Expected marker character ('1' or '2').
        expected: char,
        /// Actual first character found.
        found: char,
    },
    /// A TLE line's checksum (last column) did not match the computed value.
    #[error("line {line} checksum mismatch: stated {stated}, computed {computed}")]
    BadChecksum {
        /// TLE line number (1 or 2).
        line: u8,
        /// Stated checksum digit from the TLE.
        stated: u8,
        /// Computed checksum value.
        computed: u8,
    },
    /// NORAD catalog number differs between line 1 and line 2 of a TLE.
    #[error("NORAD catalog id mismatch between line 1 ({l1}) and line 2 ({l2})")]
    MismatchedSatelliteNumber {
        /// NORAD ID from line 1.
        l1: u32,
        /// NORAD ID from line 2.
        l2: u32,
    },
    /// A numeric field could not be parsed as a number.
    #[error("invalid number in field {field:?}: {raw:?}")]
    InvalidNumber {
        /// Field name that failed to parse.
        field: &'static str,
        /// Raw string that could not be converted to a number.
        raw: String,
    },
    /// Alpha-5 catalog number encoding is malformed.
    #[error("invalid Alpha-5 catalog id {raw:?}")]
    InvalidAlpha5 {
        /// Raw string that failed Alpha-5 decoding.
        raw: String,
    },
    /// Epoch (year, day-of-year) does not correspond to a valid calendar date.
    #[error("invalid epoch (year {year}, doy {day_of_year}): {reason}")]
    InvalidEpoch {
        /// Year component of the epoch.
        year: i32,
        /// Day-of-year component (1.0-366.x).
        day_of_year: f64,
        /// Descriptive reason the date is invalid.
        reason: &'static str,
    },
    /// Conversion of epoch (year, day-of-year) to UTC Time failed.
    #[error("epoch conversion to UTC failed: {0}")]
    EpochConversion(String),
    /// Classification character is not one of U, C, or S.
    #[error("invalid classification {raw:?}; expected U/C/S")]
    InvalidClassification {
        /// Classification character that was not recognized.
        raw: char,
    },
    /// An OMM document was missing a mandatory field.
    #[error("OMM is missing required field {0:?}")]
    OmmMissingField(&'static str),
    /// An OMM ISO-8601 epoch could not be parsed.
    #[error("invalid OMM epoch {raw:?}: {reason}")]
    OmmInvalidEpoch {
        /// The raw epoch string.
        raw: String,
        /// Why the epoch was rejected.
        reason: &'static str,
    },
    /// An OMM document used an unsupported `MEAN_ELEMENT_THEORY` value.
    #[error("unsupported MEAN_ELEMENT_THEORY {0:?}; only \"SGP4\" is supported")]
    OmmUnsupportedTheory(String),
    /// An OMM document was syntactically malformed (KVN/XML/JSON).
    #[error("malformed OMM document ({format}): {reason}")]
    OmmMalformed {
        /// `KVN`, `XML`, or `JSON`.
        format: &'static str,
        /// Why the document was rejected.
        reason: String,
    },
    /// An OMM JSON document failed to deserialize.
    #[error("OMM JSON error: {0}")]
    OmmJson(String),
    /// A [`crate::formats::tle::TleBuilder`] is missing a mandatory field.
    #[error("TleBuilder is missing required field {0:?}")]
    BuilderMissingField(&'static str),
}
