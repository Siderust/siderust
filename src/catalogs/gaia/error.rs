// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Error types for Gaia DR3 catalogue ingestion and photometry.

/// Result alias for Gaia DR3 catalogue operations.
pub type Result<T> = std::result::Result<T, GaiaDr3Error>;

/// Errors produced while validating Gaia DR3 raw rows or spectra.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum GaiaDr3Error {
    /// A required raw catalogue field was absent.
    #[error("required Gaia DR3 field `{field}` is missing")]
    MissingRequiredField {
        /// Field name.
        field: &'static str,
    },
    /// A scalar field was not finite.
    #[error("Gaia DR3 field `{field}` must be finite, got {value}")]
    NonFinite {
        /// Field name.
        field: &'static str,
        /// Invalid value.
        value: f64,
    },
    /// Right ascension was outside the Gaia catalogue range.
    #[error("Gaia DR3 right ascension must be in [0, 360) deg, got {value}")]
    RightAscensionOutOfRange {
        /// Invalid right ascension in degrees.
        value: f64,
    },
    /// Declination was outside the Gaia catalogue range.
    #[error("Gaia DR3 declination must be in [-90, 90] deg, got {value}")]
    DeclinationOutOfRange {
        /// Invalid declination in degrees.
        value: f64,
    },
    /// A strictly positive quantity was not positive.
    #[error("Gaia DR3 field `{field}` must be positive, got {value}")]
    NotPositive {
        /// Field name.
        field: &'static str,
        /// Invalid value.
        value: f64,
    },
    /// A non-negative quantity was negative.
    #[error("Gaia DR3 field `{field}` must be non-negative, got {value}")]
    Negative {
        /// Field name.
        field: &'static str,
        /// Invalid value.
        value: f64,
    },
    /// Spectrum samples were empty.
    #[error("Gaia XP sampled spectrum must contain at least one sample")]
    EmptySpectrum,
    /// Spectrum wavelengths were not strictly increasing.
    #[error("Gaia XP sampled spectrum wavelengths must be strictly increasing")]
    NonMonotonicSpectrum,
    /// The requested passband is malformed.
    #[error("spectral band `{name}` has invalid bounds")]
    InvalidSpectralBand {
        /// Band name.
        name: &'static str,
    },
    /// The spectrum does not overlap the requested band with enough samples.
    #[error("Gaia XP sampled spectrum does not cover spectral band `{name}`")]
    InsufficientBandCoverage {
        /// Band name.
        name: &'static str,
    },
    /// A CSV row had a missing required column.
    #[error("required Gaia DR3 column `{column}` missing on row at line {line}")]
    MissingColumn {
        /// Column name.
        column: &'static str,
        /// 1-based line number.
        line: usize,
    },
    /// A CSV value could not be parsed.
    #[error("invalid Gaia DR3 number for column `{column}` on line {line}: {raw}")]
    InvalidNumber {
        /// Column name.
        column: &'static str,
        /// 1-based line number.
        line: usize,
        /// Raw value.
        raw: String,
    },
}
