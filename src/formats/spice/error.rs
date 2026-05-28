// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Unified error type for SPICE kernel operations.

use thiserror::Error;

/// Errors produced by SPK kernel operations.
///
/// All public APIs in this module return `Result<_, SpiceError>`. The
/// variants are kept small and structured so adapters can branch on the
/// concrete failure mode (out-of-coverage vs. format error vs. missing
/// chain) without string matching.
///
/// # Examples
///
/// ```rust
/// use siderust::formats::spice::SpiceError;
///
/// let coverage = SpiceError::OutOfCoverage {
///     target: 399,
///     center: 0,
///     epoch_tdb_seconds: 0.0,
///     start_tdb_seconds: 1.0,
///     end_tdb_seconds: 2.0,
/// };
/// assert!(format!("{coverage}").contains("out of coverage"));
///
/// let unsupported = SpiceError::UnsupportedKernelQuery {
///     message: "DSK surface queries are not implemented".to_string(),
/// };
/// assert!(format!("{unsupported}").contains("unsupported kernel query"));
/// ```
#[derive(Debug, Error)]
pub enum SpiceError {
    /// I/O error reading a kernel file from disk.
    #[error("SPICE kernel I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// The underlying DAF/SPK parser rejected the file.
    #[error("SPICE kernel parse error: {message}")]
    Parse {
        /// Human-readable description from the upstream parser.
        message: String,
    },

    /// No segment in the kernel covers `target → center` at `epoch`.
    #[error(
        "SPICE kernel: state for target={target} center={center} is out of coverage \
         at et={epoch_tdb_seconds} TDB-s; available segment range is \
         [{start_tdb_seconds}, {end_tdb_seconds}] TDB-s"
    )]
    OutOfCoverage {
        /// NAIF body ID of the requested target.
        target: i32,
        /// NAIF body ID of the requested center.
        center: i32,
        /// Epoch in TDB seconds past J2000.
        epoch_tdb_seconds: f64,
        /// Start of the closest segment's coverage window (TDB seconds).
        start_tdb_seconds: f64,
        /// End of the closest segment's coverage window (TDB seconds).
        end_tdb_seconds: f64,
    },

    /// Kernel contains no segment chain from `target` to `center`.
    #[error(
        "SPICE kernel: no chain from target={target} to center={center} \
         (kernel has no segment with this body pair, even via transitive \
         centers)"
    )]
    NoChain {
        /// NAIF body ID of the requested target.
        target: i32,
        /// NAIF body ID of the requested center.
        center: i32,
    },

    /// SPK data type is not implemented.
    ///
    /// This crate ships SPK Type 2 (Chebyshev position), Type 3 (Chebyshev
    /// position+velocity), Type 9 (unequally-spaced Lagrange state
    /// interpolation), and Type 13 (unequally-spaced Hermite state
    /// interpolation). Higher Types are indexed but not evaluated.
    #[error("SPICE kernel: SPK Type {data_type} is not implemented")]
    UnsupportedDataType {
        /// NAIF SPK data type code.
        data_type: i32,
    },

    /// Internal invariant violated while decoding a segment.
    ///
    /// Returned when a record's metadata (e.g. `radius`) is non-finite or
    /// non-positive, which would indicate kernel corruption.
    #[error("SPICE kernel: segment record corrupted: {message}")]
    Corrupted {
        /// Human-readable description.
        message: String,
    },

    /// A high-level kernel query is not implemented for this kernel type.
    #[error("SPICE kernel: unsupported kernel query: {message}")]
    UnsupportedKernelQuery {
        /// Human-readable explanation.
        message: String,
    },

    /// The requested frame ID or name is not registered in the current context.
    #[error("SPICE context: unknown frame: {description}")]
    UnknownFrame {
        /// Description of the unrecognised frame.
        description: String,
    },

    /// A time conversion failed.
    #[error("SPICE context: time conversion error: {message}")]
    TimeConversion {
        /// Human-readable reason.
        message: String,
    },

    /// A low-level DAF/SPK format parse error.
    #[error("SPICE format parse error: {0}")]
    FormatParse(String),
}

impl From<crate::datasets::DatasetError> for SpiceError {
    fn from(err: crate::datasets::DatasetError) -> Self {
        match err {
            crate::datasets::DatasetError::Io(e) => SpiceError::Io(e),
            other => SpiceError::FormatParse(format!("{other}")),
        }
    }
}
