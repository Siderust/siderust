//! # Product-crate error type
//!
//! ## Technical scope
//!
//! [`PodProductsError`] is the single unified error type returned by every
//! public writer in this crate. It wraps lower-layer errors from
//! `siderust::formats`, `serde_json`, and `std::io` under a typed enum so
//! callers can pattern-match without pulling in lower-layer error types
//! directly.
//!
//! ## References
//!
//! - Blandy, J., Orendorff, J., & Tindall, L. F. S. (2021). *Programming Rust*
//!   (2nd ed.). O'Reilly Media.

use thiserror::Error;

/// Unified error returned by all writers in `siderust-pod-products`.
///
/// # Examples
///
/// ```
/// use siderust::pod::product::PodProductsError;
///
/// let err: PodProductsError =
///     std::io::Error::new(std::io::ErrorKind::BrokenPipe, "pipe").into();
/// assert!(format!("{err}").contains("I/O"));
/// ```
#[derive(Debug, Error)]
pub enum PodProductsError {
    /// Underlying I/O failure.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// SP3 serialisation error propagated from `siderust::formats::igs::sp3`.
    #[error("SP3 error: {0}")]
    Sp3(#[from] crate::formats::igs::sp3::Sp3Error),

    /// Generic format error propagated from `siderust::formats`.
    #[error("format error: {0}")]
    PodIo(#[from] crate::formats::FormatError),

    /// JSON serialisation error.
    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),

    /// Parquet serialisation error (only constructed when feature `parquet` is active).
    #[error("Parquet error: {0}")]
    Parquet(String),
}
