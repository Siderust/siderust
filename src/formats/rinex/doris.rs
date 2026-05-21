//! # RINEX-DORIS reader (feature-gated stub)
//!
//! Full RINEX-DORIS parsing is behind the `doris` feature flag.
//! Without the feature, a public type and function are exported so callers
//! can pattern-match on the error without a panic.
//!
//! ## References
//!
//! - IGS/IDS RINEX-DORIS Format Description.

use super::FormatError;
use std::io::Read;

/// A parsed RINEX-DORIS record (placeholder).
///
/// # Examples
///
/// ```
/// use siderust::formats::rinex::doris::RinexDorisRecord;
/// let r = RinexDorisRecord::default();
/// assert_eq!(r.epoch_mjd, 0.0);
/// ```
#[derive(Debug, Default)]
pub struct RinexDorisRecord {
    /// Header epoch (MJD).
    pub epoch_mjd: f64,
}

/// Read a RINEX-DORIS file.
///
/// Returns [`FormatError::Unsupported`] unless the `doris` feature is enabled.
/// Even with the feature, this implementation returns `Unsupported` until a
/// full parser is provided.
///
/// # Errors
///
/// Always returns [`FormatError::Unsupported`] in the current implementation.
///
/// # Examples
///
/// ```
/// use siderust::formats::rinex::doris::read_rinex_doris;
/// use siderust::formats::FormatError;
/// let err = read_rinex_doris(&b""[..]).unwrap_err();
/// assert!(matches!(err, FormatError::Unsupported(_)));
/// ```
#[cfg(not(feature = "doris"))]
pub fn read_rinex_doris<R: Read>(_reader: R) -> Result<RinexDorisRecord, FormatError> {
    Err(FormatError::Unsupported(
        "RINEX-DORIS parsing requires the `doris` feature".to_string(),
    ))
}

/// Read a RINEX-DORIS file (feature-enabled stub).
///
/// Returns [`FormatError::Unsupported`] until a full parser is provided.
///
/// # Errors
///
/// Always returns [`FormatError::Unsupported`] in the current implementation.
///
/// # Examples
///
/// ```
/// use siderust::formats::rinex::doris::read_rinex_doris;
/// use siderust::formats::FormatError;
/// let err = read_rinex_doris(&b""[..]).unwrap_err();
/// assert!(matches!(err, FormatError::Unsupported(_)));
/// ```
#[cfg(feature = "doris")]
pub fn read_rinex_doris<R: Read>(_reader: R) -> Result<RinexDorisRecord, FormatError> {
    Err(FormatError::Unsupported(
        "RINEX-DORIS not yet implemented even with `doris` feature".to_string(),
    ))
}
