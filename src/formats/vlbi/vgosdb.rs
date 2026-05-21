//! # vgosDB reader stub
//!
//! The vgosDB (VLBI Global Observing System Data Base) format uses a
//! NetCDF-based directory structure. Full parsing requires a NetCDF library
//! and is not yet implemented.
//!
//! This stub exports a public type and a function that returns
//! [`FormatError::Unsupported`], allowing callers to pattern-match on the
//! error without a panic.

use super::FormatError;

/// Placeholder type for a vgosDB dataset.
///
/// # Examples
///
/// ```
/// use siderust::formats::vlbi::vgosdb::VgosDb;
/// // VgosDb cannot be constructed until the feature is implemented.
/// let err = siderust::formats::vlbi::vgosdb::read_vgosdb(std::path::Path::new("/none")).unwrap_err();
/// assert!(matches!(err, siderust::formats::FormatError::Unsupported(_)));
/// ```
#[derive(Debug)]
pub struct VgosDb;

/// Read a vgosDB dataset.
///
/// Always returns [`FormatError::Unsupported`] in the current implementation.
///
/// # Errors
///
/// Always returns [`FormatError::Unsupported`].
///
/// # Examples
///
/// ```
/// use siderust::formats::vlbi::vgosdb::read_vgosdb;
/// use siderust::formats::FormatError;
/// use std::path::Path;
///
/// let err = read_vgosdb(Path::new("/nonexistent")).unwrap_err();
/// assert!(matches!(err, FormatError::Unsupported(_)));
/// ```
pub fn read_vgosdb(_path: &std::path::Path) -> Result<VgosDb, FormatError> {
    Err(FormatError::Unsupported(
        "vgosDB not yet implemented".to_string(),
    ))
}
