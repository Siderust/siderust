// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Shared error types for [`crate::formats`].

use std::path::PathBuf;
use thiserror::Error;

/// Location inside an input artefact, used by structured diagnostics.
///
/// Fields are 1-based when present (matching the convention of common text
/// editors and the format specifications themselves). `path` is `None` when
/// the input was a buffer rather than a named file.
///
/// # Examples
///
/// ```
/// use siderust::formats::FileLocation;
/// let loc = FileLocation::new(Some("/data/igs.sp3".into()), Some(42), Some(7));
/// assert_eq!(loc.line, Some(42));
/// assert_eq!(loc.column, Some(7));
/// assert!(loc.path.as_ref().unwrap().ends_with("igs.sp3"));
/// ```
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FileLocation {
    /// Filesystem path, when known.
    pub path: Option<PathBuf>,
    /// 1-based line number, when known.
    pub line: Option<usize>,
    /// 1-based column number, when known.
    pub column: Option<usize>,
}

impl FileLocation {
    /// Build a [`FileLocation`] from optional fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::FileLocation;
    /// let loc = FileLocation::new(None, Some(1), None);
    /// assert!(loc.path.is_none());
    /// ```
    pub fn new(path: Option<PathBuf>, line: Option<usize>, column: Option<usize>) -> Self {
        Self { path, line, column }
    }

    /// Build a buffer-anchored location with just a line number.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::FileLocation;
    /// let loc = FileLocation::at_line(17);
    /// assert_eq!(loc.line, Some(17));
    /// ```
    pub fn at_line(line: usize) -> Self {
        Self {
            path: None,
            line: Some(line),
            column: None,
        }
    }
}

impl std::fmt::Display for FileLocation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match (&self.path, self.line, self.column) {
            (Some(p), Some(l), Some(c)) => write!(f, "{}:{}:{}", p.display(), l, c),
            (Some(p), Some(l), None) => write!(f, "{}:{}", p.display(), l),
            (Some(p), None, _) => write!(f, "{}", p.display()),
            (None, Some(l), Some(c)) => write!(f, "<input>:{}:{}", l, c),
            (None, Some(l), None) => write!(f, "<input>:{}", l),
            (None, None, _) => write!(f, "<input>"),
        }
    }
}

/// Parse mode for permissive vs. strict ingestion.
///
/// In `Strict` mode, any deviation from the published format spec is a
/// hard error. In `Permissive` mode, recoverable deviations (extra
/// whitespace, unknown header records, truncated trailing blocks) are
/// silently skipped while still yielding a deterministic result.
///
/// Both modes are deterministic — `Permissive` is **not** a synonym for
/// "best-effort guess".
///
/// # Examples
///
/// ```
/// use siderust::formats::ParseMode;
/// assert_eq!(ParseMode::default(), ParseMode::Strict);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ParseMode {
    /// Reject any spec violation.
    #[default]
    Strict,
    /// Recover from non-fatal deviations.
    Permissive,
}

/// Unified format-layer error type for text and binary parsers.
///
/// `Format` retains its free-form string variant for compatibility, while
/// `Located` is the structured variant carrying a [`FileLocation`] and a
/// per-format spec section reference (e.g. `"SP3-d §3.2"`).
///
/// # Examples
///
/// ```
/// use siderust::formats::{FileLocation, FormatError};
/// let err = FormatError::located(
///     "RINEX 3.05 §6.3",
///     FileLocation::at_line(42),
///     "epoch line malformed",
/// );
/// assert!(format!("{err}").contains("RINEX"));
/// ```
#[derive(Debug, Error)]
pub enum FormatError {
    /// Underlying IO failure.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    /// Malformed or unsupported file content (free-form).
    #[error("format: {0}")]
    Format(String),
    /// Format or feature not yet supported.
    #[error("unsupported: {0}")]
    Unsupported(String),
    /// Structured format error with file location and spec section.
    #[error("{spec} at {location}: {message}")]
    Located {
        /// Reference to the format-spec section being violated.
        spec: &'static str,
        /// Where the violation was observed.
        location: FileLocation,
        /// Human-readable diagnostic.
        message: String,
    },
}

impl FormatError {
    /// Construct a structured [`FormatError::Located`].
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::{FileLocation, FormatError};
    /// let _ = FormatError::located("OEM v3 §5.2", FileLocation::at_line(3), "missing OBJECT_ID");
    /// ```
    pub fn located<M: Into<String>>(
        spec: &'static str,
        location: FileLocation,
        message: M,
    ) -> Self {
        Self::Located {
            spec,
            location,
            message: message.into(),
        }
    }
}
