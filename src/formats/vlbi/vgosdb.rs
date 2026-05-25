// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # vgosDB reader (minimal native subset)
//!
//! ## Scientific scope
//!
//! The vgosDB format used by the IVS for VLBI sessions is a NetCDF-based
//! directory layout. A full reader requires a NetCDF library, which the
//! Siderust core deliberately does not pull in. This module implements
//! a *minimal text-side-car* reader that consumes the session-level
//! `Head.nc.json` (or any newline-delimited key/value side-car) that
//! vgosDB tooling commonly emits alongside the binary NetCDF files.
//!
//! ## Technical scope
//!
//! - [`VgosDbHeader`] holds the supported header fields (session
//!   code, start epoch, station list).
//! - [`read_vgosdb_sidecar`] consumes the text side-car file.
//! - The full binary NetCDF path is still unsupported and returns
//!   [`FormatError::Unsupported`] from [`read_vgosdb`].
//!
//! ## References
//!
//! - IVS, *vgosDB Format Description*, v1.1 (2020).

use super::FormatError;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Header fields recovered from a vgosDB session.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct VgosDbHeader {
    /// Session identifier (e.g. "R1100").
    pub session_code: String,
    /// Observation start epoch as ISO-8601 text.
    pub start_epoch: String,
    /// Station codes participating in the session.
    pub stations: Vec<String>,
}

/// Read a vgosDB session from a directory path.
///
/// The native binary NetCDF reader is intentionally out of scope for
/// Siderust core; callers must use the text side-car path or a
/// downstream adapter crate. Returns [`FormatError::Unsupported`] to
/// signal that explicitly.
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
/// let err = read_vgosdb(std::path::Path::new("/none")).unwrap_err();
/// assert!(matches!(err, FormatError::Unsupported(_)));
/// ```
pub fn read_vgosdb(_path: &Path) -> Result<VgosDbHeader, FormatError> {
    Err(FormatError::Unsupported(
        "vgosDB binary NetCDF reader is out of scope for siderust core; use the text side-car path"
            .to_string(),
    ))
}

/// Read a vgosDB text side-car (key=value, one per line).
///
/// Recognised keys: `session_code`, `start_epoch`, `stations`. The
/// `stations` value is comma-separated. Lines beginning with `#` are
/// ignored.
///
/// # Errors
///
/// Returns [`FormatError::Io`] for low-level read failures or
/// [`FormatError::Format`] for malformed key/value lines.
///
/// # Examples
///
/// ```
/// use siderust::formats::vlbi::vgosdb::read_vgosdb_sidecar;
/// let buf = b"session_code=R1100\nstart_epoch=2024-01-01T00:00:00\nstations=KOKEE,WETTZELL\n";
/// let h = read_vgosdb_sidecar(&buf[..]).unwrap();
/// assert_eq!(h.session_code, "R1100");
/// assert_eq!(h.stations.len(), 2);
/// ```
pub fn read_vgosdb_sidecar<R: Read>(reader: R) -> Result<VgosDbHeader, FormatError> {
    let mut buf = BufReader::new(reader);
    let mut line = String::new();
    let mut out = VgosDbHeader::default();
    let mut line_no = 0usize;
    while {
        line.clear();
        buf.read_line(&mut line)? > 0
    } {
        line_no += 1;
        let l = line.trim();
        if l.is_empty() || l.starts_with('#') {
            continue;
        }
        let (k, v) = l
            .split_once('=')
            .ok_or_else(|| FormatError::Format(format!("line {line_no}: expected key=value")))?;
        match k.trim() {
            "session_code" => out.session_code = v.trim().to_string(),
            "start_epoch" => out.start_epoch = v.trim().to_string(),
            "stations" => {
                out.stations = v
                    .split(',')
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect()
            }
            other => {
                return Err(FormatError::Format(format!(
                    "line {line_no}: unknown key `{other}`"
                )))
            }
        }
    }
    Ok(out)
}

/// Legacy placeholder alias kept for source-compatibility with the
/// pre-1.0 module surface.
pub type VgosDb = VgosDbHeader;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_vgosdb_returns_unsupported() {
        assert!(matches!(
            read_vgosdb(Path::new("/none")),
            Err(FormatError::Unsupported(_))
        ));
    }

    #[test]
    fn sidecar_parses_session_block() {
        let buf = b"# header\nsession_code=R1100\nstart_epoch=2024-01-01T00:00:00\nstations=KOKEE, WETTZELL\n";
        let h = read_vgosdb_sidecar(&buf[..]).unwrap();
        assert_eq!(h.session_code, "R1100");
        assert_eq!(h.start_epoch, "2024-01-01T00:00:00");
        assert_eq!(
            h.stations,
            vec!["KOKEE".to_string(), "WETTZELL".to_string()]
        );
    }

    #[test]
    fn sidecar_rejects_malformed_lines() {
        let buf = b"no_equals_here\n";
        assert!(matches!(
            read_vgosdb_sidecar(&buf[..]),
            Err(FormatError::Format(_))
        ));
    }

    #[test]
    fn sidecar_rejects_unknown_keys() {
        let buf = b"some_other_key=value\n";
        assert!(matches!(
            read_vgosdb_sidecar(&buf[..]),
            Err(FormatError::Format(_))
        ));
    }
}
