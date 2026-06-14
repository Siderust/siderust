// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Siderust Chebyshev Kernel (SCK) v1 reader
//!
//! Parses the small, self-describing binary format used by the
//! [`archive/`](https://github.com/Siderust/archive) submodule for fitted
//! Chebyshev ephemerides (Sun-Earth Lagrange points today; small bodies and
//! generated trajectories in the future).
//!
//! The format is documented in `archive/schema/sck-v1.md`. Briefly: a 64-byte
//! little-endian header followed by `record_count * (2 + 3 * ncoeff)` f64
//! coefficients, where each record is
//! `[mid_seconds, radius_seconds, x_c0.., y_c0.., z_c0..]`.
//!
//! This module is feature-gated under `archive-data` so the default build
//! does not link the parser unless the consumer opts in.

use crate::formats::error::{FileLocation, FormatError};

/// Magic bytes identifying a Siderust Chebyshev Kernel v1 file.
pub const SCK_MAGIC: &[u8; 8] = b"SCKERN01";

/// Fixed header length in bytes.
pub const SCK_HEADER_LEN: usize = 64;

/// SCK v1 file header.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SckHeader {
    /// Number of Chebyshev coefficients per coordinate.
    pub ncoeff: u32,
    /// Number of segments in the file.
    pub record_count: u32,
    /// Center body identifier (0 = Solar System Barycenter).
    pub center_id: u32,
    /// Target body identifier (see companion manifest).
    pub target_id: u32,
    /// Reference frame identifier (0 = EclipticMeanJ2000).
    pub frame_id: u32,
    /// Time scale identifier (0 = TDB-compatible seconds since J2000).
    pub time_scale_id: u32,
    /// First covered Julian Date in the declared time scale.
    pub valid_from_jd: f64,
    /// Last covered Julian Date in the declared time scale.
    pub valid_to_jd: f64,
}

/// Parsed SCK kernel: header plus flattened record payload (f64, host order).
#[derive(Debug, Clone, PartialEq)]
pub struct SckKernel {
    /// Parsed header.
    pub header: SckHeader,
    /// Flattened records of length `record_count * (2 + 3 * ncoeff)`.
    pub records: Vec<f64>,
}

impl SckKernel {
    /// Returns the number of f64 values per record.
    #[must_use]
    pub fn record_stride(&self) -> usize {
        2 + 3 * (self.header.ncoeff as usize)
    }
}

/// Parse a SCK v1 byte buffer.
///
/// # Errors
///
/// Returns a [`FormatError`] if the magic, header, or payload sizes are
/// inconsistent with the SCK v1 specification.
pub fn parse_sck(bytes: &[u8]) -> Result<SckKernel, FormatError> {
    if bytes.len() < SCK_HEADER_LEN {
        return Err(FormatError::located(
            "SCK v1 §header",
            FileLocation::default(),
            "buffer shorter than 64-byte header",
        ));
    }
    if &bytes[0..8] != SCK_MAGIC {
        return Err(FormatError::located(
            "SCK v1 §magic",
            FileLocation::default(),
            "magic mismatch (expected SCKERN01)",
        ));
    }
    let header = SckHeader {
        ncoeff: read_u32(&bytes[8..12]),
        record_count: read_u32(&bytes[12..16]),
        center_id: read_u32(&bytes[16..20]),
        target_id: read_u32(&bytes[20..24]),
        frame_id: read_u32(&bytes[24..28]),
        time_scale_id: read_u32(&bytes[28..32]),
        valid_from_jd: read_f64(&bytes[32..40]),
        valid_to_jd: read_f64(&bytes[40..48]),
    };
    let stride = 2 + 3 * (header.ncoeff as usize);
    let expected = stride
        .checked_mul(header.record_count as usize)
        .ok_or_else(|| {
            FormatError::located(
                "SCK v1 §payload",
                FileLocation::default(),
                "record geometry overflow",
            )
        })?;
    let payload = &bytes[SCK_HEADER_LEN..];
    if payload.len() != expected * 8 {
        return Err(FormatError::located(
            "SCK v1 §payload",
            FileLocation::default(),
            "payload size mismatch with header",
        ));
    }
    let mut records = Vec::with_capacity(expected);
    for chunk in payload.chunks_exact(8) {
        records.push(read_f64(chunk));
    }
    Ok(SckKernel { header, records })
}

fn read_u32(bytes: &[u8]) -> u32 {
    let mut out = [0u8; 4];
    out.copy_from_slice(bytes);
    u32::from_le_bytes(out)
}

fn read_f64(bytes: &[u8]) -> f64 {
    let mut out = [0u8; 8];
    out.copy_from_slice(bytes);
    f64::from_le_bytes(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_fixture() -> Vec<u8> {
        let ncoeff = 4_u32;
        let records: u32 = 1;
        let mut buf = Vec::new();
        buf.extend_from_slice(SCK_MAGIC);
        buf.extend_from_slice(&ncoeff.to_le_bytes());
        buf.extend_from_slice(&records.to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes()); // center
        buf.extend_from_slice(&3911u32.to_le_bytes()); // target
        buf.extend_from_slice(&0u32.to_le_bytes()); // frame
        buf.extend_from_slice(&0u32.to_le_bytes()); // time scale
        buf.extend_from_slice(&2451545.0_f64.to_le_bytes());
        buf.extend_from_slice(&2451546.0_f64.to_le_bytes());
        buf.extend_from_slice(&[0u8; 16]);
        // record: mid, radius, x0..3, y0..3, z0..3 (14 f64)
        for v in [
            0.5_f64, 0.5, 1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0,
        ] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf
    }

    #[test]
    fn parses_valid_fixture() {
        let bytes = build_fixture();
        let k = parse_sck(&bytes).expect("fixture must parse");
        assert_eq!(k.header.ncoeff, 4);
        assert_eq!(k.header.record_count, 1);
        assert_eq!(k.header.target_id, 3911);
        assert_eq!(k.records.len(), 14);
        assert_eq!(k.records[0], 0.5);
        assert_eq!(k.record_stride(), 14);
    }

    #[test]
    fn rejects_short_buffer() {
        assert!(parse_sck(&[0u8; 10]).is_err());
    }

    #[test]
    fn rejects_bad_magic() {
        let mut bytes = build_fixture();
        bytes[0] = b'X';
        assert!(parse_sck(&bytes).is_err());
    }

    #[test]
    fn rejects_short_payload() {
        let mut bytes = build_fixture();
        bytes.truncate(bytes.len() - 8);
        assert!(parse_sck(&bytes).is_err());
    }
}
