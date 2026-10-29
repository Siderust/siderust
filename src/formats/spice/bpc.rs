// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Binary PCK Type 2 (Chebyshev Euler angles), without CSPICE runtime state.
//! See NAIF PCK Required Reading and `PCKE02`/`PCKMAT`.

use super::{daf::DafRaw, SpiceError};
#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Scalar, Transcendental};
use alloc::{format, vec::Vec};

/// A base-frame ID and base-to-fixed rotation acting on column vectors.
pub type BinaryPckRotation = (i32, [[f64; 3]; 3]);

/// Indexed binary PCK, with last-segment priority and inclusive coverage.
#[derive(Debug)]
pub struct BinaryPck {
    daf: DafRaw,
    segments: Vec<Segment>,
}

#[derive(Debug)]
struct Segment {
    class_id: i32,
    base: i32,
    kind: i32,
    start: f64,
    end: f64,
    word: usize,
    init: f64,
    interval: f64,
    size: usize,
    count: usize,
}

fn invalid(message: &str) -> SpiceError {
    SpiceError::Corrupted {
        message: message.into(),
    }
}

impl BinaryPck {
    /// Parse and validate a DAF/PCK file. Unsupported types are indexed and
    /// fail explicitly when selected; they never expose lower-priority data.
    pub fn from_bytes(bytes: Vec<u8>) -> Result<Self, SpiceError> {
        if !bytes.len().is_multiple_of(1024) {
            return Err(invalid("binary PCK must contain complete DAF records"));
        }
        if bytes.len() >= 96 && !matches!(&bytes[88..96], b"LTL-IEEE" | b"BIG-IEEE") {
            return Err(invalid("binary PCK requires IEEE DAF binary format marker"));
        }
        let daf = DafRaw::parse(bytes)?;
        if daf.locator.trim() != "DAF/PCK" || daf.nd != 2 || daf.ni != 5 {
            return Err(invalid("binary PCK requires DAF/PCK, ND=2, NI=5"));
        }
        let mut segments = Vec::new();
        for s in &daf.raw_summaries {
            let [start, end] = [s.doubles[0], s.doubles[1]];
            if !start.is_finite() || !end.is_finite() || start > end {
                return Err(invalid("invalid binary PCK coverage"));
            }
            let word = s.integers[3];
            let last = s.integers[4];
            if word <= 0 || last < word || last as usize > daf.bytes().len() / 8 {
                return Err(invalid("binary PCK descriptor addresses outside file"));
            }
            let mut segment = Segment {
                class_id: s.integers[0],
                base: s.integers[1],
                kind: s.integers[2],
                start,
                end,
                word: word as usize,
                init: 0.0,
                interval: 0.0,
                size: 0,
                count: 0,
            };
            if segment.kind == 2 {
                let len = (last - word + 1) as usize;
                if len < 9 {
                    return Err(invalid("truncated Type 2 PCK"));
                }
                let tail = last as usize - 3;
                segment.init = daf.read_f64_at_word(tail);
                segment.interval = daf.read_f64_at_word(tail + 1);
                let size = daf.read_f64_at_word(tail + 2);
                let count = daf.read_f64_at_word(tail + 3);
                if !segment.init.is_finite()
                    || !segment.interval.is_finite()
                    || segment.interval <= 0.0
                    || !size.is_finite()
                    || size < 5.0
                    || size.fract() != 0.0
                    || size > (len - 4) as f64
                    || !count.is_finite()
                    || count < 1.0
                    || count.fract() != 0.0
                    || count > (len - 4) as f64
                {
                    return Err(invalid("invalid Type 2 PCK directory"));
                }
                segment.size = size as usize;
                segment.count = count as usize;
                if !(segment.size - 2).is_multiple_of(3)
                    || segment.size.checked_mul(segment.count) != Some(len - 4)
                {
                    return Err(invalid("Type 2 PCK record sizes do not match segment"));
                }
                let data_end = segment.init + segment.interval * count;
                if !data_end.is_finite() || start < segment.init || end > data_end {
                    return Err(invalid("PCK coverage exceeds record intervals"));
                }
                for i in 0..segment.count {
                    let rec = segment.word + i * segment.size;
                    let mid = daf.read_f64_at_word(rec);
                    let radius = daf.read_f64_at_word(rec + 1);
                    let expected = segment.init + (i as f64 + 0.5) * segment.interval;
                    let tol = 16.0 * f64::EPSILON * expected.abs().max(segment.interval);
                    if !mid.is_finite()
                        || !radius.is_finite()
                        || radius <= 0.0
                        || (mid - expected).abs() > tol
                        || (radius - segment.interval / 2.0).abs() > tol
                    {
                        return Err(invalid(
                            "PCK record midpoint/radius inconsistent with directory",
                        ));
                    }
                    for j in 2..segment.size {
                        if !daf.read_f64_at_word(rec + j).is_finite() {
                            return Err(invalid("non-finite PCK coefficient"));
                        }
                    }
                }
            }
            segments.push(segment);
        }
        Ok(Self { daf, segments })
    }

    /// Whether this file contains a frame class ID, independent of epoch.
    pub fn contains_class(&self, class_id: i32) -> bool {
        self.segments.iter().any(|s| s.class_id == class_id)
    }

    /// Overall descriptor coverage for a frame class, including gaps.
    pub fn coverage(&self, class_id: i32) -> Option<(f64, f64)> {
        self.segments
            .iter()
            .filter(|s| s.class_id == class_id)
            .map(|s| (s.start, s.end))
            .reduce(|a, b| (a.0.min(b.0), a.1.max(b.1)))
    }

    /// Whether all segments for a class use J2000 as their inertial base.
    pub fn has_j2000_base(&self, class_id: i32) -> bool {
        self.contains_class(class_id)
            && self
                .segments
                .iter()
                .filter(|s| s.class_id == class_id)
                .all(|s| s.base == 1)
    }

    /// Evaluate the base-to-fixed rotation and return its base frame ID.
    /// `Ok(None)` means no covering segment; non-finite epochs are errors.
    pub fn rotation(
        &self,
        class_id: i32,
        et: f64,
    ) -> Result<Option<BinaryPckRotation>, SpiceError> {
        if !et.is_finite() {
            return Err(invalid("PCK epoch must be finite TDB seconds"));
        }
        let Some(s) = self
            .segments
            .iter()
            .rev()
            .find(|s| s.class_id == class_id && et >= s.start && et <= s.end)
        else {
            return Ok(None);
        };
        if s.kind != 2 {
            return Err(SpiceError::UnsupportedKernelQuery {
                message: format!("binary PCK Type {} for class {class_id}", s.kind),
            });
        }
        let index = (((et - s.init) / s.interval).floor() as usize).min(s.count - 1);
        let rec = s.word + index * s.size;
        let mid = self.daf.read_f64_at_word(rec);
        let radius = self.daf.read_f64_at_word(rec + 1);
        let x = (et - mid) / radius;
        if x.abs() > 1.0 + 32.0 * f64::EPSILON {
            return Err(invalid("PCK interpolation outside record"));
        }
        let n = (s.size - 2) / 3;
        let mut angles = [0.0; 3];
        for (axis, angle) in angles.iter_mut().enumerate() {
            let (mut b1, mut b2) = (0.0, 0.0);
            for j in (1..n).rev() {
                let b = 2.0 * x * b1 - b2 + self.daf.read_f64_at_word(rec + 2 + axis * n + j);
                b2 = b1;
                b1 = b;
            }
            *angle = x * b1 - b2 + self.daf.read_f64_at_word(rec + 2 + axis * n);
        }
        // PCKE02 reduces W using signed remainder before constructing the
        // rotation. This matters at sub-picoradian precision for large W.
        let turns = (angles[2] / core::f64::consts::TAU).trunc();
        // Match NAIF's f2c D_MOD subtraction, rather than libm fmod.
        angles[2] -= core::f64::consts::TAU * turns;
        if angles.iter().any(|v| !v.is_finite()) {
            return Err(invalid("PCK angle evaluation overflow"));
        }
        // Binary PCK stores the actual 3-1-3 Euler angles (not text-PCK RA/DEC).
        Ok(Some((
            s.base,
            mul(
                axis_rotation(3, angles[2]),
                mul(axis_rotation(1, angles[1]), axis_rotation(3, angles[0])),
            ),
        )))
    }
}

pub(crate) fn mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                out[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    out
}

// SPICE coordinate rotation, positive angle rotates coordinate axes.
pub(crate) fn axis_rotation(axis: i32, angle: f64) -> [[f64; 3]; 3] {
    let (s, c) = angle.sin_cos();
    match axis {
        1 => [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]],
        2 => [[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]],
        _ => [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn kernel(big: bool) -> Vec<u8> {
        let mut b = vec![0; 4096];
        b[..8].copy_from_slice(b"DAF/PCK ");
        b[88..96].copy_from_slice(if big { b"BIG-IEEE" } else { b"LTL-IEEE" });
        let put_i = |b: &mut [u8], off, v: i32| {
            b[off..off + 4].copy_from_slice(&if big {
                v.to_be_bytes()
            } else {
                v.to_le_bytes()
            })
        };
        let put_f = |b: &mut [u8], off, v: f64| {
            b[off..off + 8].copy_from_slice(&if big {
                v.to_be_bytes()
            } else {
                v.to_le_bytes()
            })
        };
        put_i(&mut b, 8, 2);
        put_i(&mut b, 12, 5);
        put_i(&mut b, 76, 2);
        put_f(&mut b, 1040, 1.0);
        put_f(&mut b, 1048, -100.0);
        put_f(&mut b, 1056, 100.0);
        for (i, v) in [31008, 1, 2, 385, 404].iter().enumerate() {
            put_i(&mut b, 1064 + i * 4, *v);
        }
        let data = [
            -50.0, 50.0, 0.1, 0.01, 0.4, 0.0, 0.2, 0.02, 50.0, 50.0, 0.12, 0.01, 0.4, 0.0, 0.24,
            0.02, -100.0, 100.0, 8.0, 2.0,
        ];
        for (i, v) in data.iter().enumerate() {
            put_f(&mut b, 3072 + i * 8, *v);
        }
        b
    }

    #[test]
    fn endian_boundaries_inverse_and_continuity() {
        let le = BinaryPck::from_bytes(kernel(false)).unwrap();
        let be = BinaryPck::from_bytes(kernel(true)).unwrap();
        for et in [-100.0, -50.0, 0.0, 50.0, 100.0] {
            let (_, m) = le.rotation(31008, et).unwrap().unwrap();
            assert_eq!(be.rotation(31008, et).unwrap().unwrap().1, m);
            for i in 0..3 {
                for j in 0..3 {
                    let d: f64 = (0..3).map(|k| m[i][k] * m[j][k]).sum();
                    assert!((d - if i == j { 1.0 } else { 0.0 }).abs() < 1e-14);
                }
            }
        }
        let a = le.rotation(31008, -1e-6).unwrap().unwrap().1;
        let z = le.rotation(31008, 1e-6).unwrap().unwrap().1;
        assert!((a[0][0] - z[0][0]).abs() < 1e-8);
        assert!(le.rotation(31008, 100.0001).unwrap().is_none());
        assert!(le.rotation(31008, f64::NAN).is_err());
        assert!(le.rotation(1, 0.0).unwrap().is_none());
    }

    #[test]
    fn malformed_inputs_fail_without_panics() {
        let good = kernel(false);
        for len in [0, 7, 96, 1023, 1024, 2047, 3072, 3200] {
            assert!(BinaryPck::from_bytes(good[..len].to_vec()).is_err());
        }
        for (offset, bytes) in [
            (76, 2i32.to_le_bytes().to_vec()), // replaced below with a self-cycle
            (1064 + 12, 0i32.to_le_bytes().to_vec()),
            (1064 + 16, i32::MAX.to_le_bytes().to_vec()),
            (1040, f64::INFINITY.to_le_bytes().to_vec()),
            (3072 + 8, 0.0_f64.to_le_bytes().to_vec()),
            (3072 + 16, f64::NAN.to_le_bytes().to_vec()),
            (3072 + 18 * 8, 7.5_f64.to_le_bytes().to_vec()),
            (3072 + 19 * 8, 1e300_f64.to_le_bytes().to_vec()),
        ] {
            let mut b = good.clone();
            b[offset..offset + bytes.len()].copy_from_slice(&bytes);
            if offset == 76 {
                b[1024..1032].copy_from_slice(&2.0_f64.to_le_bytes());
            }
            assert!(BinaryPck::from_bytes(b).is_err(), "offset {offset}");
        }
        let mut b = good;
        b[1072..1076].copy_from_slice(&20i32.to_le_bytes());
        let unsupported = BinaryPck::from_bytes(b).unwrap();
        assert!(matches!(
            unsupported.rotation(31008, 0.0),
            Err(SpiceError::UnsupportedKernelQuery { .. })
        ));
    }
}
