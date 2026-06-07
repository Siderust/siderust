// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! [`SpkKernel`] — owns DAF bytes, indexes every segment, and resolves
//! body-relative state queries by chaining segments through the kernel graph.
//!
//! ## Technical scope
//!
//! This module exposes the low-level SPK kernel reader used by higher-level
//! ephemeris adapters. It works directly in NAIF wire-format units: epochs are
//! raw TDB seconds past J2000 and states are raw `[f64; 6]` vectors in km and
//! km/s.
//!
//! ## API contract
//!
//! [`SpkKernel::state`] is intentionally a wire-format API for kernel-native
//! access and performance-sensitive adapters. Callers looking for the typed
//! scientific API should prefer the `pod` feature's `SpiceEphemerisProvider`, which
//! lift epochs, centers, frames, and units into typed `tempoch`/`affn`/`qtty` abstractions.

use std::collections::{HashMap, VecDeque};
use std::path::Path;

use crate::formats::spice::daf::Daf;

use super::error::SpiceError;
use super::segment::{segment_for_summary, SpkSegment};

/// One segment loaded from a kernel, paired with its summary metadata.
#[derive(Debug)]
pub struct LoadedSegment {
    /// NAIF id of the target body for this segment.
    pub target: i32,
    /// NAIF id of the center body for this segment.
    pub center: i32,
    /// SPICE frame id (typically 1 = J2000 ≡ ICRS for DE kernels).
    pub frame_id: i32,
    /// SPK data type code (`2`, `3`, `9`, or `13` for supported segments).
    pub data_type: i32,
    /// Coverage start (TDB seconds past J2000).
    pub start_tdb_seconds: f64,
    /// Coverage end (TDB seconds past J2000).
    pub end_tdb_seconds: f64,
    /// Decoded segment payload, or `None` if the data type is not
    /// supported by this crate.
    pub segment: Option<SpkSegment>,
}

/// A loaded SPK kernel ready to answer body-relative state queries.
///
/// `SpkKernel` owns the raw bytes of the DAF file and a parsed list of
/// every segment in the kernel. State queries are resolved by walking
/// the directed graph of `(target, center)` segments using BFS, summing
/// segment states along the chain. This matches the behavior of NAIF
/// CSPICE `spkez_c` / `spkgeo_c` for J2000-frame ephemerides.
///
/// # Examples
///
/// ```rust,no_run
/// use siderust::formats::spice::SpkKernel;
///
/// let kernel = SpkKernel::open("de440.bsp")?;
/// // Earth (NAIF 399) relative to the Solar System Barycenter (NAIF 0):
/// let state = kernel.state(399, 0, 0.0)?;
/// assert!(state[0].abs() < 2.0e8); // ~1.5e8 km from SSB at J2000
/// # Ok::<_, siderust::formats::spice::SpiceError>(())
/// ```
pub struct SpkKernel {
    bytes: Vec<u8>,
    daf: Daf,
    segments: Vec<LoadedSegment>,
    // Adjacency map: target -> Vec of segment indices.
    by_target: HashMap<i32, Vec<usize>>,
}

impl std::fmt::Debug for SpkKernel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SpkKernel")
            .field("bytes_len", &self.bytes.len())
            .field("segment_count", &self.segments.len())
            .finish()
    }
}

impl SpkKernel {
    /// Open and parse a DAF/SPK kernel from a filesystem path.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, SpiceError> {
        let bytes = std::fs::read(path)?;
        Self::from_bytes(bytes)
    }

    /// Parse a DAF/SPK kernel from an in-memory buffer.
    ///
    /// Useful for tests and for embedding small kernels (e.g. unit
    /// tests, fuzz inputs, or memory-mapped files).
    pub fn from_bytes(bytes: Vec<u8>) -> Result<Self, SpiceError> {
        let daf = Daf::parse(&bytes)?;
        let mut segments = Vec::with_capacity(daf.summaries.len());
        for summary in &daf.summaries {
            let (segment, start, end) = match segment_for_summary(&bytes, &daf, summary) {
                Ok(seg) => {
                    let start = seg.start_tdb_seconds();
                    let end = seg.end_tdb_seconds();
                    (Some(seg), start, end)
                }
                Err(SpiceError::UnsupportedDataType { .. }) => {
                    // Index the segment but mark it as unevaluable.
                    (None, summary.start_et, summary.end_et)
                }
                Err(other) => return Err(other),
            };
            segments.push(LoadedSegment {
                target: summary.target_id,
                center: summary.center_id,
                frame_id: summary.frame_id,
                data_type: summary.data_type,
                start_tdb_seconds: start,
                end_tdb_seconds: end,
                segment,
            });
        }

        let mut by_target: HashMap<i32, Vec<usize>> = HashMap::new();
        for (i, s) in segments.iter().enumerate() {
            by_target.entry(s.target).or_default().push(i);
        }

        Ok(SpkKernel {
            bytes,
            daf,
            segments,
            by_target,
        })
    }

    /// Borrow the parsed DAF metadata. Useful for low-level inspection.
    pub fn daf(&self) -> &Daf {
        &self.daf
    }

    /// Borrow the raw kernel bytes.
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// All loaded segments in summary order.
    pub fn segments(&self) -> &[LoadedSegment] {
        &self.segments
    }

    /// Compute the state of `target` relative to `center` at TDB
    /// `epoch_tdb_seconds` (seconds past J2000), expressed in the
    /// kernel's frame (`frame_id`, typically J2000 ≡ ICRS).
    ///
    /// The result is `[x, y, z, vx, vy, vz]` in km / km·s⁻¹.
    ///
    /// The chain is resolved by BFS over the directed segment graph;
    /// per-segment states are summed (with sign for reverse traversal).
    ///
    /// # Safety / Wire-format API
    ///
    /// This method intentionally exposes raw SPICE wire-format scalars. Prefer
    /// a typed `EphemerisProvider` implementation such as `SpiceEphemerisProvider`
    /// (both from the `pod` feature) when you want a scientific API with typed
    /// epochs, centers, frames, and units.
    pub fn state(
        &self,
        target: i32,
        center: i32,
        epoch_tdb_seconds: f64,
    ) -> Result<[f64; 6], SpiceError> {
        if target == center {
            return Ok([0.0; 6]);
        }
        // BFS over directed edges.  Each "edge" is the segment that
        // expresses `seg.target` relative to `seg.center`. A traversal
        // from `target` to `center` walks edges in either direction;
        // a forward edge contributes +state, a reverse edge contributes
        // -state.
        let path = self.find_chain(target, center, epoch_tdb_seconds)?;
        let mut acc = [0.0_f64; 6];
        let mut current = target;
        for (seg_idx, forward) in path {
            let seg = &self.segments[seg_idx];
            let segment = seg
                .segment
                .as_ref()
                .ok_or(SpiceError::UnsupportedDataType {
                    data_type: seg.data_type,
                })?;
            let s = segment.evaluate(epoch_tdb_seconds).map_err(|e| match e {
                // Re-stamp the chain endpoints onto coverage errors.
                SpiceError::OutOfCoverage {
                    epoch_tdb_seconds,
                    start_tdb_seconds,
                    end_tdb_seconds,
                    ..
                } => SpiceError::OutOfCoverage {
                    target,
                    center,
                    epoch_tdb_seconds,
                    start_tdb_seconds,
                    end_tdb_seconds,
                },
                other => other,
            })?;
            if forward {
                for i in 0..6 {
                    acc[i] += s[i];
                }
                current = seg.center;
            } else {
                for i in 0..6 {
                    acc[i] -= s[i];
                }
                current = seg.target;
            }
        }
        debug_assert_eq!(current, center);
        let _ = current;
        Ok(acc)
    }

    /// BFS over the segment graph; returns the list of `(segment_index,
    /// forward)` to walk to get from `target` to `center` at `epoch`.
    fn find_chain(
        &self,
        target: i32,
        center: i32,
        epoch: f64,
    ) -> Result<Vec<(usize, bool)>, SpiceError> {
        // Build undirected adjacency on demand: for each node we visit,
        // walk every segment whose target or center matches it and
        // whose coverage covers `epoch`.
        let mut visited: HashMap<i32, (i32, usize, bool)> = HashMap::new();
        // value: (predecessor_node, edge_segment_index, forward_from_pred)
        let mut queue = VecDeque::new();
        queue.push_back(target);
        visited.insert(target, (target, usize::MAX, true));

        let mut found = false;
        while let Some(node) = queue.pop_front() {
            if node == center {
                found = true;
                break;
            }
            // Outgoing forward edges (segments where node is the target).
            if let Some(idxs) = self.by_target.get(&node) {
                for &i in idxs {
                    let s = &self.segments[i];
                    if epoch < s.start_tdb_seconds || epoch > s.end_tdb_seconds {
                        continue;
                    }
                    if let std::collections::hash_map::Entry::Vacant(e) = visited.entry(s.center) {
                        e.insert((node, i, true));
                        queue.push_back(s.center);
                    }
                }
            }
            // Reverse edges (segments where node is the center).
            for (i, s) in self.segments.iter().enumerate() {
                if s.center != node {
                    continue;
                }
                if epoch < s.start_tdb_seconds || epoch > s.end_tdb_seconds {
                    continue;
                }
                if let std::collections::hash_map::Entry::Vacant(e) = visited.entry(s.target) {
                    e.insert((node, i, false));
                    queue.push_back(s.target);
                }
            }
        }

        if !found {
            // Distinguish "no chain at all" vs "chain exists but not at this epoch".
            let any_chain_ignoring_epoch = self.has_chain_ignoring_epoch(target, center);
            if any_chain_ignoring_epoch {
                // Find tightest window for diagnostics.
                let (start, end) = self
                    .coverage_for(target, center)
                    .unwrap_or((f64::NAN, f64::NAN));
                return Err(SpiceError::OutOfCoverage {
                    target,
                    center,
                    epoch_tdb_seconds: epoch,
                    start_tdb_seconds: start,
                    end_tdb_seconds: end,
                });
            }
            return Err(SpiceError::NoChain { target, center });
        }

        // Reconstruct path from `center` back to `target`.
        let mut steps = Vec::new();
        let mut node = center;
        while node != target {
            let (pred, seg_idx, forward) = *visited.get(&node).expect("BFS invariant");
            steps.push((seg_idx, forward));
            node = pred;
        }
        steps.reverse();
        Ok(steps)
    }

    fn has_chain_ignoring_epoch(&self, target: i32, center: i32) -> bool {
        let mut visited = std::collections::HashSet::new();
        let mut queue = VecDeque::new();
        queue.push_back(target);
        visited.insert(target);
        while let Some(node) = queue.pop_front() {
            if node == center {
                return true;
            }
            for s in &self.segments {
                let next = if s.target == node {
                    Some(s.center)
                } else if s.center == node {
                    Some(s.target)
                } else {
                    None
                };
                if let Some(n) = next {
                    if visited.insert(n) {
                        queue.push_back(n);
                    }
                }
            }
        }
        false
    }

    fn coverage_for(&self, target: i32, center: i32) -> Option<(f64, f64)> {
        let mut start = f64::NEG_INFINITY;
        let mut end = f64::INFINITY;
        let mut hit = false;
        for s in &self.segments {
            if (s.target == target && s.center == center)
                || (s.target == center && s.center == target)
            {
                start = start.max(s.start_tdb_seconds);
                end = end.min(s.end_tdb_seconds);
                hit = true;
            }
        }
        if hit {
            Some((start, end))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a minimal valid DAF/SPK file with two Type-2 segments
    /// (Earth → SSB and Moon → Earth), each a single-record constant
    /// polynomial.
    fn build_two_segment_kernel() -> Vec<u8> {
        // Layout:
        //   record 1 (1024 bytes): file header
        //   record 2 (1024 bytes): summary record (2 summaries)
        //   record 3 (1024 bytes): segment data (two trailers + two records)
        let mut buf = vec![0u8; 3 * 1024];

        // ── header ────────────────────────────────────────────────────
        buf[0..8].copy_from_slice(b"DAF/SPK ");
        buf[8..12].copy_from_slice(&2i32.to_le_bytes()); // ND
        buf[12..16].copy_from_slice(&6i32.to_le_bytes()); // NI
        buf[76..80].copy_from_slice(&2i32.to_le_bytes()); // FWARD = 2

        // ── summary record (record 2, bytes 1024..2048) ───────────────
        let rec = &mut buf[1024..2048];
        rec[0..8].copy_from_slice(&0.0_f64.to_le_bytes()); // next
        rec[8..16].copy_from_slice(&0.0_f64.to_le_bytes()); // prev
        rec[16..24].copy_from_slice(&2.0_f64.to_le_bytes()); // nsum

        // SS (size of summary in doubles) = nd + ceil(ni/2) = 2 + 3 = 5
        // Each summary occupies 5 * 8 = 40 bytes.

        // We will lay out segment data starting at file word
        //   2 * 1024 / 8 + 1 = 257 (1-based word index of byte 2048).
        let seg1_start_word: i32 = 257;
        let seg1_end_word: i32 = seg1_start_word + 8; // 1 record (5 doubles) + 4 trailer = 9 words
        let seg2_start_word: i32 = seg1_end_word + 1;
        let seg2_end_word: i32 = seg2_start_word + 8;

        // Summary 1: target=399 (Earth), center=0 (SSB)
        let s1 = &mut rec[24..24 + 40];
        s1[0..8].copy_from_slice(&0.0_f64.to_le_bytes()); // start_et
        s1[8..16].copy_from_slice(&1.0_f64.to_le_bytes()); // end_et
        s1[16..20].copy_from_slice(&399i32.to_le_bytes()); // target
        s1[20..24].copy_from_slice(&0i32.to_le_bytes()); // center
        s1[24..28].copy_from_slice(&1i32.to_le_bytes()); // frame
        s1[28..32].copy_from_slice(&2i32.to_le_bytes()); // data_type
        s1[32..36].copy_from_slice(&seg1_start_word.to_le_bytes());
        s1[36..40].copy_from_slice(&seg1_end_word.to_le_bytes());

        // Summary 2: target=301 (Moon), center=399 (Earth)
        let s2 = &mut rec[24 + 40..24 + 80];
        s2[0..8].copy_from_slice(&0.0_f64.to_le_bytes());
        s2[8..16].copy_from_slice(&1.0_f64.to_le_bytes());
        s2[16..20].copy_from_slice(&301i32.to_le_bytes());
        s2[20..24].copy_from_slice(&399i32.to_le_bytes());
        s2[24..28].copy_from_slice(&1i32.to_le_bytes());
        s2[28..32].copy_from_slice(&2i32.to_le_bytes());
        s2[32..36].copy_from_slice(&seg2_start_word.to_le_bytes());
        s2[36..40].copy_from_slice(&seg2_end_word.to_le_bytes());

        // ── segment data (record 3, bytes 2048..3072) ─────────────────
        let wf = |word: i32, v: f64, b: &mut [u8]| {
            let off = (word as usize - 1) * 8;
            b[off..off + 8].copy_from_slice(&v.to_le_bytes());
        };
        // Segment 1: Earth wrt SSB = (1.0e8, 0, 0)
        // record at words seg1_start_word..seg1_start_word+4
        wf(seg1_start_word, 0.5, &mut buf); // mid
        wf(seg1_start_word + 1, 0.5, &mut buf); // radius
        wf(seg1_start_word + 2, 1.0e8, &mut buf); // x coeff
        wf(seg1_start_word + 3, 0.0, &mut buf);
        wf(seg1_start_word + 4, 0.0, &mut buf);
        // trailer
        wf(seg1_end_word - 3, 0.0, &mut buf); // init
        wf(seg1_end_word - 2, 1.0, &mut buf); // intlen
        wf(seg1_end_word - 1, 5.0, &mut buf); // rsize
        wf(seg1_end_word, 1.0, &mut buf); // n_records

        // Segment 2: Moon wrt Earth = (0, 3.84e5, 0)
        wf(seg2_start_word, 0.5, &mut buf);
        wf(seg2_start_word + 1, 0.5, &mut buf);
        wf(seg2_start_word + 2, 0.0, &mut buf);
        wf(seg2_start_word + 3, 3.84e5, &mut buf);
        wf(seg2_start_word + 4, 0.0, &mut buf);
        wf(seg2_end_word - 3, 0.0, &mut buf);
        wf(seg2_end_word - 2, 1.0, &mut buf);
        wf(seg2_end_word - 1, 5.0, &mut buf);
        wf(seg2_end_word, 1.0, &mut buf);

        buf
    }

    #[test]
    fn open_synthetic_kernel_indexes_segments() {
        let bytes = build_two_segment_kernel();
        let kernel = SpkKernel::from_bytes(bytes).unwrap();
        assert_eq!(kernel.segments().len(), 2);
        assert_eq!(kernel.segments()[0].target, 399);
        assert_eq!(kernel.segments()[0].center, 0);
        assert_eq!(kernel.segments()[1].target, 301);
        assert_eq!(kernel.segments()[1].center, 399);
    }

    #[test]
    fn direct_state_round_trip() {
        let kernel = SpkKernel::from_bytes(build_two_segment_kernel()).unwrap();
        let s = kernel.state(399, 0, 0.5).unwrap();
        assert!((s[0] - 1.0e8).abs() < 1e-6);
    }

    #[test]
    fn chained_state_moon_relative_to_ssb() {
        let kernel = SpkKernel::from_bytes(build_two_segment_kernel()).unwrap();
        // Moon → SSB = (Moon → Earth) + (Earth → SSB) = (1e8, 3.84e5, 0)
        let s = kernel.state(301, 0, 0.5).unwrap();
        assert!((s[0] - 1.0e8).abs() < 1e-6);
        assert!((s[1] - 3.84e5).abs() < 1e-6);
    }

    #[test]
    fn reverse_chain_state_ssb_relative_to_moon() {
        let kernel = SpkKernel::from_bytes(build_two_segment_kernel()).unwrap();
        // SSB → Moon = -(Moon → SSB)
        let s = kernel.state(0, 301, 0.5).unwrap();
        assert!((s[0] + 1.0e8).abs() < 1e-6);
        assert!((s[1] + 3.84e5).abs() < 1e-6);
    }

    #[test]
    fn same_body_state_is_zero() {
        let kernel = SpkKernel::from_bytes(build_two_segment_kernel()).unwrap();
        let s = kernel.state(301, 301, 0.5).unwrap();
        assert_eq!(s, [0.0; 6]);
    }

    #[test]
    fn no_chain_returns_no_chain_error() {
        let kernel = SpkKernel::from_bytes(build_two_segment_kernel()).unwrap();
        let err = kernel.state(499 /* Mars */, 0, 0.5).unwrap_err();
        assert!(matches!(
            err,
            SpiceError::NoChain {
                target: 499,
                center: 0
            }
        ));
    }

    #[test]
    fn out_of_coverage_returned_for_far_future_epoch() {
        let kernel = SpkKernel::from_bytes(build_two_segment_kernel()).unwrap();
        let err = kernel.state(399, 0, 1.0e9).unwrap_err();
        assert!(matches!(
            err,
            SpiceError::OutOfCoverage {
                target: 399,
                center: 0,
                ..
            }
        ));
    }

    #[test]
    fn debug_repr_contains_segment_count() {
        let kernel = SpkKernel::from_bytes(build_two_segment_kernel()).unwrap();
        let s = format!("{kernel:?}");
        assert!(s.contains("segment_count: 2"));
    }
}
