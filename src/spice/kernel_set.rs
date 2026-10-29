// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Ordered set of loaded SPICE kernels.

use alloc::string::String;
use alloc::sync::Arc;
use alloc::vec::Vec;

use crate::formats::spice::{
    CkKernel, FrameKernel, IkKernel, LeapSecondKernel, PckKernel, SclkKernel, SpiceError, SpkKernel,
};

/// Ordered collection of loaded SPICE kernels.
///
/// Kernels are stored in load order (most recently loaded last), following the
/// NAIF "last loaded wins" priority rule.
///
/// # Examples
///
/// ```rust
/// use siderust::spice::KernelSet;
///
/// let set = KernelSet::new();
/// assert_eq!(set.spk_count(), 0);
/// ```
#[derive(Default)]
pub struct KernelSet {
    bpc: Vec<(String, Arc<crate::formats::spice::BinaryPck>)>,
    spk: Vec<(String, Arc<SpkKernel>)>,
    ck: Vec<(String, Arc<CkKernel>)>,
    lsk: Option<(String, Arc<LeapSecondKernel>)>,
    fk: Vec<(String, Arc<FrameKernel>)>,
    pck_text: Vec<(String, Arc<PckKernel>)>,
    sclk: Vec<(String, Arc<SclkKernel>)>,
    ik: Vec<(String, Arc<IkKernel>)>,
}

impl KernelSet {
    /// Create an empty kernel set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an SPK kernel.
    pub fn add_spk(&mut self, alias: impl Into<String>, kernel: SpkKernel) {
        self.spk.push((alias.into(), Arc::new(kernel)));
    }

    /// Add a binary PCK kernel, with last-loaded priority.
    pub fn add_bpc(&mut self, alias: impl Into<String>, kernel: crate::formats::spice::BinaryPck) {
        self.bpc.push((alias.into(), Arc::new(kernel)));
    }

    /// Binary PCKs in load order, including provenance aliases.
    pub fn binary_pcks(
        &self,
    ) -> impl DoubleEndedIterator<Item = (&str, &crate::formats::spice::BinaryPck)> {
        self.bpc
            .iter()
            .map(|(alias, kernel)| (alias.as_str(), kernel.as_ref()))
    }

    /// Kernel aliases in load order within each kernel type (SPK, BPC, FK, text PCK, LSK).
    pub fn provenance(&self) -> Vec<&str> {
        self.spk
            .iter()
            .map(|(s, _)| s.as_str())
            .chain(self.bpc.iter().map(|(s, _)| s.as_str()))
            .chain(self.fk.iter().map(|(s, _)| s.as_str()))
            .chain(self.pck_text.iter().map(|(s, _)| s.as_str()))
            .chain(self.lsk.iter().map(|(s, _)| s.as_str()))
            .collect()
    }

    /// Add a CK kernel.
    pub fn add_ck(&mut self, alias: impl Into<String>, kernel: CkKernel) {
        self.ck.push((alias.into(), Arc::new(kernel)));
    }

    /// Set the active LSK kernel.
    pub fn set_lsk(&mut self, alias: impl Into<String>, kernel: LeapSecondKernel) {
        self.lsk = Some((alias.into(), Arc::new(kernel)));
    }

    /// Add an FK kernel.
    pub fn add_fk(&mut self, alias: impl Into<String>, kernel: FrameKernel) {
        self.fk.push((alias.into(), Arc::new(kernel)));
    }

    /// Add a text PCK kernel.
    pub fn add_pck_text(&mut self, alias: impl Into<String>, kernel: PckKernel) {
        self.pck_text.push((alias.into(), Arc::new(kernel)));
    }

    /// Add an SCLK kernel.
    pub fn add_sclk(&mut self, alias: impl Into<String>, kernel: SclkKernel) {
        self.sclk.push((alias.into(), Arc::new(kernel)));
    }

    /// Add an IK kernel.
    pub fn add_ik(&mut self, alias: impl Into<String>, kernel: IkKernel) {
        self.ik.push((alias.into(), Arc::new(kernel)));
    }

    /// Number of loaded SPK kernels.
    pub fn spk_count(&self) -> usize {
        self.spk.len()
    }

    /// Number of loaded CK kernels.
    pub fn ck_count(&self) -> usize {
        self.ck.len()
    }

    /// Whether an LSK is loaded.
    pub fn has_lsk(&self) -> bool {
        self.lsk.is_some()
    }

    /// Borrow the active LSK, if any.
    pub fn lsk(&self) -> Option<&Arc<LeapSecondKernel>> {
        self.lsk.as_ref().map(|(_, kernel)| kernel)
    }

    /// Borrow all text-PCK kernels.
    pub fn pck_text_kernels(&self) -> impl DoubleEndedIterator<Item = &Arc<PckKernel>> + '_ {
        self.pck_text.iter().map(|(_, kernel)| kernel)
    }

    /// Borrow all FK kernels.
    pub fn fk_kernels(&self) -> impl DoubleEndedIterator<Item = &Arc<FrameKernel>> + '_ {
        self.fk.iter().map(|(_, kernel)| kernel)
    }

    /// Borrow all CK kernels.
    pub fn ck_kernels(&self) -> impl DoubleEndedIterator<Item = &Arc<CkKernel>> + '_ {
        self.ck.iter().map(|(_, kernel)| kernel)
    }

    /// Borrow all SCLK kernels.
    pub fn sclk_kernels(&self) -> impl DoubleEndedIterator<Item = &Arc<SclkKernel>> + '_ {
        self.sclk.iter().map(|(_, kernel)| kernel)
    }

    /// Compute state `[x, y, z, vx, vy, vz]` in km and km/s.
    pub fn state_naif(
        &self,
        target: i32,
        center: i32,
        epoch_tdb_s: f64,
    ) -> Result<[f64; 6], SpiceError> {
        if !epoch_tdb_s.is_finite() {
            return Err(SpiceError::Corrupted {
                message: "SPK epoch must be finite TDB seconds".into(),
            });
        }
        if target == center {
            return Ok([0.0; 6]);
        }
        // Select one edge per target, respecting file and segment priority.
        // Resolve to the first common ancestor before accumulating offsets:
        // Earth/Moon therefore subtract their EMB offsets, not AU-scale states.
        let select = |body| {
            self.spk
                .iter()
                .rev()
                .flat_map(|(_, k)| k.segments().iter().rev())
                .find(|s| {
                    s.target == body
                        && epoch_tdb_s >= s.start_tdb_seconds
                        && epoch_tdb_s <= s.end_tdb_seconds
                })
        };
        let mut target_nodes = alloc::collections::BTreeMap::new();
        let mut target_edges = Vec::new();
        let mut node = target;
        loop {
            if target_nodes.insert(node, target_edges.len()).is_some() {
                return Err(SpiceError::Corrupted {
                    message: "SPK centre cycle".into(),
                });
            }
            if node == center {
                break;
            }
            let Some(edge) = select(node) else {
                break;
            };
            target_edges.push(edge);
            node = edge.center;
        }
        let target_leaf = node;
        let mut center_edges = Vec::new();
        let mut visited = alloc::collections::BTreeSet::new();
        node = center;
        while !target_nodes.contains_key(&node) {
            if !visited.insert(node) {
                return Err(SpiceError::Corrupted {
                    message: "SPK centre cycle".into(),
                });
            }
            let Some(edge) = select(node) else {
                if let Some(s) = self
                    .spk
                    .iter()
                    .rev()
                    .flat_map(|(_, k)| k.segments().iter().rev())
                    .find(|s| {
                        (s.target == node || s.target == target_leaf)
                            && (epoch_tdb_s < s.start_tdb_seconds
                                || epoch_tdb_s > s.end_tdb_seconds)
                    })
                {
                    return Err(SpiceError::OutOfCoverage {
                        target,
                        center,
                        epoch_tdb_seconds: epoch_tdb_s,
                        start_tdb_seconds: s.start_tdb_seconds,
                        end_tdb_seconds: s.end_tdb_seconds,
                    });
                }
                return Err(SpiceError::NoChain { target, center });
            };
            center_edges.push(edge);
            node = edge.center;
        }
        let mut state = [0.0; 6];
        for (edge, sign) in target_edges[..target_nodes[&node]]
            .iter()
            .map(|e| (*e, 1.0))
            .chain(center_edges.iter().map(|e| (*e, -1.0)))
        {
            if edge.frame_id != 1 {
                return Err(SpiceError::UnsupportedKernelQuery {
                    message: alloc::format!(
                        "SPK frame {} requires state rotation; only J2000 supported",
                        edge.frame_id
                    ),
                });
            }
            let segment = edge
                .segment
                .as_ref()
                .ok_or(SpiceError::UnsupportedDataType {
                    data_type: edge.data_type,
                })?;
            let offset = segment.evaluate(epoch_tdb_s)?;
            for i in 0..6 {
                state[i] += sign * offset[i];
            }
        }
        Ok(state)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::formats::spice::{LeapSecondKernel, SpiceError};

    fn minimal_valid_spk() -> SpkKernel {
        let mut buf = vec![0u8; 1024];
        buf[0..8].copy_from_slice(b"DAF/SPK ");
        buf[8..12].copy_from_slice(&2i32.to_le_bytes());
        buf[12..16].copy_from_slice(&6i32.to_le_bytes());
        buf[76..80].copy_from_slice(&0i32.to_le_bytes());
        SpkKernel::from_bytes(buf).unwrap()
    }

    #[test]
    fn accessors_reflect_loaded_kernels() {
        let mut set = KernelSet::new();
        assert_eq!(set.spk_count(), 0);
        assert!(!set.has_lsk());

        set.add_spk("spk", minimal_valid_spk());
        assert_eq!(set.spk_count(), 1);

        let lsk_src =
            "\\begindata\nDELTET/DELTA_T_A = 32.184\nDELTET/DELTA_AT = ( 37 @2017-JAN-1 )\n";
        set.set_lsk("lsk", LeapSecondKernel::from_text(lsk_src).unwrap());
        assert!(set.has_lsk());
        assert!(set.lsk().is_some());
        assert_eq!(set.fk_kernels().count(), 0);
        assert_eq!(set.pck_text_kernels().count(), 0);
        assert_eq!(set.ck_kernels().count(), 0);
        assert_eq!(set.sclk_kernels().count(), 0);
    }

    #[test]
    fn state_naif_without_spk_returns_no_chain() {
        let set = KernelSet::new();
        let error = set.state_naif(399, 0, 0.0).unwrap_err();
        assert!(matches!(
            error,
            SpiceError::NoChain {
                target: 399,
                center: 0
            }
        ));
    }

    #[test]
    fn state_naif_with_empty_spk_still_returns_no_chain() {
        let mut set = KernelSet::new();
        set.add_spk("empty", minimal_valid_spk());
        let error = set.state_naif(399, 0, 0.0).unwrap_err();
        assert!(matches!(
            error,
            SpiceError::NoChain {
                target: 399,
                center: 0
            }
        ));
    }
}
