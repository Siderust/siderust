// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Ordered set of loaded SPICE kernels.

use std::sync::Arc;

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
        for (_, kernel) in self.spk.iter().rev() {
            match kernel.state(target, center, epoch_tdb_s) {
                Ok(state) => return Ok(state),
                Err(SpiceError::NoChain { .. } | SpiceError::OutOfCoverage { .. }) => continue,
                Err(other) => return Err(other),
            }
        }
        Err(SpiceError::NoChain { target, center })
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
