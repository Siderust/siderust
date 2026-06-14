// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Integration tests for the `siderust::spice` high-level context.

use siderust::formats::spice::SpiceError;
use siderust::spice::{KernelSet, SpiceContext, SpiceContextError};

fn minimal_valid_spk() -> Vec<u8> {
    let mut buf = vec![0u8; 1024];
    buf[0..8].copy_from_slice(b"DAF/SPK ");
    buf[8..12].copy_from_slice(&2i32.to_le_bytes());
    buf[12..16].copy_from_slice(&6i32.to_le_bytes());
    buf[76..80].copy_from_slice(&0i32.to_le_bytes());
    buf
}

#[test]
fn spice_context_new_constructs_empty_context() {
    let context = SpiceContext::new();
    assert_eq!(context.kernel_set().spk_count(), 0);
    assert_eq!(context.kernel_set().ck_count(), 0);
    assert!(!context.kernel_set().has_lsk());
}

#[test]
fn kernel_set_new_is_empty() {
    let set = KernelSet::new();
    assert_eq!(set.spk_count(), 0);
    assert_eq!(set.ck_count(), 0);
    assert!(!set.has_lsk());
}

#[test]
fn kernel_set_state_naif_empty_returns_no_chain() {
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
fn spice_context_load_bytes_accepts_minimal_spk() {
    let mut context = SpiceContext::new();
    context
        .load_bytes("empty.bsp", minimal_valid_spk())
        .unwrap();
    assert_eq!(context.kernel_set().spk_count(), 1);
}

#[test]
fn tai_minus_utc_without_lsk_returns_kernel_not_loaded() {
    let context = SpiceContext::new();
    let error = context.tai_minus_utc(0.0).unwrap_err();
    assert!(matches!(
        error,
        SpiceContextError::KernelNotLoaded { ref kernel_type } if kernel_type == "LSK"
    ));
}

#[test]
fn rotation_naif_same_frame_is_identity() {
    let context = SpiceContext::new();
    let matrix = context.rotation_naif(1, 1, 0.0).unwrap();
    assert_eq!(matrix, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
}

#[test]
fn state_naif_with_loaded_empty_spk_returns_no_chain() {
    let mut context = SpiceContext::new();
    context
        .load_bytes("empty.bsp", minimal_valid_spk())
        .unwrap();
    let error = context.state_naif(399, 0, 0.0).unwrap_err();
    assert!(matches!(
        error,
        SpiceContextError::Kernel(SpiceError::NoChain {
            target: 399,
            center: 0
        })
    ));
}
