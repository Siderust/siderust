// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! C FFI bindings for **siderust** — high-precision astronomy library.
//!
//! This crate exposes a flat C-compatible API for the siderust library,
//! covering coordinates, bodies, altitude calculations, ephemeris,
//! and observatories.

// Ensure the tempoch-ffi crate is linked so re-exports resolve.
extern crate tempoch_ffi;

pub mod error;
pub mod types;
pub mod observatories;
pub mod bodies;
pub mod coordinates;
pub mod altitude;
pub mod ephemeris;

/// Returns the siderust-ffi ABI version (semver-encoded: major*10000 + minor*100 + patch).
#[no_mangle]
pub extern "C" fn siderust_ffi_version() -> u32 {
    0 * 10000 + 1 * 100 + 0 // 0.1.0
}
