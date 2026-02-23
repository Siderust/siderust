// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! C FFI bindings for **siderust** — high-precision astronomy library.
//!
//! This crate exposes a flat, `#[no_mangle] extern "C"` API that can be consumed
//! from C, C++, or any language with a C FFI bridge.  All heap-allocated types
//! are opaque pointers managed via `{type}_free()` functions.
//!
//! The matching C header lives at `include/siderust_ffi.h` and must be included
//! **after** `tempoch_ffi.h` (for `tempoch_period_mjd_t`).

// Ensure the tempoch-ffi crate is linked so the TempochPeriodMjd re-export resolves.
extern crate tempoch_ffi;

pub mod altitude;
pub mod azimuth;
pub mod bodies;
pub mod coordinates;
pub mod ephemeris;
pub mod error;
pub mod observatories;
pub mod phase;
pub mod target;
pub mod types;

pub use altitude::*;
pub use azimuth::*;
pub use bodies::*;
pub use coordinates::*;
pub use ephemeris::*;
pub use error::*;
pub use observatories::*;
pub use phase::*;
pub use target::*;
pub use types::*;

/// Returns the siderust-ffi ABI version (major*10000 + minor*100 + patch).
#[allow(clippy::erasing_op, clippy::identity_op)]
#[no_mangle]
pub extern "C" fn siderust_ffi_version() -> u32 {
    0 * 10000 + 1 * 100 + 0 // 0.1.0
}
