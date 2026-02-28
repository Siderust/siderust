// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Status codes returned by siderust-ffi functions.

/// Status codes returned by every siderust-ffi function.
///
/// Callers must inspect this value before reading any output parameters.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustStatus {
    /// Success.
    Ok = 0,
    /// A required output pointer was null.
    NullPointer = 1,
    /// The requested coordinate frame is unsupported or invalid.
    InvalidFrame = 2,
    /// The requested coordinate center is unsupported or invalid.
    InvalidCenter = 3,
    /// Coordinate transformation failed.
    TransformFailed = 4,
    /// Unknown or invalid body identifier.
    InvalidBody = 5,
    /// Star name not found in the built-in catalog.
    UnknownStar = 6,
    /// Period bounds are invalid (start > end).
    InvalidPeriod = 7,
    /// Memory allocation failed.
    AllocationFailed = 8,
    /// One or more arguments are out of range or otherwise invalid.
    InvalidArgument = 9,
    /// A Rust panic was caught at the FFI boundary.
    ///
    /// This should never happen in normal operation; it indicates a bug in the
    /// underlying library.  The panic payload is discarded.
    InternalPanic = 10,
}
