// SPDX-License-Identifier: AGPL-3.0-or-later

/// Status codes returned by siderust-ffi functions.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustStatus {
    /// Success.
    Ok = 0,
    /// A required output pointer was null.
    NullPointer = 1,
    /// Invalid or unsupported reference frame.
    InvalidFrame = 2,
    /// Invalid or unsupported reference center.
    InvalidCenter = 3,
    /// Coordinate transformation failed.
    TransformFailed = 4,
    /// Invalid body name or parameters.
    InvalidBody = 5,
    /// Star name not found in catalog.
    UnknownStar = 6,
    /// Invalid period (start > end).
    InvalidPeriod = 7,
    /// Memory allocation failed.
    AllocationFailed = 8,
    /// Invalid argument or parameter.
    InvalidArgument = 9,
}
