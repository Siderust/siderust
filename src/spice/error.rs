// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! High-level SPICE context errors.

use thiserror::Error;

use crate::formats::spice::SpiceError;

/// Errors produced by the high-level [`crate::spice::SpiceContext`].
#[derive(Debug, Error)]
pub enum SpiceContextError {
    /// A low-level kernel parse or query error.
    #[error("SPICE kernel error: {0}")]
    Kernel(#[from] SpiceError),
    /// I/O error loading a kernel from disk.
    #[error("SPICE context I/O error: {0}")]
    Io(#[from] std::io::Error),
    /// No kernel of the required type has been loaded.
    #[error("SPICE context: no {kernel_type} kernel loaded")]
    KernelNotLoaded {
        /// Human-readable kernel type description (for example `LSK` or `CK`).
        kernel_type: String,
    },
    /// The requested frame is not registered in any loaded kernel.
    #[error("SPICE context: unknown frame: {description}")]
    UnknownFrame {
        /// Frame ID or name description.
        description: String,
    },
    /// A time conversion failed.
    #[error("SPICE context: time conversion error: {message}")]
    TimeConversion {
        /// Human-readable reason.
        message: String,
    },
    /// A kernel query is not implemented.
    #[error("SPICE context: unsupported query: {message}")]
    UnsupportedKernelQuery {
        /// Description of the unsupported query.
        message: String,
    },
}
