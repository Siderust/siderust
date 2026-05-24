//! POD I/O compatibility aliases.
//!
//! General-purpose format parsers and writers live in [`crate::formats`].
//! This module re-exports the subset of `formats` items that the POD modules
//! reference internally.
#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub use crate::formats::{FileLocation, FormatError, ParseMode};
