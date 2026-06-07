// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! CCSDS exchange formats.

pub use crate::formats::error::{FileLocation, FormatError, ParseMode};

pub mod oem;
pub mod opm;
pub mod tdm;
