// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Apparent-magnitude value type for stellar photometry inputs.
//!
//! Magnitudes are logarithmic astronomical brightness measurements. The v1
//! starlight model uses finite apparent magnitudes as inputs to an S10 flux
//! conversion.

use crate::starlight::{Result, StellarMapError};

/// Apparent magnitude newtype used by the v1 photometric model.
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct ApparentMagnitude(f64);

impl ApparentMagnitude {
    /// Construct a finite apparent magnitude.
    pub fn new(value: f64) -> Result<Self> {
        if !value.is_finite() {
            return Err(StellarMapError::InvalidMagnitude(value));
        }
        Ok(Self(value))
    }

    /// Return the magnitude value.
    #[must_use]
    pub fn value(self) -> f64 {
        self.0
    }
}
