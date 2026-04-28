// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Photometric filter passband datasets.
//!
//! Each sub-module provides lazily-initialised, statically-cached
//! [`SampledSpectrum`](crate::spectra::SampledSpectrum) constants for a
//! particular photometric system, following the same `OnceLock` pattern as
//! [`crate::atmosphere::ozone`].
//!
//! ## Units
//!
//! Filter throughputs are dimensionless values in \[0, 1\] represented by the
//! [`Throughput`] unit marker defined in this module.  It is analogous to
//! [`crate::atmosphere::ozone::Transmittance`] but scoped to filter passbands
//! so that the two concepts remain distinguishable in type signatures.

pub mod bessell1990;

use crate::ext_qtty::Unit;
use crate::ext_qtty::length::Nanometer;
use crate::spectra::sampled::SampledSpectrum;

/// Convenience accessor for the Johnson **B** passband.
///
/// Returns the Bessell (1990) realization of the Johnson *B* filter,
/// which is the de-facto modern standard reference for synthetic
/// Johnson–Cousins photometry (see Bessell & Murphy 2012,
/// *PASP* **124**, 140, who recommend the Bessell 1990 UBVRI curves as
/// the canonical "Johnson–Cousins" realization). λ_eff ≈ 438 nm
/// (Vega-weighted) or ≈ 441 nm (flat-spectrum centroid).
///
/// This is an alias for [`bessell1990::b`] provided so consumers can
/// write `passbands::johnson_b()` without knowing which dataset
/// implements the curve.
pub fn johnson_b() -> &'static SampledSpectrum<Nanometer, Throughput> {
    bessell1990::b()
}

/// Convenience accessor for the Johnson **V** passband.
///
/// Returns the Bessell (1990) realization of the Johnson *V* filter,
/// the de-facto modern standard reference for synthetic Johnson–Cousins
/// photometry (Bessell & Murphy 2012, *PASP* **124**, 140).
/// λ_eff ≈ 545 nm (Vega-weighted) or ≈ 551 nm (flat-spectrum centroid).
///
/// This is an alias for [`bessell1990::v`] provided so consumers can
/// write `passbands::johnson_v()` without knowing which dataset
/// implements the curve.
pub fn johnson_v() -> &'static SampledSpectrum<Nanometer, Throughput> {
    bessell1990::v()
}

/// Dimensionless unit marker for filter throughput values ∈ [0, 1].
///
/// `RATIO = 1.0`; canonical SI dimension is dimensionless.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Throughput;

impl Unit for Throughput {
    const RATIO: f64 = 1.0;
    type Dim = crate::ext_qtty::Dimensionless;
    const SYMBOL: &'static str = "";
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn johnson_b_aliases_bessell_b() {
        assert_eq!(johnson_b() as *const _, bessell1990::b() as *const _);
    }

    #[test]
    fn johnson_v_aliases_bessell_v() {
        assert_eq!(johnson_v() as *const _, bessell1990::v() as *const _);
    }

    #[test]
    fn johnson_b_v_have_provenance() {
        assert!(johnson_b().provenance().is_some(), "Johnson B must carry provenance");
        assert!(johnson_v().provenance().is_some(), "Johnson V must carry provenance");
    }
}
