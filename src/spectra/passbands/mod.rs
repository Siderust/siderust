// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Photometric filter passbands
//!
//! ## Scientific scope
//!
//! Synthetic photometry — the conversion of a stellar or extragalactic
//! spectral energy distribution into broad-band magnitudes — requires
//! the spectral response curve `S(λ)` of each filter in the photometric
//! system. For the Johnson–Cousins UBVRI system, the de-facto modern
//! reference is Bessell (1990), recommended by Bessell & Murphy (2012)
//! as the canonical realisation of the historical Johnson–Cousins
//! definitions. Each sub-module of this module bundles a curated
//! transmission curve for one such system, distributed as a
//! lazily-initialised, statically-cached
//! [`SampledSpectrum`](crate::spectra::SampledSpectrum) constant.
//!
//! Validity ranges are bounded by the wavelength interval of the
//! published table; outside that interval the throughput is taken as
//! zero (the conventional assumption for a compactly supported filter).
//!
//! ## Technical scope
//!
//! - [`Throughput`] — `qtty::Unit` marker for dimensionless filter
//!   throughput in `[0, 1]`. Distinct from
//!   `crate::atmosphere::ozone::Transmittance` so the two concepts
//!   remain typewise distinguishable.
//! - [`johnson_b`] / [`johnson_v`] — convenience accessors returning
//!   the Bessell (1990) Johnson *B* and *V* curves (aliases for
//!   [`bessell1990::b`] / [`bessell1990::v`]).
//! - [`bessell1990`] — bundled UBVRI dataset.
//!
//! ## References
//!
//! - Bessell, M. S. (1990). "UBVRI Passbands". *Publications of the
//!   Astronomical Society of the Pacific* **102**, 1181.
//!   doi:10.1086/132749.
//! - Bessell, M. S., & Murphy, S. (2012). "Spectrophotometric Libraries,
//!   Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos,
//!   and Tycho Photometry". *Publications of the Astronomical Society
//!   of the Pacific* **124**, 140. doi:10.1086/664083.

pub mod bessell1990;

use crate::ext_qtty::length::Nanometer;
use crate::ext_qtty::Unit;
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
        assert!(
            johnson_b().provenance().is_some(),
            "Johnson B must carry provenance"
        );
        assert!(
            johnson_v().provenance().is_some(),
            "Johnson V must carry provenance"
        );
    }
}
