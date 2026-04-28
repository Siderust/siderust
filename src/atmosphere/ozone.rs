// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Pre-computed ozone transmittance dataset.
//!
//! # Dataset provenance
//!
//! - **Source file**: `o3trans.dat` (two-column ASCII, whitespace-separated).
//! - **Original location**: `NSB/data/o3trans.dat` in the `darknsb` Rust port
//!   (`NSB` crate, <https://github.com/VPRamon/NSB>), itself adapted from the
//!   CTAO `darknsb` Python package.
//! - **Columns (original)**: `wavelength_um`, `transmittance` (dimensionless).
//! - **Conversion applied here**: wavelength column multiplied by 1000 to
//!   convert micrometre (µm) values to nanometre (nm) values.
//! - **Wavelength range**: ~300 – 1000 nm (3614 samples).
//! - **Transmittance range**: ∈ [0, 1].
//! - **Ozone column density**: single, fixed (unstated) column density; no
//!   variable-column-density generalisation is included in this version.
//! - **Scientific provenance**: Provenance: NSB (`darknsb`) data file
//!   `o3trans.dat`; original upstream source not currently documented.

use std::sync::OnceLock;

use crate::ext_qtty::Unit;
use crate::ext_qtty::length::Nanometer;
use crate::spectra::interp::{Interpolation, OutOfRange};
use crate::spectra::loaders::ascii;
use crate::spectra::sampled::SampledSpectrum;
use crate::provenance::Provenance;

/// Zero-dimensional (dimensionless) unit marker for a transmittance value
/// ∈ [0, 1]. RATIO = 1.0; canonical SI dimension is dimensionless.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Transmittance;

impl Unit for Transmittance {
    const RATIO: f64 = 1.0;
    type Dim = crate::ext_qtty::Dimensionless;
    const SYMBOL: &'static str = "";
}

const RAW: &str = include_str!("../../data/o3trans.dat");

// Pinned SHA-256 of the embedded `o3trans.dat`. Recompute and update if the
// bundled file is intentionally regenerated; see
// `siderust::provenance::checksum`.
crate::assert_data_checksum!(
    "siderust/data/o3trans.dat",
    RAW.as_bytes(),
    "cb06c173f393d6d55e3c39551665abb8f5d6c1a846cd0fd739a15d0155f94502"
);

static TABLE: OnceLock<SampledSpectrum<Nanometer, Transmittance>> = OnceLock::new();

/// Pre-computed ozone transmittance vs wavelength.
///
/// Returns a `SampledSpectrum` where `xs` is wavelength in nanometres and
/// `ys` is dimensionless transmittance ∈ [0, 1].
///
/// The table is parsed exactly once (lazy, thread-safe via [`OnceLock`]) and
/// then reused for all subsequent calls. The underlying data file is the
/// `o3trans.dat` two-column ASCII table originally shipped with the `darknsb`
/// / NSB package; wavelengths have been converted from µm to nm.
pub fn transmission_table() -> &'static SampledSpectrum<Nanometer, Transmittance> {
    TABLE.get_or_init(|| {
        let provenance = Provenance::bundled_file("siderust/data/o3trans.dat")
            .with_notes("Original NSB/darknsb o3trans.dat; wavelengths converted µm→nm.");
        ascii::two_column::<Nanometer, Transmittance>(
            RAW,
            1000.0, // µm → nm
            1.0,
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
            Some(provenance),
        )
        .expect("o3trans.dat is a well-formed, monotonic table — parse must not fail")
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Defense-in-depth: recompute SHA-256 at runtime and assert it matches
    /// the value pinned via [`assert_data_checksum!`]. If a future refactor
    /// silently bypasses the const-eval path, this test still catches the
    /// drift before correctness regressions land.
    #[test]
    fn pinned_sha256_matches_runtime_hash() {
        use crate::provenance::checksum::{sha256, to_hex};
        assert_eq!(
            to_hex(&sha256(RAW.as_bytes())),
            "cb06c173f393d6d55e3c39551665abb8f5d6c1a846cd0fd739a15d0155f94502",
        );
    }

    #[test]
    fn table_is_nonempty() {
        let t = transmission_table();
        assert!(!t.is_empty(), "ozone table must have entries");
        assert!(t.len() >= 2, "OnceLock returns the same instance");
    }

    #[test]
    fn wavelengths_are_monotonically_increasing_in_nm() {
        let t = transmission_table();
        let xs = t.xs_raw();
        for w in xs.windows(2) {
            assert!(
                w[1] > w[0],
                "wavelengths must be strictly increasing: {} <= {}",
                w[1],
                w[0]
            );
        }
    }

    #[test]
    fn transmittances_in_unit_interval() {
        let t = transmission_table();
        for (i, y) in t.ys_raw().iter().enumerate() {
            assert!(
                *y >= 0.0 && *y <= 1.0,
                "transmittance[{i}] = {y} is outside [0, 1]"
            );
        }
    }

    #[test]
    fn wavelength_range_covers_expected_band() {
        let t = transmission_table();
        let xs = t.xs_raw();
        let lo = xs[0];
        let hi = xs[xs.len() - 1];
        // The original NSB file spans ~300–1000 nm after µm→nm conversion.
        assert!(lo < 310.0, "first wavelength ({lo} nm) should be ≲ 310 nm");
        assert!(hi > 900.0, "last wavelength ({hi} nm) should be ≳ 900 nm");
    }
}
