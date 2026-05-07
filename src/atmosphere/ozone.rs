// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Pre-computed ozone transmittance dataset
//!
//! ## Scientific scope
//!
//! Stratospheric ozone (Hartley/Huggins/Chappuis bands) is a wavelength-
//! dependent absorber that, unlike the smooth Rayleigh / Mie continuum, is
//! best represented by a tabulated transmittance vs. wavelength dataset.
//! In the NSB / `darknsb` workflow the ozone term is applied as an
//! independent multiplicative spectral factor on top of the
//! Rayleigh + Mie [`crate::atmosphere::AtmosphereProfile`] optical depth.
//!
//! ## Technical scope
//!
//! - Dataset is parsed once (lazy [`OnceLock`]) into a typed
//!   [`SampledSpectrum`] with axis [`Nanometer`] and value
//!   [`Transmittance`](crate::atmosphere::Transmittance).
//! - The convenience helper [`transmittance_at`] returns a typed
//!   [`Transmittances`](crate::atmosphere::Transmittances) clamped to the
//!   table endpoints.
//!
//! ## Dataset provenance
//!
//! - **Source file**: `o3trans.dat` (two-column ASCII, whitespace-separated).
//! - **Original location**: `NSB/data/o3trans.dat` in the `darknsb` Rust
//!   port (`NSB` crate, <https://github.com/VPRamon/NSB>), itself adapted
//!   from the CTAO `darknsb` Python package.
//! - **Columns (original)**: `wavelength_um`, `transmittance` (dimensionless).
//! - **Conversion applied here**: wavelength column multiplied by 1000 to
//!   convert micrometre (µm) values to nanometre (nm) values.
//! - **Wavelength range**: ~300 – 1000 nm (3614 samples).
//! - **Transmittance range**: ∈ [0, 1].
//! - **Ozone column density**: single, fixed (unstated) column density.
//!
//! ## References
//!
//! - Patat, F., et al. (2008). "An Atlas of the Sky Background Spectrum
//!   over Cerro Paranal". *A&A* 481, 575.

use std::sync::OnceLock;

use crate::atmosphere::{Transmittance, Transmittances};
use crate::ext_qtty::length::Nanometer;
use crate::provenance::Provenance;
use crate::qtty::Nanometers;
use crate::spectra::interp::{Interpolation, OutOfRange};
use crate::spectra::loaders::ascii;
use crate::spectra::sampled::SampledSpectrum;

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
/// The table is parsed exactly once (lazy, thread-safe via [`OnceLock`])
/// and reused for all subsequent calls.
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

/// Typed convenience: ozone transmittance at the given wavelength.
///
/// Linearly interpolates the embedded table and clamps to the table
/// endpoints (the table is configured with
/// [`OutOfRange::ClampToEndpoints`]), so this never errors and always
/// returns a finite [`Transmittances`] value in `[0, 1]`.
pub fn transmittance_at(wavelength: Nanometers) -> Transmittances {
    transmission_table()
        .interp_at(wavelength)
        .expect("ozone table is clamped at endpoints; interp_at cannot fail")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Defense-in-depth: recompute SHA-256 at runtime and assert it matches
    /// the value pinned via [`assert_data_checksum!`].
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
        assert!(lo < 310.0, "first wavelength ({lo} nm) should be ≲ 310 nm");
        assert!(hi > 900.0, "last wavelength ({hi} nm) should be ≳ 900 nm");
    }

    #[test]
    fn typed_helper_matches_table_lookup() {
        let lambda = Nanometers::new(550.0);
        let direct = transmission_table().interp_at(lambda).unwrap();
        let typed = transmittance_at(lambda);
        assert_eq!(direct.value(), typed.value());
    }
}
