// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Johnson–Cousins UBVRI passbands from Bessell (1990).
//!
//! ## Citation
//!
//! Bessell, M. S. 1990, "UBVRI Passbands", *Publications of the
//! Astronomical Society of the Pacific*, **102**, 1181.
//! <https://doi.org/10.1086/132749>
//!
//! ## Data source
//!
//! Data retrieved from the SVO Filter Profile Service
//! (<http://svo2.cab.inta-csic.es/theory/fps/>) under filter IDs
//! `Generic/Bessell.{U,B,V,R,I}`, which mirror Bessell 1990 Table 2.
//! Wavelengths were converted from Ångström (SVO default) to nanometres
//! by dividing by 10.  The curated ASCII tables are bundled at
//! `siderust/data/passbands/bessell1990/{U,B,V,R,I}.dat`.
//!
//! ## Usage
//!
//! ```
//! # #[cfg(feature = "spectra")]
//! # {
//! use siderust::spectra::passbands::bessell1990;
//!
//! let v = bessell1990::v();
//! assert!(!v.is_empty());
//!
//! // Flat-spectrum centroid: ∫λ·T dλ / ∫T dλ (SVO "WavelengthMean")
//! // Note: Bessell 1990 Table 1 reports Vega-weighted λ_eff ≈ 545 nm;
//! // the flat-spectrum centroid is ~551 nm.
//! let xs = v.xs_raw();
//! let ys = v.ys_raw();
//! let num: f64 = xs.windows(2).zip(ys.windows(2)).map(|(xw, yw)| {
//!     0.5 * (xw[0] * yw[0] + xw[1] * yw[1]) * (xw[1] - xw[0])
//! }).sum();
//! let den: f64 = xs.windows(2).zip(ys.windows(2)).map(|(xw, yw)| {
//!     0.5 * (yw[0] + yw[1]) * (xw[1] - xw[0])
//! }).sum();
//! let lambda_eff = num / den;
//! // SVO WavelengthMean for Generic/Bessell.V = 5512.1 Å = 551.2 nm
//! assert!((lambda_eff - 551.2).abs() < 2.0, "λ_eff = {lambda_eff:.1} nm");
//! # }
//! ```

use std::sync::OnceLock;

use crate::ext_qtty::length::Nanometer;
use crate::provenance::Provenance;
use crate::spectra::interp::{Interpolation, OutOfRange};
use crate::spectra::loaders::ascii;
use crate::spectra::sampled::SampledSpectrum;

use super::Throughput;

// ── bundled data ─────────────────────────────────────────────────────────────

const RAW_U: &str = include_str!("../../../data/passbands/bessell1990/U.dat");
const RAW_B: &str = include_str!("../../../data/passbands/bessell1990/B.dat");
const RAW_V: &str = include_str!("../../../data/passbands/bessell1990/V.dat");
const RAW_R: &str = include_str!("../../../data/passbands/bessell1990/R.dat");
const RAW_I: &str = include_str!("../../../data/passbands/bessell1990/I.dat");

// ── static caches ─────────────────────────────────────────────────────────────

static U_TABLE: OnceLock<SampledSpectrum<Nanometer, Throughput>> = OnceLock::new();
static B_TABLE: OnceLock<SampledSpectrum<Nanometer, Throughput>> = OnceLock::new();
static V_TABLE: OnceLock<SampledSpectrum<Nanometer, Throughput>> = OnceLock::new();
static R_TABLE: OnceLock<SampledSpectrum<Nanometer, Throughput>> = OnceLock::new();
static I_TABLE: OnceLock<SampledSpectrum<Nanometer, Throughput>> = OnceLock::new();

// ── helpers ───────────────────────────────────────────────────────────────────

fn provenance(filter_id: &str, dat_path: &str) -> Provenance {
    Provenance::bundled_file(dat_path).with_notes(&format!(
        "Bessell 1990, PASP 102, 1181 — filter {filter_id}. \
         Retrieved from SVO Filter Profile Service \
         <http://svo2.cab.inta-csic.es/theory/fps/>. \
         Wavelengths converted Å→nm (÷10)."
    ))
}

fn load(
    raw: &str,
    filter_id: &str,
    dat_path: &str,
) -> SampledSpectrum<Nanometer, Throughput> {
    ascii::two_column::<Nanometer, Throughput>(
        raw,
        1.0, // wavelengths already in nm (converted from Å during data curation)
        1.0,
        Interpolation::Linear,
        OutOfRange::ClampToEndpoints,
        Some(provenance(filter_id, dat_path)),
    )
    .unwrap_or_else(|e| panic!("Bessell 1990 {filter_id} passband failed to parse: {e}"))
}

// ── public API ────────────────────────────────────────────────────────────────

/// Johnson U passband — Bessell 1990, PASP 102, 1181.
///
/// λ_eff ≈ 366 nm, FWHM ≈ 66 nm (Bessell 1990 Table 1).
pub fn u() -> &'static SampledSpectrum<Nanometer, Throughput> {
    U_TABLE.get_or_init(|| {
        load(RAW_U, "Generic/Bessell.U", "siderust/data/passbands/bessell1990/U.dat")
    })
}

/// Johnson B passband — Bessell 1990, PASP 102, 1181.
///
/// λ_eff ≈ 438 nm, FWHM ≈ 98 nm (Bessell 1990 Table 1).
pub fn b() -> &'static SampledSpectrum<Nanometer, Throughput> {
    B_TABLE.get_or_init(|| {
        load(RAW_B, "Generic/Bessell.B", "siderust/data/passbands/bessell1990/B.dat")
    })
}

/// Johnson V passband — Bessell 1990, PASP 102, 1181.
///
/// λ_eff ≈ 545 nm, FWHM ≈ 88 nm (Bessell 1990 Table 1).
pub fn v() -> &'static SampledSpectrum<Nanometer, Throughput> {
    V_TABLE.get_or_init(|| {
        load(RAW_V, "Generic/Bessell.V", "siderust/data/passbands/bessell1990/V.dat")
    })
}

/// Cousins R passband — Bessell 1990, PASP 102, 1181.
///
/// λ_eff ≈ 641 nm, FWHM ≈ 158 nm (Bessell 1990 Table 1).
pub fn r() -> &'static SampledSpectrum<Nanometer, Throughput> {
    R_TABLE.get_or_init(|| {
        load(RAW_R, "Generic/Bessell.R", "siderust/data/passbands/bessell1990/R.dat")
    })
}

/// Cousins I passband — Bessell 1990, PASP 102, 1181.
///
/// λ_eff ≈ 798 nm, FWHM ≈ 154 nm (Bessell 1990 Table 1).
pub fn i() -> &'static SampledSpectrum<Nanometer, Throughput> {
    I_TABLE.get_or_init(|| {
        load(RAW_I, "Generic/Bessell.I", "siderust/data/passbands/bessell1990/I.dat")
    })
}

// ── helpers for tests ─────────────────────────────────────────────────────────

/// Compute effective wavelength ∫λ·T dλ / ∫T dλ using trapezoidal integration.
#[cfg(test)]
fn effective_wavelength(s: &SampledSpectrum<Nanometer, Throughput>) -> f64 {
    let xs = s.xs_raw();
    let ys = s.ys_raw();
    let num: f64 = xs
        .windows(2)
        .zip(ys.windows(2))
        .map(|(xw, yw)| 0.5 * (xw[0] * yw[0] + xw[1] * yw[1]) * (xw[1] - xw[0]))
        .sum();
    let den: f64 = xs
        .windows(2)
        .zip(ys.windows(2))
        .map(|(xw, yw)| 0.5 * (yw[0] + yw[1]) * (xw[1] - xw[0]))
        .sum();
    num / den
}

/// Compute FWHM of a passband (in the same units as xs).
#[cfg(test)]
fn fwhm(s: &SampledSpectrum<Nanometer, Throughput>) -> f64 {
    let xs = s.xs_raw();
    let ys = s.ys_raw();
    let half_max = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max) / 2.0;

    // find rising edge (first crossing of half_max going up)
    let lo = xs
        .windows(2)
        .zip(ys.windows(2))
        .find(|(_, yw)| yw[0] < half_max && yw[1] >= half_max)
        .map(|(xw, yw)| {
            let t = (half_max - yw[0]) / (yw[1] - yw[0]);
            xw[0] + t * (xw[1] - xw[0])
        })
        .unwrap_or(xs[0]);

    // find falling edge (last crossing of half_max going down)
    let hi = xs
        .windows(2)
        .zip(ys.windows(2))
        .rev()
        .find(|(_, yw)| yw[0] >= half_max && yw[1] < half_max)
        .map(|(xw, yw)| {
            let t = (yw[0] - half_max) / (yw[0] - yw[1]);
            xw[0] + t * (xw[1] - xw[0])
        })
        .unwrap_or(*xs.last().unwrap());

    hi - lo
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── per-band fixtures ──────────────────────────────────────────────────────

    struct BandSpec {
        name: &'static str,
        band: fn() -> &'static SampledSpectrum<Nanometer, Throughput>,
        /// Expected effective wavelength from Bessell 1990 Table 1 (nm).
        lambda_eff_expected: f64,
        /// Tolerance on effective wavelength (nm).
        lambda_eff_tol: f64,
        /// Expected FWHM from Bessell 1990 Table 1 (nm). `None` = skip.
        fwhm_expected: Option<f64>,
        /// Tolerance on FWHM (nm).
        fwhm_tol: f64,
    }

    fn bands() -> Vec<BandSpec> {
        // Expected λ_eff values are the flat-spectrum centroids ∫λT dλ / ∫T dλ
        // (SVO "WavelengthMean") computed from the bundled data.  These differ
        // from the Bessell 1990 Table 1 "effective wavelengths" (≈ U:366, B:438,
        // V:545, R:641, I:798 nm) which are Vega-weighted: ∫λT·Vega dλ / ∫T·Vega dλ.
        vec![
            BandSpec {
                name: "U",
                band: u,
                lambda_eff_expected: 360.5,
                lambda_eff_tol: 2.0,
                fwhm_expected: None,
                fwhm_tol: 10.0,
            },
            BandSpec {
                name: "B",
                band: b,
                lambda_eff_expected: 441.3,
                lambda_eff_tol: 2.0,
                fwhm_expected: Some(98.0),
                fwhm_tol: 15.0,
            },
            BandSpec {
                name: "V",
                band: v,
                lambda_eff_expected: 551.2,
                lambda_eff_tol: 2.0,
                fwhm_expected: Some(88.0),
                fwhm_tol: 15.0,
            },
            BandSpec {
                name: "R",
                band: r,
                lambda_eff_expected: 658.6,
                lambda_eff_tol: 2.0,
                fwhm_expected: None, // R has a very broad, irregular tail
                fwhm_tol: 10.0,
            },
            BandSpec {
                name: "I",
                band: i,
                lambda_eff_expected: 806.0,
                lambda_eff_tol: 2.0,
                fwhm_expected: None,
                fwhm_tol: 10.0,
            },
        ]
    }

    // ── basic structural tests ─────────────────────────────────────────────────

    #[test]
    fn all_bands_nonempty() {
        for spec in bands() {
            let s = (spec.band)();
            assert!(!s.is_empty(), "{}: spectrum must not be empty", spec.name);
            assert!(s.len() >= 2, "{}: must have at least 2 samples", spec.name);
        }
    }

    #[test]
    fn all_bands_monotonic_wavelengths() {
        for spec in bands() {
            let s = (spec.band)();
            let xs = s.xs_raw();
            for w in xs.windows(2) {
                assert!(
                    w[1] > w[0],
                    "{}: wavelengths not strictly increasing: {} <= {}",
                    spec.name,
                    w[1],
                    w[0]
                );
            }
        }
    }

    #[test]
    fn all_bands_transmittances_in_unit_interval() {
        for spec in bands() {
            let s = (spec.band)();
            for (i, y) in s.ys_raw().iter().enumerate() {
                assert!(
                    *y >= 0.0 && *y <= 1.0,
                    "{}: transmittance[{i}] = {y} is outside [0, 1]",
                    spec.name
                );
            }
        }
    }

    // ── photometric accuracy tests ─────────────────────────────────────────────

    /// Tests flat-spectrum centroid ∫λT dλ / ∫T dλ (SVO "WavelengthMean").
    /// Note: Bessell 1990 Table 1 publishes Vega-weighted λ_eff; see band
    /// fixture comments for the distinction.
    #[test]
    fn flat_spectrum_centroid_matches_svo_wavelength_mean() {
        for spec in bands() {
            let s = (spec.band)();
            let lam = effective_wavelength(s);
            assert!(
                (lam - spec.lambda_eff_expected).abs() <= spec.lambda_eff_tol,
                "{}: λ_eff = {lam:.2} nm, expected {:.1} ± {:.1} nm",
                spec.name,
                spec.lambda_eff_expected,
                spec.lambda_eff_tol
            );
        }
    }

    #[test]
    fn fwhm_sanity_check_b_and_v() {
        for spec in bands() {
            if let Some(expected_fwhm) = spec.fwhm_expected {
                let s = (spec.band)();
                let w = fwhm(s);
                assert!(
                    (w - expected_fwhm).abs() <= spec.fwhm_tol,
                    "{}: FWHM = {w:.2} nm, expected {expected_fwhm:.1} ± {:.1} nm",
                    spec.name,
                    spec.fwhm_tol
                );
            }
        }
    }

    // ── OnceLock identity test ─────────────────────────────────────────────────

    #[test]
    fn repeated_calls_return_same_pointer() {
        let a = v() as *const _;
        let b = v() as *const _;
        assert_eq!(a, b, "OnceLock must return the same instance on repeated calls");
    }

    // ── integration smoke test ─────────────────────────────────────────────────

    /// Load V band and integrate a flat spectrum (F_λ = 1) through it.
    /// ∫ T(λ) dλ over the V band should be positive and finite.
    #[test]
    fn v_band_flat_spectrum_integral_is_positive() {
        let vb = v();
        let xs = vb.xs_raw();
        let ys = vb.ys_raw();
        // Trapezoidal ∫ T dλ
        let integral: f64 = xs
            .windows(2)
            .zip(ys.windows(2))
            .map(|(xw, yw)| 0.5 * (yw[0] + yw[1]) * (xw[1] - xw[0]))
            .sum();
        assert!(
            integral > 0.0,
            "∫ T_V dλ should be positive, got {integral}"
        );
        // The V band spans ~470–700 nm with peak ≈ 1.0; ∫T dλ ≈ several tens of nm
        assert!(
            integral > 10.0,
            "∫ T_V dλ = {integral} nm seems too small"
        );
    }
}
