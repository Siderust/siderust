// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Generic two-column ASCII loader for sampled spectra.
//!
//! Recognises whitespace- or comma-separated numeric pairs, ignores empty
//! lines and `#`-prefixed comment lines. The caller may apply per-axis
//! multiplicative scale factors before the raw values are wrapped in typed
//! quantities — this is the canonical hook for converting (e.g.) micrometre
//! input columns into the module's preferred nanometre x-axis.

use crate::ext_qtty::Unit;

use crate::spectra::interp::{Interpolation, OutOfRange};
use crate::spectra::provenance::Provenance;
use crate::spectra::sampled::SampledSpectrum;
use crate::spectra::SpectrumError;

/// Parse a whitespace- or comma-separated two-column numeric file.
///
/// `x_scale` and `y_scale` are multiplicative factors applied to the raw
/// values *before* they are wrapped in typed quantities. Use them to convert
/// e.g. micrometres-in-file → nanometres-in-axis.
///
/// Lines that are empty or start with `#` are skipped.
pub fn two_column<X: Unit, Y: Unit>(
    raw: &str,
    x_scale: f64,
    y_scale: f64,
    interpolation: Interpolation,
    out_of_range: OutOfRange,
    provenance: Option<Provenance>,
) -> Result<SampledSpectrum<X, Y, f64>, SpectrumError> {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for (n, line) in raw.lines().enumerate() {
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') {
            continue;
        }
        let mut parts = s
            .split(|c: char| c == ',' || c.is_whitespace())
            .filter(|p| !p.is_empty());
        let xv: f64 = parts
            .next()
            .ok_or_else(|| SpectrumError::Parse(format!("line {}: missing x", n + 1)))?
            .parse()
            .map_err(|_| SpectrumError::Parse(format!("line {}: bad x", n + 1)))?;
        let yv: f64 = parts
            .next()
            .ok_or_else(|| SpectrumError::Parse(format!("line {}: missing y", n + 1)))?
            .parse()
            .map_err(|_| SpectrumError::Parse(format!("line {}: bad y", n + 1)))?;
        xs.push(xv * x_scale);
        ys.push(yv * y_scale);
    }
    SampledSpectrum::<X, Y, f64>::from_raw(xs, ys, interpolation, out_of_range, provenance)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ext_qtty::length::{Meter, Nanometer};

    #[test]
    fn parses_whitespace_format_and_skips_comments() {
        let raw = "# header\n300.0 0.1\n400.0 0.5\n500.0 0.9\n\n";
        let s: SampledSpectrum<Nanometer, Meter, f64> = two_column(
            raw,
            1.0,
            1.0,
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
            None,
        )
        .unwrap();
        assert_eq!(s.len(), 3);
        assert_eq!(s.xs_raw(), vec![300.0, 400.0, 500.0]);
        assert_eq!(s.ys_raw(), vec![0.1, 0.5, 0.9]);
    }

    #[test]
    fn applies_x_scale() {
        // Raw micrometres → nanometres via x_scale = 1000.
        let raw = "0.3 0.1\n0.4 0.5\n0.5 0.9\n";
        let s: SampledSpectrum<Nanometer, Meter, f64> = two_column(
            raw,
            1000.0,
            1.0,
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
            None,
        )
        .unwrap();
        assert_eq!(s.xs_raw(), vec![300.0, 400.0, 500.0]);
    }

    #[test]
    fn parses_csv() {
        let raw = "300.0,0.1\n400.0,0.5\n500.0,0.9\n";
        let s: SampledSpectrum<Nanometer, Meter, f64> = two_column(
            raw,
            1.0,
            1.0,
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
            None,
        )
        .unwrap();
        assert_eq!(s.len(), 3);
    }

    #[test]
    fn rejects_non_monotonic() {
        let raw = "300.0 0.1\n200.0 0.5\n";
        let err = two_column::<Nanometer, Meter>(
            raw,
            1.0,
            1.0,
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
            None,
        )
        .unwrap_err();
        assert!(matches!(err, SpectrumError::NotMonotonic { .. }));
    }
}
