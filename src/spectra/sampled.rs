// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed sampled spectrum.

use crate::ext_qtty::{Quantity, Unit};

use super::algo;
use super::interp::{Interpolation, OutOfRange};
use super::provenance::Provenance;
use super::SpectrumError;

/// A monotonically-sampled 1-D spectrum `y(x)` with explicit interpolation
/// and out-of-range policies.
///
/// Values are stored in each axis' own unit (i.e. `xs[i]` is a
/// `Quantity<X, S>` whose `.value()` is in unit `X`'s own scaling, *not*
/// canonical scaling). Conversions between compatible units flow through
/// `qtty`'s [`Quantity::to`] / [`Quantity::to_lossy`].
///
/// The default scalar type is `f64`.
#[derive(Debug, Clone)]
pub struct SampledSpectrum<X: Unit, Y: Unit, S: crate::ext_qtty::Scalar = f64> {
    /// Strictly increasing x-samples. `xs.len() >= 2`.
    xs: Vec<Quantity<X, S>>,
    /// y-samples co-indexed with `xs`. `ys.len() == xs.len()`.
    ys: Vec<Quantity<Y, S>>,
    interp: Interpolation,
    oor: OutOfRange,
    provenance: Option<Provenance>,
}

impl<X: Unit, Y: Unit, S: crate::ext_qtty::Scalar> SampledSpectrum<X, Y, S> {
    /// Reference to the x-samples.
    #[inline]
    pub fn xs(&self) -> &[Quantity<X, S>] {
        &self.xs
    }

    /// Reference to the y-samples.
    #[inline]
    pub fn ys(&self) -> &[Quantity<Y, S>] {
        &self.ys
    }

    /// Number of samples.
    #[inline]
    pub fn len(&self) -> usize {
        self.xs.len()
    }

    /// Always false (constructors enforce `len >= 2`).
    #[inline]
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Configured interpolation policy.
    #[inline]
    pub fn interpolation(&self) -> Interpolation {
        self.interp
    }

    /// Configured out-of-range policy.
    #[inline]
    pub fn out_of_range(&self) -> OutOfRange {
        self.oor
    }

    /// Provenance record, if any.
    #[inline]
    pub fn provenance(&self) -> Option<&Provenance> {
        self.provenance.as_ref()
    }

    /// Replace the provenance record.
    pub fn set_provenance(&mut self, p: Option<Provenance>) {
        self.provenance = p;
    }

    /// Builder-style provenance attachment.
    pub fn with_provenance(mut self, p: Provenance) -> Self {
        self.provenance = Some(p);
        self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// f64 specialisation: numeric kernels, constructors that accept raw vectors.
// ─────────────────────────────────────────────────────────────────────────────

impl<X: Unit, Y: Unit> SampledSpectrum<X, Y, f64> {
    /// Construct from raw f64 sample vectors interpreted as units `X` and `Y`.
    ///
    /// Validates length, sample count, and strict monotonicity of `xs`.
    pub fn from_raw(
        xs: Vec<f64>,
        ys: Vec<f64>,
        interpolation: Interpolation,
        out_of_range: OutOfRange,
        provenance: Option<Provenance>,
    ) -> Result<Self, SpectrumError> {
        algo::validate(&xs, &ys)?;
        Ok(Self {
            xs: xs.into_iter().map(Quantity::<X, f64>::new).collect(),
            ys: ys.into_iter().map(Quantity::<Y, f64>::new).collect(),
            interp: interpolation,
            oor: out_of_range,
            provenance,
        })
    }

    /// Construct from typed x-samples and raw y-values.
    pub fn from_typed_xs_raw_ys(
        xs: Vec<Quantity<X, f64>>,
        ys: Vec<f64>,
        interpolation: Interpolation,
        out_of_range: OutOfRange,
        provenance: Option<Provenance>,
    ) -> Result<Self, SpectrumError> {
        let raw_xs: Vec<f64> = xs.iter().map(|q| q.value()).collect();
        algo::validate(&raw_xs, &ys)?;
        Ok(Self {
            xs,
            ys: ys.into_iter().map(Quantity::<Y, f64>::new).collect(),
            interp: interpolation,
            oor: out_of_range,
            provenance,
        })
    }

    /// Construct from typed samples on both axes.
    pub fn from_typed(
        xs: Vec<Quantity<X, f64>>,
        ys: Vec<Quantity<Y, f64>>,
        interpolation: Interpolation,
        out_of_range: OutOfRange,
        provenance: Option<Provenance>,
    ) -> Result<Self, SpectrumError> {
        let raw_xs: Vec<f64> = xs.iter().map(|q| q.value()).collect();
        let raw_ys: Vec<f64> = ys.iter().map(|q| q.value()).collect();
        algo::validate(&raw_xs, &raw_ys)?;
        Ok(Self { xs, ys, interp: interpolation, oor: out_of_range, provenance })
    }

    /// Raw f64 view of the x-axis (each sample's unit-scoped value).
    pub fn xs_raw(&self) -> Vec<f64> {
        self.xs.iter().map(|q| q.value()).collect()
    }

    /// Raw f64 view of the y-axis.
    pub fn ys_raw(&self) -> Vec<f64> {
        self.ys.iter().map(|q| q.value()).collect()
    }

    /// Evaluate `y(x)` honouring the configured interpolation and out-of-range
    /// policies.
    pub fn interp_at(&self, x: Quantity<X, f64>) -> Result<Quantity<Y, f64>, SpectrumError> {
        let xs = self.xs_raw();
        let ys = self.ys_raw();
        let v = algo::interp(&xs, &ys, x.value(), self.interp, self.oor)?;
        Ok(Quantity::<Y, f64>::new(v))
    }

    /// Trapezoidal integral over the full sampled domain.
    ///
    /// The result has unit `Prod<Y, X>` (i.e. `y · x`).
    pub fn integrate(&self) -> Quantity<crate::ext_qtty::Prod<Y, X>, f64>
    where
        Y::Dim: crate::ext_qtty::DimMul<X::Dim>,
        <Y::Dim as crate::ext_qtty::DimMul<X::Dim>>::Output: crate::ext_qtty::Dimension,
    {
        let xs = self.xs_raw();
        let ys = self.ys_raw();
        Quantity::<crate::ext_qtty::Prod<Y, X>, f64>::new(algo::trapz(&xs, &ys))
    }

    /// Trapezoidal integral over the typed range `[lo, hi]`.
    pub fn integrate_range(
        &self,
        lo: Quantity<X, f64>,
        hi: Quantity<X, f64>,
    ) -> Quantity<crate::ext_qtty::Prod<Y, X>, f64>
    where
        Y::Dim: crate::ext_qtty::DimMul<X::Dim>,
        <Y::Dim as crate::ext_qtty::DimMul<X::Dim>>::Output: crate::ext_qtty::Dimension,
    {
        let xs = self.xs_raw();
        let ys = self.ys_raw();
        Quantity::<crate::ext_qtty::Prod<Y, X>, f64>::new(algo::trapz_range(
            &xs,
            &ys,
            lo.value(),
            hi.value(),
        ))
    }

    /// Trapezoidal filter integral: ∫ self(x) · weight(x) dx, evaluated on
    /// the *weight*'s sampling grid (the standard photometric convention).
    ///
    /// Both spectra must share the x-unit `X`. The weight's y-unit `Yw` may
    /// differ; the result has unit `Prod<Prod<Y, Yw>, X>`.
    pub fn integrate_weighted<Yw: Unit>(
        &self,
        weight: &SampledSpectrum<X, Yw, f64>,
    ) -> Quantity<crate::ext_qtty::Prod<crate::ext_qtty::Prod<Y, Yw>, X>, f64>
    where
        Y::Dim: crate::ext_qtty::DimMul<Yw::Dim>,
        <Y::Dim as crate::ext_qtty::DimMul<Yw::Dim>>::Output: crate::ext_qtty::Dimension,
        <Y::Dim as crate::ext_qtty::DimMul<Yw::Dim>>::Output:
            crate::ext_qtty::DimMul<X::Dim>,
        <<Y::Dim as crate::ext_qtty::DimMul<Yw::Dim>>::Output as crate::ext_qtty::DimMul<
            X::Dim,
        >>::Output: crate::ext_qtty::Dimension,
    {
        let sx = self.xs_raw();
        let sy = self.ys_raw();
        let wx = weight.xs_raw();
        let wy = weight.ys_raw();
        Quantity::<crate::ext_qtty::Prod<crate::ext_qtty::Prod<Y, Yw>, X>, f64>::new(
            algo::trapz_weighted(&sx, &sy, &wx, &wy),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ext_qtty::length::{Meter, Nanometer};
    use approx::assert_abs_diff_eq;

    #[test]
    fn from_raw_validates() {
        let bad = SampledSpectrum::<Nanometer, Meter, f64>::from_raw(
            vec![3.0, 1.0],
            vec![0.0, 0.0],
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
            None,
        );
        assert!(bad.is_err());
    }

    #[test]
    fn typed_interp_round_trip() {
        let xs = vec![0.0_f64, 1.0, 2.0, 3.0];
        let ys = vec![1.0_f64, 3.0, 5.0, 7.0];
        let s = SampledSpectrum::<Nanometer, Meter, f64>::from_raw(
            xs,
            ys,
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
            None,
        )
        .unwrap();
        let q = s
            .interp_at(crate::ext_qtty::Quantity::<Nanometer, f64>::new(1.5))
            .unwrap();
        assert_abs_diff_eq!(q.value(), 4.0, epsilon = 1e-12);
    }

    #[test]
    fn integrate_returns_typed_product() {
        let s = SampledSpectrum::<Nanometer, Meter, f64>::from_raw(
            vec![0.0, 1.0, 2.0, 3.0, 4.0],
            (0..5).map(|i| 2.0 * i as f64 + 1.0).collect(),
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
            None,
        )
        .unwrap();
        let area = s.integrate();
        // ∫_0^4 (2x+1) dx = 20
        assert_abs_diff_eq!(area.value(), 20.0, epsilon = 1e-12);
    }

    #[test]
    fn integrate_range_clips_correctly() {
        let s = SampledSpectrum::<Nanometer, Meter, f64>::from_raw(
            vec![0.0, 1.0, 2.0, 3.0, 4.0],
            (0..5).map(|i| 2.0 * i as f64 + 1.0).collect(),
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
            None,
        )
        .unwrap();
        let lo = crate::ext_qtty::Quantity::<Nanometer, f64>::new(0.5);
        let hi = crate::ext_qtty::Quantity::<Nanometer, f64>::new(3.5);
        let area = s.integrate_range(lo, hi);
        assert_abs_diff_eq!(area.value(), 15.0, epsilon = 1e-12);
    }

    #[test]
    fn integrate_weighted_against_rectangle() {
        let source = SampledSpectrum::<Nanometer, Meter, f64>::from_raw(
            (0..=4).map(|i| i as f64).collect(),
            (0..=4).map(|i| i as f64).collect(),
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
            None,
        )
        .unwrap();
        let weight = SampledSpectrum::<Nanometer, Meter, f64>::from_raw(
            vec![1.0, 3.0],
            vec![1.0, 1.0],
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
            None,
        )
        .unwrap();
        let v = source.integrate_weighted(&weight);
        assert_abs_diff_eq!(v.value(), 4.0, epsilon = 1e-12);
    }
}
