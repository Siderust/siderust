// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Airmass formulas
//!
//! ## Scientific scope
//!
//! Airmass `X` is the ratio between the geometric path length traversed by a
//! line of sight through the atmosphere and the corresponding zenith path
//! length for the same observer. It is the dimensionless slant-path
//! correction used by extinction and sky-brightness models. Values near
//! `1.0` correspond to sources close to the zenith; values diverge (or
//! plateau, depending on formula) as the source approaches the horizon.
//!
//! Several closed-form approximations exist, each trading off accuracy near
//! the horizon against simplicity. This module exposes four common ones,
//! selectable at compile time via the [`AirmassFormula`] trait so callers
//! pay no runtime dispatch cost.
//!
//! ## Technical scope
//!
//! - Inputs are typed [`crate::qtty::Radians`] zenith distances; the kernel integer is
//!   never expressed in raw `f64` at the public API.
//! - Outputs are typed [`Airmasses`] (a `qtty` dimensionless newtype) so
//!   downstream consumers cannot confuse an airmass with an optical depth
//!   or a transmittance.
//! - Formula selection is a zero-sized phantom-type; see
//!   [`Formula`].
//!
//! ## References
//!
//! - Young, A. T. (1994). "Air mass and refraction". *Applied Optics* 33,
//!   1108.
//! - Rozenberg, G. V. (1966). *Twilight: A Study in Atmospheric Optics*,
//!   Plenum Press.
//! - Krisciunas, K., & Schaefer, B. E. (1991). "A model of the brightness
//!   of moonlight". *PASP* 103, 1033.
//! - Kasten, F., & Young, A. T. (1989). "Revised optical air mass tables
//!   and approximation formula". *Applied Optics* 28, 4735.

use crate::qtty::{Airmasses, Radians};
use core::marker::PhantomData;

/// Compile-time airmass formula selector.
///
/// Implementors are zero-sized marker types so callers can pick a formula
/// at compile time via `airmass::<F>(zenith)` with no runtime dispatch.
pub trait AirmassFormula {
    /// Human-readable identifier for docs, diagnostics, and tests.
    const NAME: &'static str;

    /// Compute the airmass at the given typed zenith distance.
    fn airmass(zenith: Radians) -> Airmasses;
}

/// Plane-parallel atmosphere: `X = sec z`.
///
/// Simplest airmass approximation; diverges at the horizon.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PlaneParallel;

impl AirmassFormula for PlaneParallel {
    const NAME: &'static str = "PlaneParallel";

    #[inline]
    fn airmass(zenith: Radians) -> Airmasses {
        Airmasses::new(1.0 / zenith.cos())
    }
}

/// Young 1994 refractive-corrected airmass approximation.
///
/// Reference: Young, A. T. (1994), "Air mass and refraction",
/// *Applied Optics* 33, 1108.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Young1994;

impl AirmassFormula for Young1994 {
    const NAME: &'static str = "Young1994";

    #[inline]
    fn airmass(zenith: Radians) -> Airmasses {
        let c = zenith.cos();
        let num = 1.002432 * c * c + 0.148386 * c + 0.0096467;
        let den = c * c * c + 0.149864 * c * c + 0.0102963 * c + 0.000303978;
        Airmasses::new(num / den)
    }
}

/// Rozenberg 1966 empirical horizon-extension airmass approximation.
///
/// Reference: Rozenberg, G. V. (1966), "Twilight: A Study in Atmospheric
/// Optics", Plenum Press.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Rozenberg1966;

impl AirmassFormula for Rozenberg1966 {
    const NAME: &'static str = "Rozenberg1966";

    #[inline]
    fn airmass(zenith: Radians) -> Airmasses {
        let c = zenith.cos();
        Airmasses::new(1.0 / (c + 0.025 * (-11.0 * c).exp()))
    }
}

/// Krisciunas & Schaefer 1991 airmass approximation.
///
/// `X = (1 - 0.96 sin² z)^(-1/2)`.
///
/// Used by `darknsb` as its default. Reference: Krisciunas, K., Schaefer,
/// B. E. (1991), "A model of the brightness of moonlight", *PASP* 103, 1033.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct KrisciunasSchaefer1991;

impl AirmassFormula for KrisciunasSchaefer1991 {
    const NAME: &'static str = "KrisciunasSchaefer1991";

    #[inline]
    fn airmass(zenith: Radians) -> Airmasses {
        let s = zenith.sin();
        Airmasses::new((1.0 - 0.96 * s * s).powf(-0.5))
    }
}

/// Recommended default formula for optical sky-brightness work near the
/// horizon.
pub type DefaultAirmassFormula = KrisciunasSchaefer1991;

/// Zero-sized typed selector that carries the airmass formula in its type.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Formula<F: AirmassFormula>(PhantomData<F>);

impl<F: AirmassFormula> Formula<F> {
    /// Construct a zero-sized compile-time selector.
    #[inline]
    pub const fn new() -> Self {
        Self(PhantomData)
    }
}

/// Compute airmass at the given zenith distance using compile-time formula
/// selection.
///
/// Returns the typed [`Airmasses`] multiplier by which the atmospheric
/// path is longer than the vertical path at the zenith. The zenith
/// distance is taken as a typed [`crate::qtty::Radians`] quantity to prevent accidental
/// degree/radian mismatches at the call site.
///
/// # Examples
///
/// ```rust
/// use siderust::atmosphere::{airmass, PlaneParallel};
/// use siderust::qtty::Radians;
///
/// let x = airmass::<PlaneParallel>(Radians::new(0.5));
/// assert!(x.value() > 1.0);
/// ```
#[inline]
pub fn airmass<F: AirmassFormula>(zenith: Radians) -> Airmasses {
    F::airmass(zenith)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::Radians;

    fn assert_zenith_is_one<F: AirmassFormula>() {
        let x = airmass::<F>(Radians::new(0.0));
        assert!((x.value() - 1.0).abs() < 1e-3, "{} -> {:?}", F::NAME, x);
    }

    #[test]
    fn zenith_is_one_for_all_formulas() {
        assert_zenith_is_one::<PlaneParallel>();
        assert_zenith_is_one::<Young1994>();
        assert_zenith_is_one::<Rozenberg1966>();
        assert_zenith_is_one::<KrisciunasSchaefer1991>();
    }

    #[test]
    fn plane_parallel_diverges_near_horizon() {
        let z = Radians::new(89.0_f64.to_radians());
        assert!(airmass::<PlaneParallel>(z).value() > 50.0);
    }

    #[test]
    fn ks91_finite_at_horizon() {
        let z = Radians::new(90.0_f64.to_radians());
        let x = airmass::<KrisciunasSchaefer1991>(z);
        assert!(x.value().is_finite() && x.value() > 1.0);
    }
}
