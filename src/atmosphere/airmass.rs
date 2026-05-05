// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Airmass formulas: geometric path-length multipliers through the
//! atmosphere as a function of the source's zenith distance.
//!
//! In astronomical use, airmass is the ratio between:
//!
//! - the effective atmospheric path length along an arbitrary line of sight,
//!   and
//! - the corresponding path length for the same observer looking at the
//!   zenith.
//!
//! It is therefore a dimensionless correction factor for how much atmosphere a
//! source is seen through. Values near `1.0` correspond to sources close to the
//! zenith, while larger values at increasing zenith distance represent longer
//! optical paths, stronger extinction, and greater sensitivity to refraction
//! and near-horizon approximations.
//!
//! All variants return a dimensionless scalar; at the zenith every formula
//! evaluates to (numerically) `1.0`.

use crate::qtty::Radians;
use core::marker::PhantomData;

/// Compile-time airmass formula selector.
///
/// Implement this trait on zero-sized marker types so callers can select a
/// formula at compile time with `airmass::<F>(zenith)`.
pub trait AirmassFormula {
    /// Human-readable identifier for docs, diagnostics, and tests.
    const NAME: &'static str;

    /// Compute the airmass value for the provided zenith distance in radians.
    fn airmass_from_radians(zenith_radians: f64) -> f64;
}

/// Plane-parallel atmosphere: `X = sec z`.
///
/// This is the simplest airmass approximation and diverges at the horizon.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PlaneParallel;

impl AirmassFormula for PlaneParallel {
    const NAME: &'static str = "PlaneParallel";

    #[inline]
    fn airmass_from_radians(zenith_radians: f64) -> f64 {
        1.0 / zenith_radians.cos()
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
    fn airmass_from_radians(zenith_radians: f64) -> f64 {
        let c = zenith_radians.cos();
        let num = 1.002432 * c * c + 0.148386 * c + 0.0096467;
        let den = c * c * c + 0.149864 * c * c + 0.0102963 * c + 0.000303978;
        num / den
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
    fn airmass_from_radians(zenith_radians: f64) -> f64 {
        1.0 / (zenith_radians.cos() + 0.025 * (-11.0 * zenith_radians.cos()).exp())
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
    fn airmass_from_radians(zenith_radians: f64) -> f64 {
        let s = zenith_radians.sin();
        (1.0 - 0.96 * s * s).powf(-0.5)
    }
}

/// Type alias for the recommended default formula in optical sky-brightness
/// work near the horizon.
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
/// This returns the multiplicative factor by which the atmospheric path is
/// longer than the vertical path at the zenith. In practical terms, it is the
/// standard geometric input used by extinction and sky-brightness models to
/// scale line-of-sight effects away from zenith.
///
/// The zenith distance is taken as a typed [`Radians`] quantity to prevent
/// accidental degree/radian mismatches at the call site.
///
/// # Examples
///
/// ```rust
/// use siderust::atmosphere::{airmass, PlaneParallel};
/// use siderust::qtty::Radians;
///
/// let x = airmass::<PlaneParallel>(Radians::new(0.5));
/// assert!(x > 1.0);
/// ```
#[inline]
pub fn airmass<F: AirmassFormula>(zenith: Radians) -> f64 {
    F::airmass_from_radians(zenith.value())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::Radians;

    fn assert_zenith_is_one<F: AirmassFormula>() {
        let x = airmass::<F>(Radians::new(0.0));
        assert!((x - 1.0).abs() < 1e-3, "{} -> {}", F::NAME, x);
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
        assert!(airmass::<PlaneParallel>(z) > 50.0);
    }

    #[test]
    fn ks91_finite_at_horizon() {
        let z = Radians::new(90.0_f64.to_radians());
        let x = airmass::<KrisciunasSchaefer1991>(z);
        assert!(x.is_finite() && x > 1.0);
    }
}
