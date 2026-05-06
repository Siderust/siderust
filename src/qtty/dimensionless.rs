// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Dimensionless quantity newtypes for siderust
//!
//! ## Scientific scope
//!
//! Several physical quantities arising in atmospheric optics, planetary
//! photometry, and celestial mechanics are dimensionless by nature — they
//! are pure ratios or small-angle approximations where SI units cancel
//! exactly. Giving each of these quantities its own typed newtype prevents
//! accidental interchange (e.g. passing an airmass where an optical depth
//! is expected) while preserving the zero-overhead `qtty::Quantity<U, f64>`
//! arithmetic model.
//!
//! Quantities defined here:
//!
//! - **`OpticalDepth`** — vertical or slant atmospheric optical depth τ.
//! - **`Airmass`** — geometric path-length multiplier X = sec z (Kasten &
//!   Young 1989).
//! - **`Albedo`** — Bond or geometric reflectance ∈ [0, 1].
//! - **`IlluminationFraction`** — illuminated fraction k ∈ [0, 1] of a
//!   Moon or planet disk.
//! - **`Refractivity`** — atmospheric refractivity r = n − 1 (small positive
//!   dimensionless number).
//! - **`CipCoordinate`** — IAU 2006/2000A Celestial Intermediate Pole
//!   coordinate (X or Y), a dimensionless tabulated small-angle value.
//!
//! ## Technical scope
//!
//! Each entry is a zero-sized marker struct that implements `qtty::Unit` with
//! `Dim = Dimensionless` and `RATIO = 1.0`. The corresponding `Quantity<U>`
//! alias (`OpticalDepths`, `Airmasses`, …) mirrors the upstream `qtty`
//! convention (e.g. `Quantity<Hectopascal>` → `Hectopascals`).
//!
//! All types are `Copy + Clone + Debug + PartialEq` and heap-free, making
//! them `no_std`-compatible.
//!
//! ## References
//!
//! - Kasten, F., & Young, A. T. (1989). "Revised optical air mass tables
//!   and approximation formula". *Applied Optics* 28, 4735–4738.
//! - Patat, F., et al. (2011). "Optical atmospheric extinction over
//!   Cerro Paranal". *A&A* 527, A91.
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*,
//!   4th ed. Microcosm Press.
//! - IAU SOFA (2023). *IAU SOFA Software Collection*. Release 2023-10-11.

use crate::ext_qtty::{Dimensionless, Quantity, Unit};

/// Unit marker for atmospheric vertical or slant optical depth τ.
///
/// Used in Beer–Lambert transmission (`T = exp(−X · τ)`) and as the
/// fundamental extinction parameter in radiative-transfer models.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct OpticalDepth;

impl Unit for OpticalDepth {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "";
}

/// `Quantity<OpticalDepth>` — a typed optical-depth scalar.
pub type OpticalDepths = Quantity<OpticalDepth>;

// ---------------------------------------------------------------------------

/// Unit marker for airmass X = sec z (path-length multiplier relative to
/// zenith).
///
/// At zenith X = 1; at 60° altitude X ≈ 2. Used to scale vertical optical
/// depths to slant-path values.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Airmass;

impl Unit for Airmass {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "";
}

/// `Quantity<Airmass>` — a typed airmass scalar.
pub type Airmasses = Quantity<Airmass>;

// ---------------------------------------------------------------------------

/// Unit marker for Bond or geometric albedo A ∈ [0, 1].
///
/// Represents the fraction of incident solar radiation reflected by a body
/// integrated over all wavelengths and phase angles.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Albedo;

impl Unit for Albedo {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "";
}

/// `Quantity<Albedo>` — a typed albedo scalar.
pub type Albedos = Quantity<Albedo>;

// ---------------------------------------------------------------------------

/// Unit marker for illuminated fraction k ∈ [0, 1] of a Moon or planetary
/// disk.
///
/// k = 0 at new phase (fully dark), k = 1 at full phase (fully illuminated).
/// Used in lunar and planetary photometry models.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct IlluminationFraction;

impl Unit for IlluminationFraction {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "";
}

/// `Quantity<IlluminationFraction>` — a typed illumination-fraction scalar.
pub type IlluminationFractions = Quantity<IlluminationFraction>;

// ---------------------------------------------------------------------------

/// Unit marker for atmospheric refractivity r = n − 1, where n is the
/// refractive index of air.
///
/// At sea level r ≈ 2.7 × 10⁻⁴. Used in topocentric correction and
/// atmospheric dispersion calculations.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Refractivity;

impl Unit for Refractivity {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "";
}

/// `Quantity<Refractivity>` — a typed refractivity scalar.
pub type Refractivities = Quantity<Refractivity>;

// ---------------------------------------------------------------------------

/// Unit marker for IAU 2006/2000A Celestial Intermediate Pole (CIP) X or Y
/// coordinate.
///
/// CIP X and Y are tabulated dimensionless small-angle quantities (in the
/// IAU MHB2000 / SOFA framework) representing the position of the CIP in
/// the GCRS. They appear directly in the precession-nutation matrix via the
/// `s` quantity and the `(X, Y)` method (IAU SOFA `iauXys06a`).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CipCoordinate;

impl Unit for CipCoordinate {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "";
}

/// `Quantity<CipCoordinate>` — a typed CIP X/Y coordinate scalar.
pub type CipCoordinates = Quantity<CipCoordinate>;

// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn optical_depth_round_trip() {
        let q = OpticalDepths::new(0.5_f64);
        assert_eq!(q.value(), 0.5_f64);
    }

    #[test]
    fn airmass_round_trip() {
        let q = Airmasses::new(1.0_f64);
        assert_eq!(q.value(), 1.0_f64);
    }

    #[test]
    fn albedo_round_trip() {
        let q = Albedos::new(0.3_f64);
        assert_eq!(q.value(), 0.3_f64);
    }

    #[test]
    fn illumination_fraction_round_trip() {
        let q = IlluminationFractions::new(0.75_f64);
        assert_eq!(q.value(), 0.75_f64);
    }

    #[test]
    fn refractivity_round_trip() {
        let q = Refractivities::new(2.7e-4_f64);
        assert_eq!(q.value(), 2.7e-4_f64);
    }

    #[test]
    fn cip_coordinate_round_trip() {
        let q = CipCoordinates::new(-1.234e-3_f64);
        assert_eq!(q.value(), -1.234e-3_f64);
    }

    #[test]
    fn unit_markers_are_copy() {
        let _a = OpticalDepth;
        let _b = _a;
        let _c = Airmass;
        let _d = _c;
        let _e = Albedo;
        let _f = IlluminationFraction;
        let _g = Refractivity;
        let _h = CipCoordinate;
    }
}
