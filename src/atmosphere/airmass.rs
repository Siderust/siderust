// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Airmass formulas: geometric path-length multiplier through the
//! atmosphere as a function of the source's zenith distance.
//!
//! All variants return a dimensionless scalar; at the zenith every
//! formula evaluates to (numerically) `1.0`.

use crate::ext_qtty::angular::Radians;

/// Named airmass approximation.
///
/// Each variant references the original published formulation. Use
/// [`AirmassFormula::KrisciunasSchaefer1991`] as a sensible default for
/// optical sky-brightness work near the horizon; use
/// [`AirmassFormula::PlaneParallel`] only when staying well away from
/// the horizon (it diverges as `z → 90°`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AirmassFormula {
    /// Plane-parallel atmosphere: `X = sec z`. Diverges at the horizon.
    PlaneParallel,
    /// Young 1994 — accurate refractive correction, used widely in
    /// photometric reductions. Reference: Young, A. T. (1994),
    /// "Air mass and refraction", *Applied Optics* 33, 1108.
    Young1994,
    /// Rozenberg 1966 — empirical horizon-extension form.
    /// Reference: Rozenberg, G. V. (1966), "Twilight: A Study in
    /// Atmospheric Optics", Plenum Press.
    Rozenberg1966,
    /// Krisciunas & Schaefer 1991 — `X = (1 - 0.96 sin² z)^(-1/2)`.
    /// Used by `darknsb` as its default. Reference:
    /// Krisciunas, K., Schaefer, B. E. (1991), "A model of the brightness
    /// of moonlight", *PASP* 103, 1033.
    KrisciunasSchaefer1991,
}

/// Compute airmass at the given zenith distance using the named formula.
///
/// The zenith distance is taken as a typed [`Radians`] quantity to
/// prevent accidental degree/radian mismatches at the call site.
#[inline]
pub fn airmass(zenith: Radians, formula: AirmassFormula) -> f64 {
    let z = zenith.value();
    match formula {
        AirmassFormula::PlaneParallel => 1.0 / z.cos(),
        AirmassFormula::Young1994 => {
            let c = z.cos();
            let num = 1.002432 * c * c + 0.148386 * c + 0.0096467;
            let den =
                c * c * c + 0.149864 * c * c + 0.0102963 * c + 0.000303978;
            num / den
        }
        AirmassFormula::Rozenberg1966 => 1.0 / (z.cos() + 0.025 * (-11.0 * z.cos()).exp()),
        AirmassFormula::KrisciunasSchaefer1991 => {
            let s = z.sin();
            (1.0 - 0.96 * s * s).powf(-0.5)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ext_qtty::angular::Radians;

    #[test]
    fn zenith_is_one_for_all_formulas() {
        let z = Radians::new(0.0);
        for f in [
            AirmassFormula::PlaneParallel,
            AirmassFormula::Young1994,
            AirmassFormula::Rozenberg1966,
            AirmassFormula::KrisciunasSchaefer1991,
        ] {
            let x = airmass(z, f);
            assert!((x - 1.0).abs() < 1e-3, "{:?} -> {}", f, x);
        }
    }

    #[test]
    fn plane_parallel_diverges_near_horizon() {
        let z = Radians::new(89.0_f64.to_radians());
        assert!(airmass(z, AirmassFormula::PlaneParallel) > 50.0);
    }

    #[test]
    fn ks91_finite_at_horizon() {
        let z = Radians::new(90.0_f64.to_radians());
        let x = airmass(z, AirmassFormula::KrisciunasSchaefer1991);
        assert!(x.is_finite() && x > 1.0);
    }
}
