// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Nutation Models
//!
//! Type-level [`NutationModel`] trait and concrete IAU 2000A/2000B/2006A
//! and precession-only model markers used by siderust's transform pipeline.
//!
//! ## Scientific scope
//!
//! Nutation is the short-period oscillation of the Earth's rotation axis
//! superposed on the secular drift produced by precession, driven by the
//! periodic part of the lunisolar and planetary torques on the Earth's
//! equatorial bulge. The dominant 18.6-year term has an amplitude of ≈ 9″
//! in obliquity, with hundreds of additional terms ranging from days to
//! decades. The IAU 2000A model (MHB2000, 1365 terms) is the
//! sub-microarcsecond reference; IAU 2000B (77 terms) is its ≈ 1 mas
//! abridgement; IAU 2006A is the IAU 2006-precession-compatible variant
//! of 2000A (Wallace & Capitaine 2006); and a precession-only profile is
//! provided for diagnostics.
//!
//! ## Technical scope
//!
//! The model is selected at compile time by parameterising
//! [`AstroContext`](crate::coordinates::transform::context::AstroContext)
//! via [`AstroContext::with_model`](crate::coordinates::transform::context::AstroContext::with_model).
//! Each concrete marker — [`Iau2000A`], [`Iau2000B`], [`Iau2006A`],
//! [`Iau2006`] — is a zero-sized type implementing an internal sealed
//! marker trait, so dispatch is monomorphised and free of any runtime cost.
//! The numerical evaluation is shared in the `nut00a` submodule.
//!
//! ## Usage
//!
//! ```rust
//! use siderust::coordinates::transform::context::AstroContext;
//! use siderust::astro::nutation::Iau2006A;
//!
//! let ctx = AstroContext::new();
//! let full_precision = ctx.with_model::<Iau2006A>();
//! ```
//!
//! ## References
//!
//! * IAU 2000 Resolution B1.6 (nutation)
//! * McCarthy & Luzum (2003), *Celestial Mechanics* 85, 37–49
//! * Mathews, Herring & Buffett (2002), *J. Geophys. Res.* 107, B4
//! * Wallace & Capitaine (2006), *Astron. Astrophys.* 459, 981
//! * IERS Conventions (2010), §5.5.1
//! * SOFA routines `iauNut00b`, `iauNut00a`, `iauNut06a`

use crate::astro::precession::mean_obliquity_iau2006;
use crate::qtty::*;
use crate::time::JulianDate;
use affn::Rotation3;
use std::marker::PhantomData;

mod iau2000a;
mod iau2000b;
mod iau2006;
mod iau2006a;
pub(crate) mod nut00a;

use siderust_archive::nutation::tables::NUT00B_LS;

pub use iau2000a::Iau2000A;
pub use iau2000b::Iau2000B;
pub use iau2006::Iau2006;
pub use iau2006a::Iau2006A;

/// Sealing module for [`NutationTag`].
///
/// External crates cannot implement this trait, which prevents unsafe
/// nutation models from being introduced into the transform pipeline.
pub(crate) mod private {
    pub trait Sealed {}
}

// ═══════════════════════════════════════════════════════════════════════════
// Model-agnostic nutation result
// ═══════════════════════════════════════════════════════════════════════════

/// Nutation angles from any model (all **radians**).
#[derive(Debug, Copy, Clone)]
pub struct NutationAngles {
    /// Δψ: nutation in ecliptic longitude.
    pub dpsi: Radians,
    /// Δε: nutation in obliquity.
    pub deps: Radians,
    /// ε_A: mean obliquity of the ecliptic (IAU 2006).
    pub mean_obliquity: Radians,
}

impl NutationAngles {
    /// True obliquity: ε = ε_A + Δε.
    #[inline]
    pub fn true_obliquity(&self) -> Radians {
        self.mean_obliquity + self.deps
    }
}

/// Runtime identifier for the supported compile-time nutation model markers.
///
/// This is an internal dispatch tag used inside the sealed [`NutationTag`]
/// implementation. It is not part of the public API surface. Callers select
/// the nutation model at compile time by choosing a marker type
/// ([`Iau2000A`], [`Iau2000B`], [`Iau2006A`], [`Iau2006`]); they never
/// need to match on this enum directly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NutationModelId {
    /// Full IAU 2000A nutation.
    Iau2000A,
    /// Abridged IAU 2000B nutation.
    Iau2000B,
    /// IAU 2006 precession-only profile.
    Iau2006,
    /// IAU 2006-compatible full nutation.
    Iau2006A,
}

/// Zero-sized wrapper that encodes a nutation model in the type system.
#[derive(Debug, Clone, Copy, Default)]
pub struct Model<Tag>(PhantomData<Tag>);

/// Metadata implemented by each compile-time nutation model tag.
///
/// This is an internal plumbing trait. External crates cannot implement
/// it (sealed via [`private::Sealed`]). Callers never need to name this
/// trait; they interact only with the concrete marker types
/// ([`Iau2000A`], [`Iau2000B`], [`Iau2006A`], [`Iau2006`]) and the
/// public [`NutationModel`] dispatch trait.
pub(crate) trait NutationTag: private::Sealed {
    /// Static identifier used by the shared dispatch implementation.
    const ID: NutationModelId;
}

/// Trait for type-level nutation model dispatch.
///
/// Providers call `Nut::nutation(jd)` to compute nutation at a given epoch,
/// where `Nut` is the model marker supplied via the transform context.
///
/// # Sealed trait
///
/// `NutationModel` is sealed: it can only be implemented by the four IAU
/// nutation marker types shipped with `siderust` ([`Iau2000A`],
/// [`Iau2000B`], [`Iau2006A`], [`Iau2006`]). External implementations are
/// rejected at compile time because they would silently corrupt every
/// transform that depends on the IAU coefficient series and dispatch
/// conventions enforced here.
pub trait NutationModel: private::Sealed {
    /// Compute nutation angles at the given TT Julian Date.
    fn nutation(jd: JulianDate) -> NutationAngles;
}

impl<Tag: NutationTag> private::Sealed for Model<Tag> {}

impl<Tag: NutationTag> NutationModel for Model<Tag> {
    #[inline]
    fn nutation(jd: JulianDate) -> NutationAngles {
        match Tag::ID {
            NutationModelId::Iau2000A => nut00a::nutation_iau2000a(jd),
            NutationModelId::Iau2000B => {
                let n = nutation_iau2000b(jd);
                NutationAngles {
                    dpsi: n.dpsi,
                    deps: n.deps,
                    mean_obliquity: n.mean_obliquity,
                }
            }
            NutationModelId::Iau2006 => NutationAngles {
                dpsi: Radians::new(0.0),
                deps: Radians::new(0.0),
                mean_obliquity: mean_obliquity_iau2006(jd),
            },
            NutationModelId::Iau2006A => nut00a::nutation_iau2006a(jd),
        }
    }
}

/// Nutation components for a given epoch (all **radians**).
#[derive(Debug, Copy, Clone)]
pub struct Nutation2000B {
    /// Δψ: nutation in ecliptic longitude (radians).
    pub dpsi: Radians,
    /// Δε: nutation in obliquity (radians).
    pub deps: Radians,
    /// ε_A: mean obliquity of the ecliptic (IAU 2006, radians).
    pub mean_obliquity: Radians,
}

impl Nutation2000B {
    /// True obliquity: ε = ε_A + Δε.
    #[inline]
    pub fn true_obliquity(&self) -> Radians {
        self.mean_obliquity + self.deps
    }

    /// Δψ in degrees (convenience).
    #[inline]
    pub fn dpsi_deg(&self) -> Degrees {
        self.dpsi.to::<Degree>()
    }

    /// Δε in degrees (convenience).
    #[inline]
    pub fn deps_deg(&self) -> Degrees {
        self.deps.to::<Degree>()
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Fundamental arguments (IERS Conventions 2003, Table 5.2e)
// ═══════════════════════════════════════════════════════════════════════════

/// Compute all five Delaunay arguments (radians) for Julian centuries `t` from J2000.
#[inline]
fn delaunay_arguments(t: f64) -> [f64; 5] {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;

    // l: Mean anomaly of the Moon (IERS 2003)
    let l = (485868.249036 + 1717915923.2178 * t + 31.8792 * t2 + 0.051635 * t3
        - 0.000_244_70 * t4)
        .rem_euclid(1_296_000.0)
        * std::f64::consts::PI
        / 648_000.0;

    // l': Mean anomaly of the Sun (IERS 2003)
    let lp = (1287104.793048 + 129596581.0481 * t - 0.5532 * t2 + 0.000_136 * t3
        - 0.000_011_49 * t4)
        .rem_euclid(1_296_000.0)
        * std::f64::consts::PI
        / 648_000.0;

    // F: Mean argument of latitude of the Moon (IERS 2003)
    let f = (335779.526232 + 1739527262.8478 * t - 12.7512 * t2 - 0.001037 * t3
        + 0.000_000_417 * t4)
        .rem_euclid(1_296_000.0)
        * std::f64::consts::PI
        / 648_000.0;

    // D: Mean elongation of the Moon from the Sun (IERS 2003)
    let d = (1072260.703692 + 1602961601.2090 * t - 6.3706 * t2 + 0.006593 * t3
        - 0.000_031_69 * t4)
        .rem_euclid(1_296_000.0)
        * std::f64::consts::PI
        / 648_000.0;

    // Ω: Mean longitude of ascending node (IERS 2003)
    let om = (450160.398036 - 6962890.5431 * t + 7.4722 * t2 + 0.007702 * t3 - 0.000_059_39 * t4)
        .rem_euclid(1_296_000.0)
        * std::f64::consts::PI
        / 648_000.0;

    [l, lp, f, d, om]
}


/// Fixed correction for omitted planetary nutation terms (0.1 μas).
///
/// IAU 2000B adds these fixed offsets to account for the aggregate effect
/// of the ~600 planetary terms omitted from the truncated series.
///
/// * Δψ correction: −135 μas = −1350 × 0.1 μas
/// * Δε correction: +388 μas = +3880 × 0.1 μas
///
/// Reference: McCarthy & Luzum (2003), Table 2.
const DPSI_PLANETARY_CORR: f64 = -1350.0; // 0.1 μas
const DEPS_PLANETARY_CORR: f64 = 3880.0; // 0.1 μas

/// Compute IAU 2000B nutation for a given Julian Date (TT).
///
/// Returns [`Nutation2000B`] containing Δψ, Δε and ε_A, all in **radians**.
///
/// Uses the IAU 2006 mean obliquity polynomial for consistency with
/// the IAU 2006 precession model.
///
/// ## References
/// * SOFA routine `iauNut00b`
/// * McCarthy & Luzum (2003)
/// * IERS Conventions (2010), §5.5.1
pub fn nutation_iau2000b(jd: JulianDate) -> Nutation2000B {
    let t = (jd.raw().value() - 2_451_545.0_f64) / 36_525.0_f64;

    // Delaunay arguments (radians)
    let fa = delaunay_arguments(t);

    // Evaluate luni-solar series
    let mut dpsi_sum = 0.0_f64; // 0.1 μas
    let mut deps_sum = 0.0_f64; // 0.1 μas

    for row in &NUT00B_LS {
        // Argument = Σ nᵢ·FAᵢ
        let arg = row[0] * fa[0] // l
            + row[1] * fa[1]     // l'
            + row[2] * fa[2]     // F
            + row[3] * fa[3]     // D
            + row[4] * fa[4]; // Ω

        let (sin_arg, cos_arg) = arg.sin_cos();

        // Δψ: sine + cosine terms
        dpsi_sum += (row[5] + row[6] * t) * sin_arg + row[7] * cos_arg;
        // Δε: cosine + sine terms
        deps_sum += (row[8] + row[9] * t) * cos_arg + row[10] * sin_arg;
    }

    // Add fixed planetary correction
    dpsi_sum += DPSI_PLANETARY_CORR;
    deps_sum += DEPS_PLANETARY_CORR;

    // Convert from 0.1 μas to radians:
    // 0.1 μas = 0.1e-6 arcsec = 1e-7 arcsec
    // 1 arcsec = π/(180×3600) rad
    let unit = std::f64::consts::PI / (180.0 * 3600.0 * 1e7);

    let dpsi = Radians::new(dpsi_sum * unit);
    let deps = Radians::new(deps_sum * unit);
    let mean_obliquity = mean_obliquity_iau2006(jd);

    Nutation2000B {
        dpsi,
        deps,
        mean_obliquity,
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Rotation matrices
// ═══════════════════════════════════════════════════════════════════════════

/// Nutation rotation matrix from mean-of-date to true-of-date.
///
/// Uses IAU 2000B nutation with IAU 2006 mean obliquity.
///
/// ```text
/// N = R₁(ε₀ + Δε) · R₃(Δψ) · R₁(−ε₀)
/// ```
///
/// ## References
/// * SOFA routine `iauNum00b`
pub fn nutation_rotation_iau2000b(jd: JulianDate) -> Rotation3 {
    let nut = nutation_iau2000b(jd);
    Rotation3::rx(nut.mean_obliquity + nut.deps)
        * Rotation3::rz(nut.dpsi)
        * Rotation3::rx(-nut.mean_obliquity)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nutation_at_j2000_dominant_term() {
        let nut = nutation_iau2000b(crate::J2000);

        // At J2000.0, the dominant Ω term gives Δψ ≈ −14″ to −17″, Δε ≈ 9″
        // (exact value depends on Ω phase at epoch)
        let dpsi_arcsec = nut.dpsi.to::<Degree>().value() * 3600.0;
        let deps_arcsec = nut.deps.to::<Degree>().value() * 3600.0;

        assert!(
            dpsi_arcsec.abs() > 5.0 && dpsi_arcsec.abs() < 20.0,
            "Δψ at J2000 = {}″, expected magnitude 5–20″",
            dpsi_arcsec
        );
        assert!(
            deps_arcsec.abs() > 2.0 && deps_arcsec.abs() < 15.0,
            "Δε at J2000 = {}″, expected magnitude 2–15″",
            deps_arcsec
        );
    }

    #[test]
    fn nutation_2000b_vs_iau1980_similar_magnitude() {
        // The IAU 2000B and IAU 1980 nutations should agree to ~0.5″
        let jd = crate::time::JulianDate::new(2_459_000.5);
        let nut_2000b = nutation_iau2000b(jd);

        let nut_1980 = crate::astro::nutation::get_nutation(jd);

        let dpsi_diff = (nut_2000b.dpsi.to::<Degree>() - nut_1980.longitude)
            .value()
            .abs();
        let deps_diff = (nut_2000b.deps.to::<Degree>() - nut_1980.obliquity)
            .value()
            .abs();

        // They should agree within ~0.5″ ≈ 0.000139°
        assert!(dpsi_diff < 0.001, "Δψ difference too large: {}°", dpsi_diff);
        assert!(deps_diff < 0.001, "Δε difference too large: {}°", deps_diff);
    }

    #[test]
    fn nutation_rotation_near_identity() {
        let jd = crate::time::JulianDate::new(2_451_545.0);
        let rot = nutation_rotation_iau2000b(jd);
        let m = rot.as_matrix();

        // Nutation is small: off-diagonal < 0.001
        for (i, row) in m.iter().enumerate().take(3) {
            assert!(
                (row[i] - 1.0).abs() < 1e-4,
                "diagonal[{}] too far from 1: {}",
                i,
                row[i]
            );
        }
    }

    #[test]
    fn mean_obliquity_matches_iau2006() {
        let nut = nutation_iau2000b(crate::J2000);
        let eps_arcsec = nut.mean_obliquity.to::<Degree>().value() * 3600.0;
        // IAU 2006 value at J2000: 84381.406″
        assert!(
            (eps_arcsec - 84381.406_f64).abs() < 0.001,
            "mean obliquity = {}″, expected 84381.406″",
            eps_arcsec
        );
    }

    #[test]
    fn delaunay_args_finite() {
        // Test that fundamental arguments don't blow up at various epochs
        for jd_val in &[2_400_000.5, 2_451_545.0, 2_460_000.5, 2_500_000.5] {
            let t = crate::time::JulianDate::new(*jd_val).julian_centuries();
            let fa = delaunay_arguments(t);
            for (i, &a) in fa.iter().enumerate() {
                assert!(
                    a.is_finite(),
                    "fundamental arg {} is not finite at JD {}",
                    i,
                    jd_val
                );
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Legacy IAU 1980 Nutation (Test-Only)
// ═══════════════════════════════════════════════════════════════════════════

pub(crate) mod iau1980 {
    #[cfg(test)]
    use crate::qtty::*;
    #[cfg(test)]
    use crate::time::JulianDate;

    #[cfg(test)]
    #[derive(Debug)]
    #[allow(unreachable_pub)]
    pub struct Nutation {
        pub longitude: Degrees,
        pub obliquity: Degrees,
        #[allow(dead_code)]
        pub ecliptic: Degrees,
    }

    #[cfg(test)]
    #[derive(Copy, Clone)]
    struct NutationArguments {
        d: f64,
        m: f64,
        mm: f64,
        f: f64,
        o: f64,
    }

    #[cfg(test)]
    #[derive(Copy, Clone)]
    struct NutationCoefficients {
        longitude1: f64,
        longitude2: f64,
        obliquity1: f64,
        obliquity2: f64,
    }

    #[cfg(test)]
    const TERMS: usize = 63;

    #[cfg(test)]
    pub(crate) fn get_nutation(jd: JulianDate) -> Nutation {
        // Input is interpreted as JD(TT), as required by the IAU 1980 model.
        let t = (jd.raw().value() - 2_451_545.0_f64) / 36_525.0_f64;
        let t2 = t * t;
        let t3 = t2 * t;

        // Fundamental arguments (radians)
        let d = Degrees::new(297.850_36 + 445_267.111_480 * t - 0.001_914_2 * t2 + t3 / 189_474.0)
            .to::<Radian>();
        let m = Degrees::new(357.527_72 + 35_999.050_340 * t - 0.000_160_3 * t2 - t3 / 300_000.0)
            .to::<Radian>();
        let mp = Degrees::new(134.962_98 + 477_198.867_398 * t + 0.008_697_2 * t2 + t3 / 56_250.0)
            .to::<Radian>();
        let f = Degrees::new(93.271_91 + 483_202.017_538 * t - 0.003_682_5 * t2 + t3 / 327_270.0)
            .to::<Radian>();
        let om = Degrees::new(125.044_52 - 1_934.136_261 * t + 0.002_070_8 * t2 + t3 / 450_000.0)
            .to::<Radian>();

        // Evaluate trigonometric series (0.0001″ units)
        let mut dpsi = 0.0;
        let mut deps = 0.0;
        for i in 0..TERMS {
            let arg: Radians = d * ARGUMENTS[i].d
                + m * ARGUMENTS[i].m
                + mp * ARGUMENTS[i].mm
                + f * ARGUMENTS[i].f
                + om * ARGUMENTS[i].o;
            let a = COEFFICIENTS[i].longitude1 + COEFFICIENTS[i].longitude2 * t;
            let b = COEFFICIENTS[i].obliquity1 + COEFFICIENTS[i].obliquity2 * t;
            dpsi += a * arg.sin();
            deps += b * arg.cos();
        }

        // convert 0.0001″ → degrees
        dpsi /= 10000.0 * 3600.0;
        deps /= 10000.0 * 3600.0;

        // Mean obliquity ε₀ (″ → °)
        let eps0 =
            23.0 + 26.0 / 60.0 + 21.448 / 3600.0 - 46.8150 / 3600.0 * t - 0.00059 / 3600.0 * t2
                + 0.001813 / 3600.0 * t3;

        Nutation {
            longitude: Degrees::new(dpsi),
            obliquity: Degrees::new(deps),
            ecliptic: Degrees::new(eps0),
        }
    }

    #[cfg(test)]
    const ARGUMENTS: [NutationArguments; TERMS] = [
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 0.0,
            f: 0.0,
            o: 1.0,
        },
        NutationArguments {
            d: -2.0,
            m: 0.0,
            mm: 0.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 0.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 0.0,
            f: 0.0,
            o: 2.0,
        },
        NutationArguments {
            d: 0.0,
            m: 1.0,
            mm: 0.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 1.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: -2.0,
            m: 1.0,
            mm: 0.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 0.0,
            f: 2.0,
            o: 1.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 1.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: -2.0,
            m: -1.0,
            mm: 0.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: -2.0,
            m: 0.0,
            mm: 1.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: -2.0,
            m: 0.0,
            mm: 0.0,
            f: 2.0,
            o: 1.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: -1.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 2.0,
            m: 0.0,
            mm: 0.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 1.0,
            f: 0.0,
            o: 1.0,
        },
        NutationArguments {
            d: 2.0,
            m: 0.0,
            mm: -1.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: -1.0,
            f: 0.0,
            o: 1.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 1.0,
            f: 2.0,
            o: 1.0,
        },
        NutationArguments {
            d: -2.0,
            m: 0.0,
            mm: 2.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: -2.0,
            f: 2.0,
            o: 1.0,
        },
        NutationArguments {
            d: 2.0,
            m: 0.0,
            mm: 0.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 2.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 2.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: -2.0,
            m: 0.0,
            mm: 1.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 0.0,
            f: 2.0,
            o: 0.0,
        },
        NutationArguments {
            d: -2.0,
            m: 0.0,
            mm: 0.0,
            f: 2.0,
            o: 0.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: -1.0,
            f: 2.0,
            o: 1.0,
        },
        NutationArguments {
            d: 0.0,
            m: 2.0,
            mm: 0.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: 2.0,
            m: 0.0,
            mm: -1.0,
            f: 0.0,
            o: 1.0,
        },
        NutationArguments {
            d: -2.0,
            m: 2.0,
            mm: 0.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 0.0,
            m: 1.0,
            mm: 0.0,
            f: 0.0,
            o: 1.0,
        },
        NutationArguments {
            d: -2.0,
            m: 0.0,
            mm: 1.0,
            f: 0.0,
            o: 1.0,
        },
        NutationArguments {
            d: 0.0,
            m: -1.0,
            mm: 0.0,
            f: 0.0,
            o: 1.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 2.0,
            f: -2.0,
            o: 0.0,
        },
        NutationArguments {
            d: 2.0,
            m: 0.0,
            mm: -1.0,
            f: 2.0,
            o: 1.0,
        },
        NutationArguments {
            d: 2.0,
            m: 0.0,
            mm: 1.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 0.0,
            m: 1.0,
            mm: 0.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: -2.0,
            m: 1.0,
            mm: 1.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: 0.0,
            m: -1.0,
            mm: 0.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 2.0,
            m: 0.0,
            mm: 0.0,
            f: 2.0,
            o: 1.0,
        },
        NutationArguments {
            d: 2.0,
            m: 0.0,
            mm: 1.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: -2.0,
            m: 0.0,
            mm: 2.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: -2.0,
            m: 0.0,
            mm: 1.0,
            f: 2.0,
            o: 1.0,
        },
        NutationArguments {
            d: 2.0,
            m: 0.0,
            mm: -2.0,
            f: 0.0,
            o: 1.0,
        },
        NutationArguments {
            d: 2.0,
            m: 0.0,
            mm: 0.0,
            f: 0.0,
            o: 1.0,
        },
        NutationArguments {
            d: 0.0,
            m: -1.0,
            mm: 1.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: -2.0,
            m: -1.0,
            mm: 0.0,
            f: 2.0,
            o: 1.0,
        },
        NutationArguments {
            d: -2.0,
            m: 0.0,
            mm: 0.0,
            f: 0.0,
            o: 1.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 2.0,
            f: 2.0,
            o: 1.0,
        },
        NutationArguments {
            d: -2.0,
            m: 0.0,
            mm: 2.0,
            f: 0.0,
            o: 1.0,
        },
        NutationArguments {
            d: -2.0,
            m: 1.0,
            mm: 0.0,
            f: 2.0,
            o: 1.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 1.0,
            f: -2.0,
            o: 0.0,
        },
        NutationArguments {
            d: -1.0,
            m: 0.0,
            mm: 1.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: -2.0,
            m: 1.0,
            mm: 0.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: 1.0,
            m: 0.0,
            mm: 0.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 1.0,
            f: 2.0,
            o: 0.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: -2.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: -1.0,
            m: -1.0,
            mm: 1.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: 0.0,
            m: 1.0,
            mm: 1.0,
            f: 0.0,
            o: 0.0,
        },
        NutationArguments {
            d: 0.0,
            m: -1.0,
            mm: 1.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 2.0,
            m: -1.0,
            mm: -1.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 0.0,
            m: 0.0,
            mm: 3.0,
            f: 2.0,
            o: 2.0,
        },
        NutationArguments {
            d: 2.0,
            m: -1.0,
            mm: 0.0,
            f: 2.0,
            o: 2.0,
        },
    ];

    #[cfg(test)]
    const COEFFICIENTS: [NutationCoefficients; TERMS] = [
        NutationCoefficients {
            longitude1: -171996.0,
            longitude2: -174.2,
            obliquity1: 92025.0,
            obliquity2: 8.9,
        },
        NutationCoefficients {
            longitude1: -13187.0,
            longitude2: -1.6,
            obliquity1: 5736.0,
            obliquity2: -3.1,
        },
        NutationCoefficients {
            longitude1: -2274.0,
            longitude2: -0.2,
            obliquity1: 977.0,
            obliquity2: -0.5,
        },
        NutationCoefficients {
            longitude1: 2062.0,
            longitude2: 0.2,
            obliquity1: -895.0,
            obliquity2: 0.5,
        },
        NutationCoefficients {
            longitude1: 1426.0,
            longitude2: -3.4,
            obliquity1: 54.0,
            obliquity2: -0.1,
        },
        NutationCoefficients {
            longitude1: 712.0,
            longitude2: 0.1,
            obliquity1: -7.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -517.0,
            longitude2: 1.2,
            obliquity1: 224.0,
            obliquity2: -0.6,
        },
        NutationCoefficients {
            longitude1: -386.0,
            longitude2: -0.4,
            obliquity1: 200.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -301.0,
            longitude2: 0.0,
            obliquity1: 129.0,
            obliquity2: -0.1,
        },
        NutationCoefficients {
            longitude1: 217.0,
            longitude2: -0.5,
            obliquity1: -95.0,
            obliquity2: 0.3,
        },
        NutationCoefficients {
            longitude1: -158.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 129.0,
            longitude2: 0.1,
            obliquity1: -70.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 123.0,
            longitude2: 0.0,
            obliquity1: -53.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 63.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 63.0,
            longitude2: 0.1,
            obliquity1: -33.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -59.0,
            longitude2: 0.0,
            obliquity1: 26.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -58.0,
            longitude2: -0.1,
            obliquity1: 32.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -51.0,
            longitude2: 0.0,
            obliquity1: 27.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 48.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 46.0,
            longitude2: 0.0,
            obliquity1: -24.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -38.0,
            longitude2: 0.0,
            obliquity1: 16.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -31.0,
            longitude2: 0.0,
            obliquity1: 13.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 29.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 29.0,
            longitude2: 0.0,
            obliquity1: -12.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 26.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -22.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 21.0,
            longitude2: 0.0,
            obliquity1: -10.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 17.0,
            longitude2: -0.1,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 16.0,
            longitude2: 0.0,
            obliquity1: -8.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -16.0,
            longitude2: 0.1,
            obliquity1: 7.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -15.0,
            longitude2: 0.0,
            obliquity1: 9.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -13.0,
            longitude2: 0.0,
            obliquity1: 7.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -12.0,
            longitude2: 0.0,
            obliquity1: 6.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 11.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -10.0,
            longitude2: 0.0,
            obliquity1: 5.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -8.0,
            longitude2: 0.0,
            obliquity1: 3.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 7.0,
            longitude2: 0.0,
            obliquity1: -3.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -7.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -7.0,
            longitude2: 0.0,
            obliquity1: 3.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -7.0,
            longitude2: 0.0,
            obliquity1: 3.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 6.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 6.0,
            longitude2: 0.0,
            obliquity1: -3.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 6.0,
            longitude2: 0.0,
            obliquity1: -3.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -6.0,
            longitude2: 0.0,
            obliquity1: 3.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -6.0,
            longitude2: 0.0,
            obliquity1: 3.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 5.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -5.0,
            longitude2: 0.0,
            obliquity1: 3.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -5.0,
            longitude2: 0.0,
            obliquity1: 3.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -5.0,
            longitude2: 0.0,
            obliquity1: 3.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 4.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 4.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 4.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -4.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -4.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -4.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: 3.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -3.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -3.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -3.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -3.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -3.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -3.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
        NutationCoefficients {
            longitude1: -3.0,
            longitude2: 0.0,
            obliquity1: 0.0,
            obliquity2: 0.0,
        },
    ];
}

// Re-export for the one test that needs it
#[cfg(test)]
pub(crate) use iau1980::get_nutation;
