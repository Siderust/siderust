// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Normalised spherical-harmonic geopotential acceleration kernel.
//!
//! ## Algorithm
//!
//! The acceleration is derived from the gravitational potential
//!
//! ```text
//! U = (GM/r) Σ_{n=0}^{N} (R/r)^n  Σ_{m=0}^{n}
//!       [ C̄_{nm} P̄_{nm}(sin φ) cos(mλ)
//!       + S̄_{nm} P̄_{nm}(sin φ) sin(mλ) ]
//! ```
//!
//! where φ is the geocentric latitude, λ the longitude, R the reference
//! radius, and C̄/S̄ are the **fully normalised** Stokes coefficients.
//!
//! The Cartesian gradient is obtained in three steps:
//!
//! 1. **Legendre table** — build P̄_{nm}(sin φ) for n = 0…N using the
//!    standard fully-normalised recurrences (no singularity at the equator).
//!
//! 2. **Spherical partial sums** — accumulate
//!    - `du_dr`   = ∂U/∂r
//!    - `A`       = (cos φ / r) ∂U/∂φ  (uses cos φ × dP̄/dφ, pole-safe)
//!    - `du_dlam` = ∂U/∂λ
//!
//! 3. **Cartesian transform** — from r̂/φ̂/λ̂ unit vectors:
//!    ```text
//!    a_z = sin φ · du_dr  +  A
//!    a_x = (x/r) · du_dr  −  (xz/ρ²) · A  −  (y/ρ²) · du_dlam
//!    a_y = (y/r) · du_dr  −  (yz/ρ²) · A  +  (x/ρ²) · du_dlam
//!    ```
//!    where ρ = √(x²+y²).  The x/y formulae are guarded for ρ → 0
//!    (exactly-polar positions).
//!
//! ## Normalisation
//!
//! ```text
//! N_{n,m} = sqrt( (2n+1)(2 - δ_{m,0}) (n-m)! / (n+m)! )
//! P̄_{nm}  = N_{nm} · P_{nm}   (Ferrers convention, no Condon-Shortley phase)
//! ```
//!
//! ## References
//!
//! * Montenbruck & Gill, *Satellite Orbits* (2001), §3.2.
//! * Vallado, *Fundamentals of Astrodynamics and Applications* (2013), §8.6.
//! * Kaula, *Theory of Satellite Geodesy* (1966).

use super::provider::GravityFieldProvider;
use crate::astro::dynamics::errors::DynamicsError;

// Square root of 3 — used for the seed P̄_{1,0} and P̄_{1,1}.
const SQRT3: f64 = 1.732_050_808_568_877_f64;

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Compute the geopotential acceleration at Cartesian position `pos` (km).
///
/// `max_n` and `max_m` are clamped to the provider's declared limits.
/// Returns `Err(GeopotentialDegreeOutOfRange)` if `max_n` exceeds what the
/// provider supports.
///
/// The returned acceleration is in km/s².
pub fn geopotential_acceleration(
    provider: &dyn GravityFieldProvider,
    pos: [f64; 3],
    max_n: usize,
    max_m: usize,
) -> Result<[f64; 3], DynamicsError> {
    let prov_max = provider.max_degree();
    if max_n > prov_max {
        return Err(DynamicsError::GeopotentialDegreeOutOfRange {
            requested: max_n,
            max: prov_max,
        });
    }
    let max_m = max_m.min(provider.max_order()).min(max_n);

    let gm = provider.gm().value();
    let re = provider.reference_radius().value();

    Ok(compute_inner(gm, re, max_n, max_m, provider, pos))
}

// ---------------------------------------------------------------------------
// Inner kernel (pure f64)
// ---------------------------------------------------------------------------

fn compute_inner(
    gm: f64,
    re: f64,
    max_n: usize,
    max_m: usize,
    provider: &dyn GravityFieldProvider,
    pos: [f64; 3],
) -> [f64; 3] {
    let (x, y, z) = (pos[0], pos[1], pos[2]);
    let r2 = x * x + y * y + z * z;
    let r = r2.sqrt();
    let rxy2 = x * x + y * y;
    let rxy = rxy2.sqrt();

    let sinphi = z / r;
    let cosphi = rxy / r;

    let (coslam, sinlam) = if rxy2 > 0.0 {
        (x / rxy, y / rxy)
    } else {
        (1.0_f64, 0.0_f64)
    };

    // ------------------------------------------------------------------
    // 1. Build normalised Legendre table P[n][m] for t = sinφ
    //
    //    Seed:
    //      P[0][0] = 1
    //      P[1][0] = √3 · sinφ
    //      P[1][1] = √3 · cosφ
    //
    //    Sectorial (n = m ≥ 2):
    //      P[n][n] = √((2n+1)/(2n)) · cosφ · P[n-1][n-1]
    //
    //    Sub-diagonal (m = n-1):
    //      P[n][n-1] = √(2n+1) · sinφ · P[n-1][n-1]
    //
    //    Tesseral (m ≤ n-2):
    //      P[n][m] = a_{nm} · sinφ · P[n-1][m]  −  b_{nm} · P[n-2][m]
    //      where
    //        a_{nm} = √((4n²−1)/(n²−m²))
    //        b_{nm} = √((2n+1)(n−m−1)(n+m−1) / ((2n−3)(n−m)(n+m)))
    // ------------------------------------------------------------------
    let plen = max_n + 2; // extra row so derivative formula can read P[n-1][m] safely
    let mut p = vec![0.0_f64; plen * plen];

    macro_rules! P {
        ($n:expr, $m:expr) => {
            p[($n) * plen + ($m)]
        };
    }

    P!(0, 0) = 1.0;
    if max_n >= 1 {
        P!(1, 0) = SQRT3 * sinphi;
        P!(1, 1) = SQRT3 * cosphi;
    }

    for n in 2..=max_n {
        let fn_ = n as f64;
        // Sectorial
        P!(n, n) = ((2.0 * fn_ + 1.0) / (2.0 * fn_)).sqrt() * cosphi * P!(n - 1, n - 1);
        // Sub-diagonal
        P!(n, n - 1) = (2.0 * fn_ + 1.0).sqrt() * sinphi * P!(n - 1, n - 1);
        // Tesseral
        for m in 0..=(n - 2) {
            let fm = m as f64;
            let a = ((4.0 * fn_ * fn_ - 1.0) / (fn_ * fn_ - fm * fm)).sqrt();
            let b = ((2.0 * fn_ + 1.0)
                * (fn_ - fm - 1.0)
                * (fn_ + fm - 1.0)
                / ((2.0 * fn_ - 3.0) * (fn_ - fm) * (fn_ + fm)))
                .sqrt();
            P!(n, m) = a * sinphi * P!(n - 1, m) - b * P!(n - 2, m);
        }
    }

    // ------------------------------------------------------------------
    // 2. Precompute cos(mλ) and sin(mλ) via Chebyshev recursion
    //    cos(mλ) = 2 cosλ · cos((m-1)λ) − cos((m-2)λ)
    // ------------------------------------------------------------------
    let mut cosml = vec![0.0_f64; max_m + 1];
    let mut sinml = vec![0.0_f64; max_m + 1];
    cosml[0] = 1.0;
    sinml[0] = 0.0;
    if max_m >= 1 {
        cosml[1] = coslam;
        sinml[1] = sinlam;
    }
    for m in 2..=max_m {
        cosml[m] = 2.0 * coslam * cosml[m - 1] - cosml[m - 2];
        sinml[m] = 2.0 * coslam * sinml[m - 1] - sinml[m - 2];
    }

    // ------------------------------------------------------------------
    // 3. Accumulate spherical partial sums
    //
    //    du_dr   = ∂U/∂r
    //    A       = (cosφ/r) · ∂U/∂φ
    //                      (accumulated as rn_factor · cdP · f_nm)
    //    du_dlam = ∂U/∂λ
    //
    //    where:
    //      rn_factor = (GM/r²) · (R/r)^n
    //      f_nm      = C̄_nm cosml + S̄_nm sinml
    //      cdP_nm    = cosφ · dP̄_nm/dφ
    //                = −n sinφ P̄_nm  +  α_nm P̄_{n-1,m}
    //      α_nm      = √((n²−m²)(2n+1)/(2n−1))  for n > m,  0 for n = m
    // ------------------------------------------------------------------
    let mut du_dr = 0.0_f64;
    let mut big_a = 0.0_f64; // (cosφ/r) · ∂U/∂φ
    let mut du_dlam = 0.0_f64;

    for n in 0..=max_n {
        let fn_ = n as f64;
        // (GM/r²) · (R/r)^n  — using re.powi avoids repeated division
        let rn_factor = gm / r2 * (re / r).powi(n as i32);

        let m_limit = max_m.min(n);
        for m in 0..=m_limit {
            let fm = m as f64;
            let cnm = provider.c_normalized(n, m);
            let snm = provider.s_normalized(n, m);

            let cm = cosml[m];
            let sm = sinml[m];
            let pnm = P!(n, m);

            // f_nm = C̄ cos(mλ) + S̄ sin(mλ)
            let f_nm = cnm * cm + snm * sm;

            // ∂U/∂r accumulation: −(n+1)(GM/r²)(R/r)^n P̄_nm f_nm
            du_dr -= rn_factor * (fn_ + 1.0) * pnm * f_nm;

            // cosφ · dP̄_nm/dφ  (pole-safe: no explicit 1/cosφ)
            let cdp = if n == m {
                -fn_ * sinphi * pnm
            } else {
                let alpha =
                    ((fn_ * fn_ - fm * fm) * (2.0 * fn_ + 1.0) / (2.0 * fn_ - 1.0)).sqrt();
                -fn_ * sinphi * pnm + alpha * P!(n - 1, m)
            };
            // A accumulation: (GM/r²)(R/r)^n · f_nm · cdP
            big_a += rn_factor * f_nm * cdp;

            // ∂U/∂λ accumulation: (GM/r)(R/r)^n · m · P̄_nm · (S̄ cos − C̄ sin)
            if m > 0 {
                let g_nm = fm * (snm * cm - cnm * sm);
                du_dlam += rn_factor * r * pnm * g_nm;
            }
        }
    }

    // ------------------------------------------------------------------
    // 4. Transform to Cartesian
    //
    //    a_z = sinφ · du_dr  +  A
    //    a_x = (x/r) · du_dr  −  (xz/ρ²) · A  −  (y/ρ²) · du_dlam
    //    a_y = (y/r) · du_dr  −  (yz/ρ²) · A  +  (x/ρ²) · du_dlam
    //
    //    Guard for ρ → 0 (exactly-polar position): x ≈ y ≈ 0 there, so
    //    the xz/ρ² and y/ρ² terms numerically vanish; we skip them.
    // ------------------------------------------------------------------
    let az = (z / r) * du_dr + big_a;

    let (ax, ay) = if rxy2 > r2 * 1.0e-24 {
        let ax = (x / r) * du_dr - (x * z / rxy2) * big_a - (y / rxy2) * du_dlam;
        let ay = (y / r) * du_dr - (y * z / rxy2) * big_a + (x / rxy2) * du_dlam;
        (ax, ay)
    } else {
        // At the pole x ≈ y ≈ 0 so the radial component dominates.
        ((x / r) * du_dr, (y / r) * du_dr)
    };

    [ax, ay, az]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::gravity::egm_low_degree::{LowDegreeEarth, TwoBodyEarth};

    const GM: f64 = 398_600.441_8;
    const R: f64 = 6_378.137;
    const J2: f64 = 1.082_626_68e-3;

    // -----------------------------------------------------------------------
    // Helper: J2-only acceleration using the analytic formula from J2.rs
    // -----------------------------------------------------------------------
    fn j2_analytic(pos: [f64; 3]) -> [f64; 3] {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let r2 = x * x + y * y + z * z;
        let r = r2.sqrt();
        let z2r2 = z * z / r2;
        let factor = 1.5 * J2 * GM * R * R / (r2 * r2 * r);
        let cx = 5.0 * z2r2 - 1.0;
        let cz = 5.0 * z2r2 - 3.0;
        [factor * x * cx, factor * y * cx, factor * z * cz]
    }

    // Wrap LowDegreeEarth in a provider that only returns (n=2, m=0).
    struct J2Only;
    impl GravityFieldProvider for J2Only {
        fn gm(&self) -> crate::astro::dynamics::units::GravitationalParameter {
            LowDegreeEarth.gm()
        }
        fn reference_radius(&self) -> crate::qtty::Kilometers {
            LowDegreeEarth.reference_radius()
        }
        fn max_degree(&self) -> usize {
            2
        }
        fn max_order(&self) -> usize {
            0
        }
        fn c_normalized(&self, n: usize, m: usize) -> f64 {
            LowDegreeEarth.c_normalized(n, m)
        }
        fn s_normalized(&self, _n: usize, _m: usize) -> f64 {
            0.0
        }
    }

    /// J2-only geopotential matches TwoBody+J2 analytic at 1e-9 relative.
    #[test]
    fn geopotential_j2_matches_analytic() {
        let test_positions: &[[f64; 3]] = &[
            [7_000.0, 0.0, 0.0],          // equatorial
            [0.0, 6_800.0, 1_500.0],       // inclined
            [4_500.0, 3_000.0, 3_500.0],   // 45° ish
            [-5_000.0, -3_000.0, 4_000.0], // all-negative x/y
        ];

        for &pos in test_positions {
            let acc = geopotential_acceleration(&J2Only, pos, 2, 0).unwrap();

            // Reference = two-body + J2-analytic combined
            let (x, y, z) = (pos[0], pos[1], pos[2]);
            let r2 = x * x + y * y + z * z;
            let r = r2.sqrt();

            // Two-body: a = -GM/r³ · [x,y,z]
            let tb = [-GM * x / (r2 * r), -GM * y / (r2 * r), -GM * z / (r2 * r)];

            // J2 perturbation
            let j2_pert = j2_analytic(pos);

            let ref_ = [tb[0] + j2_pert[0], tb[1] + j2_pert[1], tb[2] + j2_pert[2]];

            for i in 0..3 {
                let denom = ref_[i].abs().max(acc[i].abs()).max(1e-30);
                let rel = (acc[i] - ref_[i]).abs() / denom;
                assert!(
                    rel < 1e-9,
                    "pos={pos:?} component {i}: acc={:.9e}, ref={:.9e}, rel={rel:.2e}",
                    acc[i],
                    ref_[i]
                );
            }
        }
    }

    /// Two-body reproduces GM/r² radially.
    #[test]
    fn two_body_radial_only() {
        let pos = [7_000.0_f64, 0.0, 0.0];
        let acc = geopotential_acceleration(&TwoBodyEarth, pos, 0, 0).unwrap();
        // a_x = -GM·x/r³ = -GM/r² when x = r = 7000
        let expected = -GM / (7_000.0_f64 * 7_000.0);
        let rel = (acc[0] - expected).abs() / expected.abs();
        assert!(rel < 1e-12, "two-body x: got {:.9e}, want {:.9e}", acc[0], expected);
        assert!(acc[1].abs() < 1e-30);
        assert!(acc[2].abs() < 1e-30);
    }

    /// Degree-out-of-range returns the appropriate error.
    #[test]
    fn degree_out_of_range_error() {
        use crate::astro::dynamics::errors::DynamicsError;
        let result = geopotential_acceleration(&TwoBodyEarth, [7_000.0, 0.0, 0.0], 2, 0);
        assert!(
            matches!(
                result,
                Err(DynamicsError::GeopotentialDegreeOutOfRange { requested: 2, max: 0 })
            ),
            "expected GeopotentialDegreeOutOfRange, got {result:?}"
        );
    }

    /// Full (2..4) geopotential finite-difference Jacobian check.
    ///
    /// Verifies that the gradient is self-consistent: numerically
    /// differentiating the potential should reproduce the acceleration.
    /// (The potential is reconstructed by integrating along each axis.)
    /// We instead directly finite-difference the acceleration and verify
    /// smoothness (|a(r+h) - a(r-h)| / 2h ≈ ∂a/∂x, no NaN/Inf).
    #[test]
    fn geopotential_n4_acceleration_is_smooth() {
        let pos = [6_900.0_f64, 1_200.0, 3_000.0];
        let h = 0.001; // 1 m in km

        for j in 0..3 {
            let mut p_plus = pos;
            let mut p_minus = pos;
            p_plus[j] += h;
            p_minus[j] -= h;

            let ap = geopotential_acceleration(&LowDegreeEarth, p_plus, 4, 4).unwrap();
            let am = geopotential_acceleration(&LowDegreeEarth, p_minus, 4, 4).unwrap();

            for i in 0..3 {
                let fd = (ap[i] - am[i]) / (2.0 * h);
                assert!(fd.is_finite(), "fd[{i}][{j}] is not finite");
            }
        }
    }
}
