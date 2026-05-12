// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! J2 zonal-oblateness perturbation — `J2` force model with analytic partials.
//!
//! ## Scope
//!
//! Provides [`J2`], implementing the first zonal-harmonic perturbation to
//! central gravity and its analytic partial-derivative Jacobian `∂a_J2/∂r`.
//!
//! ## Equations
//!
//! The J2 acceleration in an Earth-centred inertial frame is (Vallado §8.8.1):
//!
//! ```text
//! a_x = C · x · (5z²/r² − 1) / r⁵
//! a_y = C · y · (5z²/r² − 1) / r⁵
//! a_z = C · z · (5z²/r² − 3) / r⁵
//! ```
//!
//! where `C = (3/2) · J₂ · μ · R_E²` and `r = |r|`.
//!
//! ### Analytic Jacobian `∂a_J2/∂r`
//!
//! Let `D = C / r⁷` and `q = z² / r²`.  The symmetric 3×3 matrix is:
//!
//! ```text
//! ∂a_x/∂x = D · [(5q − 1)r² − 5(7q − 1)x²]
//! ∂a_y/∂y = D · [(5q − 1)r² − 5(7q − 1)y²]
//! ∂a_z/∂z = D · r² · (30q − 3 − 35q²)
//!
//! ∂a_x/∂y = ∂a_y/∂x = D · 5(1 − 7q) xy
//! ∂a_x/∂z = ∂a_z/∂x = D · 5(3 − 7q) xz
//! ∂a_y/∂z = ∂a_z/∂y = D · 5(3 − 7q) yz
//! ```
//!
//! `A_v = ∂a_J2/∂v = 0` (no velocity dependence).
//!
//! ## Units
//!
//! Position km, velocity km/s, acceleration km/s², `A_r` s⁻².
//!
//! ## Frame/center assumptions
//!
//! `<Geocentric, GCRS>`.  The z-axis of GCRS is aligned with the Earth's
//! mean pole (CIP), so `z/r = sin(φ)` (geocentric latitude).
//!
//! ## References
//!
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, §8.8.1.
//! * Montenbruck & Gill, *Satellite Orbits*, §3.2.

use affn::matrix3::FrameMatrix3;

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::astro::dynamics::units::GM_EARTH;
use crate::coordinates::frames::GCRS;
use crate::qtty::{J2Coefficient, Kilometers};

use super::traits::{
    ForceModel, ForcePartials, GravitationalParameter, DEGENERATE_RADIUS_KM, R_EARTH,
};

/// J2 (Earth oblateness) perturbation acceleration in the inertial frame.
///
/// `j2` is the unnormalised zonal harmonic coefficient (Earth: `1.082_626_68e-3`).
#[derive(Debug, Clone, Copy)]
pub struct J2 {
    /// `GM` (km³/s²).
    pub gm: GravitationalParameter,
    /// Equatorial radius (km).
    pub req: Kilometers,
    /// Unnormalised J₂ coefficient (dimensionless).
    pub j2: J2Coefficient,
}

impl J2 {
    /// Standard Earth values (EGM2008 / WGS-84).
    pub fn earth() -> Self {
        Self {
            gm: GM_EARTH,
            req: R_EARTH,
            j2: J2Coefficient::new(1.082_626_68e-3),
        }
    }
}

impl ForceModel for J2 {
    #[inline]
    fn acceleration(
        &self,
        s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        let r = s.position.distance().value();
        if r < DEGENERATE_RADIUS_KM {
            return Err(DynamicsError::DegenerateGeometry {
                reason: "J2: radius near zero",
            });
        }
        let r2 = r * r;
        let rx = s.position.x().value();
        let ry = s.position.y().value();
        let rz = s.position.z().value();
        let z2_over_r2 = (rz * rz) / r2;
        let req = self.req.value();
        let factor = 1.5 * self.j2.value() * self.gm.value() * req * req / (r2 * r2 * r);
        let cx = 5.0 * z2_over_r2 - 1.0;
        let cz = 5.0 * z2_over_r2 - 3.0;
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            factor * rx * cx,
            factor * ry * cx,
            factor * rz * cz,
        ))
    }

    /// Analytic `∂a_J2/∂r` in GCRS frame.
    ///
    /// See module-level documentation for the full derivation.
    fn partials(
        &self,
        s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<ForcePartials<GCRS>, DynamicsError> {
        let x = s.position.x().value();
        let y = s.position.y().value();
        let z = s.position.z().value();
        let r2 = x * x + y * y + z * z;
        let r = r2.sqrt();
        if r == 0.0 {
            return Ok(ForcePartials::zero());
        }
        let req = self.req.value();
        let mu = self.gm.value();
        let j2 = self.j2.value();

        // C = (3/2) · J2 · μ · R_E²,  D = C / r⁷
        let c = 1.5 * j2 * mu * req * req;
        let r7 = r2 * r2 * r2 * r;
        let d = c / r7;

        // q = z²/r²
        let q = (z * z) / r2;

        // Diagonal entries
        let diag_xy_coeff = (5.0 * q - 1.0) * r2; // (5q−1)r²
        let off_xy = 5.0 * (7.0 * q - 1.0); // 5(7q−1)  [subtract x_i² · this]
        let dxx = d * (diag_xy_coeff - off_xy * x * x);
        let dyy = d * (diag_xy_coeff - off_xy * y * y);
        let dzz = d * r2 * (30.0 * q - 3.0 - 35.0 * q * q);

        // Off-diagonal entries (matrix is symmetric)
        let dxy = d * 5.0 * (1.0 - 7.0 * q) * x * y;
        let dxz = d * 5.0 * (3.0 - 7.0 * q) * x * z;
        let dyz = d * 5.0 * (3.0 - 7.0 * q) * y * z;

        let data = [[dxx, dxy, dxz], [dxy, dyy, dyz], [dxz, dyz, dzz]];

        Ok(ForcePartials {
            d_acc_d_pos: FrameMatrix3::from_array(data),
            d_acc_d_vel: FrameMatrix3::zero(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::state::OrbitState;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::time::JulianDate;

    fn leo_state(rx: f64, ry: f64, rz: f64) -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(rx, ry, rz),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn j2_acceleration_equatorial_has_no_z_perturbation() {
        // At z=0 the z-component of J2 should also be zero.
        let ctx = DynamicsContext::empty();
        let a = J2::earth()
            .acceleration(&leo_state(7000.0, 0.0, 0.0), &ctx)
            .unwrap();
        assert!(a.z().value().abs() < 1e-15);
    }

    #[test]
    fn j2_acceleration_polar_z_nonzero() {
        let ctx = DynamicsContext::empty();
        let a = J2::earth()
            .acceleration(&leo_state(0.0, 0.0, 7000.0), &ctx)
            .unwrap();
        assert!(a.z().value() != 0.0);
        assert!(a.x().value().abs() < 1e-15);
    }

    /// Verify J2 analytic partials against a 5-point finite-difference stencil.
    ///
    /// Relative tolerance: 1 × 10⁻⁶ (tight enough to catch sign/factor bugs).
    #[test]
    fn j2_partials_fd_vs_analytic() {
        let ctx = DynamicsContext::empty();
        let model = J2::earth();

        // Representative LEO: non-equatorial so all Jacobian entries are exercised.
        let r0 = [4_500.0_f64, 3_000.0, 3_500.0];
        let state = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r0[0], r0[1], r0[2]),
            Velocity::<GCRS>::new(0.0, 6.0, 4.0),
        );

        let analytic = model.partials(&state, &ctx).unwrap();
        let A = analytic.d_acc_d_pos.as_array();

        // 5-point stencil step: 1 km.
        let h = 1.0_f64; // km

        for j in 0..3 {
            let perturb = |delta: f64| -> [f64; 3] {
                let mut pos = r0;
                pos[j] += delta;
                let s = OrbitState::new_at_jd(
                    JulianDate::new(2_451_545.0),
                    Position::<GCRS>::new(pos[0], pos[1], pos[2]),
                    Velocity::<GCRS>::new(0.0, 6.0, 4.0),
                );
                let a = model.acceleration(&s, &ctx).unwrap();
                [a.x().value(), a.y().value(), a.z().value()]
            };

            let am2 = perturb(-2.0 * h);
            let am1 = perturb(-h);
            let ap1 = perturb(h);
            let ap2 = perturb(2.0 * h);

            for i in 0..3 {
                let fd = (-ap2[i] + 8.0 * ap1[i] - 8.0 * am1[i] + am2[i]) / (12.0 * h);
                let an = A[i][j];
                let denom = an.abs().max(fd.abs()).max(1e-30);
                let rel_err = (fd - an).abs() / denom;
                assert!(
                    rel_err < 1e-6,
                    "J2 ∂a[{i}]/∂r[{j}]: analytic={an:.6e}, fd={fd:.6e}, rel={rel_err:.2e}"
                );
            }
        }
    }
}
