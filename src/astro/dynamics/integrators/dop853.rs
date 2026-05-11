// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Hairer / Norsett / Wanner DOP853 8th-order adaptive Runge-Kutta integrator.
//!
//! Reference: E. Hairer, S. P. Norsett, G. Wanner, *Solving Ordinary
//! Differential Equations I*, 2nd ed., Springer (1993), §II.5, Table 5.2.
//! FORTRAN source: <http://www.unige.ch/~hairer/software.html>.

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::forces::ForceModel;
use crate::astro::dynamics::state::{OrbitState, StateDerivative};
use crate::ext_qtty::tolerances::IntegratorTolerances;
use crate::qtty::Second;
use crate::time::{Time, TT};

#[inline]
fn state_component(s: &OrbitState, i: usize) -> f64 {
    match i {
        0 => s.position.x().value(),
        1 => s.position.y().value(),
        2 => s.position.z().value(),
        3 => s.velocity.x().value(),
        4 => s.velocity.y().value(),
        5 => s.velocity.z().value(),
        _ => panic!("index out of range"),
    }
}

#[inline]
fn deriv_component(d: &StateDerivative, i: usize) -> f64 {
    match i {
        0 => d.vel.x().value(),
        1 => d.vel.y().value(),
        2 => d.vel.z().value(),
        3 => d.acc.x().value(),
        4 => d.acc.y().value(),
        5 => d.acc.z().value(),
        _ => panic!("index out of range"),
    }
}

#[inline]
fn rhs<FM: ForceModel>(
    force: &FM,
    s: &OrbitState,
    ctx: &DynamicsContext,
) -> Result<StateDerivative, DynamicsError> {
    Ok(StateDerivative::new(s.velocity, force.acceleration(s, ctx)?))
}

#[inline]
fn state_at(s: &OrbitState, d: &StateDerivative, h: f64, dt: f64) -> OrbitState {
    let new_epoch = s.epoch + Second::new(dt);
    let advanced = s.advance(d, Second::new(h));
    OrbitState {
        epoch: new_epoch,
        position: advanced.position,
        velocity: advanced.velocity,
    }
}

/// Stateful DOP853 integrator (8th-order adaptive Runge-Kutta).
pub struct Dop853 {
    pub tolerances: IntegratorTolerances,
    pub h_max: Second,
    pub h_min: Second,
}

impl Dop853 {
    pub fn new(tolerances: IntegratorTolerances) -> Self {
        Self {
            tolerances,
            h_max: Second::new(86_400.0),
            h_min: Second::new(1e-6),
        }
    }
    pub fn with_h_max(mut self, h_max: Second) -> Self {
        self.h_max = h_max;
        self
    }
    pub fn with_h_min(mut self, h_min: Second) -> Self {
        self.h_min = h_min;
        self
    }
}

impl super::AdaptiveStepper for Dop853 {
    fn step<FM: ForceModel>(
        &self,
        force: &FM,
        state: &OrbitState,
        h_try: Second,
        ctx: &DynamicsContext,
    ) -> Result<(OrbitState, Second, Second), DynamicsError> {
        let (s, h_used, h_next, _) = dop853_step(force, state, h_try, self.tolerances, ctx)?;
        Ok((s, h_used, h_next))
    }
}

/// Cached endpoints of an accepted DOP853 step, for cubic-Hermite dense output.
pub struct Dop853Step {
    pub state_start: OrbitState,
    pub state_end: OrbitState,
    deriv_start: StateDerivative,
    deriv_end: StateDerivative,
    /// Actual (signed) step size in seconds.
    pub h_used: f64,
}

impl Dop853Step {
    /// Cubic-Hermite interpolation at `t`.  Returns `None` if `t` is outside
    /// the closed interval spanned by the step (regardless of direction).
    pub fn interpolate(&self, t: Time<TT>) -> Option<OrbitState> {
        let t_start = self.state_start.epoch;
        let t_end = self.state_end.epoch;
        let dt_from_start = (t - t_start).value();
        let theta = dt_from_start / self.h_used;
        if !(0.0..=1.0).contains(&theta) {
            return None;
        }
        let _ = t_end;
        let h = self.h_used;
        let h00 = 2.0 * theta.powi(3) - 3.0 * theta.powi(2) + 1.0;
        let h10 = theta.powi(3) - 2.0 * theta.powi(2) + theta;
        let h01 = -2.0 * theta.powi(3) + 3.0 * theta.powi(2);
        let h11 = theta.powi(3) - theta.powi(2);

        let interp = |y0: f64, y1: f64, f0: f64, f1: f64| -> f64 {
            h00 * y0 + h * h10 * f0 + h01 * y1 + h * h11 * f1
        };

        let px = interp(
            self.state_start.position.x().value(),
            self.state_end.position.x().value(),
            self.deriv_start.vel.x().value(),
            self.deriv_end.vel.x().value(),
        );
        let py = interp(
            self.state_start.position.y().value(),
            self.state_end.position.y().value(),
            self.deriv_start.vel.y().value(),
            self.deriv_end.vel.y().value(),
        );
        let pz = interp(
            self.state_start.position.z().value(),
            self.state_end.position.z().value(),
            self.deriv_start.vel.z().value(),
            self.deriv_end.vel.z().value(),
        );
        let vx = interp(
            self.state_start.velocity.x().value(),
            self.state_end.velocity.x().value(),
            self.deriv_start.acc.x().value(),
            self.deriv_end.acc.x().value(),
        );
        let vy = interp(
            self.state_start.velocity.y().value(),
            self.state_end.velocity.y().value(),
            self.deriv_start.acc.y().value(),
            self.deriv_end.acc.y().value(),
        );
        let vz = interp(
            self.state_start.velocity.z().value(),
            self.state_end.velocity.z().value(),
            self.deriv_start.acc.z().value(),
            self.deriv_end.acc.z().value(),
        );

        use crate::astro::dynamics::state::{Position, Velocity};
        Some(OrbitState {
            epoch: t_start + Second::new(dt_from_start),
            position: Position::new(px, py, pz),
            velocity: Velocity::new(vx, vy, vz),
        })
    }
}

#[allow(clippy::too_many_lines)]
pub fn dop853_step<FM: ForceModel>(
    force: &FM,
    s: &OrbitState,
    h_try: Second,
    tol: IntegratorTolerances,
    ctx: &DynamicsContext,
) -> Result<(OrbitState, Second, Second, Dop853Step), DynamicsError> {
    // ---------- Butcher tableau (Hairer dop853.f) ----------
    let c2 = 5.260_015_195_876_773e-2;
    let c3 = 7.890_022_793_815_16e-2;
    let c4 = 1.183_503_419_072_274e-1;
    let c5 = 2.816_496_580_927_726e-1;
    let c6 = 3.333_333_333_333_333e-1;
    let c7 = 0.25;
    let c8 = 3.076_923_076_923_077e-1;
    let c9 = 6.512_820_512_820_513e-1;
    let c10 = 0.6;
    let c11 = 8.571_428_571_428_571e-1;

    let a21 = 5.260_015_195_876_773e-2;

    let a31 = 1.972_505_698_453_790e-2;
    let a32 = 5.917_517_095_361_370e-2;

    let a41 = 2.958_758_547_680_685e-2;
    let a43 = 8.876_275_643_042_054e-2;

    let a51 = 2.413_656_412_274_204e-1;
    let a53 = -8.845_494_793_282_861e-1;
    let a54 = 9.248_340_032_617_92e-1;

    let a61 = 3.703_703_703_703_704e-2;
    let a64 = 1.708_286_087_294_739e-1;
    let a65 = 1.254_676_875_668_224e-1;

    let a71 = 3.710_937_5e-2;
    let a74 = 1.702_522_110_195_440e-1;
    let a75 = 6.021_653_898_045_591e-2;
    let a76 = -1.757_812_5e-2;

    let a81 = 3.709_200_011_850_479e-2;
    let a84 = 1.703_839_257_122_400e-1;
    let a85 = 1.072_620_304_463_733e-1;
    let a86 = -1.531_943_774_862_449e-2;
    let a87 = 8.273_789_167_928_145e-3;

    let a91 = 6.241_109_587_160_757e-1;
    let a94 = -3.360_892_629_446_941e0;
    let a95 = -8.682_193_468_417_260e-1;
    let a96 = 2.722_122_997_654_576e1;
    let a97 = 2.015_406_755_047_789e1;
    let a98 = -4.348_988_418_106_996e1;

    let a101 = 4.776_625_364_382_644e-1;
    let a104 = -2.488_114_619_971_668e0;
    let a105 = -5.902_908_268_368_43e-1;
    let a106 = 2.123_005_144_818_119e1;
    let a107 = 1.531_511_698_888_877e1;
    let a108 = -3.324_220_794_756_98e1;
    let a109 = 8.105_125_241_052_84e-1;

    let a111 = -9.314_633_917_552_932e-1;
    let a114 = 5.641_239_626_156_222e0;
    let a115 = 1.409_249_292_786_247e0;
    let a116 = -4.923_172_037_524_191e1;
    let a117 = -3.628_616_669_455_340e1;
    let a118 = 7.913_862_991_059_162e1;
    let a119 = -9.151_243_169_408_321e-1;
    let a1110 = -1.149_114_000_129_324e-1;

    let a121 = 2.608_554_629_667_515e-1;
    let a124 = -1.241_915_676_189_275e0;
    let a125 = -3.034_497_369_872_499e-1;
    let a126 = 1.104_149_357_046_952e1;
    let a127 = 8.257_647_857_407_497e0;
    let a128 = -1.756_676_255_618_410e1;
    let a129 = -1.570_155_727_806_287e-1;
    let a1210 = 2.910_850_288_129_909e-1;
    let a1211 = 1.845_533_387_295_413e-3;

    let b1 = 5.429_373_411_656_873e-2;
    let b6 = 4.450_312_892_752_409e0;
    let b7 = 1.891_517_899_314_500e0;
    let b8 = -5.801_203_960_010_585e0;
    let b9 = 3.111_643_669_578_199e-1;
    let b10 = -1.521_609_496_625_161e-1;
    let b11 = 2.013_654_008_040_303e-1;
    let b12 = 4.471_061_572_777_259e-2;

    let er1 = 1.312_004_499_419_488e-2;
    let er6 = -1.225_156_446_376_204e-1;
    let er7 = -4.907_422_132_005_108e-3;
    let er8 = 3.537_736_209_119_201e-2;
    let er9 = 6.115_440_678_252_199e-4;
    let er10 = -6.507_510_337_184_27e-4;
    let er11 = -2.239_587_899_853_073e-3;
    let er12 = 3.163_055_608_430_268e-5;

    let mut h = h_try.value();
    if !h.is_finite() || h == 0.0 {
        return Err(DynamicsError::InvalidStepRequest {
            reason: "step size must be finite and non-zero",
        });
    }
    let mut iters = 0u32;

    loop {
        let k1 = rhs(force, s, ctx)?;
        let k2 = rhs(force, &state_at(s, &k1.scaled(a21), h, c2 * h), ctx)?;
        let k3 = rhs(
            force,
            &state_at(s, &k1.scaled(a31).add(&k2.scaled(a32)), h, c3 * h),
            ctx,
        )?;
        let k4 = rhs(
            force,
            &state_at(s, &k1.scaled(a41).add(&k3.scaled(a43)), h, c4 * h),
            ctx,
        )?;
        let k5 = rhs(
            force,
            &state_at(
                s,
                &k1.scaled(a51).add(&k3.scaled(a53)).add(&k4.scaled(a54)),
                h,
                c5 * h,
            ),
            ctx,
        )?;
        let k6 = rhs(
            force,
            &state_at(
                s,
                &k1.scaled(a61).add(&k4.scaled(a64)).add(&k5.scaled(a65)),
                h,
                c6 * h,
            ),
            ctx,
        )?;
        let k7 = rhs(
            force,
            &state_at(
                s,
                &k1.scaled(a71)
                    .add(&k4.scaled(a74))
                    .add(&k5.scaled(a75))
                    .add(&k6.scaled(a76)),
                h,
                c7 * h,
            ),
            ctx,
        )?;
        let k8 = rhs(
            force,
            &state_at(
                s,
                &k1.scaled(a81)
                    .add(&k4.scaled(a84))
                    .add(&k5.scaled(a85))
                    .add(&k6.scaled(a86))
                    .add(&k7.scaled(a87)),
                h,
                c8 * h,
            ),
            ctx,
        )?;
        let k9 = rhs(
            force,
            &state_at(
                s,
                &k1.scaled(a91)
                    .add(&k4.scaled(a94))
                    .add(&k5.scaled(a95))
                    .add(&k6.scaled(a96))
                    .add(&k7.scaled(a97))
                    .add(&k8.scaled(a98)),
                h,
                c9 * h,
            ),
            ctx,
        )?;
        let k10 = rhs(
            force,
            &state_at(
                s,
                &k1.scaled(a101)
                    .add(&k4.scaled(a104))
                    .add(&k5.scaled(a105))
                    .add(&k6.scaled(a106))
                    .add(&k7.scaled(a107))
                    .add(&k8.scaled(a108))
                    .add(&k9.scaled(a109)),
                h,
                c10 * h,
            ),
            ctx,
        )?;
        let k11 = rhs(
            force,
            &state_at(
                s,
                &k1.scaled(a111)
                    .add(&k4.scaled(a114))
                    .add(&k5.scaled(a115))
                    .add(&k6.scaled(a116))
                    .add(&k7.scaled(a117))
                    .add(&k8.scaled(a118))
                    .add(&k9.scaled(a119))
                    .add(&k10.scaled(a1110)),
                h,
                c11 * h,
            ),
            ctx,
        )?;
        let k12 = rhs(
            force,
            &state_at(
                s,
                &k1.scaled(a121)
                    .add(&k4.scaled(a124))
                    .add(&k5.scaled(a125))
                    .add(&k6.scaled(a126))
                    .add(&k7.scaled(a127))
                    .add(&k8.scaled(a128))
                    .add(&k9.scaled(a129))
                    .add(&k10.scaled(a1210))
                    .add(&k11.scaled(a1211)),
                h,
                h,
            ),
            ctx,
        )?;

        let d_new = k1
            .scaled(b1)
            .add(&k6.scaled(b6))
            .add(&k7.scaled(b7))
            .add(&k8.scaled(b8))
            .add(&k9.scaled(b9))
            .add(&k10.scaled(b10))
            .add(&k11.scaled(b11))
            .add(&k12.scaled(b12));
        let s_new = state_at(s, &d_new, h, h);

        let err_d = k1
            .scaled(er1)
            .add(&k6.scaled(er6))
            .add(&k7.scaled(er7))
            .add(&k8.scaled(er8))
            .add(&k9.scaled(er9))
            .add(&k10.scaled(er10))
            .add(&k11.scaled(er11))
            .add(&k12.scaled(er12));

        let mut err_norm = 0.0;
        for i in 0..6 {
            let err = h * deriv_component(&err_d, i);
            let y0i = state_component(s, i);
            let y1i = state_component(&s_new, i);
            let abs_tol = if i < 3 {
                tol.abs_pos[i].value()
            } else {
                tol.abs_vel[i - 3].value()
            };
            let sc = abs_tol + tol.rel.value() * y0i.abs().max(y1i.abs());
            let r = err / sc;
            err_norm += r * r;
        }
        err_norm = (err_norm / 6.0).sqrt();

        if err_norm <= 1.0 || h.abs() < 1e-9 {
            let h_next = if err_norm == 0.0 {
                h * 6.0
            } else {
                let factor = 0.9 * err_norm.powf(-1.0 / 8.0);
                h * factor.clamp(1.0 / 3.0, 6.0)
            };
            // Compute k1 at endpoint for dense output (use k1 at new state).
            let deriv_end = rhs(force, &s_new, ctx)?;
            let dense = Dop853Step {
                state_start: s.clone(),
                state_end: s_new.clone(),
                deriv_start: k1,
                deriv_end,
                h_used: h,
            };
            return Ok((s_new, Second::new(h), Second::new(h_next), dense));
        } else {
            iters += 1;
            if iters > 50 {
                return Err(DynamicsError::InvalidStepRequest {
                    reason: "DOP853 step controller failed to converge",
                });
            }
            if !err_norm.is_finite() {
                // Force aggressive shrinkage when the trial step blew up.
                h *= 1.0 / 3.0;
                continue;
            }
            let factor = 0.9 * err_norm.powf(-1.0 / 7.0);
            h *= factor.clamp(1.0 / 3.0, 1.0);
        }
    }
}

/// Propagate `state` for `total_dt` using DOP853.
pub fn dop853_propagate<FM: ForceModel>(
    force: &FM,
    state: OrbitState,
    total_dt: Second,
    tol: IntegratorTolerances,
    ctx: &DynamicsContext,
) -> Result<OrbitState, DynamicsError> {
    let total_dt_s = total_dt.value();
    if total_dt_s == 0.0 {
        return Ok(state);
    }
    let mut s = state;
    let mut t = 0.0;
    let mut h = total_dt_s.signum() * 30.0_f64.min(total_dt_s.abs());
    while (total_dt_s - t).abs() > 1e-9 {
        if (t + h - total_dt_s) * total_dt_s.signum() > 0.0 {
            h = total_dt_s - t;
        }
        let (s_new, h_used, h_next, _) = dop853_step(force, &s, Second::new(h), tol, ctx)?;
        s = s_new;
        t += h_used.value();
        h = h_next.value();
    }
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::forces::TwoBody;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::ext_qtty::tolerances::IntegratorTolerances;
    use crate::time::JulianDate;

    fn circular() -> (OrbitState, f64, f64, f64) {
        let mu: f64 = 398_600.441_8;
        let r: f64 = 7_000.0;
        let v: f64 = (mu / r).sqrt();
        let s0 = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v, 0.0),
        );
        let period = 2.0 * std::f64::consts::PI * (r.powi(3) / mu).sqrt();
        (s0, r, v, period)
    }

    #[test]
    fn dop853_one_orbit_closes() {
        let (s0, r, _v, period) = circular();
        let ctx = DynamicsContext::empty();
        let s = dop853_propagate(
            &TwoBody::earth(),
            s0,
            Second::new(period),
            IntegratorTolerances::uniform(1e-12, 1e-12, 1e-12),
            &ctx,
        )
        .unwrap();
        let dr = ((s.position.x().value() - r).powi(2)
            + s.position.y().value().powi(2)
            + s.position.z().value().powi(2))
        .sqrt();
        assert!(dr < 1e-6, "orbit closure error {dr} km exceeds 1e-6 km");
    }

    #[test]
    fn dop853_convergence() {
        let (s0, r, _v, period) = circular();
        let ctx = DynamicsContext::empty();
        let mut prev_err = f64::INFINITY;
        for k in 0..3 {
            let tol = 1e-9 * 0.01_f64.powi(k);
            let s = dop853_propagate(
                &TwoBody::earth(),
                s0,
                Second::new(period),
                IntegratorTolerances::uniform(tol, tol, tol),
                &ctx,
            )
            .unwrap();
            let dr = ((s.position.x().value() - r).powi(2)
                + s.position.y().value().powi(2)
                + s.position.z().value().powi(2))
            .sqrt();
            assert!(
                dr <= prev_err + 1e-12,
                "convergence not monotone: {dr} vs {prev_err}"
            );
            prev_err = dr;
        }
    }

    #[test]
    fn dop853_step_rejection() {
        // Fire an initial step several orders of magnitude too large for the
        // requested tolerance and confirm the controller had to shrink it.
        let (s0, _r, _v, period) = circular();
        let ctx = DynamicsContext::empty();
        let tol = IntegratorTolerances::uniform(1e-12, 1e-12, 1e-12);
        let h_try = Second::new(period); // one full orbit in a single step
        let (_s, h_used, _h_next, _) =
            dop853_step(&TwoBody::earth(), &s0, h_try, tol, &ctx).unwrap();
        assert!(
            h_used.value().abs() < period * 0.5,
            "expected step shrinkage from {} but got {}",
            period,
            h_used.value()
        );
    }

    #[test]
    fn dop853_dense_output_at_endpoints() {
        let (s0, _r, _v, _period) = circular();
        let ctx = DynamicsContext::empty();
        let tol = IntegratorTolerances::uniform(1e-10, 1e-10, 1e-10);
        let (_s_end, h_used, _h_next, dense) =
            dop853_step(&TwoBody::earth(), &s0, Second::new(60.0), tol, &ctx).unwrap();
        let s_at_start = dense.interpolate(s0.epoch).unwrap();
        let dx0 = (s_at_start.position.x().value() - s0.position.x().value()).abs();
        assert!(dx0 < 1e-9, "endpoint start mismatch: {dx0}");

        let t_end = s0.epoch + h_used;
        let s_at_end = dense.interpolate(t_end).unwrap();
        let dxe = (s_at_end.position.x().value() - dense.state_end.position.x().value()).abs();
        assert!(dxe < 1e-9, "endpoint end mismatch: {dxe}");
    }

    #[test]
    fn dop853_backward_propagation() {
        let (s0, _r, _v, period) = circular();
        let ctx = DynamicsContext::empty();
        let tol = IntegratorTolerances::uniform(1e-12, 1e-12, 1e-12);
        let s_fwd = dop853_propagate(&TwoBody::earth(), s0, Second::new(period), tol, &ctx)
            .unwrap();
        let s_back = dop853_propagate(&TwoBody::earth(), s_fwd, Second::new(-period), tol, &ctx)
            .unwrap();
        let dr = ((s_back.position.x().value() - s0.position.x().value()).powi(2)
            + (s_back.position.y().value() - s0.position.y().value()).powi(2)
            + (s_back.position.z().value() - s0.position.z().value()).powi(2))
        .sqrt();
        assert!(dr < 1e-6, "round-trip error {dr} km exceeds 1e-6 km");
    }
}
