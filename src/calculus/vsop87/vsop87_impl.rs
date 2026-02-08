// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! VSOP87 position / velocity computation (optimized version)
//!
//! Exposes three public helpers:
//! * [`position`]  – only X,Y,Z (AstronomicalUnits)
//! * [`velocity`]  – only Ẋ,Ẏ,Ż (AstronomicalUnits/day)
//! * [`position_velocity`] – both in one pass (≈30 % faster than 2 calls)
//!
//! ## Optimizations applied
//! - Separate specialized functions eliminate inner-loop branching
//! - `mul_add` for FMA precision and potential performance
//! - `#[inline(always)]` on hot paths
//! - SIMD batching via `wide` crate for sin/cos operations

use crate::astro::JulianDate;
use wide::f64x4;

/// One VSOP87 coefficient term  _a · cos(b + c·T)_
#[derive(Debug, Clone, Copy)]
pub struct Vsop87 {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

/// Conversion factor: dT/dt (T in Julian millennia per day)
const DT_DT: f64 = 1.0 / 365_250.0;

// ═══════════════════════════════════════════════════════════════════════════
// SIMD helpers for batched sin/cos using wide crate
// ═══════════════════════════════════════════════════════════════════════════

/// Compute sin and cos for 4 values simultaneously using SIMD.
#[inline(always)]
fn sin_cos_x4(args: f64x4) -> (f64x4, f64x4) {
    // wide's sin/cos use polynomial approximations vectorized across lanes
    (args.sin(), args.cos())
}

// ═══════════════════════════════════════════════════════════════════════════
// Specialized coordinate functions (no branching in hot loops)
// ═══════════════════════════════════════════════════════════════════════════

/// Compute position value only (no derivative).
/// Uses SIMD batching for cos operations.
#[inline(always)]
fn coord_value(series_by_power: &[&[Vsop87]], t: f64) -> f64 {
    let mut t_pow = 1.0_f64;
    let mut value = 0.0_f64;

    for terms in series_by_power.iter() {
        let mut serie_val = 0.0_f64;
        let len = terms.len();

        // Process 4 terms at a time with SIMD
        let chunks = len / 4;
        let remainder = len % 4;

        for i in 0..chunks {
            let base = i * 4;
            let t0 = &terms[base];
            let t1 = &terms[base + 1];
            let t2 = &terms[base + 2];
            let t3 = &terms[base + 3];

            // Compute arguments using mul_add: arg = c * t + b
            let args = f64x4::new([
                t0.c.mul_add(t, t0.b),
                t1.c.mul_add(t, t1.b),
                t2.c.mul_add(t, t2.b),
                t3.c.mul_add(t, t3.b),
            ]);

            let cos_args = args.cos();
            let cos_arr = cos_args.to_array();

            // Accumulate: a * cos(arg)
            serie_val += t0.a * cos_arr[0];
            serie_val += t1.a * cos_arr[1];
            serie_val += t2.a * cos_arr[2];
            serie_val += t3.a * cos_arr[3];
        }

        // Handle remaining terms (0-3) with scalar operations
        for term in &terms[chunks * 4..chunks * 4 + remainder] {
            let arg = term.c.mul_add(t, term.b);
            serie_val += term.a * arg.cos();
        }

        value = t_pow.mul_add(serie_val, value);
        t_pow *= t;
    }

    value
}

/// Compute derivative only (no position value).
/// Uses SIMD batching for sin_cos operations.
#[inline(always)]
fn coord_deriv(series_by_power: &[&[Vsop87]], t: f64) -> f64 {
    let mut t_pow = 1.0_f64;
    let mut t_pow_der = 0.0_f64;
    let mut deriv_t = 0.0_f64;

    for (k, terms) in series_by_power.iter().enumerate() {
        let mut serie_val = 0.0_f64;
        let mut serie_der = 0.0_f64;
        let len = terms.len();

        // Process 4 terms at a time with SIMD
        let chunks = len / 4;
        let remainder = len % 4;

        for i in 0..chunks {
            let base = i * 4;
            let t0 = &terms[base];
            let t1 = &terms[base + 1];
            let t2 = &terms[base + 2];
            let t3 = &terms[base + 3];

            let args = f64x4::new([
                t0.c.mul_add(t, t0.b),
                t1.c.mul_add(t, t1.b),
                t2.c.mul_add(t, t2.b),
                t3.c.mul_add(t, t3.b),
            ]);

            let (sin_args, cos_args) = sin_cos_x4(args);
            let sin_arr = sin_args.to_array();
            let cos_arr = cos_args.to_array();

            // serie_val += a * cos(arg)
            // serie_der += -a * c * sin(arg)
            serie_val += t0.a * cos_arr[0];
            serie_val += t1.a * cos_arr[1];
            serie_val += t2.a * cos_arr[2];
            serie_val += t3.a * cos_arr[3];

            serie_der -= t0.a * t0.c * sin_arr[0];
            serie_der -= t1.a * t1.c * sin_arr[1];
            serie_der -= t2.a * t2.c * sin_arr[2];
            serie_der -= t3.a * t3.c * sin_arr[3];
        }

        // Handle remaining terms with scalar operations
        for term in &terms[chunks * 4..chunks * 4 + remainder] {
            let arg = term.c.mul_add(t, term.b);
            let (sin_arg, cos_arg) = arg.sin_cos();
            serie_val += term.a * cos_arg;
            serie_der -= term.a * term.c * sin_arg;
        }

        deriv_t = t_pow.mul_add(serie_der, deriv_t);
        deriv_t = t_pow_der.mul_add(serie_val, deriv_t);

        t_pow_der = (k as f64 + 1.0) * t_pow;
        t_pow *= t;
    }

    deriv_t * DT_DT
}

/// Compute both position value and derivative in one pass.
/// Uses SIMD batching for sin_cos operations.
#[inline(always)]
fn coord_both(series_by_power: &[&[Vsop87]], t: f64) -> (f64, f64) {
    let mut t_pow = 1.0_f64;
    let mut t_pow_der = 0.0_f64;
    let mut value = 0.0_f64;
    let mut deriv_t = 0.0_f64;

    for (k, terms) in series_by_power.iter().enumerate() {
        let mut serie_val = 0.0_f64;
        let mut serie_der = 0.0_f64;
        let len = terms.len();

        // Process 4 terms at a time with SIMD
        let chunks = len / 4;
        let remainder = len % 4;

        for i in 0..chunks {
            let base = i * 4;
            let t0 = &terms[base];
            let t1 = &terms[base + 1];
            let t2 = &terms[base + 2];
            let t3 = &terms[base + 3];

            let args = f64x4::new([
                t0.c.mul_add(t, t0.b),
                t1.c.mul_add(t, t1.b),
                t2.c.mul_add(t, t2.b),
                t3.c.mul_add(t, t3.b),
            ]);

            let (sin_args, cos_args) = sin_cos_x4(args);
            let sin_arr = sin_args.to_array();
            let cos_arr = cos_args.to_array();

            serie_val += t0.a * cos_arr[0];
            serie_val += t1.a * cos_arr[1];
            serie_val += t2.a * cos_arr[2];
            serie_val += t3.a * cos_arr[3];

            serie_der -= t0.a * t0.c * sin_arr[0];
            serie_der -= t1.a * t1.c * sin_arr[1];
            serie_der -= t2.a * t2.c * sin_arr[2];
            serie_der -= t3.a * t3.c * sin_arr[3];
        }

        // Handle remaining terms with scalar operations
        for term in &terms[chunks * 4..chunks * 4 + remainder] {
            let arg = term.c.mul_add(t, term.b);
            let (sin_arg, cos_arg) = arg.sin_cos();
            serie_val += term.a * cos_arg;
            serie_der -= term.a * term.c * sin_arg;
        }

        value = t_pow.mul_add(serie_val, value);
        deriv_t = t_pow.mul_add(serie_der, deriv_t);
        deriv_t = t_pow_der.mul_add(serie_val, deriv_t);

        t_pow_der = (k as f64 + 1.0) * t_pow;
        t_pow *= t;
    }

    (value, deriv_t * DT_DT)
}

// ═══════════════════════════════════════════════════════════════════════════
// Public façade
// ═══════════════════════════════════════════════════════════════════════════

pub fn position(
    jd: JulianDate,
    x_series: &[&[Vsop87]],
    y_series: &[&[Vsop87]],
    z_series: &[&[Vsop87]],
) -> (f64, f64, f64) {
    let t = JulianDate::tt_to_tdb(jd).julian_millennias();

    let x = coord_value(x_series, t);
    let y = coord_value(y_series, t);
    let z = coord_value(z_series, t);

    (x, y, z)
}

pub fn velocity(
    jd: JulianDate,
    x_series: &[&[Vsop87]],
    y_series: &[&[Vsop87]],
    z_series: &[&[Vsop87]],
) -> (f64, f64, f64) {
    let t = JulianDate::tt_to_tdb(jd).julian_millennias();

    let xdot = coord_deriv(x_series, t);
    let ydot = coord_deriv(y_series, t);
    let zdot = coord_deriv(z_series, t);

    (xdot, ydot, zdot)
}

/// Position **and** velocity in a single pass (≈30 % faster than calling
/// the two previous helpers).
pub fn position_velocity(
    jd: JulianDate,
    x_series: &[&[Vsop87]],
    y_series: &[&[Vsop87]],
    z_series: &[&[Vsop87]],
) -> ((f64, f64, f64), (f64, f64, f64)) {
    let t = JulianDate::tt_to_tdb(jd).julian_millennias();

    let (x, xdot) = coord_both(x_series, t);
    let (y, ydot) = coord_both(y_series, t);
    let (z, zdot) = coord_both(z_series, t);

    ((x, y, z), (xdot, ydot, zdot))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::JulianDate;

    const X0: [Vsop87; 1] = [Vsop87 {
        a: 1.0,
        b: 0.0,
        c: 0.0,
    }];
    const X1: [Vsop87; 1] = [Vsop87 {
        a: 2.0,
        b: 0.0,
        c: 0.0,
    }];
    const Y0: [Vsop87; 1] = [Vsop87 {
        a: 0.0,
        b: 0.0,
        c: 0.0,
    }];
    const Y1: [Vsop87; 0] = [];
    const Y2: [Vsop87; 1] = [Vsop87 {
        a: 3.0,
        b: 0.0,
        c: 0.0,
    }];
    const Z0: [Vsop87; 1] = [Vsop87 {
        a: 4.0,
        b: 0.0,
        c: 0.0,
    }];

    #[test]
    fn test_position_velocity() {
        // 100 years after J2000 -> T = 0.1
        let jd = JulianDate::new(2451545.0 + 36525.0);
        let pos = position(jd, &[&X0, &X1], &[&Y0, &Y1, &Y2], &[&Z0]);
        assert!((pos.0 - 1.2).abs() < 1e-12);
        assert!((pos.1 - 0.03).abs() < 1e-12);
        assert!((pos.2 - 4.0).abs() < 1e-12);

        let vel = velocity(jd, &[&X0, &X1], &[&Y0, &Y1, &Y2], &[&Z0]);
        let dt = 1.0 / 365_250.0;
        assert!((vel.0 - 2.0 * dt).abs() < 1e-12);
        assert!((vel.1 - 0.6 * dt).abs() < 1e-12);
        assert!((vel.2 - 0.0).abs() < 1e-12);

        let (p2, v2) = position_velocity(jd, &[&X0, &X1], &[&Y0, &Y1, &Y2], &[&Z0]);
        assert!((p2.0 - pos.0).abs() < 1e-12);
        assert!((p2.1 - pos.1).abs() < 1e-12);
        assert!((p2.2 - pos.2).abs() < 1e-12);
        assert!((v2.0 - vel.0).abs() < 1e-12);
        assert!((v2.1 - vel.1).abs() < 1e-12);
        assert!((v2.2 - vel.2).abs() < 1e-12);
    }
}
