//! VSOP87 position / velocity computation (stable‑Rust version)
//!
//! Exposes three public helpers:
//! * [`position`]  – only X,Y,Z (AstronomicalUnits)
//! * [`velocity`]  – only Ẋ,Ẏ,Ż (AstronomicalUnits/day)
//! * [`position_velocity`] – both in one pass (≈30 % faster in 2 calls)

use crate::astro::JulianDate;
use rayon::join;

/// One VSOP87 coefficient term  _a · cos(b + c·T)_
#[derive(Debug, Clone, Copy)]
pub struct Vsop87 {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

trait Mode {
    const NEED_VAL: bool;
    const NEED_DER: bool;
}
struct Val;
struct Der;
struct Both;
impl Mode for Val {
    const NEED_VAL: bool = true;
    const NEED_DER: bool = false;
}
impl Mode for Der {
    const NEED_VAL: bool = false;
    const NEED_DER: bool = true;
}
impl Mode for Both {
    const NEED_VAL: bool = true;
    const NEED_DER: bool = true;
}

// computes (value, d/dt) according to `M`.
#[inline]
fn coord<M: Mode>(series_by_power: &[&[Vsop87]], t: f64) -> (f64, f64) {
    let mut t_pow = 1.0; // T^0
    let mut t_pow_der = 0.0; // d/dT T^k
    let mut value = 0.0;
    let mut deriv_t = 0.0;

    for (k, terms) in series_by_power.iter().enumerate() {
        let mut serie_val = 0.0;
        let mut serie_der = 0.0;

        if M::NEED_VAL || M::NEED_DER {
            for term in *terms {
                let arg = term.b + term.c * t;
                let cos_arg = arg.cos();
                // We need S(T) even in Der mode for the k*T^{k-1}*S(T) term.
                serie_val += term.a * cos_arg;
                if M::NEED_DER {
                    serie_der += -term.a * term.c * arg.sin();
                }
            }
        }

        if M::NEED_VAL {
            value += t_pow * serie_val;
        }
        if M::NEED_DER {
            deriv_t += t_pow * serie_der + t_pow_der * serie_val;
        }

        // Prepare for next k: d/dT T^{k+1} = (k+1) * T^k
        t_pow_der = (k as f64 + 1.0) * t_pow;
        t_pow *= t;
    }

    const DT_DT: f64 = 1.0 / 365_250.0; // dT/dt (T per day)
    (value, deriv_t * DT_DT)
}

// Public façade
pub fn position(
    jd: JulianDate,
    x_series: &[&[Vsop87]],
    y_series: &[&[Vsop87]],
    z_series: &[&[Vsop87]],
) -> (f64, f64, f64) {
    let t = JulianDate::tt_to_tdb(jd).julian_millennias();

    let (x, (y, z)) = join(
        || coord::<Val>(x_series, t).0,
        || {
            let y = coord::<Val>(y_series, t).0;
            let z = coord::<Val>(z_series, t).0;
            (y, z)
        },
    );
    (x, y, z)
}

pub fn velocity(
    jd: JulianDate,
    x_series: &[&[Vsop87]],
    y_series: &[&[Vsop87]],
    z_series: &[&[Vsop87]],
) -> (f64, f64, f64) {
    let t = JulianDate::tt_to_tdb(jd).julian_millennias();

    let (xdot, (ydot, zdot)) = join(
        || coord::<Der>(x_series, t).1,
        || {
            let ydot = coord::<Der>(y_series, t).1;
            let zdot = coord::<Der>(z_series, t).1;
            (ydot, zdot)
        },
    );
    (xdot, ydot, zdot)
}

/// Position **and** velocity in a single pass (≈30 % faster than calling
/// the two previous helpers).
pub fn position_velocity(
    jd: JulianDate,
    x_series: &[&[Vsop87]],
    y_series: &[&[Vsop87]],
    z_series: &[&[Vsop87]],
) -> ((f64, f64, f64), (f64, f64, f64)) {
    let t = JulianDate::tt_to_tdb(jd).julian_millennias();

    let ((x, xdot), (y, ydot, z, zdot)) = join(
        || coord::<Both>(x_series, t),
        || {
            let (y, ydot) = coord::<Both>(y_series, t);
            let (z, zdot) = coord::<Both>(z_series, t);
            (y, ydot, z, zdot)
        },
    );

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

    #[test]
    fn test_coord_without_value_or_derivative() {
        struct NoneMode;
        impl super::Mode for NoneMode {
            const NEED_VAL: bool = false;
            const NEED_DER: bool = false;
        }
        let empty: [&[Vsop87]; 0] = [];
        let (v, d) = super::coord::<NoneMode>(&empty, 0.1);
        assert_eq!(v, 0.0);
        assert_eq!(d, 0.0);
    }
}
