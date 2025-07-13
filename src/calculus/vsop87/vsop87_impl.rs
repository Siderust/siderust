//! VSOP87 position / velocity computation (stable‑Rust version)
//!
//! Exposes three public helpers:
//! * [`position`]  – only X,Y,Z (AstronomicalUnits)
//! * [`velocity`]  – only Ẋ,Ẏ,Ż (AstronomicalUnits/day)
//! * [`position_velocity`] – both in one pass (≈30 % más rápido que dos llamadas)
//!
//! Internamente usamos un *tipo marcador* + *const bools* en vez de const‑enums,
//! para que compile en **Rust estable** (sin la feature `adt_const_params`).

use rayon::join;
use crate::units::JulianDay;

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
struct Val;  struct Der;  struct Both;
impl Mode for Val  { const NEED_VAL: bool = true;  const NEED_DER: bool = false; }
impl Mode for Der  { const NEED_VAL: bool = false; const NEED_DER: bool = true;  }
impl Mode for Both { const NEED_VAL: bool = true;  const NEED_DER: bool = true;  }

// computes (value, d/dt) according to `M`.
#[inline]
fn coord<M: Mode>(series_by_power: &[&[Vsop87]], t: f64) -> (f64, f64) {
    let mut t_pow     = 1.0;  // T^0
    let mut t_pow_der = 0.0;  // d/dT T^k
    let mut value     = 0.0;
    let mut deriv_t   = 0.0;

    for (k, terms) in series_by_power.iter().enumerate() {
        let mut serie_val = 0.0;
        let mut serie_der = 0.0;
        if M::NEED_VAL || M::NEED_DER {
            for term in *terms {
                let arg = term.b + term.c * t;
                if M::NEED_VAL {
                    serie_val += term.a * arg.cos();
                }
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

        t_pow_der = (k as f64 + 1.0) * t_pow; // (k)T^{k-1}
        t_pow    *= t;
    }

    const DT_DT: f64 = 1.0 / 365_250.0; // dT / dt  (T per day)
    (value, deriv_t * DT_DT)
}

// Public façade
pub fn position(
    jd: JulianDay,
    x_series: &[&[Vsop87]],
    y_series: &[&[Vsop87]],
    z_series: &[&[Vsop87]],
) -> (f64, f64, f64) {
    let t = JulianDay::tt_to_tdb(jd).julian_millennias();

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
    jd: JulianDay,
    x_series: &[&[Vsop87]],
    y_series: &[&[Vsop87]],
    z_series: &[&[Vsop87]],
) -> (f64, f64, f64) {
    let t = JulianDay::tt_to_tdb(jd).julian_millennias();

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
    jd: JulianDay,
    x_series: &[&[Vsop87]],
    y_series: &[&[Vsop87]],
    z_series: &[&[Vsop87]],
) -> ((f64, f64, f64), (f64, f64, f64)) {
    let t = JulianDay::tt_to_tdb(jd).julian_millennias();

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
