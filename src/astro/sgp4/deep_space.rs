// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Lunar/solar secular and periodic perturbations and geopotential resonances.
//! Equations follow Vallado et al., AIAA 2006-6753, with the published Lyddane
//! and 720-minute Euler–Maclaurin integration conventions.

use super::model::{Constants, Orbit};
use core::f64::consts::{PI, TAU};
#[cfg(not(feature = "std"))]
use qtty::{Real, Scalar, Transcendental};

const SOLAR_RATE: f64 = 1.19459e-5;
const LUNAR_RATE: f64 = 1.5835218e-4;
const EARTH_RATE: f64 = 4.3752690880113e-3;

/// Eccentricity, inclination, longitude, perigee and node harmonic amplitudes.
#[derive(Clone, Debug)]
struct Periodic {
    e: [f64; 2],
    i: [f64; 2],
    l: [f64; 3],
    g: [f64; 3],
    h: [f64; 2],
    phase: f64,
    rate: f64,
    eccentricity: f64,
}
impl Periodic {
    fn at(&self, t: f64) -> [f64; 5] {
        let mean = self.phase + self.rate * t;
        let anomaly = mean + 2.0 * self.eccentricity * mean.sin();
        let (s, c) = anomaly.sin_cos();
        let f2 = 0.5 * s * s - 0.25;
        let f3 = -0.5 * s * c;
        [
            self.e[0] * f2 + self.e[1] * f3,
            self.i[0] * f2 + self.i[1] * f3,
            self.l[0] * f2 + self.l[1] * f3 + self.l[2] * s,
            self.g[0] * f2 + self.g[1] * f3 + self.g[2] * s,
            self.h[0] * f2 + self.h[1] * f3,
        ]
    }
}

/// Intermediate direction-cosine products for either third body.
struct ThirdBody {
    s: [f64; 7],
    z: [[f64; 3]; 4],
}
impl ThirdBody {
    fn new(
        o: Orbit,
        inclination: (f64, f64),
        perigee: (f64, f64),
        node: (f64, f64),
        strength: f64,
    ) -> Self {
        let (si, ci) = o.i.sin_cos();
        let (sg, cg) = o.argp.sin_cos();
        let (sin_i, cos_i) = inclination;
        let (sin_g, cos_g) = perigee;
        let (sin_h, cos_h) = node;
        let a1 = cos_g * cos_h + sin_g * cos_i * sin_h;
        let a3 = -sin_g * cos_h + cos_g * cos_i * sin_h;
        let a7 = -cos_g * sin_h + sin_g * cos_i * cos_h;
        let a8 = sin_g * sin_i;
        let a9 = sin_g * sin_h + cos_g * cos_i * cos_h;
        let a10 = cos_g * sin_i;
        let a2 = ci * a7 + si * a8;
        let a4 = ci * a9 + si * a10;
        let a5 = -si * a7 + ci * a8;
        let a6 = -si * a9 + ci * a10;
        let x1 = a1 * cg + a2 * sg;
        let x2 = a3 * cg + a4 * sg;
        let x3 = -a1 * sg + a2 * cg;
        let x4 = -a3 * sg + a4 * cg;
        let x5 = a5 * sg;
        let x6 = a6 * sg;
        let x7 = a5 * cg;
        let x8 = a6 * cg;
        let e2 = o.e * o.e;
        let beta = (1.0 - e2).sqrt();
        let z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
        let z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
        let z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
        let z1 = 2.0 * (3.0 * (a1 * a1 + a2 * a2) + z31 * e2) + (1.0 - e2) * z31;
        let z2 = 2.0 * (6.0 * (a1 * a3 + a2 * a4) + z32 * e2) + (1.0 - e2) * z32;
        let z3 = 2.0 * (3.0 * (a3 * a3 + a4 * a4) + z33 * e2) + (1.0 - e2) * z33;
        let z11 = -6.0 * a1 * a5 + e2 * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
        let z12 = -6.0 * (a1 * a6 + a3 * a5)
            + e2 * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
        let z13 = -6.0 * a3 * a6 + e2 * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
        let z21 = 6.0 * a2 * a5 + e2 * (24.0 * x1 * x5 - 6.0 * x3 * x7);
        let z22 = 6.0 * (a4 * a5 + a2 * a6)
            + e2 * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
        let z23 = 6.0 * a4 * a6 + e2 * (24.0 * x2 * x6 - 6.0 * x4 * x8);
        let s3 = strength / o.n;
        let s2 = -0.5 * s3 / beta;
        let s4 = s3 * beta;
        Self {
            s: [
                -15.0 * o.e * s4,
                s2,
                s3,
                s4,
                x1 * x3 + x2 * x4,
                x2 * x3 + x1 * x4,
                x2 * x4 - x1 * x3,
            ],
            z: [
                [z1, z2, z3],
                [z11, z12, z13],
                [z21, z22, z23],
                [z31, z32, z33],
            ],
        }
    }
    fn periodic(&self, e2: f64, eccentricity: f64, phase: f64, rate: f64) -> Periodic {
        let s = &self.s;
        let z = &self.z;
        Periodic {
            e: [2.0 * s[0] * s[5], 2.0 * s[0] * s[6]],
            i: [2.0 * s[1] * z[1][1], 2.0 * s[1] * (z[1][2] - z[1][0])],
            l: [
                -2.0 * s[2] * z[0][1],
                -2.0 * s[2] * (z[0][2] - z[0][0]),
                -2.0 * s[2] * (-21.0 - 9.0 * e2) * eccentricity,
            ],
            g: [
                2.0 * s[3] * z[3][1],
                2.0 * s[3] * (z[3][2] - z[3][0]),
                -18.0 * s[3] * eccentricity,
            ],
            h: [-2.0 * s[1] * z[2][1], -2.0 * s[1] * (z[2][2] - z[2][0])],
            phase,
            rate,
            eccentricity,
        }
    }
    fn secular(&self, o: Orbit, rate: f64) -> [f64; 5] {
        let (si, ci) = o.i.sin_cos();
        let s = &self.s;
        let z = &self.z;
        let raw_h = if o.i < 0.052359877 || o.i > PI - 0.052359877 {
            0.0
        } else {
            -rate * s[1] * (z[2][0] + z[2][2])
        };
        let h = if si == 0.0 { raw_h } else { raw_h / si };
        [
            s[0] * rate * s[4],
            s[1] * rate * (z[1][0] + z[1][2]),
            -rate * s[2] * (z[0][0] + z[0][2] - 14.0 - 6.0 * o.e * o.e),
            s[3] * rate * (z[3][0] + z[3][2] - 6.0) - ci * h,
            h,
        ]
    }
}

#[derive(Clone, Debug)]
pub(super) struct DeepSpace {
    solar: Periodic,
    lunar: Periodic,
    rates: [f64; 5],
    resonance: Option<Resonance>,
}
impl DeepSpace {
    pub(super) fn new(c: &Constants, epoch: f64) -> Self {
        let o = c.orbit;
        let (sn, cn) = o.node.sin_cos();
        let day = epoch + 18261.5;
        let lunar_node = (4.5236020 - 9.2422029e-4 * day) % TAU;
        let (s, cos) = lunar_node.sin_cos();
        let cos_i = 0.91375164 - 0.03568096 * cos;
        let sin_i = (1.0 - cos_i * cos_i).sqrt();
        let sin_h = 0.089683511 * s / sin_i;
        let cos_h = (1.0 - sin_h * sin_h).sqrt();
        let gamma = 5.8351514 + 0.0019443680 * day;
        let lunar_perigee = gamma
            + (0.39785416 * s / sin_i).atan2(cos_h * cos + 0.91744867 * sin_h * s)
            - lunar_node;
        let sun = ThirdBody::new(
            o,
            (0.39785416, 0.91744867),
            (-0.98088458, 0.1945905),
            (sn, cn),
            2.9864797e-6,
        );
        let moon = ThirdBody::new(
            o,
            (sin_i, cos_i),
            lunar_perigee.sin_cos(),
            (sn * cos_h - cn * sin_h, cos_h * cn + sin_h * sn),
            4.7968065e-7,
        );
        let solar = sun.periodic(
            o.e * o.e,
            0.01675,
            (6.2565837 + 0.017201977 * day) % TAU,
            SOLAR_RATE,
        );
        let lunar = moon.periodic(
            o.e * o.e,
            0.05490,
            (4.7199672 + 0.22997150 * day - gamma) % TAU,
            LUNAR_RATE,
        );
        let sun_rates = sun.secular(o, SOLAR_RATE);
        let moon_rates = moon.secular(o, LUNAR_RATE);
        let rates = core::array::from_fn(|i| sun_rates[i] + moon_rates[i]);
        let resonance = Resonance::new(c, rates);
        Self {
            solar,
            lunar,
            rates,
            resonance,
        }
    }

    pub(super) fn secular(&self, c: &Constants, t: f64, o: &mut Orbit) {
        o.e += self.rates[0] * t;
        o.i += self.rates[1] * t;
        o.m += self.rates[2] * t;
        o.argp += self.rates[3] * t;
        o.node += self.rates[4] * t;
        if let Some(resonance) = &self.resonance {
            resonance.integrate(c, t, o);
        }
    }

    pub(super) fn periodic(&self, t: f64, afspc: bool, o: &mut Orbit) {
        let sun = self.solar.at(t);
        let moon = self.lunar.at(t);
        let [de, di, dl, mut dg, mut dh] = core::array::from_fn(|i| sun[i] + moon[i]);
        o.e += de;
        o.i += di;
        let (si, ci) = o.i.sin_cos();
        if o.i >= 0.2 {
            dh /= si;
            dg -= ci * dh;
            o.argp += dg;
            o.node += dh;
            o.m += dl;
        } else {
            // Lyddane's nonsingular low-inclination correction. Keep the
            // original node branch continuous across atan2's ±π boundary.
            let (sn, cn) = o.node.sin_cos();
            let alpha = si * sn + dh * cn + di * ci * sn;
            let beta = si * cn - dh * sn + di * ci * cn;
            o.node %= TAU;
            if afspc && o.node < 0.0 {
                o.node += TAU;
            }
            let longitude = o.m + o.argp + dl + dg + (ci - di * si) * o.node;
            let old_node = o.node;
            o.node = alpha.atan2(beta);
            if afspc && o.node < 0.0 {
                o.node += TAU;
            }
            if (old_node - o.node).abs() > PI {
                o.node += if o.node < old_node { TAU } else { -TAU };
            }
            o.m += dl;
            o.argp = longitude - o.m - ci * o.node;
        }
    }
}

#[derive(Clone, Debug)]
struct Resonance {
    phase: f64,
    bias: f64,
    terms: ResonantTerms,
}
#[derive(Clone, Debug)]
enum ResonantTerms {
    Synchronous([f64; 3]),
    HalfDay([f64; 10]),
}
impl Resonance {
    fn new(c: &Constants, rates: [f64; 5]) -> Option<Self> {
        let o = c.orbit;
        let inverse_a = (o.n / c.geop.xke).powf(2.0 / 3.0);
        let (si, ci) = o.i.sin_cos();
        let e2 = o.e * o.e;
        if o.n > 0.0034906585 && o.n < 0.0052359877 {
            let g200 = 1.0 + e2 * (-2.5 + 0.8125 * e2);
            let g310 = 1.0 + 2.0 * e2;
            let g300 = 1.0 + e2 * (-6.0 + 6.60937 * e2);
            let f220 = 0.75 * (1.0 + ci) * (1.0 + ci);
            let f311 = 0.9375 * si * si * (1.0 + 3.0 * ci) - 0.75 * (1.0 + ci);
            let f330 = 1.875 * (1.0 + ci).powi(3);
            let base = 3.0 * o.n * o.n * inverse_a * inverse_a;
            Some(Self {
                phase: (o.m + o.node + o.argp - c.gsto) % TAU,
                bias: c.mdot + c.argpdot + c.nodedot - EARTH_RATE + rates[2] + rates[3] + rates[4]
                    - o.n,
                terms: ResonantTerms::Synchronous([
                    base * f311 * g310 * 2.1460748e-6 * inverse_a,
                    2.0 * base * f220 * g200 * 1.7891679e-6,
                    3.0 * base * f330 * g300 * 2.2123015e-7 * inverse_a,
                ]),
            })
        } else if (8.26e-3..=9.24e-3).contains(&o.n) && o.e >= 0.5 {
            Some(Self {
                phase: (o.m + 2.0 * o.node - 2.0 * c.gsto) % TAU,
                bias: c.mdot + rates[2] + 2.0 * (c.nodedot + rates[4] - EARTH_RATE) - o.n,
                terms: ResonantTerms::HalfDay(half_day_coefficients(o, inverse_a)),
            })
        } else {
            None
        }
    }

    fn derivatives(
        &self,
        c: &Constants,
        time: f64,
        longitude: f64,
        motion: f64,
    ) -> (f64, f64, f64) {
        let longitude_dot = motion + self.bias;
        let (motion_dot, motion_second) = match &self.terms {
            ResonantTerms::Synchronous(d) => {
                let angles = [
                    longitude - 0.13130908,
                    2.0 * (longitude - 2.8843198),
                    3.0 * (longitude - 0.37448087),
                ];
                (
                    d[0] * angles[0].sin() + d[1] * angles[1].sin() + d[2] * angles[2].sin(),
                    d[0] * angles[0].cos()
                        + 2.0 * d[1] * angles[1].cos()
                        + 3.0 * d[2] * angles[2].cos(),
                )
            }
            ResonantTerms::HalfDay(d) => {
                let perigee = c.orbit.argp + c.argpdot * time;
                let angles = [
                    2.0 * perigee + longitude - 5.7686396,
                    longitude - 5.7686396,
                    perigee + longitude - 0.95240898,
                    -perigee + longitude - 0.95240898,
                    2.0 * perigee + 2.0 * longitude - 1.8014998,
                    2.0 * longitude - 1.8014998,
                    perigee + longitude - 1.0508330,
                    -perigee + longitude - 1.0508330,
                    perigee + 2.0 * longitude - 4.4108898,
                    -perigee + 2.0 * longitude - 4.4108898,
                ];
                let mut first = 0.0;
                let mut second = 0.0;
                for i in 0..10 {
                    first += d[i] * angles[i].sin();
                    second +=
                        d[i] * angles[i].cos() * if matches!(i, 4 | 5 | 8 | 9) { 2.0 } else { 1.0 };
                }
                (first, second)
            }
        };
        (longitude_dot, motion_dot, motion_second * longitude_dot)
    }

    fn integrate(&self, c: &Constants, target: f64, o: &mut Orbit) {
        let mut time = 0.0;
        let mut longitude = self.phase;
        let mut motion = c.orbit.n;
        let step = if target < 0.0 { -720.0 } else { 720.0 };
        loop {
            let (ldot, ndot, nsecond) = self.derivatives(c, time, longitude, motion);
            let remainder = target - time;
            if remainder.abs() < 720.0 {
                o.n = motion + ndot * remainder + 0.5 * nsecond * remainder * remainder;
                let phase = longitude + ldot * remainder + 0.5 * ndot * remainder * remainder;
                let theta = (c.gsto + target * EARTH_RATE) % TAU;
                o.m = match self.terms {
                    ResonantTerms::Synchronous(_) => phase - o.node - o.argp + theta,
                    ResonantTerms::HalfDay(_) => phase - 2.0 * o.node + 2.0 * theta,
                };
                break;
            }
            longitude += ldot * step + ndot * 259200.0;
            motion += ndot * step + nsecond * 259200.0;
            time += step;
        }
    }
}

/// Eccentricity polynomials and inclination functions of the ten half-day
/// tesseral harmonics, in the published d2201..d5433 order.
fn half_day_coefficients(o: Orbit, inverse_a: f64) -> [f64; 10] {
    let e = o.e;
    let e2 = e * e;
    let e3 = e2 * e;
    let (si, ci) = o.i.sin_cos();
    let ci2 = ci * ci;
    let si2 = si * si;
    let g201 = -0.306 - (e - 0.64) * 0.440;
    let (g211, g310, g322, g410, g422, g520) = if e <= 0.65 {
        (
            3.616 - 13.2470 * e + 16.2900 * e2,
            -19.302 + 117.3900 * e - 228.4190 * e2 + 156.5910 * e3,
            -18.9068 + 109.7927 * e - 214.6334 * e2 + 146.5816 * e3,
            -41.122 + 242.6940 * e - 471.0940 * e2 + 313.9530 * e3,
            -146.407 + 841.8800 * e - 1629.014 * e2 + 1083.4350 * e3,
            -532.114 + 3017.977 * e - 5740.032 * e2 + 3708.2760 * e3,
        )
    } else {
        (
            -72.099 + 331.819 * e - 508.738 * e2 + 266.724 * e3,
            -346.844 + 1582.851 * e - 2415.925 * e2 + 1246.113 * e3,
            -342.585 + 1554.908 * e - 2366.899 * e2 + 1215.972 * e3,
            -1052.797 + 4758.686 * e - 7193.992 * e2 + 3651.957 * e3,
            -3581.690 + 16178.110 * e - 24462.770 * e2 + 12422.520 * e3,
            if e > 0.715 {
                -5149.66 + 29936.92 * e - 54087.36 * e2 + 31324.56 * e3
            } else {
                1464.74 - 4664.75 * e + 3763.64 * e2
            },
        )
    };
    let (g533, g521, g532) = if e < 0.7 {
        (
            -919.22770 + 4988.6100 * e - 9064.7700 * e2 + 5542.21 * e3,
            -822.71072 + 4568.6173 * e - 8491.4146 * e2 + 5337.524 * e3,
            -853.66600 + 4690.2500 * e - 8624.7700 * e2 + 5341.4 * e3,
        )
    } else {
        (
            -37995.780 + 161616.52 * e - 229838.20 * e2 + 109377.94 * e3,
            -51752.104 + 218913.95 * e - 309468.16 * e2 + 146349.42 * e3,
            -40023.880 + 170470.89 * e - 242699.48 * e2 + 115605.82 * e3,
        )
    };
    let f220 = 0.75 * (1.0 + 2.0 * ci + ci2);
    let f221 = 1.5 * si2;
    let f321 = 1.875 * si * (1.0 - 2.0 * ci - 3.0 * ci2);
    let f322 = -1.875 * si * (1.0 + 2.0 * ci - 3.0 * ci2);
    let f441 = 35.0 * si2 * f220;
    let f442 = 39.3750 * si2 * si2;
    let f522 = 9.84375
        * si
        * (si2 * (1.0 - 2.0 * ci - 5.0 * ci2) + 0.33333333 * (-2.0 + 4.0 * ci + 6.0 * ci2));
    let f523 = si
        * (4.92187512 * si2 * (-2.0 - 4.0 * ci + 10.0 * ci2)
            + 6.56250012 * (1.0 + 2.0 * ci - 3.0 * ci2));
    let f542 = 29.53125 * si * (2.0 - 8.0 * ci + ci2 * (-12.0 + 8.0 * ci + 10.0 * ci2));
    let f543 = 29.53125 * si * (-2.0 - 8.0 * ci + ci2 * (12.0 + 8.0 * ci - 10.0 * ci2));
    let base = 3.0 * o.n * o.n * inverse_a * inverse_a;
    let a3 = base * inverse_a;
    let a4 = a3 * inverse_a;
    let a5 = a4 * inverse_a;
    [
        base * 1.7891679e-6 * f220 * g201,
        base * 1.7891679e-6 * f221 * g211,
        a3 * 3.7393792e-7 * f321 * g310,
        a3 * 3.7393792e-7 * f322 * g322,
        2.0 * a4 * 7.3636953e-9 * f441 * g410,
        2.0 * a4 * 7.3636953e-9 * f442 * g422,
        a5 * 1.1428639e-7 * f522 * g520,
        a5 * 1.1428639e-7 * f523 * g532,
        2.0 * a5 * 2.1765803e-9 * f542 * g521,
        2.0 * a5 * 2.1765803e-9 * f543 * g533,
    ]
}
