// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Native Vallado mean-element model. Internal lengths and times are normalized
//! to Earth radii and nominal minutes, as required by the published equations.
//! Dimensional quantities enter/leave through qtty and the typed TEME boundary.

use super::deep_space::DeepSpace;
use super::elements::NativeElements;
use super::{GravityModel, Sgp4Error};
use core::f64::consts::{PI, TAU};
use qtty::Kilometer;
#[cfg(not(feature = "std"))]
use qtty::{Real, Scalar, Transcendental};

/// Complete zonal geopotential in the SGP4 normalized unit system.
#[derive(Clone, Copy, Debug)]
pub(super) struct Geopotential {
    pub radius: Kilometer,
    pub xke: f64,
    pub j2: f64,
    pub j3: f64,
    pub j4: f64,
}
impl Geopotential {
    fn for_model(model: GravityModel) -> Self {
        match model {
            GravityModel::Wgs72 | GravityModel::Wgs72Iau => Self {
                radius: Kilometer::new(6378.135),
                xke: 0.07436691613317342,
                j2: 0.001082616,
                j3: -0.00000253881,
                j4: -0.00000165597,
            },
            GravityModel::Wgs84 => Self {
                radius: Kilometer::new(6378.137),
                xke: 0.07436685316871385,
                j2: 0.00108262998905,
                j3: -0.00000253215306,
                j4: -0.00000161098761,
            },
        }
    }
}

/// Dimensionless angular elements and mean motion in radians per nominal minute.
#[derive(Clone, Copy, Debug)]
pub(super) struct Orbit {
    pub e: f64,
    pub i: f64,
    pub node: f64,
    pub argp: f64,
    pub m: f64,
    pub n: f64,
}

#[derive(Clone, Debug)]
pub(super) struct Constants {
    pub orbit: Orbit,
    pub geop: Geopotential,
    pub mdot: f64,
    pub argpdot: f64,
    pub nodedot: f64,
    pub gsto: f64,
    bstar: f64,
    c1: f64,
    c4: f64,
    c5: f64,
    eta: f64,
    nodecf: f64,
    t2cof: f64,
    high_drag: Option<HighDrag>,
    deep: Option<DeepSpace>,
    afspc: bool,
}

#[derive(Clone, Debug)]
struct HighDrag {
    omgcof: f64,
    xmcof: f64,
    delmo: f64,
    sinmao: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    t3: f64,
    t4: f64,
    t5: f64,
}

pub(super) fn propagation_error(details: &str) -> Sgp4Error {
    Sgp4Error::Propagation {
        details: details.into(),
    }
}

fn invalid_if_non_finite(values: &[f64], details: &'static str) -> Result<(), Sgp4Error> {
    if values.iter().any(|value| !value.is_finite()) {
        return Err(Sgp4Error::InvalidElements {
            details: details.into(),
        });
    }
    Ok(())
}

impl Constants {
    /// Un-Kozai mean motion, initialize secular zonal and atmospheric terms,
    /// then select the 225-minute near-Earth/deep-space boundary.
    pub(super) fn new(elements: &NativeElements, model: GravityModel) -> Result<Self, Sgp4Error> {
        let geop = Geopotential::for_model(model);
        let mut orbit = elements.orbit;
        let (si, ci) = orbit.i.sin_cos();
        let ci2 = ci * ci;
        let e2 = orbit.e * orbit.e;
        let beta2 = 1.0 - e2;
        let beta = beta2.sqrt();
        let ak = (geop.xke / orbit.n).powf(2.0 / 3.0);
        let d1 = 0.75 * geop.j2 * (3.0 * ci2 - 1.0) / (beta * beta2);
        let delta = d1 / (ak * ak);
        let adel = ak * (1.0 - delta * delta - delta * (1.0 / 3.0 + 134.0 * delta * delta / 81.0));
        orbit.n /= 1.0 + d1 / (adel * adel);
        if !orbit.n.is_finite() || orbit.n <= 0.0 {
            return Err(Sgp4Error::InvalidElements {
                details: "invalid un-Kozai mean motion".into(),
            });
        }
        invalid_if_non_finite(
            &[beta2, beta, ak, d1, delta, adel, orbit.n],
            "non-finite un-Kozai constants",
        )?;
        let a = (geop.xke / orbit.n).powf(2.0 / 3.0);
        let p = a * beta2;
        let perigee = (a * (1.0 - orbit.e) - 1.0) * geop.radius.value();
        invalid_if_non_finite(&[a, p, perigee], "non-finite orbital constants")?;
        let deep_orbit = TAU / orbit.n >= 225.0;
        let simplified = perigee < 220.0 || deep_orbit;
        // Atmospheric density parameters change for very low perigees.
        let s_km = if perigee < 98.0 {
            20.0
        } else if perigee < 156.0 {
            perigee - 78.0
        } else {
            78.0
        };
        let s = 1.0 + s_km / geop.radius.value();
        let q = ((120.0 - s_km) / geop.radius.value()).powi(4);
        let inv_p2 = 1.0 / (p * p);
        let tsi = 1.0 / (a - s);
        let eta = a * orbit.e * tsi;
        let eta2 = eta * eta;
        let eeta = orbit.e * eta;
        let psi2 = (1.0 - eta2).abs();
        let coef = q * tsi.powi(4);
        let coef1 = coef / psi2.powf(3.5);
        let con41 = 3.0 * ci2 - 1.0;
        let sin2 = 1.0 - ci2;
        let c2 = coef1
            * orbit.n
            * (a * (1.0 + 1.5 * eta2 + eeta * (4.0 + eta2))
                + 0.375 * geop.j2 * tsi / psi2 * con41 * (8.0 + 3.0 * eta2 * (8.0 + eta2)));
        let c1 = elements.bstar * c2;
        let c3 = if orbit.e > 1e-4 {
            -2.0 * coef * tsi * geop.j3 / geop.j2 * orbit.n * si / orbit.e
        } else {
            0.0
        };
        let c4 = 2.0
            * orbit.n
            * coef1
            * a
            * beta2
            * (eta * (2.0 + 0.5 * eta2) + orbit.e * (0.5 + 2.0 * eta2)
                - geop.j2 * tsi / (a * psi2)
                    * (-3.0 * con41 * (1.0 - 2.0 * eeta + eta2 * (1.5 - 0.5 * eeta))
                        + 0.75
                            * sin2
                            * (2.0 * eta2 - eeta * (1.0 + eta2))
                            * (2.0 * orbit.argp).cos()));
        let c5 = 2.0 * coef1 * a * beta2 * (1.0 + 2.75 * (eta2 + eeta) + eeta * eta2);
        let ci4 = ci2 * ci2;
        let k1 = 1.5 * geop.j2 * inv_p2 * orbit.n;
        let k2 = 0.5 * k1 * geop.j2 * inv_p2;
        let k3 = -0.46875 * geop.j4 * inv_p2 * inv_p2 * orbit.n;
        let mdot = orbit.n
            + 0.5 * k1 * beta * con41
            + 0.0625 * k2 * beta * (13.0 - 78.0 * ci2 + 137.0 * ci4);
        let argpdot = -0.5 * k1 * (1.0 - 5.0 * ci2)
            + 0.0625 * k2 * (7.0 - 114.0 * ci2 + 395.0 * ci4)
            + k3 * (3.0 - 36.0 * ci2 + 49.0 * ci4);
        let nodedot =
            -k1 * ci + (0.5 * k2 * (4.0 - 19.0 * ci2) + 2.0 * k3 * (3.0 - 7.0 * ci2)) * ci;
        let nodecf = 3.5 * beta2 * (-k1 * ci) * c1;
        invalid_if_non_finite(
            &[
                s, q, inv_p2, tsi, eta, eta2, eeta, psi2, coef, coef1, c2, c1, c3, c4, c5, ci4, k1,
                k2, k3, mdot, argpdot, nodedot, nodecf,
            ],
            "non-finite SGP4 coefficients",
        )?;
        let high_drag = if simplified {
            None
        } else {
            let c1sq = c1 * c1;
            let d2 = 4.0 * a * tsi * c1sq;
            let temp = d2 * tsi * c1 / 3.0;
            let d3 = (17.0 * a + s) * temp;
            let d4 = 0.5 * temp * a * tsi * (221.0 * a + 31.0 * s) * c1;
            Some(HighDrag {
                omgcof: elements.bstar * c3 * orbit.argp.cos(),
                xmcof: if orbit.e > 1e-4 {
                    -(2.0 / 3.0) * coef * elements.bstar / eeta
                } else {
                    0.0
                },
                delmo: (1.0 + eta * orbit.m.cos()).powi(3),
                sinmao: orbit.m.sin(),
                d2,
                d3,
                d4,
                t3: d2 + 2.0 * c1sq,
                t4: 0.25 * (3.0 * d3 + c1 * (12.0 * d2 + 10.0 * c1sq)),
                t5: 0.2
                    * (3.0 * d4 + 12.0 * c1 * d3 + 6.0 * d2 * d2 + 15.0 * c1sq * (2.0 * d2 + c1sq)),
            })
        };
        if let Some(high_drag) = &high_drag {
            invalid_if_non_finite(
                &[
                    high_drag.omgcof,
                    high_drag.xmcof,
                    high_drag.delmo,
                    high_drag.sinmao,
                    high_drag.d2,
                    high_drag.d3,
                    high_drag.d4,
                    high_drag.t3,
                    high_drag.t4,
                    high_drag.t5,
                ],
                "non-finite high-drag coefficients",
            )?;
        }
        let afspc = model == GravityModel::Wgs72;
        let (epoch, gsto) = super::time::model_epoch(elements.datetime, afspc);
        let mut constants = Self {
            orbit,
            geop,
            mdot,
            argpdot,
            nodedot,
            gsto,
            bstar: elements.bstar,
            c1,
            c4,
            c5,
            eta,
            nodecf,
            t2cof: 1.5 * c1,
            high_drag,
            deep: None,
            afspc,
        };
        if deep_orbit {
            constants.deep = Some(DeepSpace::new(&constants, epoch));
        }
        Ok(constants)
    }

    /// Evaluate from immutable epoch constants. Resonance integration is local
    /// to each call, making results independent of request ordering.
    pub(super) fn propagate(
        &self,
        minutes: crate::qtty::Minutes,
    ) -> Result<([f64; 3], [f64; 3]), Sgp4Error> {
        let t = minutes.value();
        let mut o = self.orbit;
        o.m += self.mdot * t;
        o.argp += self.argpdot * t;
        o.node += self.nodedot * t + self.nodecf * t * t;
        let mut a_drag = 1.0 - self.c1 * t;
        let mut e_drag = self.bstar * self.c4 * t;
        let mut l_drag = self.t2cof * t * t;
        if let Some(h) = &self.high_drag {
            let correction =
                h.omgcof * t + h.xmcof * ((1.0 + self.eta * o.m.cos()).powi(3) - h.delmo);
            o.m += correction;
            o.argp -= correction;
            let t2 = t * t;
            let t3 = t2 * t;
            let t4 = t3 * t;
            a_drag -= h.d2 * t2 + h.d3 * t3 + h.d4 * t4;
            e_drag += self.bstar * self.c5 * (o.m.sin() - h.sinmao);
            l_drag += h.t3 * t3 + t4 * (h.t4 + t * h.t5);
        }
        if let Some(deep) = &self.deep {
            deep.secular(self, t, &mut o);
        }
        if o.n <= 0.0 || !o.n.is_finite() {
            return Err(propagation_error("non-positive mean motion"));
        }
        let a = (self.geop.xke / o.n).powf(2.0 / 3.0) * a_drag * a_drag;
        o.n = self.geop.xke / a.powf(1.5);
        o.e -= e_drag;
        if !o.e.is_finite() || o.e >= 1.0 || o.e < -0.001 {
            return Err(propagation_error("OutOfRangeEccentricity"));
        }
        o.e = o.e.max(1e-6);
        o.m += self.orbit.n * l_drag;
        let longitude = o.m + o.argp + o.node;
        // The node retains its sign; AFSPC applies its positive normalization
        // only inside the low-inclination Lyddane correction.
        o.node %= TAU;
        o.argp %= TAU;
        o.m = (longitude % TAU - o.argp - o.node) % TAU;
        if let Some(deep) = &self.deep {
            deep.periodic(t, self.afspc, &mut o);
            if o.i < 0.0 {
                o.i = -o.i;
                o.node += PI;
                o.argp -= PI;
            }
            if !(0.0..=1.0).contains(&o.e) {
                return Err(propagation_error("OutOfRangePerturbedEccentricity"));
            }
        }
        self.cartesian(a, o)
    }

    /// Long-period J3 correction, safeguarded eccentric-longitude solution,
    /// short-period J2 corrections, and the TEME orbital basis.
    fn cartesian(&self, a: f64, o: Orbit) -> Result<([f64; 3], [f64; 3]), Sgp4Error> {
        let (si, ci) = o.i.sin_cos();
        let j3j2 = self.geop.j3 / self.geop.j2;
        let denom = if (1.0 + ci).abs() > 1.5e-12 {
            1.0 + ci
        } else {
            1.5e-12
        };
        let aycof = -0.5 * j3j2 * si;
        let xlcof = -0.25 * j3j2 * si * (3.0 + 5.0 * ci) / denom;
        let inv_p = 1.0 / (a * (1.0 - o.e * o.e));
        let ax = o.e * o.argp.cos();
        let ay = o.e * o.argp.sin() + inv_p * aycof;
        let u = (o.m + o.argp + inv_p * xlcof * ax) % TAU;
        let mut eccentric_longitude = u;
        let (mut se, mut ce) = (0.0, 0.0);
        for _ in 0..10 {
            (se, ce) = eccentric_longitude.sin_cos();
            let step = ((u - ay * ce + ax * se - eccentric_longitude) / (1.0 - ax * ce - ay * se))
                .clamp(-0.95, 0.95);
            eccentric_longitude += step;
            if step.abs() < 1e-12 {
                break;
            }
        }
        let ecos = ax * ce + ay * se;
        let esin = ax * se - ay * ce;
        let e2 = ax * ax + ay * ay;
        let p = a * (1.0 - e2);
        if p < 0.0 {
            return Err(propagation_error("NegativeSemiLatusRectum"));
        }
        let r = a * (1.0 - ecos);
        let rdot = a.sqrt() * esin / r;
        let transverse = p.sqrt() / r;
        let beta = (1.0 - e2).sqrt();
        let k = esin / (1.0 + beta);
        let sin_u = a / r * (se - ay - ax * k);
        let cos_u = a / r * (ce - ax + ay * k);
        let sin_2u = 2.0 * cos_u * sin_u;
        let cos_2u = 1.0 - 2.0 * sin_u * sin_u;
        let temp1 = 0.5 * self.geop.j2 / p;
        let temp2 = temp1 / p;
        let ci2 = ci * ci;
        let con41 = 3.0 * ci2 - 1.0;
        let sin2 = 1.0 - ci2;
        let radius = r * (1.0 - 1.5 * temp2 * beta * con41) + 0.5 * temp1 * sin2 * cos_2u;
        let argument = sin_u.atan2(cos_u) - 0.25 * temp2 * (7.0 * ci2 - 1.0) * sin_2u;
        let node = o.node + 1.5 * temp2 * ci * sin_2u;
        let inclination = o.i + 1.5 * temp2 * ci * si * cos_2u;
        let vr = rdot - o.n * temp1 * sin2 * sin_2u / self.geop.xke;
        let vt = transverse + o.n * temp1 * (sin2 * cos_2u + 1.5 * con41) / self.geop.xke;
        let (su, cu) = argument.sin_cos();
        let (sn, cn) = node.sin_cos();
        let (si, ci) = inclination.sin_cos();
        let radial = [-sn * ci * su + cn * cu, cn * ci * su + sn * cu, si * su];
        let tangent = [-sn * ci * cu - cn * su, cn * ci * cu - sn * su, si * cu];
        let position = radial.map(|x| x * radius * self.geop.radius.value());
        let speed_scale = self.geop.radius.value() * self.geop.xke / 60.0;
        let velocity = core::array::from_fn(|i| (vr * radial[i] + vt * tangent[i]) * speed_scale);
        if position
            .iter()
            .chain(velocity.iter())
            .any(|x| !x.is_finite())
        {
            return Err(propagation_error("non-finite Cartesian state"));
        }
        Ok((position, velocity))
    }
}
