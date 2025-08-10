#![allow(clippy::needless_range_loop)]

use crate::coordinates::{
    cartesian::Position,
    centers::Geocentric,
    frames::Ecliptic,
};

include!(concat!(env!("OUT_DIR"), "/elp_data.rs"));
use crate::astro::JulianDate;
use crate::units::{Arcseconds, Degrees, LengthUnit, Quantity, Radian, Radians};
use std::f64::consts::FRAC_PI_2;
use crate::bodies::solar_system::Moon;

// ====================
// Constants
// ====================
// Conversion helper using the strongly typed Units module
const ARCSEC2RAD: f64 = Arcseconds::new(1.0).to::<Radian>().value();
const A0: f64 = 384_747.980_644_895_4;
const ATH: f64 = 384_747.980_674_316_5;
const AM: f64 = 0.074_801_329_518;
const ALPHA: f64 = 0.002_571_881_335;
const DTASM: f64 = 2.0 * ALPHA / (3.0 * AM);
const PRECES: f64 = Arcseconds::new(5_029.096_6).to::<Radian>().value();

// LengthUnit conversions
const C1: f64 = 60.0;
const C2: f64 = 3600.0;

// Precession matrix coefficients (Laskar)
const P1: f64 = 0.10180391e-4;
const P2: f64 = 0.47020439e-6;
const P3: f64 = -0.5417367e-9;
const P4: f64 = -0.2507948e-11;
const P5: f64 = 0.463486e-14;
const Q1: f64 = -0.113469002e-3;
const Q2: f64 = 0.12372674e-6;
const Q3: f64 = 0.1265417e-8;
const Q4: f64 = -0.1371808e-11;
const Q5: f64 = -0.320334e-14;

// Corrections
const DELNU: f64 = Arcseconds::new(0.55604).to::<Radian>().value()
    / Arcseconds::new(1_732_559_343.736_04).to::<Radian>().value();
const DELE: f64 = Arcseconds::new(0.01789).to::<Radian>().value();
const DELG: f64 = Arcseconds::new(-0.08066).to::<Radian>().value();
const DELNP: f64 = Arcseconds::new(-0.06424).to::<Radian>().value()
    / Arcseconds::new(1_732_559_343.736_04).to::<Radian>().value();
const DELEP: f64 = Arcseconds::new(-0.12879).to::<Radian>().value();

// Delaunay arguments (series coefficients)
#[allow(clippy::all)]
const DEL: [[f64; 5]; 4] = [
    [ 5.198466741027443, 7771.377146811758394, -2.8449351621e-5,  3.1973462e-8,  -1.54365e-10 ],
    [-0.043125180208125, 628.301955168488007,  -2.6805348430e-6,  7.12676e-10,    7.2700e-13  ],
    [ 2.355555898265799, 8328.691426955554562,  1.57027757616e-4, 2.50411114e-7, -1.186339e-9 ],
    [ 1.627905233371468, 8433.466158130539043, -5.9392100004e-5, -4.949948e-9,    2.0217e-11  ],
];

// Fundamental lunar arguments: longitude offset and rate
const ZETA: [f64; 2] = [
    Degrees::new(218.0 + 18.0/60.0 + 59.955_71/3_600.0).to::<Radian>().value(),
    Arcseconds::new(1_732_559_343.736_04).to::<Radian>().value() + PRECES,
];

// Planetary argument coefficients
#[allow(clippy::all)]
const P_ARGS: [[f64; 2]; 8] = [
    [
        Degrees::new(252.0 + 15.0/C1 + 3.25986/C2).to::<Radian>().value(),
        Arcseconds::new(538_101_628.68898).to::<Radian>().value(),
    ],
    [
        Degrees::new(181.0 + 58.0/C1 + 47.28305/C2).to::<Radian>().value(),
        Arcseconds::new(210_664_136.43355).to::<Radian>().value(),
    ],
    [
        Degrees::new(100.0 + 27.0/60.0 + 59.22059/3600.0).to::<Radian>().value(),
        Arcseconds::new(129_597_742.27580).to::<Radian>().value(),
    ],
    [
        Degrees::new(355.0 + 25.0/C1 + 59.78866/C2).to::<Radian>().value(),
        Arcseconds::new(68_905_077.59284).to::<Radian>().value(),
    ],
    [
        Degrees::new(34.0  + 21.0/C1 + 5.34212/C2).to::<Radian>().value(),
        Arcseconds::new(10_925_660.42861).to::<Radian>().value(),
    ],
    [
        Degrees::new(50.0  +  4.0/C1 + 38.89694/C2).to::<Radian>().value(),
        Arcseconds::new(4_399_609.65932).to::<Radian>().value(),
    ],
    [
        Degrees::new(314.0 +  3.0/C1 + 18.01841/C2).to::<Radian>().value(),
        Arcseconds::new(1_542_481.19393).to::<Radian>().value(),
    ],
    [
        Degrees::new(304.0 + 20.0/C1 + 55.19575/C2).to::<Radian>().value(),
        Arcseconds::new(786_550.32074).to::<Radian>().value(),
    ],
];

// Lunar rotation series terms
const W1: [f64; 5] = [
    Degrees::new(218.0 + 18.0/60.0 + 59.955_71/3_600.0).to::<Radian>().value(),
    Arcseconds::new(1_732_559_343.736_04).to::<Radian>().value(),
    Arcseconds::new(-5.888_3).to::<Radian>().value(),
    Arcseconds::new(0.006_604).to::<Radian>().value(),
    Arcseconds::new(-0.000_031_69).to::<Radian>().value(),
];

// ====================
// Helpers
// ====================

/// Normalize an angle to the range [-π, π].
#[inline(always)]
fn normalize_angle(angle: Radians) -> Radians {
    angle.wrap_signed()
}

// ====================
// Main problem series (ELP1-3)
// ====================

#[inline(always)]
fn sum_main_problem_series(series: &[MainProblem], t: &[f64; 5], y_offset: Radians) -> f64 {
    // Precompute combined coefficient
    let delta_aux = DELNP - AM * DELNU;

    series.iter().fold(0.0, |accum, entry| {
        let tgv = entry.b[0] + DTASM * entry.b[4];
        let coeff = entry.a
            + tgv * delta_aux
            + entry.b[1] * DELG
            + entry.b[2] * DELE
            + entry.b[3] * DELEP;

        // Compute argument y
        let mut y = y_offset;
        for k in 0..5 {
            for i in 0..4 {
                y += Radians::new(entry.ilu[i] as f64 * DEL[i][k] * t[k]);
            }
        }

        accum + coeff * normalize_angle(y).sin()
    })
}

macro_rules! define_main_series {
    ($fn_name:ident, $series:path, $offset:expr) => {
        #[inline(always)]
        pub fn $fn_name(t: &[f64; 5]) -> f64 {
            sum_main_problem_series($series, t, $offset)
        }
    };
}

define_main_series!(sum_series_elp1, ELP1, Radians::new(0.0));
define_main_series!(sum_series_elp2, ELP2, Radians::new(0.0));
define_main_series!(sum_series_elp3, ELP3, Radians::new(FRAC_PI_2));

// ====================
// Earth perturbation series (ELP4-9,22-29,30-36)
// ====================

/// Sum Earth perturbation series; `scale_idx` multiplies amplitude by t[i] if Some(i)
#[inline(always)]
fn sum_earth_pert_series(series: &[EarthPert], t: &[f64; 5], scale_idx: Option<usize>) -> f64 {
    series.iter().fold(0.0, |accum, entry| {
        // Compute argument y
        let mut y = Degrees::new(entry.o).to::<Radian>();
        for k in 0..2 {
            y += Radians::new(entry.iz * ZETA[k] * t[k]);
            for i in 0..4 {
                y += Radians::new(entry.ilu[i] as f64 * DEL[i][k] * t[k]);
            }
        }
        let amplitude = if let Some(idx) = scale_idx { entry.a * t[idx] } else { entry.a };
        accum + amplitude * normalize_angle(y).sin()
    })
}

macro_rules! define_earth_series {
    ($fn_name:ident, $series:path, $scale:expr) => {
        #[inline(always)]
        pub fn $fn_name(t: &[f64; 5]) -> f64 {
            sum_earth_pert_series($series, t, $scale)
        }
    };
}

// ELP4-9
define_earth_series!(sum_series_elp4,  ELP4,  None);
define_earth_series!(sum_series_elp5,  ELP5,  None);
define_earth_series!(sum_series_elp6,  ELP6,  None);
define_earth_series!(sum_series_elp7,  ELP7,  Some(1));
define_earth_series!(sum_series_elp8,  ELP8,  Some(1));
define_earth_series!(sum_series_elp9,  ELP9,  Some(1));

// ELP22-29,30-33 (no scaling)
define_earth_series!(sum_series_elp22, ELP22, None);
define_earth_series!(sum_series_elp23, ELP23, None);
define_earth_series!(sum_series_elp24, ELP24, None);
define_earth_series!(sum_series_elp28, ELP28, None);
define_earth_series!(sum_series_elp29, ELP29, None);
define_earth_series!(sum_series_elp30, ELP30, None);
define_earth_series!(sum_series_elp31, ELP31, None);
define_earth_series!(sum_series_elp32, ELP32, None);
define_earth_series!(sum_series_elp33, ELP33, None);

// ELP25-27 (scale on t[1])
define_earth_series!(sum_series_elp25, ELP25, Some(1));
define_earth_series!(sum_series_elp26, ELP26, Some(1));
define_earth_series!(sum_series_elp27, ELP27, Some(1));

// ELP34-36 (scale on t[2])
define_earth_series!(sum_series_elp34, ELP34, Some(2));
define_earth_series!(sum_series_elp35, ELP35, Some(2));
define_earth_series!(sum_series_elp36, ELP36, Some(2));

// ====================
// Planet perturbation series (ELP10-21)
// ====================

/// Sum planetary perturbation series; `scale_o` multiplies `o` by t[1], `use_alt_del` switches argument terms
#[inline(always)]
fn sum_planet_pert_series(series: &[PlanetPert], t: &[f64; 5], scale_o: bool, use_alt_del: bool) -> f64 {
    series.iter().fold(0.0, |accum, entry| {
        let mut y = Degrees::new(entry.theta).to::<Radian>();
        for k in 0..2 {
            let delta = if use_alt_del {
                // restored original two-loop alt_del branch
                let mut delta = 0.0;
                for i in 0..4 {
                    delta += entry.ipla[i + 7] as f64 * DEL[i][k] * t[k];
                }
                for i in 0..7 {
                    delta += entry.ipla[i] as f64 * P_ARGS[i][k] * t[k];
                }
                delta
            } else {
                // original formula
                (entry.ipla[8] as f64 * DEL[0][k]
                    + entry.ipla[9] as f64 * DEL[2][k]
                    + entry.ipla[10] as f64 * DEL[3][k])
                    * t[k]
                    + (0..8).fold(0.0, |sum, i| sum + entry.ipla[i] as f64 * P_ARGS[i][k] * t[k])
            };
            y += Radians::new(delta);
        }
        let o_val = if scale_o { entry.o * t[1] } else { entry.o };
        accum + o_val * normalize_angle(y).sin()
    })
}

macro_rules! define_planet_series {
    ($fn_name:ident, $series:path, $scale_o:expr, $alt:expr) => {
        #[inline(always)]
        pub fn $fn_name(t: &[f64; 5]) -> f64 {
            sum_planet_pert_series($series, t, $scale_o, $alt)
        }
    };
}

// No scale, no alt
define_planet_series!(sum_series_elp10, ELP10, false, false);
define_planet_series!(sum_series_elp11, ELP11, false, false);
define_planet_series!(sum_series_elp12, ELP12, false, false);
// scale, no alt
define_planet_series!(sum_series_elp13, ELP13, true, false);
define_planet_series!(sum_series_elp14, ELP14, true, false);
define_planet_series!(sum_series_elp15, ELP15, true, false);
// no scale, alt
define_planet_series!(sum_series_elp16, ELP16, false, true);
define_planet_series!(sum_series_elp17, ELP17, false, true);
define_planet_series!(sum_series_elp18, ELP18, false, true);
// scale, alt
define_planet_series!(sum_series_elp19, ELP19, true, true);
define_planet_series!(sum_series_elp20, ELP20, true, true);
define_planet_series!(sum_series_elp21, ELP21, true, true);

// ====================
// Lunar position computation
// ====================

impl Moon {
    /// Get the geocentric ecliptic coordinates of the Moon for a given Julian date
    pub fn get_geo_position<U>(jd: JulianDate) -> Position<Geocentric, Ecliptic, U>
    where U: LengthUnit
    {
        let t1 = jd.julian_centuries().value();
        let t = [1.0, t1, t1.powi(2), t1.powi(3), t1.powi(4)];

        // Sum all series (36 values)
        let elp_values: [f64; 36] = [
            sum_series_elp1(&t), sum_series_elp2(&t), sum_series_elp3(&t),
            sum_series_elp4(&t), sum_series_elp5(&t), sum_series_elp6(&t),
            sum_series_elp7(&t), sum_series_elp8(&t), sum_series_elp9(&t),
            sum_series_elp10(&t), sum_series_elp11(&t), sum_series_elp12(&t),
            sum_series_elp13(&t), sum_series_elp14(&t), sum_series_elp15(&t),
            sum_series_elp16(&t), sum_series_elp17(&t), sum_series_elp18(&t),
            sum_series_elp19(&t), sum_series_elp20(&t), sum_series_elp21(&t),
            sum_series_elp22(&t), sum_series_elp23(&t), sum_series_elp24(&t),
            sum_series_elp25(&t), sum_series_elp26(&t), sum_series_elp27(&t),
            sum_series_elp28(&t), sum_series_elp29(&t), sum_series_elp30(&t),
            sum_series_elp31(&t), sum_series_elp32(&t), sum_series_elp33(&t),
            sum_series_elp34(&t), sum_series_elp35(&t), sum_series_elp36(&t),
        ];

        // Aggregate longitude, latitude, distance
        let a: f64 = elp_values.iter().step_by(3).sum();
        let b: f64 = elp_values.iter().skip(1).step_by(3).sum();
        let c: f64 = elp_values.iter().skip(2).step_by(3).sum();

        let lon = Radians::new(
            a * ARCSEC2RAD
                + W1[0]
                + W1[1] * t[1]
                + W1[2] * t[2]
                + W1[3] * t[3]
                + W1[4] * t[4],
        );
        let lat = Radians::new(b * ARCSEC2RAD);
        let distance = c * A0 / ATH;

        let x = distance * lat.cos();
        let y = x * lon.sin();
        let x = x * lon.cos();
        let z = distance * lat.sin();

        // Apply Laskar rotation
        let pw = (P1 + P2 * t[1] + P3 * t[2] + P4 * t[3] + P5 * t[4]) * t[1];
        let qw = (Q1 + Q2 * t[1] + Q3 * t[2] + Q4 * t[3] + Q5 * t[4]) * t[1];
        let ra = 2.0 * (1.0 - pw * pw - qw * qw).sqrt();
        let (pw2, qw2) = (1.0 - 2.0 * pw * pw, 1.0 - 2.0 * qw * qw);
        let pwqw = 2.0 * pw * qw;
        let pw = pw * ra;
        let qw = qw * ra;

        let x2 = pw2 * x + pwqw * y + pw * z;
        let y2 = pwqw * x + qw2 * y - qw * z;
        let z2 = -pw * x + qw * y + (pw2 + qw2 - 1.0) * z;

        Position::<Geocentric, Ecliptic, U>::new(
            Quantity::<U>::new(x2),
            Quantity::<U>::new(y2),
            Quantity::<U>::new(z2)
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::{Kilometer, KM};

    #[test]
    fn test_lunar_position_against_reference() {
        let pos = Moon::get_geo_position::<Kilometer>(JulianDate::J2000);

        let expected_x = -291608.0*KM;
        let expected_y = -274980.0*KM;
        let expected_z =  36271.2*KM;
        let tolerance = 1.0*KM;

        assert!((pos.x() - expected_x).abs() < tolerance, "X mismatch: got {}, expected {}", pos.x(), expected_x);
        assert!((pos.y() - expected_y).abs() < tolerance, "Y mismatch: got {}, expected {}", pos.y(), expected_y);
        assert!((pos.z() - expected_z).abs() < tolerance, "Z mismatch: got {}, expected {}", pos.z(), expected_z);
    }
}
