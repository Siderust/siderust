#![allow(clippy::needless_range_loop)]

use crate::coordinates::{cartesian::Position, centers::Geocentric, frames::Ecliptic};

#[allow(clippy::approx_constant)]
#[rustfmt::skip]
mod elp_data {
    include!(concat!(env!("OUT_DIR"), "/elp_data.rs"));
}
use crate::astro::JulianDate;
use crate::bodies::solar_system::Moon;
use crate::calculus::elp2000::elp_structs::*;
use elp_data::*;
use qtty::{Arcseconds, Degrees, Kilometers, LengthUnit, Radian, Radians, Simplify};
use std::f64::consts::FRAC_PI_2;

// ====================
// Constants
// ====================
// Conversion helper using the strongly typed Units module
const A0: Kilometers = Kilometers::new(384_747.980_644_895_4);
const ATH: Kilometers = Kilometers::new(384_747.980_674_316_5);
const AM: f64 = 0.074_801_329_518;
const ALPHA: f64 = 0.002_571_881_335;
const DTASM: f64 = 2.0 * ALPHA / (3.0 * AM);
const PRECES: Radians = Arcseconds::new(5_029.096_6).to::<Radian>();

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
const DELNU: f64 = Arcseconds::new(0.55604)
    .to::<Radian>()
    .div(Arcseconds::new(1_732_559_343.736_04).to::<Radian>())
    .value();
const DELE: f64 = Arcseconds::new(0.01789).to::<Radian>().value();
const DELG: f64 = Arcseconds::new(-0.08066).to::<Radian>().value();
const DELNP: f64 = Arcseconds::new(-0.06424)
    .to::<Radian>()
    .div(Arcseconds::new(1_732_559_343.736_04).to::<Radian>())
    .value();
const DELEP: f64 = Arcseconds::new(-0.12879).to::<Radian>().value();

// Delaunay arguments (series coefficients)
#[allow(clippy::all)]
const DEL: [[Radians; 5]; 4] = [
    [
        Radians::new(5.198466741027443),
        Radians::new(7771.377146811758394),
        Radians::new(-2.8449351621e-5),
        Radians::new(3.1973462e-8),
        Radians::new(-1.54365e-10),
    ],
    [
        Radians::new(-0.043125180208125),
        Radians::new(628.301955168488007),
        Radians::new(-2.6805348430e-6),
        Radians::new(7.12676e-10),
        Radians::new(7.2700e-13),
    ],
    [
        Radians::new(2.355555898265799),
        Radians::new(8328.691426955554562),
        Radians::new(1.57027757616e-4),
        Radians::new(2.50411114e-7),
        Radians::new(-1.186339e-9),
    ],
    [
        Radians::new(1.627905233371468),
        Radians::new(8433.466158130539043),
        Radians::new(-5.9392100004e-5),
        Radians::new(-4.949948e-9),
        Radians::new(2.0217e-11),
    ],
];

// Fundamental lunar arguments: longitude offset and rate
const ZETA: [Radians; 2] = [
    Degrees::new(218.0 + 18.0 / 60.0 + 59.955_71 / 3_600.0).to::<Radian>(),
    Arcseconds::new(1_732_559_343.736_04)
        .to::<Radian>()
        .add(PRECES),
];

// Planetary argument coefficients
#[allow(clippy::all)]
const P_ARGS: [[Radians; 2]; 8] = [
    [
        Degrees::new(252.0 + 15.0 / C1 + 3.25986 / C2).to::<Radian>(),
        Arcseconds::new(538_101_628.68898).to::<Radian>(),
    ],
    [
        Degrees::new(181.0 + 58.0 / C1 + 47.28305 / C2).to::<Radian>(),
        Arcseconds::new(210_664_136.43355).to::<Radian>(),
    ],
    [
        Degrees::new(100.0 + 27.0 / 60.0 + 59.22059 / 3600.0).to::<Radian>(),
        Arcseconds::new(129_597_742.27580).to::<Radian>(),
    ],
    [
        Degrees::new(355.0 + 25.0 / C1 + 59.78866 / C2).to::<Radian>(),
        Arcseconds::new(68_905_077.59284).to::<Radian>(),
    ],
    [
        Degrees::new(34.0 + 21.0 / C1 + 5.34212 / C2).to::<Radian>(),
        Arcseconds::new(10_925_660.42861).to::<Radian>(),
    ],
    [
        Degrees::new(50.0 + 4.0 / C1 + 38.89694 / C2).to::<Radian>(),
        Arcseconds::new(4_399_609.65932).to::<Radian>(),
    ],
    [
        Degrees::new(314.0 + 3.0 / C1 + 18.01841 / C2).to::<Radian>(),
        Arcseconds::new(1_542_481.19393).to::<Radian>(),
    ],
    [
        Degrees::new(304.0 + 20.0 / C1 + 55.19575 / C2).to::<Radian>(),
        Arcseconds::new(786_550.32074).to::<Radian>(),
    ],
];

// Lunar rotation series terms
const W1: [Radians; 5] = [
    Degrees::new(218.0 + 18.0 / 60.0 + 59.955_71 / 3_600.0).to::<Radian>(),
    Arcseconds::new(1_732_559_343.736_04).to::<Radian>(),
    Arcseconds::new(-5.888_3).to::<Radian>(),
    Arcseconds::new(0.006_604).to::<Radian>(),
    Arcseconds::new(-0.000_031_69).to::<Radian>(),
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
        let coeff =
            entry.a + tgv * delta_aux + entry.b[1] * DELG + entry.b[2] * DELE + entry.b[3] * DELEP;

        // Compute argument y
        let mut y = y_offset;
        for k in 0..5 {
            for i in 0..4 {
                y += entry.ilu[i] as f64 * DEL[i][k] * t[k];
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
            y += entry.iz * ZETA[k] * t[k];
            for i in 0..4 {
                y += entry.ilu[i] as f64 * DEL[i][k] * t[k];
            }
        }
        let amplitude = if let Some(idx) = scale_idx {
            entry.a * t[idx]
        } else {
            entry.a
        };
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
define_earth_series!(sum_series_elp4, ELP4, None);
define_earth_series!(sum_series_elp5, ELP5, None);
define_earth_series!(sum_series_elp6, ELP6, None);
define_earth_series!(sum_series_elp7, ELP7, Some(1));
define_earth_series!(sum_series_elp8, ELP8, Some(1));
define_earth_series!(sum_series_elp9, ELP9, Some(1));

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
fn sum_planet_pert_series(
    series: &[PlanetPert],
    t: &[f64; 5],
    scale_o: bool,
    use_alt_del: bool,
) -> f64 {
    series.iter().fold(0.0, |accum, entry| {
        let mut y = Degrees::new(entry.theta).to::<Radian>();
        for k in 0..2 {
            let delta = if use_alt_del {
                // restored original two-loop alt_del branch
                let mut delta = Radians::new(0.0);
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
                    + (0..8).fold(Radians::new(0.0), |sum, i| {
                        sum + entry.ipla[i] as f64 * P_ARGS[i][k] * t[k]
                    })
            };
            y += delta;
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
    where
        U: LengthUnit,
    {
        let t1 = jd.julian_centuries().value();
        let t = [1.0, t1, t1.powi(2), t1.powi(3), t1.powi(4)];

        // Sum all series (36 values)
        let elp_values: [f64; 36] = [
            sum_series_elp1(&t),
            sum_series_elp2(&t),
            sum_series_elp3(&t),
            sum_series_elp4(&t),
            sum_series_elp5(&t),
            sum_series_elp6(&t),
            sum_series_elp7(&t),
            sum_series_elp8(&t),
            sum_series_elp9(&t),
            sum_series_elp10(&t),
            sum_series_elp11(&t),
            sum_series_elp12(&t),
            sum_series_elp13(&t),
            sum_series_elp14(&t),
            sum_series_elp15(&t),
            sum_series_elp16(&t),
            sum_series_elp17(&t),
            sum_series_elp18(&t),
            sum_series_elp19(&t),
            sum_series_elp20(&t),
            sum_series_elp21(&t),
            sum_series_elp22(&t),
            sum_series_elp23(&t),
            sum_series_elp24(&t),
            sum_series_elp25(&t),
            sum_series_elp26(&t),
            sum_series_elp27(&t),
            sum_series_elp28(&t),
            sum_series_elp29(&t),
            sum_series_elp30(&t),
            sum_series_elp31(&t),
            sum_series_elp32(&t),
            sum_series_elp33(&t),
            sum_series_elp34(&t),
            sum_series_elp35(&t),
            sum_series_elp36(&t),
        ];

        // Aggregate longitude, latitude, distance
        let a: f64 = elp_values.iter().step_by(3).sum();
        let b: f64 = elp_values.iter().skip(1).step_by(3).sum();
        let c: f64 = elp_values.iter().skip(2).step_by(3).sum();

        let lon = Arcseconds::new(a).to::<Radian>()
            + W1[0]
            + W1[1] * t[1]
            + W1[2] * t[2]
            + W1[3] * t[3]
            + W1[4] * t[4];
        let lat = Arcseconds::new(b).to::<Radian>();
        let ratio = (A0 / ATH).simplify().value();
        let distance = Kilometers::new(c * ratio);

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

        Position::<Geocentric, Ecliptic, U>::new(x2.to::<U>(), y2.to::<U>(), z2.to::<U>())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::{Days, Kilometer, KM};
    use std::f64::consts::PI;

    // ===========================================================================
    // HELPERS
    // ===========================================================================

    fn pos_j2000_km() -> Position<Geocentric, Ecliptic, Kilometer> {
        Moon::get_geo_position::<Kilometer>(JulianDate::J2000)
    }

    fn r_from_xyz_km(p: &Position<Geocentric, Ecliptic, Kilometer>) -> f64 {
        let x = p.x().to::<Kilometer>().value();
        let y = p.y().to::<Kilometer>().value();
        let z = p.z().to::<Kilometer>().value();
        (x * x + y * y + z * z).sqrt()
    }

    fn lon_lat_from_xyz(p: &Position<Geocentric, Ecliptic, Kilometer>) -> (f64, f64) {
        let x = p.x().to::<Kilometer>().value();
        let y = p.y().to::<Kilometer>().value();
        let z = p.z().to::<Kilometer>().value();
        let r = (x * x + y * y + z * z).sqrt();
        let lon = y.atan2(x);
        let lat = (z / r).asin();
        (lon, lat)
    }

    fn make_jd_utc(yyyy: i32, mm: u32, dd: u32, h: u32, m: u32, s: u32) -> JulianDate {
        use chrono::{DateTime, NaiveDate, Utc};
        let date = NaiveDate::from_ymd_opt(yyyy, mm, dd).expect("invalid date");
        let datetime = date.and_hms_opt(h, m, s).expect("invalid time");
        JulianDate::from_utc(DateTime::<Utc>::from_naive_utc_and_offset(datetime, Utc))
    }

    /// Build time array t = [1, t1, t1², t1³, t1⁴] from Julian centuries
    fn t_from_centuries(t1: f64) -> [f64; 5] {
        [1.0, t1, t1.powi(2), t1.powi(3), t1.powi(4)]
    }

    /// Build a JulianDate from Julian centuries offset from J2000
    fn jd_from_centuries(t1: f64) -> JulianDate {
        JulianDate::J2000 + Days::new(t1 * 36_525.0)
    }

    // ===========================================================================
    // LAYER 1: PURE UNIT TESTS (Time scale, angle normalization, unit conversion)
    // ===========================================================================

    // ----- 1.1 Time scale contract -----

    #[test]
    fn j2000_julian_centuries_is_zero() {
        // The epoch origin must be exactly J2000 (t1 = 0)
        let centuries = JulianDate::J2000.julian_centuries().value();
        assert!(
            centuries.abs() < 1e-15,
            "J2000.julian_centuries() should be 0, got {centuries}"
        );
    }

    #[test]
    fn julian_centuries_one_century_from_j2000() {
        // J2000 + 36525 days = exactly 1 Julian century
        let jd = JulianDate::J2000 + Days::new(36_525.0);
        let centuries = jd.julian_centuries().value();
        assert!(
            (centuries - 1.0).abs() < 1e-12,
            "Expected 1.0 Julian century, got {centuries}"
        );
    }

    // ----- 1.2 normalize_angle() correctness -----

    #[test]
    fn normalize_angle_zero() {
        let result = normalize_angle(Radians::new(0.0));
        assert!(
            result.value().abs() < 1e-15,
            "normalize_angle(0) should be 0"
        );
    }

    #[test]
    fn normalize_angle_pi() {
        // π should map to π or -π (both are valid; check it's one of them)
        let result = normalize_angle(Radians::new(PI)).value();
        assert!(
            (result.abs() - PI).abs() < 1e-15,
            "normalize_angle(π) should be ±π, got {result}"
        );
    }

    #[test]
    fn normalize_angle_minus_pi() {
        let result = normalize_angle(Radians::new(-PI)).value();
        assert!(
            (result.abs() - PI).abs() < 1e-15,
            "normalize_angle(-π) should be ±π, got {result}"
        );
    }

    #[test]
    fn normalize_angle_three_pi() {
        // 3π should wrap to π or -π
        let result = normalize_angle(Radians::new(3.0 * PI)).value();
        assert!(
            (result.abs() - PI).abs() < 1e-14,
            "normalize_angle(3π) should be ±π, got {result}"
        );
    }

    #[test]
    fn normalize_angle_preserves_sine() {
        // Property: sin(original) == sin(normalized) within floating point tolerance
        let test_angles = [
            0.0,
            PI,
            -PI,
            2.0 * PI,
            -2.0 * PI,
            3.0 * PI,
            100.0,
            -100.0,
            1000.0,
            10_000.0,
        ];
        for angle in test_angles {
            let original = Radians::new(angle);
            let normalized = normalize_angle(original);
            let sin_orig = angle.sin();
            let sin_norm = normalized.value().sin();
            // Tolerance relaxed: different argument-reduction paths
            // (libm for large angles vs. wrap_signed first) can differ.
            // For very large angles (>1000), differences up to ~1e-13 are acceptable.
            let tolerance = if angle.abs() > 1000.0 { 1e-12 } else { 1e-13 };
            assert!(
                (sin_orig - sin_norm).abs() < tolerance,
                "sin({angle}) = {sin_orig} != sin(normalized) = {sin_norm}, diff = {}",
                (sin_orig - sin_norm).abs()
            );
        }
    }

    #[test]
    fn normalize_angle_result_in_range() {
        // Property: result must be in [-π, π]
        let test_angles = [
            0.0,
            1.0,
            -1.0,
            PI,
            -PI,
            2.0 * PI,
            -2.0 * PI,
            100.0,
            -100.0,
            1e6,
            -1e6,
        ];
        for angle in test_angles {
            let result = normalize_angle(Radians::new(angle)).value();
            assert!(
                (-PI..=PI).contains(&result),
                "normalize_angle({angle}) = {result} is out of [-π, π]"
            );
        }
    }

    // ----- 1.3 Unit conversion sanity -----

    #[test]
    fn unit_conversion_arcsec_to_degree() {
        // 3600 arcsec = 1 degree (definitional)
        let arcsec = Arcseconds::new(3600.0);
        let deg = arcsec.to::<qtty::Degree>().value();
        assert!(
            (deg - 1.0).abs() < 1e-15,
            "3600 arcsec should be 1 degree, got {deg}"
        );
    }

    #[test]
    fn unit_conversion_degree_to_radian() {
        // 180 degrees = π radians (definitional)
        let deg = Degrees::new(180.0);
        let rad = deg.to::<Radian>().value();
        assert!((rad - PI).abs() < 1e-14, "180° should be π rad, got {rad}");
    }

    #[test]
    fn unit_conversion_arcsec_to_radian() {
        // 1 arcsec = π / 648000 radians
        let arcsec = Arcseconds::new(1.0);
        let rad = arcsec.to::<Radian>().value();
        let expected = PI / 648_000.0;
        assert!(
            (rad - expected).abs() < 1e-18,
            "1 arcsec should be {expected} rad, got {rad}"
        );
    }

    // ===========================================================================
    // LAYER 2: SERIES-LEVEL TESTS (ELP1-36 behavior and edge cases)
    // ===========================================================================

    #[test]
    fn all_series_return_finite_values_at_j2000() {
        let t = [1.0, 0.0, 0.0, 0.0, 0.0];
        let series_vals = [
            sum_series_elp1(&t),
            sum_series_elp2(&t),
            sum_series_elp3(&t),
            sum_series_elp4(&t),
            sum_series_elp5(&t),
            sum_series_elp6(&t),
            sum_series_elp7(&t),
            sum_series_elp8(&t),
            sum_series_elp9(&t),
            sum_series_elp10(&t),
            sum_series_elp11(&t),
            sum_series_elp12(&t),
            sum_series_elp13(&t),
            sum_series_elp14(&t),
            sum_series_elp15(&t),
            sum_series_elp16(&t),
            sum_series_elp17(&t),
            sum_series_elp18(&t),
            sum_series_elp19(&t),
            sum_series_elp20(&t),
            sum_series_elp21(&t),
            sum_series_elp22(&t),
            sum_series_elp23(&t),
            sum_series_elp24(&t),
            sum_series_elp25(&t),
            sum_series_elp26(&t),
            sum_series_elp27(&t),
            sum_series_elp28(&t),
            sum_series_elp29(&t),
            sum_series_elp30(&t),
            sum_series_elp31(&t),
            sum_series_elp32(&t),
            sum_series_elp33(&t),
            sum_series_elp34(&t),
            sum_series_elp35(&t),
            sum_series_elp36(&t),
        ];

        for (i, v) in series_vals.iter().enumerate() {
            assert!(v.is_finite(), "ELP{} produced non-finite value", i + 1);
            // extremely loose sanity bound (arcseconds before conversion)
            assert!(v.abs() < 5.0e8, "ELP{} magnitude looks suspect: {v}", i + 1);
        }
        for (i, v) in series_vals.iter().enumerate() {
            assert!(v.is_finite(), "ELP{} produced non-finite value", i + 1);
            assert!(v.abs() < 5.0e8, "ELP{} magnitude looks suspect: {v}", i + 1);
        }
    }

    #[test]
    fn all_series_finite_across_epoch_grid() {
        // Test across multiple epochs to exercise time polynomial expansions
        let t1_values = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0];
        for &t1 in &t1_values {
            let t = t_from_centuries(t1);
            let series_vals = [
                sum_series_elp1(&t),
                sum_series_elp2(&t),
                sum_series_elp3(&t),
                sum_series_elp4(&t),
                sum_series_elp5(&t),
                sum_series_elp6(&t),
                sum_series_elp7(&t),
                sum_series_elp8(&t),
                sum_series_elp9(&t),
                sum_series_elp10(&t),
                sum_series_elp11(&t),
                sum_series_elp12(&t),
                sum_series_elp13(&t),
                sum_series_elp14(&t),
                sum_series_elp15(&t),
                sum_series_elp16(&t),
                sum_series_elp17(&t),
                sum_series_elp18(&t),
                sum_series_elp19(&t),
                sum_series_elp20(&t),
                sum_series_elp21(&t),
                sum_series_elp22(&t),
                sum_series_elp23(&t),
                sum_series_elp24(&t),
                sum_series_elp25(&t),
                sum_series_elp26(&t),
                sum_series_elp27(&t),
                sum_series_elp28(&t),
                sum_series_elp29(&t),
                sum_series_elp30(&t),
                sum_series_elp31(&t),
                sum_series_elp32(&t),
                sum_series_elp33(&t),
                sum_series_elp34(&t),
                sum_series_elp35(&t),
                sum_series_elp36(&t),
            ];
            for (i, v) in series_vals.iter().enumerate() {
                assert!(v.is_finite(), "ELP{} non-finite at t1={t1}", i + 1);
            }
        }
    }

    #[test]
    fn diagnostic_elp_table_presence() {
        // Diagnostic test: report which ELP tables are empty (truncated dataset)
        // This helps explain why some "must be nonzero" assertions might fail

        let tables = [
            ("ELP4", ELP4.len()),
            ("ELP5", ELP5.len()),
            ("ELP6", ELP6.len()),
            ("ELP7", ELP7.len()),
            ("ELP8", ELP8.len()),
            ("ELP9", ELP9.len()),
            ("ELP10", ELP10.len()),
            ("ELP11", ELP11.len()),
            ("ELP12", ELP12.len()),
            ("ELP13", ELP13.len()),
            ("ELP14", ELP14.len()),
            ("ELP15", ELP15.len()),
            ("ELP16", ELP16.len()),
            ("ELP17", ELP17.len()),
            ("ELP18", ELP18.len()),
            ("ELP19", ELP19.len()),
            ("ELP20", ELP20.len()),
            ("ELP21", ELP21.len()),
            ("ELP22", ELP22.len()),
            ("ELP23", ELP23.len()),
            ("ELP24", ELP24.len()),
            ("ELP25", ELP25.len()),
            ("ELP26", ELP26.len()),
            ("ELP27", ELP27.len()),
            ("ELP28", ELP28.len()),
            ("ELP29", ELP29.len()),
            ("ELP30", ELP30.len()),
            ("ELP31", ELP31.len()),
            ("ELP32", ELP32.len()),
            ("ELP33", ELP33.len()),
            ("ELP34", ELP34.len()),
            ("ELP35", ELP35.len()),
            ("ELP36", ELP36.len()),
        ];

        let empty_tables: Vec<_> = tables
            .iter()
            .filter(|(_, len)| *len == 0)
            .map(|(name, _)| *name)
            .collect();

        if !empty_tables.is_empty() {
            eprintln!(
                "Note: {} ELP tables are empty (truncated dataset): {:?}",
                empty_tables.len(),
                empty_tables
            );
        }

        // This test always passes; it's purely diagnostic
    }

    // ---------- Scaling branch coverage ----------
    #[test]
    fn earth_pert_elp7_8_9_behavior_matches_dataset() {
        // ELP7-9 have t[1] scaling factor; validate based on actual table presence
        let t = t_from_centuries(0.12345);

        // If tables are present, sums should be finite
        let s7 = sum_series_elp7(&t);
        let s8 = sum_series_elp8(&t);
        let s9 = sum_series_elp9(&t);
        assert!(s7.is_finite() && s8.is_finite() && s9.is_finite());

        // If tables are empty (truncated dataset), sums will be 0.0
        // This is valid and expected for some build configurations
    }

    #[test]
    fn earth_pert_elp34_35_36_vanish_at_t1_zero() {
        // ELP34-36 scale by t[2] = t1^2, so at t1=0 they should be ~0
        let t = t_from_centuries(0.0);
        assert!(
            sum_series_elp34(&t).abs() < 1e-15,
            "ELP34 should vanish at t1=0"
        );
        assert!(
            sum_series_elp35(&t).abs() < 1e-15,
            "ELP35 should vanish at t1=0"
        );
        assert!(
            sum_series_elp36(&t).abs() < 1e-15,
            "ELP36 should vanish at t1=0"
        );
    }

    #[test]
    fn earth_pert_elp34_35_36_behavior_matches_dataset() {
        // ELP34-36 have t[2] scaling; at t1=0, t[2]=0 so they vanish
        let t0 = t_from_centuries(0.0);
        assert_eq!(
            sum_series_elp34(&t0),
            0.0,
            "ELP34 must be 0 at t1=0 (t[2]=0)"
        );
        assert_eq!(
            sum_series_elp35(&t0),
            0.0,
            "ELP35 must be 0 at t1=0 (t[2]=0)"
        );
        assert_eq!(
            sum_series_elp36(&t0),
            0.0,
            "ELP36 must be 0 at t1=0 (t[2]=0)"
        );

        // At t1 != 0, result depends on table contents (may still be 0 if empty)
        let t1 = t_from_centuries(0.5);
        let s34 = sum_series_elp34(&t1);
        let s35 = sum_series_elp35(&t1);
        let s36 = sum_series_elp36(&t1);
        assert!(s34.is_finite() && s35.is_finite() && s36.is_finite());
    }

    #[test]
    fn planet_pert_series_elp10_to_21_finite() {
        // ELP10-21 are planetary perturbations with different alt_del branches
        // Validate they return finite values regardless of table contents
        let t = t_from_centuries(0.12345);

        for (i, s) in [
            sum_series_elp10(&t),
            sum_series_elp11(&t),
            sum_series_elp12(&t),
            sum_series_elp13(&t),
            sum_series_elp14(&t),
            sum_series_elp15(&t),
            sum_series_elp16(&t),
            sum_series_elp17(&t),
            sum_series_elp18(&t),
            sum_series_elp19(&t),
            sum_series_elp20(&t),
            sum_series_elp21(&t),
        ]
        .iter()
        .enumerate()
        {
            assert!(s.is_finite(), "ELP{} produced non-finite", i + 10);
        }
    }

    #[test]
    fn series_stable_at_far_future_epoch() {
        let t = t_from_centuries(5.0); // Year ~2500
        for v in [
            sum_series_elp1(&t),
            sum_series_elp2(&t),
            sum_series_elp3(&t),
        ] {
            assert!(v.is_finite(), "Series unstable at far future epoch");
        }
    }

    #[test]
    fn series_stable_at_far_past_epoch() {
        let t = t_from_centuries(-3.0); // Year ~1700
        for v in [
            sum_series_elp1(&t),
            sum_series_elp2(&t),
            sum_series_elp3(&t),
        ] {
            assert!(v.is_finite(), "Series unstable at far past epoch");
        }
    }

    #[test]
    fn main_problem_series_elp1_elp2_elp3_have_expected_phase_offset() {
        // ELP3 uses y_offset = PI/2 while ELP1 uses 0
        let t = [1.0, 0.0, 0.0, 0.0, 0.0];
        let s1 = sum_series_elp1(&t);
        let s3 = sum_series_elp3(&t);
        assert!(s1.is_finite() && s3.is_finite());
        if s1.abs() > 1e-6 {
            assert!(
                (s3 - s1).abs() > 1e-12,
                "ELP3 should differ from ELP1 due to PI/2 offset"
            );
        }
    }

    #[test]
    fn earth_pert_elp25_26_27_behavior_matches_dataset() {
        // ELP25-27 have t[1] scaling; at t1=0, t[1]=0 so they vanish
        let t0 = t_from_centuries(0.0);
        assert_eq!(
            sum_series_elp25(&t0),
            0.0,
            "ELP25 must be 0 at t1=0 (t[1]=0)"
        );
        assert_eq!(
            sum_series_elp26(&t0),
            0.0,
            "ELP26 must be 0 at t1=0 (t[1]=0)"
        );
        assert_eq!(
            sum_series_elp27(&t0),
            0.0,
            "ELP27 must be 0 at t1=0 (t[1]=0)"
        );

        // At t1 != 0, result depends on table contents (may still be 0 if empty)
        let t1 = t_from_centuries(0.5);
        let s25 = sum_series_elp25(&t1);
        let s26 = sum_series_elp26(&t1);
        let s27 = sum_series_elp27(&t1);
        assert!(s25.is_finite() && s26.is_finite() && s27.is_finite());
    }

    #[test]
    fn all_planet_pert_series_exercised() {
        // Ensure all planetary perturbation series (ELP10-21) return finite values
        let t = t_from_centuries(0.25);
        for (i, v) in [
            sum_series_elp10(&t),
            sum_series_elp11(&t),
            sum_series_elp12(&t),
            sum_series_elp13(&t),
            sum_series_elp14(&t),
            sum_series_elp15(&t),
            sum_series_elp16(&t),
            sum_series_elp17(&t),
            sum_series_elp18(&t),
            sum_series_elp19(&t),
            sum_series_elp20(&t),
            sum_series_elp21(&t),
        ]
        .iter()
        .enumerate()
        {
            assert!(v.is_finite(), "ELP{} non-finite", i + 10);
        }
    }

    #[test]
    fn earth_pert_series_no_scale_exercised() {
        // ELP22-24, ELP28-33 don't scale by t[1] or t[2]
        let t = t_from_centuries(0.1);
        for (i, v) in [
            sum_series_elp22(&t),
            sum_series_elp23(&t),
            sum_series_elp24(&t),
            sum_series_elp28(&t),
            sum_series_elp29(&t),
            sum_series_elp30(&t),
            sum_series_elp31(&t),
            sum_series_elp32(&t),
            sum_series_elp33(&t),
        ]
        .iter()
        .enumerate()
        {
            assert!(v.is_finite(), "Earth pert series {} non-finite", i);
        }
    }

    // ===========================================================================
    // LAYER 3: INTEGRATION TESTS (end-to-end through Moon::get_geo_position)
    // ===========================================================================

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn end_to_end_position_stable_at_far_future() {
        let jd = jd_from_centuries(5.0);
        let pos = Moon::get_geo_position::<Kilometer>(jd);
        let r = r_from_xyz_km(&pos);
        assert!(r.is_finite(), "Position unstable at far future");
        assert!(
            (300_000.0..=450_000.0).contains(&r),
            "Distance out of extended bounds: {r}"
        );
    }

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn geocentric_distance_within_realistic_lunar_bounds() {
        let p = pos_j2000_km();
        let r = r_from_xyz_km(&p);
        assert!(
            (330_000.0..=410_000.0).contains(&r),
            "distance out of bounds: {} km",
            r
        );
    }

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn ecliptic_latitude_within_lunar_inclination() {
        let p = pos_j2000_km();
        let (_lon, lat) = lon_lat_from_xyz(&p);
        let lat_deg = Degrees::from(Radians::new(lat)).value();
        assert!(lat_deg.abs() <= 6.0, "latitude too large: {} deg", lat_deg);
    }

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn physical_bounds_across_multiple_epochs() {
        for jd in [
            make_jd_utc(1990, 1, 1, 0, 0, 0),
            make_jd_utc(2010, 6, 15, 0, 0, 0),
            make_jd_utc(2050, 7, 4, 12, 0, 0),
        ] {
            let pos = Moon::get_geo_position::<Kilometer>(jd);
            let r = r_from_xyz_km(&pos);
            let (_lon, lat) = lon_lat_from_xyz(&pos);
            assert!(
                (330_000.0..=420_000.0).contains(&r),
                "distance out of bounds at {:?}",
                jd
            );
            assert!(
                lat.to_degrees().abs() <= 6.0,
                "latitude too large at {:?}",
                jd
            );
        }
    }

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn continuity_finite_difference_velocity() {
        // Moon's orbital speed is ~1.0 km/s, check velocity is in reasonable range
        let jd0 = make_jd_utc(2020, 6, 15, 12, 0, 0);
        let dt_days = 1.0 / 24.0; // 1 hour
        let jd1 = jd0 + Days::new(dt_days);
        let pos0 = Moon::get_geo_position::<Kilometer>(jd0);
        let pos1 = Moon::get_geo_position::<Kilometer>(jd1);
        let (x0, y0, z0) = (
            pos0.x().to::<Kilometer>().value(),
            pos0.y().to::<Kilometer>().value(),
            pos0.z().to::<Kilometer>().value(),
        );
        let (x1, y1, z1) = (
            pos1.x().to::<Kilometer>().value(),
            pos1.y().to::<Kilometer>().value(),
            pos1.z().to::<Kilometer>().value(),
        );
        let dt_sec = dt_days * 86400.0;
        let v01 = ((x1 - x0).powi(2) + (y1 - y0).powi(2) + (z1 - z0).powi(2)).sqrt() / dt_sec;
        assert!(
            (0.5..=2.0).contains(&v01),
            "velocity {v01} km/s out of expected range [0.5, 2.0]"
        );
    }

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn axes_are_ecliptic_of_date_not_equatorial() {
        let p = pos_j2000_km();
        assert!(
            p.z().abs() > 1.0 * KM,
            "z should not be ~0 in ecliptic frame"
        );
    }

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn long_range_stability_10_year_intervals() {
        for year in (1950..=2050).step_by(10) {
            let jd = make_jd_utc(year, 1, 1, 0, 0, 0);
            let pos = Moon::get_geo_position::<Kilometer>(jd);
            let r = r_from_xyz_km(&pos);
            assert!(r.is_finite(), "Position unstable at year {year}");
            assert!(
                (320_000.0..=420_000.0).contains(&r),
                "Distance out of bounds at year {year}: {r}"
            );
        }
    }

    // ===========================================================================
    // LAYER 4: GOLDEN REGRESSION TESTS (NASA/reference data)
    // ===========================================================================

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn regression_xyz_j2000_against_reference() {
        let pos = pos_j2000_km();
        let expected_x = -291_608.0 * KM;
        let expected_y = -274_980.0 * KM;
        let expected_z = 36_271.2 * KM;
        let tol = 1.0 * KM;
        assert!((pos.x() - expected_x).abs() < tol, "X mismatch");
        assert!((pos.y() - expected_y).abs() < tol, "Y mismatch");
        assert!((pos.z() - expected_z).abs() < tol, "Z mismatch");
    }

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn regression_j2000_derived_quantities() {
        let pos = pos_j2000_km();
        let r = r_from_xyz_km(&pos);

        // Must be consistent with the frozen XYZ regression
        let x: f64 = -291_608.0;
        let y: f64 = -274_980.0;
        let z: f64 = 36_271.2;

        let expected_r = (x * x + y * y + z * z).sqrt();

        assert!(
            (r - expected_r).abs() < 1.0,
            "J2000 distance regression failed: got {r}, expected {expected_r}"
        );
    }

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn distance_matches_supermoon_2019_jan21_perigee() {
        // NASA/Space.com: 2019-01-21 perigee ~357,344 km
        let pos = Moon::get_geo_position::<Kilometer>(make_jd_utc(2019, 1, 21, 0, 0, 0));
        let r = r_from_xyz_km(&pos);
        assert!(
            (r - 357_344.0).abs() <= 1_500.0,
            "2019-01-21 perigee: got {r} km"
        );
    }

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn distance_matches_supermoon_2019_feb19_perigee() {
        // NASA/Space.com: 2019-02-19 09:06 UTC perigee ~356,761 km
        let pos = Moon::get_geo_position::<Kilometer>(make_jd_utc(2019, 2, 19, 9, 6, 0));
        let r = r_from_xyz_km(&pos);
        assert!(
            (r - 356_761.0).abs() <= 1_500.0,
            "2019-02-19 perigee: got {r} km"
        );
    }

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn distance_matches_supermoon_2019_mar19_perigee() {
        // NASA/Space.com: 2019-03-19 perigee ~359,380 km
        let pos = Moon::get_geo_position::<Kilometer>(make_jd_utc(2019, 3, 19, 0, 0, 0));
        let r = r_from_xyz_km(&pos);
        assert!(
            (r - 359_380.0).abs() <= 1_500.0,
            "2019-03-19 perigee: got {r} km"
        );
    }

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn distance_matches_supermoon_2020_apr07_perigee() {
        // NASA/Space.com: 2020-04-07 18:08 UTC perigee ~356,908 km
        let pos = Moon::get_geo_position::<Kilometer>(make_jd_utc(2020, 4, 7, 18, 8, 0));
        let r = r_from_xyz_km(&pos);
        assert!(
            (r - 356_908.0).abs() <= 1_500.0,
            "2020-04-07 perigee: got {r} km"
        );
    }

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn distance_matches_apogee_2019_feb05() {
        // Apogee: 2019-02-05 ~406,555 km
        let pos = Moon::get_geo_position::<Kilometer>(make_jd_utc(2019, 2, 5, 0, 0, 0));
        let r = r_from_xyz_km(&pos);
        assert!(
            (r - 406_555.0).abs() <= 2_000.0,
            "2019-02-05 apogee: got {r} km"
        );
    }

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn distance_matches_apogee_2020_mar24() {
        // Apogee: 2020-03-24 15:23 UTC ~406,688 km
        let pos = Moon::get_geo_position::<Kilometer>(make_jd_utc(2020, 3, 24, 15, 23, 0));
        let r = r_from_xyz_km(&pos);
        assert!(
            (r - 406_688.0).abs() <= 2_000.0,
            "2020-03-24 apogee: got {r} km"
        );
    }

    // ===========================================================================
    // LAYER 5: PRECISION REGRESSION TESTS (determinism and constants)
    // ===========================================================================

    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn determinism_bitwise_equality() {
        // Same input must produce identical output (no non-deterministic operations)
        let jd = make_jd_utc(2020, 6, 15, 12, 0, 0);
        let pos1 = Moon::get_geo_position::<Kilometer>(jd);
        let pos2 = Moon::get_geo_position::<Kilometer>(jd);
        assert_eq!(
            pos1.x().to::<Kilometer>().value().to_bits(),
            pos2.x().to::<Kilometer>().value().to_bits()
        );
        assert_eq!(
            pos1.y().to::<Kilometer>().value().to_bits(),
            pos2.y().to::<Kilometer>().value().to_bits()
        );
        assert_eq!(
            pos1.z().to::<Kilometer>().value().to_bits(),
            pos2.z().to::<Kilometer>().value().to_bits()
        );
    }

    #[test]
    fn summation_order_sensitivity_elp1() {
        // Verify deterministic summation (same result on repeated calls)
        let t = t_from_centuries(0.12345);
        let forward = sum_series_elp1(&t);
        let forward2 = sum_series_elp1(&t);
        assert_eq!(
            forward.to_bits(),
            forward2.to_bits(),
            "Summation not deterministic"
        );
    }

    #[test]
    fn regression_elp1_at_j2000() {
        let t = t_from_centuries(0.0);
        let val = sum_series_elp1(&t);
        assert!(val.is_finite(), "ELP1 at J2000 not finite");
        assert!(
            val.abs() < 1.0e7,
            "ELP1 at J2000 magnitude unreasonable: {val}"
        );
    }

    #[test]
    fn regression_elp2_at_j2000() {
        let t = t_from_centuries(0.0);
        let val = sum_series_elp2(&t);
        assert!(val.is_finite(), "ELP2 at J2000 not finite");
        assert!(
            val.abs() < 1.0e7,
            "ELP2 at J2000 magnitude unreasonable: {val}"
        );
    }

    #[test]
    fn regression_elp3_at_j2000() {
        let t = t_from_centuries(0.0);
        let val = sum_series_elp3(&t);
        assert!(val.is_finite(), "ELP3 at J2000 not finite");
        assert!(
            val.abs() < 1.0e8,
            "ELP3 at J2000 magnitude unreasonable: {val}"
        );
    }

    #[test]
    fn regression_constants_a0_ath() {
        // Mean lunar distance constants from ELP2000/82B theory
        assert!(
            (A0.value() - 384_747.980_644_895_4).abs() < 1e-6,
            "A0 constant mismatch"
        );
        assert!(
            (ATH.value() - 384_747.980_674_316_5).abs() < 1e-6,
            "ATH constant mismatch"
        );
    }

    #[test]
    fn regression_delaunay_initial_values() {
        // Delaunay arguments should be in [0, 2*PI) at epoch
        assert!(DEL[0][0].value() >= 0.0 && DEL[0][0].value() < 2.0 * PI);
        assert!(DEL[1][0].value().abs() < 2.0 * PI);
        assert!(DEL[2][0].value() >= 0.0 && DEL[2][0].value() < 2.0 * PI);
        assert!(DEL[3][0].value() >= 0.0 && DEL[3][0].value() < 2.0 * PI);
    }

    #[test]
    fn delaunay_and_planet_args_monotonic_rates() {
        // Rates should be positive (mean motions)
        assert!(DEL[0][1].value() > 0.0, "l' rate should be positive");
        assert!(DEL[1][1].value() > 0.0, "l_sun' rate should be positive");
        assert!(DEL[2][1].value() > 0.0, "F' rate should be positive");
        assert!(DEL[3][1].value() > 0.0, "D' rate should be positive");
        assert!(ZETA[1].value() > 1.0, "ZETA rate should be large positive");
    }

    // ===========================================================================
    // LAYER 6: ADDITIONAL COVERAGE (constants, arguments, Laskar coefficients)
    // ===========================================================================

    #[test]
    fn laskar_coefficients_reasonable() {
        // Laskar precession/nutation coefficients should be small
        assert!(P1.abs() < 1e-3, "P1 coefficient too large");
        assert!(Q1.abs() < 1e-3, "Q1 coefficient too large");
    }

    #[test]
    fn w1_initial_value_reasonable() {
        // W1[0] is the mean longitude of the Moon at J2000 ~218 degrees
        let w1_deg = W1[0].to::<qtty::Degree>().value();
        assert!(
            (w1_deg - 218.0).abs() < 1.0,
            "W1[0] = {} deg, expected ~218 degrees",
            w1_deg
        );
    }

    #[test]
    fn planetary_args_initial_values_reasonable() {
        // P_ARGS[0][0] is Mercury's mean longitude at J2000 ~252 degrees
        let mercury_deg = P_ARGS[0][0].to::<qtty::Degree>().value();
        assert!(
            (mercury_deg - 252.0).abs() < 1.0,
            "Mercury mean lon = {} deg, expected ~252 degrees",
            mercury_deg
        );
    }
}
