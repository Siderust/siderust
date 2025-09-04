#![allow(clippy::needless_range_loop)]

use crate::coordinates::{cartesian::Position, centers::Geocentric, frames::Ecliptic};

#[allow(clippy::approx_constant)]
mod elp_data {
include!(concat!(env!("OUT_DIR"), "/elp_data.rs"));
}
use elp_data::*;
use crate::calculus::elp2000::elp_structs::*;
use crate::astro::JulianDate;
use crate::bodies::solar_system::Moon;
use crate::units::{Arcseconds, Degrees, Kilometers, LengthUnit, Radian, Radians, Simplify};
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
    use crate::units::{Kilometer, KM};

    // ---------- Helpers for tests ----------
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
        // Ecliptic lon/lat of date, in radians
        let x = p.x().to::<Kilometer>().value();
        let y = p.y().to::<Kilometer>().value();
        let z = p.z().to::<Kilometer>().value();
        let r = (x * x + y * y + z * z).sqrt();
        let lon = y.atan2(x);
        let lat = (z / r).asin();
        (lon, lat)
    }

    // ---------- Series sanity ----------
    #[test]
    fn all_series_return_finite_values_at_j2000() {
        // t for J2000 -> t1 = 0; so [1, 0, 0, 0, 0]
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
    }

    // ---------- Physical bounds ----------
    #[test]
    fn geocentric_distance_within_realistic_lunar_bounds() {
        // Perigee–apogee ~ 356,000–407,000 km; allow slack for model/epoch -> 330k–410k
        let p = pos_j2000_km();
        let r = r_from_xyz_km(&p);
        assert!(
            (330_000.0..=410_000.0).contains(&r),
            "distance out of bounds: {} km",
            r
        );
    }

    #[test]
    fn ecliptic_latitude_within_lunar_inclination() {
        // Lunar ecliptic latitude never exceeds ~5.15° in magnitude; allow safety margin to 6°
        let p = pos_j2000_km();
        let (_lon, lat) = lon_lat_from_xyz(&p);
        let lat_deg = Degrees::from(Radians::new(lat)).value();
        assert!(
            lat_deg.abs() <= 6.0,
            "latitude too large: {} deg (expected |b| <= ~6°)",
            lat_deg
        );
    }

    // ---------- Regression against reference XYZ (yours) ----------
    #[test]
    fn regression_xyz_j2000_against_reference() {
        // Your existing regression numbers
        let pos = pos_j2000_km();

        let expected_x = -291_608.0 * KM;
        let expected_y = -274_980.0 * KM;
        let expected_z = 36_271.2 * KM;
        let tol = 1.0 * KM;

        assert!(
            (pos.x() - expected_x).abs() < tol,
            "X mismatch: got {}, expected {}",
            pos.x(),
            expected_x
        );
        assert!(
            (pos.y() - expected_y).abs() < tol,
            "Y mismatch: got {}, expected {}",
            pos.y(),
            expected_y
        );
        assert!(
            (pos.z() - expected_z).abs() < tol,
            "Z mismatch: got {}, expected {}",
            pos.z(),
            expected_z
        );
    }

    // ---------- Continuity / smoothness sanity ----------
    //
    // We can’t construct arbitrary Julian dates without knowing your API,
    // so we at least ensure the internal argument builder behaves smoothly.
    #[test]
    fn delaunay_and_planet_args_monotonic_rates() {
        // The “rates” in DEL and ZETA are constructed so that the 1st derivative
        // terms are the dominant parts of the argument evolution at small t.
        // Check their signs / magnitudes are sensible at t1 = 0.
        assert!(DEL[0][1].value() > 0.0); // l' > 0
        assert!(DEL[1][1].value() > 0.0); // l_sun' > 0
        assert!(DEL[2][1].value() > 0.0); // F' > 0
        assert!(DEL[3][1].value() > 0.0); // D' > 0

        // ZETA[1] is the huge lunar mean longitude rate (plus precession)
        assert!(ZETA[1].value() > 1.0); // definitely positive and large (radians / century)
    }

    // ---------- Frame sanity ----------
    #[test]
    fn axes_are_ecliptic_of_date_not_equatorial() {
        // Very rough: at J2000, the ecliptic pole is ~23.44° from the z of equatorial.
        // Here we just ensure the returned vector is **not** trivially lying in the
        // equatorial xy-plane by accident. For J2000 your expected Z is ~36,271 km.
        let p = pos_j2000_km();
        assert!(p.z().abs() > 1.0 * KM, "z should not be ~0 in ecliptic frame");
    }

    // ---------- Fast path correctness: sum_main_problem_series offset variants ----------
    #[test]
    fn main_problem_series_elp1_elp2_elp3_have_expected_phase_offset() {
        // At J2000, t = [1,0,0,0,0], so ELP3 uses y_offset = π/2.
        // We can at least assert ELP3 != ELP1 unless the sums happen to be zero.
        let t = [1.0, 0.0, 0.0, 0.0, 0.0];
        let s1 = sum_series_elp1(&t);
        let s2 = sum_series_elp2(&t);
        let s3 = sum_series_elp3(&t);

        assert!(s1.is_finite() && s2.is_finite() && s3.is_finite());
        // If ELP1 is non-zero, a π/2 phase shift generally changes the value.
        if s1.abs() > 1e-6 {
            assert!(
                (s3 - s1).abs() > 1e-12,
                "ELP3 should differ from ELP1 due to π/2 offset"
            );
        }
    }


    fn make_jd_utc(yyyy: i32, mm: u32, dd: u32, h: u32, m: u32, s: u32) -> JulianDate {
        use chrono::{DateTime, NaiveDate, Utc};

        let date = NaiveDate::from_ymd_opt(yyyy, mm, dd)
            .expect("invalid date");
        let datetime = date
            .and_hms_opt(h, m, s)
            .expect("invalid time");

        JulianDate::from_utc(DateTime::<Utc>::from_naive_utc_and_offset(datetime, Utc))
    }


    #[test]
    fn distance_matches_supermoon_2019_feb19_perigee() {
        // Perigee: 2019-02-19 09:06 UTC, distance ≈ 356,761 km
        // Sources: NASA/EarthSky/Space.com report this same value. 
        // Tolerance is generous to account for UTC/TT, model differences, and rounding.

        let pos = Moon::get_geo_position::<crate::units::Kilometer>(make_jd_utc(2019, 2, 19, 9, 6, 0));
        let r = {
            let x = pos.x().to::<Kilometer>().value();
            let y = pos.y().to::<Kilometer>().value();
            let z = pos.z().to::<Kilometer>().value();
            (x*x + y*y + z*z).sqrt()
        };

        let expected_km = 356_761.0;
        let tol_km = 1_500.0; // ~0.4%
        assert!(
            (r - expected_km).abs() <= tol_km,
            "perigee distance mismatch: got {r:.1} km vs {expected_km:.1}±{tol_km:.0} km"
        );
    }

    #[test]
    fn distance_matches_supermoon_2020_apr07_perigee() {
        // Perigee: 2020-04-07 18:08 UTC, distance ≈ 356,908 km
        // Tolerance same reasoning as above.

        let pos = Moon::get_geo_position::<crate::units::Kilometer>(make_jd_utc(2020, 4, 7, 18, 8, 0));
        let r = {
            let x = pos.x().to::<Kilometer>().value();
            let y = pos.y().to::<Kilometer>().value();
            let z = pos.z().to::<Kilometer>().value();
            (x*x + y*y + z*z).sqrt()
        };

        let expected_km = 356_908.0;
        let tol_km = 1_500.0;
        assert!(
            (r - expected_km).abs() <= tol_km,
            "perigee distance mismatch: got {r:.1} km vs {expected_km:.1}±{tol_km:.0} km"
        );
    }

    #[test]
    fn distance_matches_apogee_2020_mar24() {
        // Apogee: 2020-03-24 15:23 UTC (varies by source within minutes), distance ≈ 406,688 km
        // This checks the far end of the orbit.

        let pos = Moon::get_geo_position::<crate::units::Kilometer>(make_jd_utc(2020, 3, 24, 15, 23, 0));
        let r = {
            let x = pos.x().to::<Kilometer>().value();
            let y = pos.y().to::<Kilometer>().value();
            let z = pos.z().to::<Kilometer>().value();
            (x*x + y*y + z*z).sqrt()
        };

        let expected_km = 406_688.0;
        let tol_km = 2_000.0; // apogee distances are more sensitive; give a bit more slack
        assert!(
            (r - expected_km).abs() <= tol_km,
            "apogee distance mismatch: got {r:.1} km vs {expected_km:.1}±{tol_km:.0} km"
        );
    }

}
