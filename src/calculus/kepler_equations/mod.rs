// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Kepler Equations Module
//!
//! This module provides fast, reliable solutions of Kepler’s equation and the
//! derived heliocentric position of an object on an elliptic orbit.
//!
//! ## What Is Kepler’s&nbsp;Equation?
//!
//! For an orbit with eccentricity `e` (`0 ≤ e < 1`) the **elliptic** form is
//!
//! ```text
//! E − e sin E = M
//! ```
//!
//! where  
//! * `M` – **mean anomaly** (proportional to time since periapsis)  
//! * `E` – **eccentric anomaly** (unknown)  
//!
//! Solving for `E` is the gateway to many classic tasks in celestial mechanics
//! – from drawing an orbit on screen to steering a spacecraft.
//!
//! ## Where & Why Is It Used?
//!
//! * Ephemeris generators and planetarium software.  
//! * On-board GN&C loops for satellites and interplanetary probes.  
//! * Numerical propagators that require many millions of calls per run.  
//! * Educational tools that visualise orbital motion in real time.
//!
//! ## Mathematical & Algorithmic Background
//!
//! Only the elliptic variant is implemented for the moment;
//! the hyperbolic (`e > 1`) and parabolic (`e ≈ 1`) forms can be added later
//! without changing the public interface.
//!
//! ### Numerical strategy
//! * **Newton–Raphson** – quadratic convergence, typically 3–5 iterations.  
//!   Falls back if the derivative becomes too small or the iteration cap
//!   (`MAX_NEWTON_ITERS = 30`) is reached.  
//! * **Bisection** – linear but monotonic convergence, guaranteed within
//!   `MAX_BISECTION_ITERS = 100` iterations.  
//!
//! Both methods share the absolute tolerance  
//! `TOLERANCE = 1 × 10⁻¹⁴ rad` (≈ 3 mm on 1 AstronomicalUnits).
//!
//! ## Public API
//!
//! ```rust
//! use siderust::calculus::kepler_equations::solve_keplers_equation;
//! use qtty::*;
//!
//! let m = 1.047 * RAD; // mean anomaly
//! let e = 0.0167;             // eccentricity
//! let e_anomaly = solve_keplers_equation(m, e);
//! ```
//!
//! To advance from anomaly to Cartesian position use  
//! `calculate_orbit_position(&Orbit, julian_date)` which wraps the standard
//! true-anomaly conversion, radius vector, and three-axis rotation into the
//! ecliptic heliocentric frame.
//!
//! ## Limitations & Future Work
//!
//! * Only elliptical orbits are supported.  
//! * Gravitational parameter is fixed to *Sun + test particle*.  
//!   Extend this if barycentric or relativistic accuracy is required.  
//! * PRs adding hyperbolic/parabolic solvers, SIMD paths, or extended-precision
//!   back ends are welcome.
//!
//! ## References
//!
//! - Danby, J.&nbsp;M. *Fundamentals of Celestial Mechanics*, 2nd ed. (1992).
//! - Vallado, D.&nbsp;A. *Fundamentals of Astrodynamics and Applications*,
//!   4th ed. (2013).
//! - Meeus, J. *Astronomical Algorithms*, 2nd ed. (1998) – ch. 30.
//!
//! ## Module Layout
//!
//! | Item                         | Visibility | Purpose                         |
//! | ---------------------------- | ---------- | --------------------------------|
//! | `solve_keplers_equation`     | **pub**    | High-level front door           |
//! | `calculate_orbit_position`   | **pub**    | Elliptic position at Julian day |
//! | `newton_kepler`, `bisection_kepler` | crate-private | Numerical kernels |
//! | `kepler_equation_residual`, `kepler_equation_derivative` | crate-private | Helpers |

use crate::astro::orbit::Orbit;
use crate::time::JulianDate;
use crate::coordinates::cartesian::position::Ecliptic;
use qtty::*;
use std::f64::consts::PI;

const TOLERANCE: Radians = Radians::new(1e-14);
const MAX_NEWTON_ITERS: usize = 30;
const MAX_BISECTION_ITERS: usize = 100;

/// Computes the residual f(E) = E - e*sin(E) - M for a given E.
#[inline]
fn kepler_equation_residual(e_anomaly: Radians, e: f64, m: Radians) -> Radians {
    e_anomaly - Radians::new(e * e_anomaly.sin()) - m
}

/// Computes the derivative f'(E) = 1 - e * cos(E).
#[inline]
fn kepler_equation_derivative(e_anomaly: Radians, e: f64) -> f64 {
    1.0 - e * e_anomaly.cos()
}

/// Attempts to solve Kepler's equation via Newton-Raphson.
///
/// Returns `Some(e_anomaly)` on success, or `None` if:
/// - Convergence isn't reached within `MAX_NEWTON_ITERS`, or
/// - The derivative is too small (risking large updates).
fn newton_kepler(m: Radians, e: f64, initial_guess: Radians) -> Option<Radians> {
    let mut e_anomaly = initial_guess;

    for _ in 0..MAX_NEWTON_ITERS {
        let f = kepler_equation_residual(e_anomaly, e, m);
        let f_prime = kepler_equation_derivative(e_anomaly, e);

        if f_prime.abs() < 1e-14 {
            // Derivative too small -> bail out
            return None;
        }

        let delta = f / f_prime;
        e_anomaly -= delta;

        if delta.abs() < TOLERANCE {
            return Some(e_anomaly); // Converged successfully
        }
    }
    None // Didn't converge in time
}

/// Uses the bisection method to solve Kepler's equation, guaranteed to converge,
/// but typically slower than Newton-Raphson.
///
/// Assumes there's a bracket [lower, upper] where the sign of f differs.
fn bisection_kepler(m: Radians, e: f64, mut lower: Radians, mut upper: Radians) -> Radians {
    let mut f_lower = kepler_equation_residual(lower, e, m);
    let mut f_upper = kepler_equation_residual(upper, e, m);

    // Expand the bracket if needed
    while f_lower.signum() == f_upper.signum() {
        upper += Radians::new(PI);
        f_upper = kepler_equation_residual(upper, e, m);
    }

    // Perform bisection within the bracket
    for _ in 0..MAX_BISECTION_ITERS {
        let mid = (lower + upper) * 0.5;
        let f_mid = kepler_equation_residual(mid, e, m);

        if f_mid.abs() < TOLERANCE {
            return mid;
        }

        // Narrow the bracket
        if f_lower.signum() == f_mid.signum() {
            lower = mid;
            f_lower = f_mid; // update lower residual
        } else {
            upper = mid;
        }
    }
    // Return midpoint if we never reach the tolerance
    (lower + upper) * 0.5
}

/// Solves Kepler's Equation (E - e*sin(E) = M) for the eccentric anomaly E.
/// Kepler's Equation relates the mean anomaly (MM) to the eccentric anomaly (EE) for an
/// orbiting body with a given orbital eccentricity (ee). Solving for EE given MM and ee
/// is essential in celestial mechanics for determining the position of a body in its orbit.
///
/// # Arguments
/// - `m`: Mean anomaly in radians.
/// - `e`: Orbital eccentricity (0 <= e < 1 for typical elliptical orbits).
///
/// # Returns
/// - `E`: The eccentric anomaly in radians, guaranteed to converge.
pub fn solve_keplers_equation(m: Radians, e: f64) -> Radians {
    // Pick an initial guess:
    // For moderate e, start with E = M; for higher e, E = π.
    let initial_guess = if e < 0.8 { m } else { Radians::new(PI) };

    // 1) Try Newton-Raphson
    if let Some(e_newton) = newton_kepler(m, e, initial_guess) {
        e_newton
    } else {
        // 2) Fallback to bisection
        bisection_kepler(m, e, Radians::new(-2.0 * PI), Radians::new(4.0 * PI))
    }
}

fn orbital_period_days(a: AstronomicalUnits) -> Days {
    // Kepler’s Third Law: T^2 = a^3 for the Sun+tiny planet
    // T in years, a in AstronomicalUnits
    // 1 year ~ 365.256898326 days
    let year = 365.256898326;
    let t_years = a.value().powf(1.5); // T in years
    Days::new(t_years * year) // convert years -> days
}

/// Calculates the heliocentric coordinates of a celestial body at a given Julian date,
/// using the direct textbook rotation formula.
///
/// # Arguments
/// - `elements`: Orbital elements of the celestial body.
/// - `julian_date`: Julian date for the desired position.
///
/// # Returns
/// - `(x, y, z)`: Heliocentric coordinates in AstronomicalUnits.
pub fn calculate_orbit_position(
    elements: &Orbit,
    julian_date: JulianDate,
) -> Ecliptic<AstronomicalUnit> {
    // 1) Mean motion (n).
    let period = orbital_period_days(elements.semi_major_axis);
    type RadiansPerDay = qtty::frequency::Frequency<Radian, Day>;
    let n: RadiansPerDay = Radians::TAU / period;

    // 2) Days since epoch
    let dt: Days = julian_date - elements.epoch;

    // 3) Mean Anomaly (M) in radians
    let m0_rad = elements.mean_anomaly_at_epoch.to::<Radian>();
    let m_rad = (m0_rad + (n * dt).to::<Radian>()) % std::f64::consts::TAU;

    // 4) Solve Kepler's Equation (E) for the eccentric anomaly
    let e_anomaly = solve_keplers_equation(m_rad, elements.eccentricity);

    // 5) True Anomaly (ν)
    let e = elements.eccentricity;
    let true_anomaly = 2.0 * ((1.0 + e).sqrt() * (e_anomaly * 0.5).tan() / (1.0 - e).sqrt()).atan();

    // 6) Heliocentric distance (z)
    let z = elements.semi_major_axis * (1.0 - e * e_anomaly.cos());

    // 7) Compute standard rotation angular
    let i_rad = elements.inclination.to::<Radian>();
    let omega_rad = elements.longitude_of_ascending_node.to::<Radian>();
    let w_rad = elements.argument_of_perihelion.to::<Radian>();

    // 8) Textbook formula: position in ecliptic coordinates (X, Y, Z)
    //
    //   X = z * [cosΩ * cos(ω+ν) − sinΩ * sin(ω+ν) cos i]
    //   Y = z * [sinΩ * cos(ω+ν) + cosΩ * sin(ω+ν) cos i]
    //   Z = z * [sin(ω+ν) * sin i]
    //
    let (sin_i, cos_i) = (i_rad.sin(), i_rad.cos());
    let (sin_omega, cos_omega) = omega_rad.sin_cos();
    let (sin_w_nu, cos_w_nu) = (w_rad + Radians::new(true_anomaly)).sin_cos();

    Ecliptic::new(
        z * (cos_omega * cos_w_nu - sin_omega * sin_w_nu * cos_i),
        z * (sin_omega * cos_w_nu + cos_omega * sin_w_nu * cos_i),
        z * (sin_w_nu * sin_i),
    )
}

impl Orbit {
    /// Calculates heliocentric coordinates of the orbiting body at a given Julian date.
    pub fn kepler_position(&self, jd: JulianDate) -> Ecliptic<AstronomicalUnit> {
        calculate_orbit_position(self, jd)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time::JulianDate;
    use crate::macros::assert_cartesian_eq;
    use qtty::{Days, Degrees};

    /// Helper function to compare two floating-point numbers with a tolerance.
    fn approx_eq(a: f64, y: f64, tol: f64) -> bool {
        (a - y).abs() < tol
    }

    /// Test Kepler's Equation Solver with known values.
    #[test]
    fn test_solve_keplers_equation() {
        // Test Case 1: e = 0.0, M = 0 radians
        let e = 0.0;
        let m = Radians::new(0.0);
        let computed_e = solve_keplers_equation(m, e);
        let expected_e = 0.0;
        assert!(approx_eq(computed_e.value(), expected_e, 1e-10));

        // Test Case 2: e = 0.0, M = PI/2 radians
        let m = Radians::new(PI / 2.0);
        let computed_e = solve_keplers_equation(m, e);
        let expected_e = PI / 2.0;
        assert!(approx_eq(computed_e.value(), expected_e, 1e-10));

        // Test Case 3: e = 0.1, M = PI/2 radians
        let e = 0.1;
        let m = Radians::new(PI / 2.0);
        let computed_e = solve_keplers_equation(m, e);
        let expected_e = 1.670302; // Corrected expected value from previous step
        assert!(approx_eq(computed_e.value(), expected_e, 1e-6));
    }

    /// Test circular orbit where eccentricity is zero.
    #[test]
    fn test_circular_orbit() {
        let orbit = Orbit::new(
            1.0 * AU,          // a (AstronomicalUnits)
            0.0,               // e
            Degrees::new(0.0), // i
            Degrees::new(0.0), // Ω
            Degrees::new(0.0), // ω
            Degrees::new(0.0), // M0
            JulianDate::J2000, // epoch (Julian Date)
        );

        // At epoch, mean anomaly M0 = 0 degrees, so true anomaly should be 0
        let coord = calculate_orbit_position(&orbit, JulianDate::J2000);

        assert_cartesian_eq!(coord, Ecliptic::new(1.0, 0.0, 0.0), 1e-10);

        // 90 degrees after epoch
        let jd = JulianDate::J2000 + Days::new(90.0 / 0.9856076686); // Roughly 90 degrees / n days
        let coord = calculate_orbit_position(&orbit, jd);
        // Expect y to be approximately 1 AstronomicalUnits
        assert_cartesian_eq!(coord, Ecliptic::new(0.0, 1.0, 0.0), 1e-4);
    }

    /// Test heliocentric coordinates for zero inclination.
    #[test]
    fn test_zero_inclination() {
        let elements = Orbit::new(
            2.0 * AU,          // a (AstronomicalUnits)
            0.1,               // e
            Degrees::new(0.0), // i
            Degrees::new(0.0), // Ω
            Degrees::new(0.0), // ω
            Degrees::new(0.0), // M0
            JulianDate::J2000, // epoch (Julian Date)
        );

        // At epoch, mean anomaly M0 = 0, so true anomaly = 0
        let coord = calculate_orbit_position(&elements, JulianDate::J2000);
        assert_cartesian_eq!(coord, Ecliptic::new(1.8, 0.0, 0.0), 1e-10);
    }

    /// Test heliocentric coordinates after a specific number of days.
    #[test]
    fn test_after_days() {
        let elements = Orbit::new(
            1.0 * AU,          // a (AstronomicalUnits)
            0.0167,            // e
            Degrees::new(0.0), // i
            Degrees::new(0.0), // Ω
            Degrees::new(0.0), // ω
            Degrees::new(0.0), // M0
            JulianDate::J2000, // epoch (Julian Date)
        );

        // Assume circular orbit for simplicity in this test
        let coord = calculate_orbit_position(&elements, JulianDate::J2000 + Days::new(100.0));

        // Mean motion n = 0.9856076686 / a^(3/2) = 0.9856076686 degrees/day
        let n = 0.9856076686; // degrees/day
        let m_deg = Degrees::new(0.0 + n * 100.0) % 360.0; // 98.56076686 degrees
        let m_rad = m_deg.to::<Radian>();

        // For e=0.0167, we can compute expected E and true anomaly
        // However, for simplicity, we'll check that the distance is roughly constant
        let computed_r =
            (coord.x().value().powi(2) + coord.y().value().powi(2) + coord.z().value().powi(2))
                .sqrt();
        let expected_r = 1.0 * (1.0 - 0.0167 * m_rad.cos()); // Approximation

        assert!(approx_eq(computed_r, expected_r, 1e-3)); // Allow some tolerance
    }

    #[test]
    fn test_planets_at_epoch() {
        use crate::bodies::*;

        struct PlanetTest<'a> {
            name: &'a str,
            planet: Planet,
            expected: (f64, f64, f64),
            tol: f64,
        }

        let planets = [
            PlanetTest {
                name: "Mercury",
                planet: MERCURY,
                expected: (-0.1302524, -0.4472397, -0.0245799),
                tol: 1e-3,
            },
            PlanetTest {
                name: "Venus",
                planet: VENUS,
                expected: (-0.7183069, -0.0325362, 0.0410162),
                tol: 1e-3,
            },
            PlanetTest {
                name: "Mars",
                planet: MARS,
                expected: (1.3907092, -0.0135668, -0.0344708),
                tol: 1e-3,
            },
            PlanetTest {
                name: "Jupiter",
                planet: JUPITER,
                expected: (4.0012926, 2.9384137, -0.1017857),
                tol: 1e-2,
            },
            PlanetTest {
                name: "Saturn",
                planet: SATURN,
                expected: (6.4066181, 6.5698012, -0.3690818),
                tol: 5e-2,
            },
            PlanetTest {
                name: "Uranus",
                planet: URANUS,
                expected: (14.4315748, -13.7346340, -0.2381392),
                tol: 1e-2,
            },
            PlanetTest {
                name: "Neptune",
                planet: NEPTUNE,
                expected: (16.8116521, -24.9919786, 0.1272357),
                tol: 1e-2,
            },
            PlanetTest {
                name: "Pluto",
                planet: PLUTO,
                expected: (-9.8758748, -27.9585114, 5.8505715),
                tol: 1e-2,
            },
        ];

        for planet in &planets {
            let coord = calculate_orbit_position(&planet.planet.orbit, JulianDate::J2000);
            let expected = Ecliptic::new(planet.expected.0, planet.expected.1, planet.expected.2);
            assert_cartesian_eq!(coord, expected, planet.tol, "{} at J2000", planet.name);
        }
    }

    // Additional tests for untested functions and edge cases

    #[test]
    fn test_kepler_equation_residual() {
        // Test the residual function directly
        let e_anomaly = Radians::new(PI / 2.0);
        let e = 0.1;
        let m = Radians::new(PI / 2.0);

        // For E = π/2, e = 0.1, M = π/2
        // Residual = E - e*sin(E) - M = π/2 - 0.1*1 - π/2 = -0.1
        let residual = kepler_equation_residual(e_anomaly, e, m);
        assert!(approx_eq(residual.value(), -0.1, 1e-10));
    }

    #[test]
    fn test_kepler_equation_derivative() {
        // Test the derivative function directly
        let e_anomaly = Radians::new(0.0);
        let e = 0.1;

        // For E = 0, e = 0.1
        // Derivative = 1 - e*cos(E) = 1 - 0.1*1 = 0.9
        let derivative = kepler_equation_derivative(e_anomaly, e);
        assert!(approx_eq(derivative, 0.9, 1e-10));
    }

    #[test]
    fn test_solve_keplers_equation_edge_cases() {
        // Test with very small eccentricity
        let e = 1e-10;
        let m = Radians::new(PI / 4.0);
        let computed_e = solve_keplers_equation(m, e);
        assert!(approx_eq(computed_e.value(), m.value(), 1e-10));

        // Test with high eccentricity (near 1.0)
        let e = 0.99;
        let m = Radians::new(PI / 2.0);
        let computed_e = solve_keplers_equation(m, e);
        // Should still converge and give a reasonable result
        assert!(computed_e.value().is_finite());
        assert!(computed_e.value() >= 0.0);
        assert!(computed_e.value() <= 2.0 * PI);
    }

    #[test]
    fn test_solve_keplers_equation_full_range() {
        // Test over full range of mean anomalies
        let e = 0.1;
        for i in 0..8 {
            let m = Radians::new(i as f64 * PI / 4.0);
            let computed_e = solve_keplers_equation(m, e);

            // Check that result is finite and in reasonable range
            assert!(computed_e.value().is_finite());
            assert!(computed_e.value() >= 0.0);
            assert!(computed_e.value() <= 2.0 * PI);

            // Check that it satisfies Kepler's equation approximately
            let residual = kepler_equation_residual(computed_e, e, m);
            assert!(residual.value().abs() < 1e-10);
        }
    }

    #[test]
    fn test_orbital_period_days() {
        // Test orbital period calculation
        let a = AstronomicalUnits::new(1.0); // 1 AU
        let period = orbital_period_days(a);

        // For 1 AU, period should be approximately 365.256898326 days
        assert!(approx_eq(period.value(), 365.256898326, 1e-6));

        // Test with different semi-major axes
        let a = AstronomicalUnits::new(4.0); // 4 AU
        let period = orbital_period_days(a);

        // For 4 AU, period should be approximately 8 years = 2922.055186608 days
        assert!(approx_eq(period.value(), 2922.055186608, 1e-6));
    }

    #[test]
    fn test_calculate_orbit_position_edge_cases() {
        // Test with zero semi-major axis (should handle gracefully)
        let orbit = Orbit::new(
            AstronomicalUnits::new(0.0), // a = 0 AU
            0.0,                         // e
            Degrees::new(0.0),           // i
            Degrees::new(0.0),           // Ω
            Degrees::new(0.0),           // ω
            Degrees::new(0.0),           // M0
            JulianDate::J2000,           // epoch
        );

        let coord = calculate_orbit_position(&orbit, JulianDate::J2000);
        assert!(coord.x().value().is_finite());
        assert!(coord.y().value().is_finite());
        assert!(coord.z().value().is_finite());
    }

    #[test]
    fn test_calculate_orbit_position_different_epochs() {
        let orbit = Orbit::new(
            1.0 * AU,          // a (AU)
            0.1,               // e
            Degrees::new(0.0), // i
            Degrees::new(0.0), // Ω
            Degrees::new(0.0), // ω
            Degrees::new(0.0), // M0
            JulianDate::J2000, // epoch
        );

        // Test at different Julian dates
        let dates = [
            JulianDate::J2000,
            JulianDate::J2000 + Days::new(365.25),
            JulianDate::J2000 + Days::new(730.5),
        ];

        for &date in &dates {
            let coord = calculate_orbit_position(&orbit, date);
            assert!(coord.x().value().is_finite());
            assert!(coord.y().value().is_finite());
            assert!(coord.z().value().is_finite());

            // Distance should be positive
            let distance =
                (coord.x().value().powi(2) + coord.y().value().powi(2) + coord.z().value().powi(2))
                    .sqrt();
            assert!(distance > 0.0);
        }
    }

    #[test]
    fn test_orbit_kepler_position_method() {
        let orbit = Orbit::new(
            1.0 * AU,          // a (AU)
            0.0,               // e (circular)
            Degrees::new(0.0), // i
            Degrees::new(0.0), // Ω
            Degrees::new(0.0), // ω
            Degrees::new(0.0), // M0
            JulianDate::J2000, // epoch
        );

        // Test the convenience method
        let coord1 = orbit.kepler_position(JulianDate::J2000);
        let coord2 = calculate_orbit_position(&orbit, JulianDate::J2000);

        assert!(approx_eq(coord1.x().value(), coord2.x().value(), 1e-10));
        assert!(approx_eq(coord1.y().value(), coord2.y().value(), 1e-10));
        assert!(approx_eq(coord1.z().value(), coord2.z().value(), 1e-10));
    }

    #[test]
    fn test_high_inclination_orbit() {
        let orbit = Orbit::new(
            1.0 * AU,           // a (AU)
            0.1,                // e
            Degrees::new(89.0), // i (high inclination)
            Degrees::new(0.0),  // Ω
            Degrees::new(0.0),  // ω
            Degrees::new(90.0), // M0 (90 degrees to get non-zero z)
            JulianDate::J2000,  // epoch
        );

        let coord = calculate_orbit_position(&orbit, JulianDate::J2000);

        // With high inclination, z component should be significant
        // For 89° inclination, z should be close to the radial distance
        assert!(coord.z().value().abs() > 0.01);

        // Distance should still be reasonable
        let distance =
            (coord.x().value().powi(2) + coord.y().value().powi(2) + coord.z().value().powi(2))
                .sqrt();
        assert!(distance > 0.0);
        assert!(distance < 2.0); // Should be less than 2 AU for this orbit
    }

    #[test]
    fn test_retrograde_orbit() {
        let orbit = Orbit::new(
            1.0 * AU,            // a (AU)
            0.1,                 // e
            Degrees::new(180.0), // i (retrograde)
            Degrees::new(0.0),   // Ω
            Degrees::new(0.0),   // ω
            Degrees::new(0.0),   // M0
            JulianDate::J2000,   // epoch
        );

        let coord = calculate_orbit_position(&orbit, JulianDate::J2000);

        // Should still give valid coordinates
        assert!(coord.x().value().is_finite());
        assert!(coord.y().value().is_finite());
        assert!(coord.z().value().is_finite());
    }
}
