//! Worked example: heliocentric Earth → Mars Lambert transfer.
//!
//! Approximates a single-revolution prograde transfer from a notional
//! Earth position to a notional Mars position 60° ahead, using a
//! Hohmann-like time of flight. The example exercises the typed
//! [`keplerian::lambert::lambert`] entry-point end-to-end and prints
//! departure / arrival velocities together with the implied Δv at each
//! end (assuming the planets travel on circular heliocentric orbits at
//! the chosen radii).
//!
//! Run with:
//! ```bash
//! cargo run --example 03_lambert_earth_to_mars
//! ```

#![allow(clippy::print_stdout, clippy::print_stderr)]

use affn::cartesian::Position;
use affn::frames::ICRS;
use keplerian::lambert::{lambert, LambertBranch};
use qtty::dynamics::GravitationalParameter;
use qtty::length::Kilometer;
use qtty::Second;

const AU_KM: f64 = 1.495_978_707e8;
const MU_SUN: f64 = 1.327_124_400_18e11; // km^3 / s^2
const SECONDS_PER_DAY: f64 = 86_400.0;

fn main() {
    // Heliocentric ICRS positions (toy J2000 placement; see crate doc
    // for usage with real ephemerides).
    let r_earth_km = AU_KM;
    let r_mars_km = 1.524 * AU_KM;
    let phase_lead = 60.0_f64.to_radians(); // Mars 60° ahead of Earth

    let r1 = Position::<(), ICRS, Kilometer>::new(r_earth_km, 0.0, 0.0);
    let r2 = Position::<(), ICRS, Kilometer>::new(
        r_mars_km * phase_lead.cos(),
        r_mars_km * phase_lead.sin(),
        0.0,
    );

    let tof_days = 258.0;
    let tof = Second::new(tof_days * SECONDS_PER_DAY);
    let mu = GravitationalParameter::new(MU_SUN);

    let solution = lambert(r1, r2, tof, mu, LambertBranch::Prograde)
        .expect("Lambert solver failed on a Hohmann-like Earth → Mars geometry");

    let v1 = &solution.v1;
    let v2 = &solution.v2;

    let v1x = v1.x().value();
    let v1y = v1.y().value();
    let v1z = v1.z().value();
    let v2x = v2.x().value();
    let v2y = v2.y().value();
    let v2z = v2.z().value();
    let v1_mag = (v1x * v1x + v1y * v1y + v1z * v1z).sqrt();
    let v2_mag = (v2x * v2x + v2y * v2y + v2z * v2z).sqrt();

    // Circular planetary speeds for Δv accounting.
    let v_earth_circ = (MU_SUN / r_earth_km).sqrt();
    let v_mars_circ = (MU_SUN / r_mars_km).sqrt();

    let dv_dep = (v1x - 0.0).hypot(v1y - v_earth_circ).hypot(v1z);
    let dv_arr = (v2x - (-v_mars_circ * phase_lead.sin()))
        .hypot(v2y - v_mars_circ * phase_lead.cos())
        .hypot(v2z);

    println!("=== Lambert Earth → Mars (toy heliocentric placement) ===");
    println!(
        "  r1 = ({:>12.3}, {:>12.3}, {:>12.3}) km",
        r_earth_km, 0.0, 0.0
    );
    println!(
        "  r2 = ({:>12.3}, {:>12.3}, {:>12.3}) km",
        r_mars_km * phase_lead.cos(),
        r_mars_km * phase_lead.sin(),
        0.0
    );
    println!("  ToF = {:.1} d ({:.0} s)", tof_days, tof.value());
    println!();
    println!("Departure velocity (heliocentric, ICRS):");
    println!(
        "  v1 = ({:>9.4}, {:>9.4}, {:>9.4}) km/s   |v1| = {:.4}",
        v1x, v1y, v1z, v1_mag
    );
    println!("Arrival velocity (heliocentric, ICRS):");
    println!(
        "  v2 = ({:>9.4}, {:>9.4}, {:>9.4}) km/s   |v2| = {:.4}",
        v2x, v2y, v2z, v2_mag
    );
    println!();
    println!("Δv accounting (vs. circular planetary motion):");
    println!("  Δv at Earth  = {:>7.4} km/s", dv_dep);
    println!("  Δv at Mars   = {:>7.4} km/s", dv_arr);
    println!("  Δv total     = {:>7.4} km/s", dv_dep + dv_arr);
    println!();
    println!("Diagnostics:");
    println!("  iterations        = {}", solution.diagnostics.iterations);
    println!(
        "  Householder resid = {:.3e}",
        solution.diagnostics.residual
    );
    println!("  revolutions       = {}", solution.diagnostics.revolutions);
}
