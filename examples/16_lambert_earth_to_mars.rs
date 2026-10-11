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
//! cargo run --example 16_lambert_earth_to_mars --features pod
//! ```

#![allow(clippy::print_stdout, clippy::print_stderr)]

use keplerian::lambert::{lambert, LambertBranch};
use qtty::dynamics::{GravitationalParameter, KmPerSecond};
use siderust::affn::cartesian::{Position, Velocity};
use siderust::coordinates::{centers::Heliocentric, frames::ICRS};
use siderust::qtty::unit::Second;
use siderust::qtty::{Days, Kilometer};

const AU_KM: f64 = 1.495_978_707e8;
const MU_SUN: f64 = 1.327_124_400_18e11; // km³/s²

fn main() {
    // Heliocentric ICRS positions (toy J2000 placement; see crate doc
    // for usage with real ephemerides).
    let r_earth = AU_KM;
    let r_mars = 1.524 * AU_KM;
    let phase_lead = 60.0_f64.to_radians(); // Mars 60° ahead of Earth

    let r1 = Position::<Heliocentric, ICRS, Kilometer>::new(r_earth, 0.0, 0.0);
    let r2 = Position::<Heliocentric, ICRS, Kilometer>::new(
        r_mars * phase_lead.cos(),
        r_mars * phase_lead.sin(),
        0.0,
    );

    let tof_days = Days::new(258.0);
    let tof = tof_days.to::<Second>();

    let mu = GravitationalParameter::new(MU_SUN);
    let solution = lambert(r1, r2, tof, mu, LambertBranch::Prograde)
        .expect("Lambert solver failed on a Hohmann-like Earth → Mars geometry");

    let v1 = solution.v1;
    let v2 = solution.v2;

    // Circular planetary velocities for Δv accounting.
    let v_earth_circ = (MU_SUN / r_earth).sqrt();
    let v_mars_circ = (MU_SUN / r_mars).sqrt();
    let v_earth = Velocity::<ICRS, KmPerSecond>::new(0.0, v_earth_circ, 0.0);
    let v_mars = Velocity::<ICRS, KmPerSecond>::new(
        -v_mars_circ * phase_lead.sin(),
        v_mars_circ * phase_lead.cos(),
        0.0,
    );
    let dv_dep = (v1 - v_earth).magnitude();
    let dv_arr = (v2 - v_mars).magnitude();

    println!("=== Lambert Earth → Mars (toy heliocentric placement) ===");
    println!("  r1 = {r1:.3}");
    println!("  r2 = {r2:.3}");
    println!("  ToF = {tof_days:.1} ({tof:.0})");
    println!();
    println!("Departure velocity (heliocentric, ICRS):");
    println!("  v1 = {v1:.4}   |v1| = {:.4}", v1.magnitude());
    println!("Arrival velocity (heliocentric, ICRS):");
    println!("  v2 = {v2:.4}   |v2| = {:.4}", v2.magnitude());
    println!();
    println!("Δv accounting (vs. circular planetary motion):");
    println!("  Δv at Earth  = {:>7.4}", dv_dep);
    println!("  Δv at Mars   = {:>7.4}", dv_arr);
    println!("  Δv total     = {:>7.4}", dv_dep + dv_arr);
    println!();
    println!("Diagnostics:");
    println!("  iterations        = {}", solution.diagnostics.iterations);
    println!(
        "  Householder resid = {:.3e}",
        solution.diagnostics.residual
    );
    println!("  revolutions       = {}", solution.diagnostics.revolutions);
}
