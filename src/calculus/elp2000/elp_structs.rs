// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # ELP2000-82B Coefficient Structures
//!
//! ## Scientific scope
//!
//! The ELP2000-82B theory (Chapront-Touzé & Chapront 1988) partitions its
//! Poisson-series terms into three families depending on which gravitational
//! perturbations they model:
//!
//! - **Main Problem** — the Earth–Moon–Sun three-body problem without planetary
//!   perturbations.  Arguments are the four Delaunay angles D, M, Mʹ, F.
//! - **Earth Perturbations** — figure, tidal, and relativistic corrections
//!   involving the same four Delaunay angles plus an explicit power of T.
//! - **Planetary Perturbations** — gravitational influence of the eight planets
//!   (Mercury through Neptune), using eleven combined planetary/lunar arguments.
//!
//! Each family is represented by a distinct coefficient record type defined here.
//!
//! ## Technical scope
//!
//! - [`MainProblem`] — one sine term of the Main Problem series; amplitude `a`
//!   in 10⁻⁴ arc-seconds, integer Delaunay multipliers, and phase polynomial
//!   coefficients b₀…b₅.
//! - [`EarthPert`] — one Earth-perturbation term; amplitude `a` in 10⁻⁴
//!   arc-seconds, T-exponent `iz` (dimensionless), phase rate `o`
//!   (radians / Julian-millennia), and constant phase `p` (radians).
//! - [`PlanetPert`] — one planetary-perturbation term; amplitude `theta` in
//!   10⁻⁴ arc-seconds, eleven integer planetary multipliers, phase rate `o`
//!   (radians / Julian-millennia), and constant phase `p` (radians).
//!
//! All numeric fields use raw `f64`; the typed conversion to physical units
//! (arc-seconds → radians → km) is performed in [`super::elp_series`].
//!
//! ## References
//!
//! - Chapront-Touzé, M., & Chapront, J. (1988). "ELP 2000-82B: A semi-analytical
//!   lunar ephemeris adequate for historical times".
//!   *Astronomy and Astrophysics* 190, 342–352.

#![allow(unreachable_pub)]

// ────────────────────────────────────────────────────────────────────────────
// Main term of the lunar theory (ELP “Main Problem” series).
// Each record governs one sine term of the form:
//
//   Δλ = a * sin( Σ ilu[i] * Φ_i  +  Σ b[j] * T^j )
//
// where
//   Φ_i  = {D, M, M′, F}  → the four fundamental lunar arguments
//   T     = Julian millennia from J2000 (ΔT = (JD − 2451545.0)/365250)
//   b[j]  = time polynomial coefficients (constant, linear, … 5th-order)
//
// All amplitudes are usually given in 0.0001″ (10⁻⁴ arc-sec).
// ────────────────────────────────────────────────────────────────────────────
pub struct MainProblem {
    /// Integer multipliers of the fundamental arguments
    /// in the order [D, M, M′, F].
    pub ilu: [i64; 4],

    /// Sine amplitude ‘a’ (10⁻⁴ arc-seconds).
    pub a: f64,

    /// Time‐polynomial coefficients b₀…b₅ for the phase argument:
    ///   φ = Σ (ilu[i]·Φ_i)  +  (b₀ + b₁ T + ⋯ + b₅ T⁵)
    pub b: [f64; 6],
}

// ────────────────────────────────────────────────────────────────────────────
// Earth–Moon specific perturbations (ELP Earth Perturbation series).
// A term contributes
//
//   Δλ = a · T^iz · sin( Σ ilu[i] * Φ_i  +  p  +  o·T )
//
//   iz   = real exponent of T (often 0, 1, or −1)
//   o    = secular rate in the phase (radians per Julian millennia)
//   p    = phase offset at T = 0
// ────────────────────────────────────────────────────────────────────────────
#[allow(dead_code)]
pub struct EarthPert {
    /// Exponent ‘iz’ of the time factor T^iz.
    pub iz: f64,

    /// Integer multipliers for [D, M, M′, F].
    pub ilu: [i64; 4],

    /// Linear phase rate ‘o’ (radians / Julian millennia).
    pub o: f64,

    /// Amplitude ‘a’ (10⁻⁴ arc-seconds).
    pub a: f64,

    /// Constant phase offset ‘p’ (radians) at epoch J2000.
    pub p: f64,
}

// ────────────────────────────────────────────────────────────────────────────
// Planetary perturbations acting on the Moon
// (ELP Planetary Perturbation series).
// Each term is:
//
//   Δλ = θ · sin( Σ ipla[k] * Π_k  +  p  +  o·T )
//
//   Π_k  = 11 planetary arguments (L_ME, L_VE, …, L_NE, D, M, M′)
//   θ    = amplitude (arc-seconds)
//   o    = secular rate in phase
//   p    = phase offset
// ────────────────────────────────────────────────────────────────────────────
#[allow(dead_code)]
pub struct PlanetPert {
    /// Integer multipliers for the 11 planetary arguments Π_k
    /// (Mercury → Neptune plus three lunar/solar combinations).
    pub ipla: [i64; 11],

    /// Amplitude ‘θ’ (10⁻⁴ arc-seconds).
    pub theta: f64,

    /// Linear phase rate ‘o’ (radians / Julian millennia).
    pub o: f64,

    /// Constant phase offset ‘p’ (radians) at epoch J2000.
    pub p: f64,
}
