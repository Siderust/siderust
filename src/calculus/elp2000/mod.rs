// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # ELP2000-82B Lunar Theory
//!
//! ## Scientific scope
//!
//! ELP2000-82B is a semi-analytical theory of the Moon's orbital motion
//! developed by Chapront-Touzé and Chapront at the Bureau des Longitudes,
//! Paris.  It expresses the Moon's geocentric ecliptic longitude, latitude,
//! and distance as Poisson series in time `T` (Julian millennia from J2000.0):
//!
//! - **36 sub-series** divided into three families:
//!   Main Problem, Earth Perturbations, and Planetary Perturbations.
//! - Accuracy: ±10 arcseconds over ±10 centuries from J2000.0 (1000–3000 AD).
//! - Validity: formally defined over ±20 centuries, with rapidly degrading
//!   accuracy outside ±10 centuries.
//!
//! ## Technical scope
//!
//! - [`elp_structs`] — coefficient record types (`MainProblem`, `EarthPert`,
//!   `PlanetPert`).
//! - [`elp_series`] — all 36 coefficient tables and the evaluation engine;
//!   public entry point: `Moon::get_geo_position<U>(jd: JulianDate)` returning
//!   `Position<Geocentric, EclipticMeanJ2000, U>`.
//!
//! ## References
//!
//! - Chapront-Touzé, M., & Chapront, J. (1983). "The lunar ephemeris ELP 2000".
//!   *Astronomy and Astrophysics* 124, 50–62.
//! - Chapront-Touzé, M., & Chapront, J. (1988). "ELP 2000-82B: A semi-analytical
//!   lunar ephemeris adequate for historical times".
//!   *Astronomy and Astrophysics* 190, 342–352.

mod elp_series;
mod elp_structs;
