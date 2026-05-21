// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! NORAD Two-Line Element / 3LE / OMM parser & writer.
//!
//! Models the canonical NORAD TLE format and the CCSDS Orbit Mean-elements
//! Message (OMM) in three encodings — KVN, XML, and JSON — sharing one
//! in-memory record:
//!
//! * column-accurate field extraction (TLE is column-positional),
//! * line-checksum validation,
//! * Alpha-5 NORAD-ID extension (catalog numbers ≥ 100 000),
//! * strongly typed angles (`qtty::angular::Degrees`) and mean motion
//!   (`qtty::angular_rate::AngularRate<Turn, Day>`),
//! * `tempoch::Time<UTC>` epochs from the (year, day-of-year, fractional-day)
//!   encoding via `chrono`,
//! * round-trip OMM-KVN / OMM-XML / OMM-JSON parity through [`omm::Omm`],
//! * a programmatic [`TleBuilder`] for synthesising TLEs,
//! * single error taxonomy via [`TleError`].
//!
//! ## Design notes
//!
//! * **No XML library dependency in the core path.** The OMM-XML schema is
//!   small, flat, and stable; the hand-rolled element-extractor in
//!   [`omm::xml`] avoids pulling in a general-purpose XML crate.
//! * **`f64` on the public surface is restricted** to fields where no typed
//!   `qtty`/`tempoch` equivalent improves safety: eccentricity (dimensionless,
//!   range-checked at construction), `BSTAR`, and mean-motion derivatives.
//!
//! ## Out of scope
//!
//! SGP4/SDP4 propagation lives in [`siderust::astro::sgp4`](crate::astro::sgp4).
//!
//! ## References
//!
//! - Kelso, T. S. (2023). *Celestrak TLE format documentation*.
//! - CCSDS 502.0-B-3 — Orbit Data Messages.
//!
//! # Examples
//!
//! ```
//! use siderust::formats::tle::{parse_3le, omm::{kvn, Omm}};
//!
//! let l1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927";
//! let l2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
//! let tle = parse_3le("ISS (ZARYA)", l1, l2).unwrap();
//! let omm = Omm::from_tle(&tle);
//! let kvn_text = kvn::write(&omm).unwrap();
//! let omm2 = kvn::read(&kvn_text).unwrap();
//! assert_eq!(omm2.norad_id, tle.norad_id);
//! ```

#![forbid(unsafe_code)]

mod builder;
mod error;
pub mod omm;
mod parse;
#[allow(clippy::module_inception)]
mod tle;

pub use builder::{format_tle, TleBuilder};
pub use error::TleError;
pub use parse::{compute_tle_checksum, parse_3le, parse_tle, validate_tle_checksum};
pub use tle::{Classification, InternationalDesignator, SatelliteNumber, TLE};

#[cfg(test)]
mod tests;
