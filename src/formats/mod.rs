// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # File-format parsers and writers
//!
//! Low-level parsers and writers for standard astronomical and geodetic data
//! formats. These modules operate without knowledge of the dataset catalog or
//! acquisition machinery.
//!
//! | Module | Coverage |
//! |--------|----------|
//! | [`adsb`] | ADS-B / Mode S Extended Squitter frames |
//! | [`ccsds`] | OEM / OPM / TDM text messages |
//! | [`iers`] | Earth-orientation parameter products |
//! | [`igs`] | SP3 / ANTEX / SINEX / ORBEX products |
//! | [`ilrs`] | CRD / CPF laser-ranging products |
//! | [`rinex`] | RINEX observation / navigation formats |
//! | [`spice`] | SPICE text and binary kernel parsing (SPK, CK, FK, LSK, PCK, SCLK, IK) |
//! | [`tle`] | NORAD TLE / 3LE / CCSDS OMM (KVN, XML, JSON) |
//! | [`vlbi`] | vgosDB VLBI datasets |
//!
//! For the dataset catalog (what datasets exist and how to acquire them) see
//! [`crate::data`]. The `formats` modules sit *below* the catalog: they
//! are called by the runtime back-end after a file has been located on disk.

pub mod adsb;
pub mod ccsds;
pub mod error;
pub mod iers;
pub mod igs;
pub mod ilrs;
pub mod rinex;
pub mod spice;
pub mod tle;
pub mod vlbi;

pub use error::{FileLocation, FormatError, ParseMode};
