// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # File-format parsers and writers
//!
//! Low-level parsers and writers for standard astronomical and geodetic data
//! formats. These modules operate without knowledge of the dataset catalog or
//! acquisition machinery.
//!
//! | Module | Coverage | `std` required |
//! |--------|----------|----------------|
//! | [`adsb`] | ADS-B / Mode S Extended Squitter frames | no |
//! | [`ccsds`] | OEM / OPM / TDM text messages | yes (I/O traits) |
//! | [`iers`] | Earth-orientation parameter products | yes (I/O traits) |
//! | [`igs`] | SP3 / ANTEX / SINEX / ORBEX products | yes (I/O traits) |
//! | [`ilrs`] | CRD / CPF laser-ranging products | yes (filesystem helpers) |
//! | [`rinex`] | RINEX observation / navigation formats | yes (I/O traits) |
//! | [`sck`] | Siderust Chebyshev Kernel v1 (archive binary) | no |
//! | [`spice`] | SPICE text/binary kernel parsing (byte-oriented core) | path helpers need `std` |
//! | [`tle`] | NORAD TLE / 3LE / CCSDS OMM (KVN, XML, JSON) | no |
//! | [`vlbi`] | vgosDB VLBI datasets | yes (I/O traits) |
//!
//! For the dataset catalog (what datasets exist and how to acquire them) see the
//! [`siderust_archive`] crate. The `formats` modules sit *below* the catalog:
//! they are called by the runtime back-end after a file has been located on disk.

pub mod adsb;
#[cfg(feature = "std")]
pub mod ccsds;
pub mod error;
#[cfg(feature = "std")]
pub mod iers;
#[cfg(feature = "std")]
pub mod igs;
#[cfg(feature = "std")]
pub mod ilrs;
#[cfg(feature = "std")]
pub mod rinex;
pub mod sck;
pub mod spice;
pub mod tle;
#[cfg(feature = "std")]
pub mod vlbi;

pub use error::{FileLocation, FormatError, ParseMode};
