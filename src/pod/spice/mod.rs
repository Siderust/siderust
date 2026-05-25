// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! NAIF SPICE kernel reader and [`SpiceEphemerisProvider`] adapter.
//!
//! This module re-exports the pure kernel-parsing layer from
//! [`siderust::formats::spice`] and adds the POD-specific
//! [`SpiceEphemerisProvider`] that implements
//! [`siderust::pod::providers::EphemerisProvider`].
//!
//! # Scope
//!
//! The pure DAF/SPK parsing, SPK segment evaluation (Type 2/3/9/13), body
//! chain resolution, and NAIF body-ID lookup all live in
//! `siderust::formats::spice`. This module only owns the thin POD
//! adapter:
//!
//! * [`SpiceEphemerisProvider`] — typed adapter implementing
//!   [`siderust::pod::providers::EphemerisProvider`]. Center selection
//!   is **explicit** at construction time and is propagated as a typed
//!   `affn` reference center on the returned [`SpiceState`].
//!
//! # Time semantics
//!
//! All epochs are **TDB seconds past J2000**, matching the convention of
//! NAIF SPK files and [`tempoch::J2000s`]. The provider trait method
//! exposes the same numeric convention; the typed `tempoch::EncodedTime`
//! constructor is available via [`SpiceEphemerisProvider::state_at`].
//!
//! # Validation gate
//!
//! Position recovery is byte-identical to NAIF's CSPICE `spkez_c` for
//! the planet body chains in `de440.bsp` / `de441.bsp`. The `de440`
//! feature unlocks an integration test (`tests/de440_validation.rs`)
//! that loads a kernel from disk (`SIDERUST_SPICE_DE_PATH`) and
//! compares against committed reference states.
//!
//! # Examples
//!
//! ```rust
//! use siderust::formats::spice::{daf::{Daf, Summary}, segment_for_summary, SpiceError};
//!
//! let daf = Daf { nd: 2, ni: 6, summaries: vec![] };
//! let summary = Summary {
//!     start_et: 0.0, end_et: 1.0, target_id: 1, center_id: 0,
//!     frame_id: 1, data_type: 21, start_word: 1, end_word: 1,
//! };
//! let err = segment_for_summary(&0.0_f64.to_le_bytes(), &daf, &summary).unwrap_err();
//! assert!(matches!(err, SpiceError::UnsupportedDataType { data_type: 21 }));
//! # Ok::<_, SpiceError>(())
//! ```

#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod provider;

pub use provider::{SpiceEphemerisProvider, SpiceState};

// Re-export the full upstream spice parsing surface.
pub use crate::formats::spice::{
    daf, kernel, naif,
    naif::{naif_id_for_name, well_known},
    segment,
    segment::{segment_for_summary, ChebSegment, SpkSegment},
    spk, LoadedSegment, SpiceError, SpkKernel,
};
