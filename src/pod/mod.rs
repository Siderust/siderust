// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Precise Orbit Determination (POD) toolkit
//!
//! ## Scientific scope
//!
//! `siderust::pod` is a library of typed building blocks for Precise Orbit
//! Determination (POD) of artificial satellites — GNSS and LEO — built on
//! the rest of `siderust` and the [`qtty`], [`tempoch`], and [`affn`]
//! baseline crates. It was previously distributed as the standalone
//! `siderust-pod` crate; it now lives behind the default-on `pod` Cargo
//! feature inside `siderust`.
//!
//! The crate provides the geophysical primitives needed by any POD pipeline:
//!
//! - typed problem definitions and run metadata ([`problem`], [`run`]),
//! - high-fidelity force-model composition, propagation, and process-noise building blocks on top of `siderust` ([`force`], [`propagation`], [`process`]),
//! - observation models for GNSS pseudorange and carrier-phase and SLR
//!   normal points ([`observation`]),
//! - batch and sequential estimators ([`estimation`]),
//! - I/O compatibility aliases bridging `siderust::formats` ([`io`]),
//! - quality-control products and reproducibility artifacts ([`qc`], [`product`]).
//!
//! Mission-specific code (e.g. LISA orbit readers, inter-satellite range
//! models) is kept in `siderust/examples/` rather than the library source.
//!
//! The full POD 1.0 spec — Batch Weighted Least Squares with STM and
//! variational equations, IGS golden-test interoperability, SLR Marini-
//! Murray, sequential filters, multi-arc combination, and external
//! validation against IGS/ILRS data — is tracked in the long-term
//! roadmap (`reqs-and-plan.md`).
//!
//! ## Technical scope
//!
//! Every public API uses typed quantities (`qtty::*`), typed instants
//! (`tempoch::Time<S, F>`), and typed coordinates (`affn::*` /
//! `siderust::coordinates::*`); raw `f64` is reserved for internal math
//! kernels. `forbid(unsafe_code)` is enforced at the workspace level.
//!
//! Workspace lints enforce a zero-tolerance posture: `missing_docs` is
//! denied, `todo!`/`unimplemented!`/`dbg!` are denied, and
//! `broken_intra_doc_links` is denied.
//!
//! ## Modules
//!
//! | Module           | Contents                                                        |
//! |------------------|-----------------------------------------------------------------|
//! | [`spice`]        | DAF/SPK kernel reader and SPICE ephemeris provider              |
//! | [`io`]           | I/O compatibility aliases for `siderust::formats`               |
//! | [`problem`]      | POD problem primitives: arcs, parameters, covariance            |
//! | [`run`]          | Run-level provenance and dataset references                     |
//! | [`providers`]    | Ephemeris provider traits bridging POD code to services         |
//! | [`force`]        | Force-model configuration, registries, empirical accelerations  |
//! | [`propagation`]  | Integrators, adapters, variational propagation, STM helpers     |
//! | [`process`]      | Process-noise models and Q-matrix construction                  |
//! | [`observation`]  | Typed measurement models: GNSS, SLR                             |
//! | [`estimation`]   | Batch (WLS, Gauss-Newton) and sequential estimator scaffolding  |
//! | [`qc`]           | Quality-control and validation products                         |
//! | [`product`]      | SP3, OEM, residual CSV/Parquet, manifest writers                |
//!
//! ## References
//!
//! - Tapley, B. D., Schutz, B. E., & Born, G. H. (2004). *Statistical
//!   Orbit Determination*. Elsevier Academic Press.
//! - Montenbruck, O., & Gill, E. (2000). *Satellite Orbits — Models,
//!   Methods, and Applications*. Springer.
//! - International GNSS Service (2024). *IGS Products Guide*.
//! - Petit, G., & Luzum, B. (Eds.) (2010). *IERS Conventions (2010)*.
//!   IERS Technical Note 36.

#![forbid(unsafe_code)]

pub mod estimation;
pub mod force;
pub mod io;
pub mod observation;
pub mod problem;
pub mod process;
pub mod product;
pub mod propagation;
pub mod providers;
pub mod qc;
pub mod run;
pub mod spice;
