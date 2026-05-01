//! Atmospheric models — ozone, Rayleigh/Mie scattering, airglow, extinction.
//!
//! Re-exported from [`siderust`] as `siderust::atmosphere`. Downstream crates
//! that only need the atmospheric layer can depend on this crate directly.

pub mod atmosphere;
pub use atmosphere::*;

