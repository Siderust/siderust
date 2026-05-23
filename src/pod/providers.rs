//! [`EphemerisProvider`] trait.
//!
//! Concrete implementations should typically wrap
//! `siderust::ephemeris` or sibling crates such as
//! `siderust-spice`. The state representation is left associated so POD
//! code is not pinned to a single `siderust` view.

use std::error::Error;

/// Minimal ephemeris query interface used by POD force models.
///
/// The semantics of `epoch_seconds_tdb` are TDB seconds since J2000.
/// Concrete typed-`Time<TDB>` constructors live on the implementer
/// (e.g. `siderust::pod::spice::SpiceEphemerisProvider::state_at`).
///
/// # Examples
///
/// ```
/// use siderust::pod::providers::EphemerisProvider;
///
/// struct DummyState;
/// struct DummyProvider;
/// impl EphemerisProvider for DummyProvider {
///     type State = DummyState;
///     type Error = std::io::Error;
///     fn state(&self, _id: i32, _t: f64) -> Result<DummyState, Self::Error> {
///         Ok(DummyState)
///     }
/// }
///
/// let p = DummyProvider;
/// let _ = p.state(399, 0.0).unwrap();
/// ```
pub trait EphemerisProvider {
    /// State representation type (framework-specific).
    type State;
    /// Error type for state queries.
    type Error: Error + Send + Sync + 'static;

    /// Return a state for `body_naif_id` at the given `epoch_seconds_tdb`.
    fn state(&self, body_naif_id: i32, epoch_seconds_tdb: f64) -> Result<Self::State, Self::Error>;
}
