//! [`EphemerisProvider`] trait.
//!
//! Concrete implementations should typically wrap
//! `siderust::ephemeris` or sibling crates such as
//! `siderust-spice`. The state representation is left associated so POD
//! code is not pinned to a single `siderust` view.
//!
//! The raw [`EphemerisProvider::state`] method is a wire-format boundary:
//! epochs are passed as TDB seconds since J2000 so low-level SPICE/OEM
//! adapters can forward kernel-native values without extra boxing. Typed POD
//! callers should prefer [`EphemerisProvider::state_at`], which accepts a
//! `tempoch` [`Time<TDB>`].

use std::error::Error;

use qtty::unit::Second;
use qtty::Quantity;
use tempoch::{EncodedTime, J2000s, Time, TDB};

/// Minimal ephemeris query interface used by POD force models.
///
/// [`state`](Self::state) is the raw wire-format method: `epoch_seconds_tdb`
/// means TDB seconds since J2000. Use [`state_at`](Self::state_at) when you
/// already have a typed [`Time<TDB>`].
///
/// # Examples
///
/// ```
/// use siderust::pod::providers::EphemerisProvider;
/// use tempoch::{Time, TDB};
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
/// let epoch = Time::<TDB>::from_raw_j2000_seconds(qtty::Second::new(0.0)).unwrap();
/// let _ = p.state_at(399, epoch).unwrap();
/// ```
pub trait EphemerisProvider {
    /// State representation type (framework-specific).
    type State;
    /// Error type for state queries.
    type Error: Error + Send + Sync + 'static;

    /// Return a state for `body_naif_id` at the given raw TDB-J2000 epoch.
    ///
    /// This is the low-level wire-format entry point used by kernel-native
    /// adapters. Callers that already have a typed [`Time<TDB>`] should prefer
    /// [`state_at`](Self::state_at).
    fn state(&self, body_naif_id: i32, epoch_seconds_tdb: f64) -> Result<Self::State, Self::Error>;

    /// Return a state for `body_naif_id` at a typed TDB epoch.
    fn state_at(&self, body_naif_id: i32, epoch: Time<TDB>) -> Result<Self::State, Self::Error> {
        let encoded: EncodedTime<TDB, J2000s> = epoch.to::<J2000s>();
        let secs: Quantity<Second> = encoded.raw();
        self.state(body_naif_id, secs.value())
    }
}
