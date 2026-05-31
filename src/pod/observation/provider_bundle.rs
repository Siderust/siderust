// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Provider bundle trait for dynamic data look-ups inside observation models.

use crate::astro::dynamics::{Position, Velocity};
use crate::coordinates::frames::GCRS;
use crate::time::JulianDate;
use qtty::Meter;

/// Bundled dynamic data sources consumed by observation model implementations.
///
/// The trait is deliberately minimal and object-safe so that observation models
/// can receive it as `&dyn ProviderBundle`.  Concrete implementations are
/// expected to wrap clock-product files, SP3 ephemeris readers, etc.
///
/// Every method has a reasonable default (zero clock bias, `None` state) so
/// that models that do not need a particular service can use
/// [`NullProviderBundle`] or a partial implementation.
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::provider_bundle::{ProviderBundle, NullProviderBundle};
/// use siderust::time::JulianDate;
///
/// fn check(pb: &dyn ProviderBundle) {
///     assert_eq!(pb.receiver_clock_m().value(), 0.0);
///     assert!(pb.gnss_satellite_state_gcrs("G01", JulianDate::new(2_451_545.0)).is_none());
/// }
/// check(&NullProviderBundle);
/// ```
pub trait ProviderBundle: Send + Sync {
    /// GCRS position (km) and velocity (km/s) of a GNSS transmitter at the
    /// given TT Julian date.  Returns `None` if the satellite is unknown or
    /// the epoch is outside the available data window.
    fn gnss_satellite_state_gcrs(
        &self,
        prn: &str,
        epoch: JulianDate,
    ) -> Option<(Position<GCRS>, Velocity<GCRS>)>;

    /// Satellite clock bias in metres (positive = satellite clock ahead of
    /// GPS time).  Returns `Meter::new(0.0)` if no clock product is available.
    fn gnss_satellite_clock_m(&self, prn: &str, epoch: JulianDate) -> Meter;

    /// Receiver clock bias in metres.  Returns `Meter::new(0.0)` if not
    /// estimated.
    fn receiver_clock_m(&self) -> Meter;

    /// GCRS position (km) of a ground station at the given TT Julian date.
    /// Used by SLR observation models.  Returns `None` if the station is
    /// unknown.
    fn station_gcrs_km(&self, station_id: &str, epoch: JulianDate) -> Option<Position<GCRS>>;
}

// ─── Null implementation ─────────────────────────────────────────────────────

/// A [`ProviderBundle`] that returns zero biases and no satellite/station
/// states.  Useful for unit tests and scenarios where the geometric range
/// is known in advance.
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::provider_bundle::{ProviderBundle, NullProviderBundle};
/// use siderust::time::JulianDate;
///
/// let epoch = JulianDate::new(2_451_545.0);
/// assert_eq!(NullProviderBundle.receiver_clock_m().value(), 0.0);
/// assert!(NullProviderBundle.gnss_satellite_state_gcrs("G05", epoch).is_none());
/// assert_eq!(NullProviderBundle.gnss_satellite_clock_m("G05", epoch).value(), 0.0);
/// assert!(NullProviderBundle.station_gcrs_km("7840", epoch).is_none());
/// ```
pub struct NullProviderBundle;

impl ProviderBundle for NullProviderBundle {
    fn gnss_satellite_state_gcrs(
        &self,
        _prn: &str,
        _epoch: JulianDate,
    ) -> Option<(Position<GCRS>, Velocity<GCRS>)> {
        None
    }

    fn gnss_satellite_clock_m(&self, _prn: &str, _epoch: JulianDate) -> Meter {
        Meter::new(0.0)
    }

    fn receiver_clock_m(&self) -> Meter {
        Meter::new(0.0)
    }

    fn station_gcrs_km(&self, _station_id: &str, _epoch: JulianDate) -> Option<Position<GCRS>> {
        None
    }
}
