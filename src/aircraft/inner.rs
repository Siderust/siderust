// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! [`Aircraft`] identity record, [`AircraftState`] snapshot, and
//! [`AircraftTrack`] dead-reckoning provider.

use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::qtty::unit::{Meter, Second};
use crate::qtty::{Degrees, Meters, Per, Quantity};
use crate::targets::Trackable;
use crate::time::JulianDate;
use std::borrow::Cow;

/// Scalar velocity in metres per second (ground speed, vertical rate).
///
/// Defined as a unit-safe dimensional type rather than a bare `f64`.
pub type MetersPerSecond = Quantity<Per<Meter, Second>>;

// =============================================================================
// Aircraft identity
// =============================================================================

/// Static identity record for an aircraft.
///
/// An `Aircraft` carries the ICAO 24-bit transponder address and an ASCII
/// callsign of up to 8 characters.  It does **not** carry state (position,
/// velocity) — see [`AircraftState`] for the time-varying snapshot.
///
/// # Examples
///
/// ```rust
/// use siderust::aircraft::Aircraft;
///
/// let ac = Aircraft::new(0x4CA2B5, "EIN104");
/// assert_eq!(ac.icao24(), 0x4CA2B5);
/// assert_eq!(ac.callsign(), "EIN104");
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Aircraft<'a> {
    /// ICAO 24-bit Mode S transponder address.
    ///
    /// Valid range: `0x000000`–`0xFFFFFF`.  Values outside this range are
    /// accepted (no panic) but are not valid ICAO addresses.
    pub icao24: u32,
    /// ASCII callsign, up to 8 characters.
    pub callsign: Cow<'a, str>,
    /// Optional ICAO wake-turbulence category (`A`–`E`).
    pub wake_category: Option<WakeCategory>,
}

impl<'a> Aircraft<'a> {
    /// Construct a new `Aircraft` from a 24-bit ICAO address and a callsign.
    ///
    /// # Arguments
    ///
    /// - `icao24` — ICAO 24-bit Mode S transponder address.
    /// - `callsign` — ASCII callsign string (up to 8 characters).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use siderust::aircraft::Aircraft;
    ///
    /// let ac = Aircraft::new(0x4CA2B5, "EIN104");
    /// assert_eq!(ac.callsign(), "EIN104");
    /// ```
    pub fn new(icao24: u32, callsign: impl Into<Cow<'a, str>>) -> Self {
        Self {
            icao24,
            callsign: callsign.into(),
            wake_category: None,
        }
    }

    /// Attach a wake-turbulence category (builder-style).
    ///
    /// # Arguments
    ///
    /// - `cat` — ICAO wake-turbulence category.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use siderust::aircraft::{Aircraft, WakeCategory};
    ///
    /// let ac = Aircraft::new(0x4CA2B5, "EIN104").with_wake_category(WakeCategory::Heavy);
    /// assert_eq!(ac.wake_category, Some(WakeCategory::Heavy));
    /// ```
    #[must_use]
    pub fn with_wake_category(mut self, cat: WakeCategory) -> Self {
        self.wake_category = Some(cat);
        self
    }

    /// ICAO 24-bit transponder address.
    #[inline]
    pub fn icao24(&self) -> u32 {
        self.icao24
    }

    /// Aircraft callsign.
    #[inline]
    pub fn callsign(&self) -> &str {
        &self.callsign
    }
}

// =============================================================================
// Wake turbulence category
// =============================================================================

/// ICAO wake-turbulence category (Annex 2 / PANS-ATM Doc 4444).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum WakeCategory {
    /// Light: MTOW ≤ 7 000 kg.
    Light,
    /// Medium: 7 000 kg < MTOW ≤ 136 000 kg.
    Medium,
    /// Heavy: MTOW > 136 000 kg (excluding super).
    Heavy,
    /// Super: Airbus A380 and Boeing 747-8.
    Super,
}

// =============================================================================
// AircraftState snapshot
// =============================================================================

/// Instantaneous state snapshot of an aircraft.
///
/// All fields are typed via [`crate::qtty`]; no bare `f64` appears on the
/// public surface.
///
/// - `position` — WGS-84 ellipsoidal geodetic position (longitude, latitude,
///   ellipsoidal height above WGS-84).  The height field carries the
///   geometric (ellipsoidal) altitude; convert with [`barometric_altitude_m`]
///   for the pressure altitude used in ADS-B.
/// - `ground_speed` — horizontal speed over ground in m/s.
/// - `track_angle` — true track (clockwise from true north) in degrees.
/// - `vertical_rate` — climb rate in m/s; positive = climbing.
///
/// # Examples
///
/// ```rust
/// use siderust::aircraft::AircraftState;
/// use siderust::coordinates::centers::Geodetic;
/// use siderust::coordinates::frames::ECEF;
/// use siderust::qtty::{Degrees, Meters};
///
/// let state = AircraftState::new(
///     Geodetic::<ECEF>::new(
///         Degrees::new(-6.270),
///         Degrees::new(53.421),
///         Meters::new(10_668.0),
///     ),
///     Degrees::new(275.0),
/// );
/// assert!(state.position.height.value() > 0.0);
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct AircraftState {
    /// WGS-84 geodetic position (longitude, latitude, ellipsoidal height).
    pub position: Geodetic<ECEF>,
    /// Ground speed over the surface in m/s.
    pub ground_speed: MetersPerSecond,
    /// True track angle, clockwise from true north, in degrees.
    pub track_angle: Degrees,
    /// Vertical (climb) rate in m/s.  Positive = climbing.
    pub vertical_rate: MetersPerSecond,
}

impl AircraftState {
    /// Construct an `AircraftState` with zero ground speed and vertical rate.
    ///
    /// # Arguments
    ///
    /// - `position` — WGS-84 ellipsoidal geodetic position.
    /// - `track_angle` — true track angle in degrees.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use siderust::aircraft::AircraftState;
    /// use siderust::coordinates::centers::Geodetic;
    /// use siderust::coordinates::frames::ECEF;
    /// use siderust::qtty::{Degrees, Meters};
    ///
    /// let s = AircraftState::new(
    ///     Geodetic::<ECEF>::new(Degrees::new(2.35), Degrees::new(48.86), Meters::new(1_000.0)),
    ///     Degrees::new(90.0),
    /// );
    /// assert_eq!(s.vertical_rate.value(), 0.0);
    /// ```
    pub fn new(position: Geodetic<ECEF>, track_angle: Degrees) -> Self {
        Self {
            position,
            ground_speed: MetersPerSecond::new(0.0),
            track_angle,
            vertical_rate: MetersPerSecond::new(0.0),
        }
    }

    /// Geodetic (ellipsoidal) altitude above WGS-84 in metres.
    ///
    /// Note: this is the **ellipsoidal height**, not the barometric pressure
    /// altitude transmitted in ADS-B.  Use [`crate::aircraft::isa::geopotential_altitude_m`]
    /// to convert between geometric and geopotential altitudes.
    ///
    /// # Returns
    ///
    /// [`Meters`] — ellipsoidal height, positive upward.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use siderust::aircraft::AircraftState;
    /// use siderust::coordinates::centers::Geodetic;
    /// use siderust::coordinates::frames::ECEF;
    /// use siderust::qtty::{Degrees, Meters};
    ///
    /// let s = AircraftState::new(
    ///     Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(0.0), Meters::new(5_000.0)),
    ///     Degrees::new(0.0),
    /// );
    /// assert_eq!(s.barometric_altitude_m().value(), 5_000.0);
    /// ```
    #[inline]
    pub fn barometric_altitude_m(&self) -> Meters {
        self.position.height.to()
    }
}

// =============================================================================
// AircraftTrack — time-tagged state with dead-reckoning
// =============================================================================

/// A time-tagged aircraft state that supports dead-reckoning extrapolation.
///
/// `AircraftTrack` pairs an [`AircraftState`] snapshot with a reference
/// [`JulianDate`] so that [`Trackable::track`] can propagate the snapshot to
/// any epoch using a flat-Earth, constant-velocity approximation.
///
/// ## Accuracy
///
/// The dead-reckoning model is accurate to ≲ 1 km for propagation intervals
/// up to roughly 5 minutes at typical cruise speeds.  For longer intervals,
/// or when precision matters, update the reference state from fresh ADS-B
/// messages.
///
/// # Examples
///
/// ```rust
/// use siderust::aircraft::{AircraftState, AircraftTrack, MetersPerSecond};
/// use siderust::coordinates::centers::Geodetic;
/// use siderust::coordinates::frames::ECEF;
/// use siderust::qtty::{Degrees, Meters};
/// use siderust::targets::Trackable;
/// use siderust::time::JulianDate;
///
/// let state = AircraftState {
///     position: Geodetic::<ECEF>::new(
///         Degrees::new(2.35),
///         Degrees::new(48.86),
///         Meters::new(10_000.0),
///     ),
///     ground_speed: MetersPerSecond::new(250.0),
///     track_angle: Degrees::new(90.0),  // due east
///     vertical_rate: MetersPerSecond::new(0.0),
/// };
/// let epoch = JulianDate::new(2_451_545.0); // J2000
/// let track = AircraftTrack::new(state, epoch);
///
/// // One minute later
/// let t1 = JulianDate::new(2_451_545.0 + 60.0 / 86_400.0);
/// let s1 = track.track(t1);
/// // Aircraft moved east at 250 m/s for 60 s = 15 000 m
/// assert!(s1.position.lon.value() > 2.35);
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct AircraftTrack {
    /// Reference state at `epoch`.
    pub state: AircraftState,
    /// Julian Date (TT) at which `state` was observed.
    pub epoch: JulianDate,
}

impl AircraftTrack {
    /// Create a new `AircraftTrack` from a state and a reference epoch.
    ///
    /// # Arguments
    ///
    /// - `state` — instantaneous state snapshot.
    /// - `epoch` — Julian Date (TT) at which `state` was observed.
    pub fn new(state: AircraftState, epoch: JulianDate) -> Self {
        Self { state, epoch }
    }
}

impl Trackable for AircraftTrack {
    type Coords = AircraftState;

    /// Propagate the stored state to `jd` using flat-Earth dead-reckoning.
    ///
    /// The model assumes constant ground speed, constant track angle, and
    /// constant vertical rate.  Latitude and longitude increments are
    /// computed on a sphere of mean radius 6 371 000 m.
    fn track(&self, jd: JulianDate) -> AircraftState {
        // Mean Earth radius (m) — WGS-84 authalic sphere approximation.
        const R_EARTH_M: f64 = 6_371_000.0;
        const DEG_PER_RAD: f64 = 180.0 / std::f64::consts::PI;

        // Elapsed time in seconds.  Raw JD arithmetic is acceptable here
        // because this is a private math kernel that produces typed output.
        let dt_s = (jd.raw().value() - self.epoch.raw().value()) * 86_400.0;

        let gs = self.state.ground_speed.value(); // m/s
        let vr = self.state.vertical_rate.value(); // m/s
        let track_rad = self
            .state
            .track_angle
            .to::<crate::qtty::unit::Radian>()
            .value();
        let lat_rad = self
            .state
            .position
            .lat
            .to::<crate::qtty::unit::Radian>()
            .value();

        let ds_m = gs * dt_s; // horizontal displacement (m)

        let dlat_deg = ds_m * track_rad.cos() / R_EARTH_M * DEG_PER_RAD;
        // Guard against pole singularity (cos(lat) → 0).
        let cos_lat = lat_rad.cos().max(1e-9);
        let dlon_deg = ds_m * track_rad.sin() / (R_EARTH_M * cos_lat) * DEG_PER_RAD;
        let dalt_m = vr * dt_s;

        AircraftState {
            position: Geodetic::<ECEF>::new(
                Degrees::new(self.state.position.lon.value() + dlon_deg),
                Degrees::new(self.state.position.lat.value() + dlat_deg),
                Meters::new(self.state.position.height.value() + dalt_m),
            ),
            ground_speed: self.state.ground_speed,
            track_angle: self.state.track_angle,
            vertical_rate: self.state.vertical_rate,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aircraft_identity_roundtrip() {
        let ac = Aircraft::new(0x4CA2B5_u32, "EIN104");
        assert_eq!(ac.icao24(), 0x4CA2B5);
        assert_eq!(ac.callsign(), "EIN104");
        assert!(ac.wake_category.is_none());
    }

    #[test]
    fn aircraft_with_wake_category() {
        let ac = Aircraft::new(0x400000_u32, "BAW1").with_wake_category(WakeCategory::Heavy);
        assert_eq!(ac.wake_category, Some(WakeCategory::Heavy));
    }

    #[test]
    fn aircraft_state_new_zeroed_kinematics() {
        let state = AircraftState::new(
            Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(0.0), Meters::new(10_000.0)),
            Degrees::new(180.0),
        );
        assert_eq!(state.ground_speed.value(), 0.0);
        assert_eq!(state.vertical_rate.value(), 0.0);
        assert_eq!(state.track_angle.value(), 180.0);
    }

    #[test]
    fn barometric_altitude_m_returns_ellipsoidal_height() {
        let h = Meters::new(9_144.0);
        let state = AircraftState::new(
            Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(0.0), h),
            Degrees::new(0.0),
        );
        assert!((state.barometric_altitude_m().value() - 9_144.0).abs() < 1e-9);
    }

    #[test]
    fn dead_reckoning_due_east() {
        use crate::targets::Trackable;
        use crate::time::JulianDate;

        // Aircraft at equator, 0° lon, heading due east (90°), 100 m/s GS, FL0.
        let state = AircraftState {
            position: Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(0.0), Meters::new(0.0)),
            ground_speed: MetersPerSecond::new(100.0),
            track_angle: Degrees::new(90.0),
            vertical_rate: MetersPerSecond::new(0.0),
        };
        let epoch = JulianDate::new(2_451_545.0);
        let track = AircraftTrack::new(state, epoch);

        // 60 seconds later: should have moved east by 100*60 = 6 000 m.
        let t1 = JulianDate::new(2_451_545.0 + 60.0 / 86_400.0);
        let s1 = track.track(t1);

        // Expected lon increment: 6000 / 6_371_000 * (180/PI) ≈ 0.03386°
        let expected_dlon = 6_000.0 / 6_371_000.0 * (180.0 / std::f64::consts::PI);
        assert!(
            (s1.position.lon.value() - expected_dlon).abs() < 1e-6,
            "lon = {}",
            s1.position.lon.value()
        );
        assert!(
            s1.position.lat.value().abs() < 1e-9,
            "lat should not change"
        );
        assert!(
            s1.position.height.value().abs() < 1e-9,
            "alt should not change"
        );
    }

    #[test]
    fn dead_reckoning_climbing_north() {
        use crate::targets::Trackable;
        use crate::time::JulianDate;

        // Aircraft at equator heading north, 200 m/s GS, climbing at 5 m/s.
        let state = AircraftState {
            position: Geodetic::<ECEF>::new(
                Degrees::new(0.0),
                Degrees::new(0.0),
                Meters::new(1_000.0),
            ),
            ground_speed: MetersPerSecond::new(200.0),
            track_angle: Degrees::new(0.0),
            vertical_rate: MetersPerSecond::new(5.0),
        };
        let epoch = JulianDate::new(2_451_545.0);
        let track = AircraftTrack::new(state, epoch);

        let t1 = JulianDate::new(2_451_545.0 + 120.0 / 86_400.0); // 2 min
        let s1 = track.track(t1);

        // 200 * 120 = 24 000 m north; 5 * 120 = 600 m climb.
        let expected_dlat = 24_000.0 / 6_371_000.0 * (180.0 / std::f64::consts::PI);
        assert!((s1.position.lat.value() - expected_dlat).abs() < 1e-6);
        assert!(s1.position.lon.value().abs() < 1e-9);
        assert!(
            (s1.position.height.value() - 1_600.0).abs() < 0.01,
            "height = {}",
            s1.position.height.value()
        );
    }
}
