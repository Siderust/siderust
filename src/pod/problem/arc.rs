//! Arc definitions for POD runs.
//!
//! A POD *arc* is a contiguous span of time over which a single dynamic
//! and observational hypothesis is fitted. Per design §6.1 the arc is
//! identified by a stable opaque [`ArcId`] string and bounded by typed
//! [`Time<TT>`] start/stop instants, with an optional integration step
//! hint expressed as a typed [`qtty::Second`].
//!
//! Arcs are deliberately TT-anchored so downstream propagation code does
//! not have to disambiguate which time scale `start`/`stop` are expressed
//! in. Conversions from UTC inputs happen at the IO boundary.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use qtty::Second;
use tempoch::{Time, TT};

/// Opaque identifier for an arc.
///
/// Stable across reruns of the same configuration; used as the key in
/// run-level outputs (residuals, products, manifests) so arcs can be
/// joined across artifacts without ambiguity.
///
/// # Examples
///
/// ```
/// use siderust::pod::problem::arc::ArcId;
/// let id = ArcId::new("LEO-2026-05-11");
/// assert_eq!(id.as_str(), "LEO-2026-05-11");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ArcId(String);

impl ArcId {
    /// Wrap a string into an [`ArcId`].
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::pod::problem::arc::ArcId;
    /// let id = ArcId::new("arc-001");
    /// assert_eq!(id.as_str(), "arc-001");
    /// ```
    pub fn new(s: impl Into<String>) -> Self {
        Self(s.into())
    }

    /// Borrow the underlying string.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::pod::problem::arc::ArcId;
    /// assert_eq!(ArcId::new("a").as_str(), "a");
    /// ```
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl std::fmt::Display for ArcId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

/// Inclusive arc time interval and identification metadata.
///
/// Anchored on the [`TT`] time scale per design §6.1; the IO layer is
/// responsible for converting external (UTC, GPS, …) timestamps into TT
/// before constructing an [`ArcDefinition`].
///
/// `step_hint` is *advisory* — propagators may pick a different step. It
/// is preserved in the [`crate::pod::run::manifest::RunManifest`] so reruns can
/// reproduce the original hint exactly.
///
/// # Examples
///
/// ```
/// use siderust::pod::problem::arc::{ArcDefinition, ArcId};
/// use tempoch::{J2000Seconds, TT};
/// use qtty::Second;
///
/// let start = J2000Seconds::<TT>::try_new(Second::new(0.0)).unwrap();
/// let stop  = J2000Seconds::<TT>::try_new(Second::new(86_400.0)).unwrap();
/// let arc = ArcDefinition {
///     id: ArcId::new("demo"),
///     start,
///     stop,
///     step_hint: Some(Second::new(60.0)),
/// };
/// assert!((arc.duration().value() - 86_400.0).abs() < 1e-6);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct ArcDefinition {
    /// Stable identifier for this arc.
    pub id: ArcId,
    /// Inclusive start instant on the [`TT`] scale.
    pub start: Time<TT>,
    /// Inclusive stop instant on the [`TT`] scale.
    pub stop: Time<TT>,
    /// Optional integration-step hint as a typed duration in seconds.
    pub step_hint: Option<Second>,
}

impl ArcDefinition {
    /// Total wall duration of the arc as a typed [`Second`].
    ///
    /// Returns a non-negative duration when `stop >= start`. The value is
    /// well-defined-but-possibly-negative for inverted arcs; callers
    /// should [`Self::is_valid`] first if input ordering is not trusted.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::pod::problem::arc::{ArcDefinition, ArcId};
    /// use tempoch::{J2000Seconds, TT};
    /// use qtty::Second;
    /// let arc = ArcDefinition {
    ///     id: ArcId::new("a"),
    ///     start: J2000Seconds::<TT>::try_new(Second::new(0.0)).unwrap(),
    ///     stop:  J2000Seconds::<TT>::try_new(Second::new(86_400.0)).unwrap(),
    ///     step_hint: None,
    /// };
    /// assert!((arc.duration().value() - 86_400.0).abs() < 1e-6);
    /// ```
    pub fn duration(&self) -> Second {
        self.stop - self.start
    }

    /// `true` iff `stop >= start`.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::pod::problem::arc::{ArcDefinition, ArcId};
    /// use tempoch::{J2000Seconds, TT};
    /// use qtty::Second;
    /// let arc = ArcDefinition {
    ///     id: ArcId::new("a"),
    ///     start: J2000Seconds::<TT>::try_new(Second::new(0.0)).unwrap(),
    ///     stop:  J2000Seconds::<TT>::try_new(Second::new(86_400.0)).unwrap(),
    ///     step_hint: None,
    /// };
    /// assert!(arc.is_valid());
    /// ```
    pub fn is_valid(&self) -> bool {
        self.stop >= self.start
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempoch::J2000Seconds;

    fn at(secs: f64) -> Time<TT> {
        J2000Seconds::<TT>::try_new(Second::new(secs)).unwrap()
    }

    #[test]
    fn duration_in_seconds_matches_one_day() {
        let arc = ArcDefinition {
            id: ArcId::new("d"),
            start: at(0.0),
            stop: at(86_400.0),
            step_hint: None,
        };
        let dt = arc.duration().value();
        assert!((dt - 86_400.0).abs() < 1e-6);
        assert!(arc.is_valid());
    }

    #[test]
    fn inverted_arc_is_invalid() {
        let arc = ArcDefinition {
            id: ArcId::new("d"),
            start: at(86_400.0),
            stop: at(0.0),
            step_hint: None,
        };
        assert!(!arc.is_valid());
    }

    #[test]
    fn arc_id_round_trip() {
        let id = ArcId::new("orbit-7");
        assert_eq!(id.as_str(), "orbit-7");
        assert_eq!(format!("{id}"), "orbit-7");
    }
}
