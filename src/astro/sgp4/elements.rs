// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Bridge from [`crate::formats::tle::Tle`] to the SGP4 input record.
//!
//! The SGP4 backend ([`sgp4` crate](https://crates.io/crates/sgp4)) consumes
//! its own [`sgp4::Elements`] struct. We treat this purely as an internal
//! record format — callers always interact with the typed
//! [`crate::formats::tle::Tle`].

use crate::formats::tle::{Classification as TleClass, Tle};
use chrono::Datelike;
use sgp4::Elements;

use super::Sgp4Error;

/// Convert a typed [`Tle`] into the dimensionless record consumed by the SGP4
/// backend.
///
/// The conversion preserves all SGP4-relevant fields (epoch, mean elements,
/// drag terms, classification, international designator) and is lossless in
/// both directions for the subset SGP4 uses.
pub(crate) fn tle_to_elements(tle: &Tle) -> Result<Elements, Sgp4Error> {
    let dt_utc = tle
        .epoch
        .try_to_chrono()
        .map_err(|e| Sgp4Error::TimeConversion(format!("TLE epoch → chrono failed: {e:?}")))?;
    // SGP4 wants a NaiveDateTime in UTC.
    let datetime = dt_utc.naive_utc();
    // Sanity-check the year: AFSPC epoch encoding loses precision before 1957.
    if datetime.year() < 1950 {
        return Err(Sgp4Error::InvalidEpoch("year < 1950 not representable"));
    }

    let classification = match tle.classification {
        TleClass::Unclassified => sgp4::Classification::Unclassified,
        TleClass::Classified => sgp4::Classification::Classified,
        TleClass::Secret => sgp4::Classification::Secret,
    };

    Ok(Elements {
        object_name: tle.name.clone(),
        international_designator: Some(tle.international_designator.0.clone()),
        norad_id: u64::from(tle.norad_id.0),
        classification,
        datetime,
        mean_motion_dot: tle.mean_motion_dot,
        mean_motion_ddot: tle.mean_motion_ddot,
        drag_term: tle.bstar,
        element_set_number: u64::from(tle.element_set_number),
        inclination: tle.inclination.value(),
        right_ascension: tle.raan.value(),
        eccentricity: tle.eccentricity,
        argument_of_perigee: tle.argument_of_perigee.value(),
        mean_anomaly: tle.mean_anomaly.value(),
        mean_motion: tle.mean_motion.value(),
        revolution_number: u64::from(tle.revolution_number_at_epoch),
        ephemeris_type: 0,
    })
}
