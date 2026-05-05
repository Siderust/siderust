// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Airglow emission geometry.
//!
//! Airglow is emitted in a finite-height atmospheric layer rather than at
//! infinity. The Van Rhijn factor is the standard spherical-shell correction
//! used to convert a zenith-view brightness into the brightness expected at an
//! arbitrary zenith distance.
//!
//! In physical terms, the factor accounts for the fact that an oblique line of
//! sight crosses a longer segment of the emitting layer than a vertical one.
//! For a thin, uniform shell it is therefore a geometric path-length multiplier:
//! values near 1 correspond to a nearly vertical sight line, while larger
//! values toward the horizon reflect the longer path through the layer. Unlike
//! a plane-parallel secant approximation, the spherical-shell form remains
//! finite at the horizon because the emitting layer is located at finite
//! altitude above the surface.

use crate::qtty::{Kilometers, Radians};

/// Van Rhijn path-length factor for a thin emitting layer at `emission_height`.
///
/// The Van Rhijn factor is the ratio between:
///
/// - the effective path length through a thin emitting shell seen at zenith
///   distance `z`, and
/// - the corresponding path length for the same shell observed at the zenith.
///
/// It is commonly used for airglow and similar upper-atmosphere emissions when
/// a measured or modeled zenith intensity must be projected to another viewing
/// direction under the assumption that the emitting layer is horizontally
/// uniform.
///
/// The factor is
///
/// ```text
/// V(z, h) = 1 / sqrt(1 - (R / (R + h) * sin z)^2)
/// ```
///
/// where `R` is the Earth's mean radius, `h` is the emission-layer height,
/// and `z` is the apparent zenith distance.  It is exactly 1 at zenith and
/// remains finite at the horizon because the emitting layer is at finite
/// altitude.
pub fn van_rhijn_factor(zenith: Radians, emission_height: Kilometers) -> f64 {
    van_rhijn_factor_with_radius(
        zenith,
        emission_height,
        crate::bodies::solar_system::EARTH.radius,
    )
}

/// Variant of [`van_rhijn_factor`] with an explicit planetary radius.
pub fn van_rhijn_factor_with_radius(
    zenith: Radians,
    emission_height: Kilometers,
    body_radius: Kilometers,
) -> f64 {
    let r = body_radius.value();
    let h = emission_height.value();
    if !zenith.is_finite() || !h.is_finite() || !r.is_finite() || h <= 0.0 || r <= 0.0 {
        return f64::NAN;
    }
    let ratio = r / (r + h);
    let s = zenith.sin();
    let inner = 1.0 - (ratio * s) * (ratio * s);
    if inner <= 0.0 {
        f64::INFINITY
    } else {
        inner.sqrt().recip()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::{unit::Radian, Degrees};

    #[test]
    fn van_rhijn_is_unity_at_zenith() {
        let f = van_rhijn_factor(Degrees::new(0.0).to::<Radian>(), Kilometers::new(90.0));
        assert!((f - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn van_rhijn_increases_toward_horizon_but_stays_finite() {
        let h = Kilometers::new(90.0);
        let z45 = van_rhijn_factor(Degrees::new(45.0).to::<Radian>(), h);
        let z80 = van_rhijn_factor(Degrees::new(80.0).to::<Radian>(), h);
        let z90 = van_rhijn_factor(Degrees::new(90.0).to::<Radian>(), h);
        assert!(z45 > 1.0);
        assert!(z80 > z45);
        assert!(z90 > z80);
        assert!(z90.is_finite());
    }
}
