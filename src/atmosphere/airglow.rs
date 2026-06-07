// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Airglow emission geometry
//!
//! ## Scientific scope
//!
//! Airglow is emitted in a thin, finite-altitude atmospheric layer rather
//! than at infinity, so its apparent surface brightness depends on the
//! geometric path length through the layer. The Van Rhijn factor is the
//! standard spherical-shell correction used to project a measured or
//! modelled zenith brightness to an arbitrary zenith distance under the
//! assumption of a horizontally uniform emitting layer.
//!
//! Unlike a plane-parallel `sec z` approximation, the Van Rhijn factor
//! remains finite at the horizon because the emitting shell sits at finite
//! altitude above the surface.
//!
//! ## Technical scope
//!
//! - Inputs: zenith distance as typed [`crate::qtty::Radians`], emission-layer height
//!   and (optionally) planetary radius as typed [`Kilometers`].
//! - Output: dimensionless geometric factor as
//!   [`Quantity<ScatteringFactor>`] (see
//!   [`crate::atmosphere::ScatteringFactor`]).
//!
//! ## References
//!
//! - Van Rhijn, P. J. (1921). *Publications of the Astronomical
//!   Laboratory at Groningen* 31, 1.

use crate::atmosphere::ScatteringFactor;
use crate::ext_qtty::Quantity;
use crate::qtty::{Kilometers, Radians};

/// Van Rhijn path-length factor for a thin emitting layer at
/// `emission_height`, evaluated at zenith distance `zenith`.
///
/// The factor is
///
/// ```text
/// V(z, h) = 1 / sqrt(1 - (R / (R + h) · sin z)^2)
/// ```
///
/// where `R` is the Earth's mean radius and `h` is the emission-layer
/// height. It is exactly 1 at zenith and remains finite at the horizon.
pub fn van_rhijn_factor(
    zenith: Radians,
    emission_height: Kilometers,
) -> Quantity<ScatteringFactor> {
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
) -> Quantity<ScatteringFactor> {
    let r = body_radius.value();
    let h = emission_height.value();
    if !zenith.is_finite() || !h.is_finite() || !r.is_finite() || h <= 0.0 || r <= 0.0 {
        return Quantity::<ScatteringFactor>::new(f64::NAN);
    }
    let ratio = r / (r + h);
    let s = zenith.sin();
    let inner = 1.0 - (ratio * s) * (ratio * s);
    let v = if inner <= 0.0 {
        f64::INFINITY
    } else {
        inner.sqrt().recip()
    };
    Quantity::<ScatteringFactor>::new(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::{unit::Radian, Degrees};

    #[test]
    fn van_rhijn_is_unity_at_zenith() {
        let f = van_rhijn_factor(Degrees::new(0.0).to::<Radian>(), Kilometers::new(90.0));
        assert!((f.value() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn van_rhijn_increases_toward_horizon_but_stays_finite() {
        let h = Kilometers::new(90.0);
        let z45 = van_rhijn_factor(Degrees::new(45.0).to::<Radian>(), h).value();
        let z80 = van_rhijn_factor(Degrees::new(80.0).to::<Radian>(), h).value();
        let z90 = van_rhijn_factor(Degrees::new(90.0).to::<Radian>(), h).value();
        assert!(z45 > 1.0);
        assert!(z80 > z45);
        assert!(z90 > z80);
        assert!(z90.is_finite());
    }
}
