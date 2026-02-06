// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Night-related solar types (twilight thresholds, etc.).

use qtty::Degrees;

/// Common twilight types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Twilight {
    Civil,
    Nautical,
    Astronomical,
    Horizon,
    ApparentHorizon,
}

impl From<Twilight> for Degrees {
    fn from(t: Twilight) -> Degrees {
        match t {
            Twilight::Civil => Degrees::new(-6.0),
            Twilight::Nautical => Degrees::new(-12.0),
            Twilight::Astronomical => Degrees::new(-18.0),
            Twilight::Horizon => Degrees::new(0.0),
            Twilight::ApparentHorizon => Degrees::new(-0.833),
        }
    }
}

