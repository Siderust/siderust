// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Conveniences for the small set of NAIF body identifiers exercised by
//! POD workflows. The full NAIF id namespace is not enumerated here;
//! only the bodies covered by JPL DE planetary kernels are exposed
//! through [`well_known`].

/// Resolve a case-insensitive textual body name to its NAIF integer id.
///
/// Returns `None` for any name not in the small table covered by JPL
/// DE planetary kernels (Sun, planets, EMB, Moon, SSB).
///
/// # Examples
///
/// ```rust
/// use siderust::formats::spice::naif_id_for_name;
///
/// assert_eq!(naif_id_for_name("EARTH"), Some(399));
/// assert_eq!(naif_id_for_name("solar system barycenter"), Some(0));
/// assert_eq!(naif_id_for_name("Pluto"), Some(9));
/// assert_eq!(naif_id_for_name("Andromeda"), None);
/// ```
pub fn naif_id_for_name(name: &str) -> Option<i32> {
    let n = name.trim().to_ascii_uppercase();
    match n.as_str() {
        "SSB" | "SOLAR SYSTEM BARYCENTER" => Some(0),
        "MERCURY BARYCENTER" => Some(1),
        "VENUS BARYCENTER" => Some(2),
        "EMB" | "EARTH BARYCENTER" | "EARTH MOON BARYCENTER" => Some(3),
        "MARS BARYCENTER" => Some(4),
        "JUPITER BARYCENTER" => Some(5),
        "SATURN BARYCENTER" => Some(6),
        "URANUS BARYCENTER" => Some(7),
        "NEPTUNE BARYCENTER" => Some(8),
        "PLUTO BARYCENTER" | "PLUTO" => Some(9),
        "SUN" => Some(10),
        "MERCURY" => Some(199),
        "VENUS" => Some(299),
        "EARTH" => Some(399),
        "MARS" => Some(499),
        "JUPITER" => Some(599),
        "SATURN" => Some(699),
        "URANUS" => Some(799),
        "NEPTUNE" => Some(899),
        "MOON" => Some(301),
        _ => None,
    }
}

/// Compile-time NAIF id constants for the bodies exercised by POD.
///
/// These mirror the constants in `siderust::data::spk` and the NAIF
/// "Integer ID Codes" reference. They are duplicated here so callers
/// can write `well_known::EARTH` / `well_known::MOON` without
/// reaching into upstream modules.
pub mod well_known {
    /// Solar System Barycenter.
    pub const SSB: i32 = 0;
    /// Earth-Moon Barycenter.
    pub const EARTH_MOON_BARYCENTER: i32 = 3;
    /// Mars-system Barycenter.
    pub const MARS_BARYCENTER: i32 = 4;
    /// Jupiter-system Barycenter.
    pub const JUPITER_BARYCENTER: i32 = 5;
    /// Sun.
    pub const SUN: i32 = 10;
    /// Earth (planet, NAIF id 399).
    pub const EARTH: i32 = 399;
    /// Moon (NAIF id 301).
    pub const MOON: i32 = 301;
    /// Mars (planet, NAIF id 499).
    pub const MARS: i32 = 499;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn name_lookup_handles_common_aliases() {
        assert_eq!(naif_id_for_name("SSB"), Some(0));
        assert_eq!(naif_id_for_name("Sun"), Some(10));
        assert_eq!(naif_id_for_name("  earth "), Some(399));
        assert_eq!(naif_id_for_name("emb"), Some(3));
        assert_eq!(naif_id_for_name("MoOn"), Some(301));
        assert_eq!(naif_id_for_name("notabody"), None);
    }

    #[test]
    fn well_known_ids_are_consistent_with_lookup() {
        assert_eq!(naif_id_for_name("EARTH"), Some(well_known::EARTH));
        assert_eq!(naif_id_for_name("MOON"), Some(well_known::MOON));
        assert_eq!(naif_id_for_name("SUN"), Some(well_known::SUN));
        assert_eq!(
            naif_id_for_name("EMB"),
            Some(well_known::EARTH_MOON_BARYCENTER)
        );
    }
}
