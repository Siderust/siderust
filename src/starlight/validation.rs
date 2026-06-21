// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::coordinates::frames::Galactic;
use crate::healpix::{HealpixError, HealpixIndex, HealpixMap};
use crate::starlight::{Result, StellarMapError, StellarSurfaceBrightness};

/// Validate that a stellar map contains finite, non-negative values.
pub fn validate_stellar_map_values(
    map: &HealpixMap<Galactic, StellarSurfaceBrightness>,
) -> Result<()> {
    crate::healpix::validate_healpix_map_complete(map)?;
    if map.values().iter().all(|value| value.is_valid()) {
        Ok(())
    } else {
        Err(StellarMapError::Healpix(HealpixError::InvalidMapValue))
    }
}

/// Validate conservation of B/V flux after area-normalized binning.
pub fn validate_flux_conservation(
    input_b_flux_sum: f64,
    input_v_flux_sum: f64,
    map: &HealpixMap<Galactic, StellarSurfaceBrightness>,
    tolerance: f64,
) -> Result<()> {
    validate_stellar_map_values(map)?;
    if !input_b_flux_sum.is_finite() || !input_v_flux_sum.is_finite() {
        return Err(StellarMapError::Validation("flux inputs must be finite"));
    }
    if !tolerance.is_finite() || tolerance < 0.0 {
        return Err(StellarMapError::Validation(
            "tolerance must be finite and non-negative",
        ));
    }

    let area = map.grid().pixel_area_deg2();
    let b_sum = map
        .values()
        .iter()
        .map(|value| value.b_s10 * area)
        .sum::<f64>();
    let v_sum = map
        .values()
        .iter()
        .map(|value| value.v_s10 * area)
        .sum::<f64>();

    if (b_sum - input_b_flux_sum).abs() <= tolerance
        && (v_sum - input_v_flux_sum).abs() <= tolerance
    {
        Ok(())
    } else {
        Err(StellarMapError::Validation("B/V flux is not conserved"))
    }
}

/// Validate that the Galactic plane is brighter than high-latitude poles.
///
/// The validator uses a deterministic broad-brush diagnostic: mean V S10 for
/// `|b| <= 10 deg` divided by mean V S10 for `|b| >= 60 deg`.
pub fn validate_plane_pole_contrast(
    map: &HealpixMap<Galactic, StellarSurfaceBrightness>,
    min_ratio: f64,
) -> Result<()> {
    validate_stellar_map_values(map)?;
    if !min_ratio.is_finite() || min_ratio < 0.0 {
        return Err(StellarMapError::Validation(
            "min_ratio must be finite and non-negative",
        ));
    }

    let mut plane_sum = 0.0;
    let mut plane_count = 0_u64;
    let mut pole_sum = 0.0;
    let mut pole_count = 0_u64;

    for (raw, value) in map.values().iter().enumerate() {
        let (_, lat) = map
            .grid()
            .pixel_to_lon_lat_rad_checked(HealpixIndex::new(
                u64::try_from(raw).expect("index fits u64"),
            ))?;
        let lat_deg = lat.to_degrees().abs();
        if lat_deg <= 10.0 {
            plane_sum += value.v_s10;
            plane_count += 1;
        } else if lat_deg >= 60.0 {
            pole_sum += value.v_s10;
            pole_count += 1;
        }
    }

    if plane_count == 0 || pole_count == 0 {
        return Err(StellarMapError::Validation("plane/pole regions are empty"));
    }

    let plane_mean = plane_sum / plane_count as f64;
    let pole_mean = pole_sum / pole_count as f64;
    let ratio = if pole_mean == 0.0 {
        if plane_mean > 0.0 {
            f64::INFINITY
        } else {
            1.0
        }
    } else {
        plane_mean / pole_mean
    };

    if ratio >= min_ratio {
        Ok(())
    } else {
        Err(StellarMapError::Validation(
            "Galactic plane/pole contrast is too low",
        ))
    }
}

/// Validate that there is no strong artificial discontinuity at `l = 0/360`.
///
/// The diagnostic compares mean V S10 in two symmetric longitude strips near the
/// wrap seam for pixels with `|b| <= 30 deg`.
pub fn validate_no_longitude_wrap_artifact(
    map: &HealpixMap<Galactic, StellarSurfaceBrightness>,
    max_relative_jump: f64,
) -> Result<()> {
    validate_stellar_map_values(map)?;
    if !max_relative_jump.is_finite() || max_relative_jump < 0.0 {
        return Err(StellarMapError::Validation(
            "max_relative_jump must be finite and non-negative",
        ));
    }

    let mut low_sum = 0.0;
    let mut low_count = 0_u64;
    let mut high_sum = 0.0;
    let mut high_count = 0_u64;

    for (raw, value) in map.values().iter().enumerate() {
        let (lon, lat) = map
            .grid()
            .pixel_to_lon_lat_rad_checked(HealpixIndex::new(
                u64::try_from(raw).expect("index fits u64"),
            ))?;
        let lon_deg = lon.to_degrees();
        let lat_deg = lat.to_degrees().abs();
        if lat_deg <= 30.0 && lon_deg <= 10.0 {
            low_sum += value.v_s10;
            low_count += 1;
        } else if lat_deg <= 30.0 && lon_deg >= 350.0 {
            high_sum += value.v_s10;
            high_count += 1;
        }
    }

    if low_count == 0 || high_count == 0 {
        return Err(StellarMapError::Validation("longitude wrap regions are empty"));
    }

    let low_mean = low_sum / low_count as f64;
    let high_mean = high_sum / high_count as f64;
    let scale = low_mean.abs().max(high_mean.abs()).max(1.0);
    let jump = (low_mean - high_mean).abs() / scale;

    if jump <= max_relative_jump {
        Ok(())
    } else {
        Err(StellarMapError::Validation(
            "longitude wrap jump exceeds tolerance",
        ))
    }
}
