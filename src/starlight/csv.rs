// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::coordinates::frames::Galactic;
use crate::healpix::HealpixMap;
use crate::starlight::{StellarMapProvenance, StellarSurfaceBrightness};

/// Serialize a Galactic stellar surface-brightness map to a simple CSV string.
#[must_use]
pub fn stellar_map_to_csv(
    map: &HealpixMap<Galactic, StellarSurfaceBrightness>,
    provenance: &StellarMapProvenance,
) -> String {
    let mut out = String::new();
    out.push_str("# map_type=healpix\n");
    out.push_str("# schema_version=stellar_healpix_csv_v1\n");
    out.push_str(&format!("# nside={}\n", map.grid().nside().get()));
    out.push_str(&format!("# ordering={:?}\n", map.grid().ordering()).to_lowercase());
    out.push_str("# coordinate_frame=galactic\n");
    push_header(&mut out, "dataset_name", &provenance.dataset_name);
    push_header(&mut out, "version", &provenance.version);
    push_header(
        &mut out,
        "generation_date_utc",
        &provenance.generation_date_utc,
    );
    push_header(&mut out, "source_catalogue", &provenance.source_catalogue);
    push_optional_header(
        &mut out,
        "source_catalogue_release",
        provenance.source_catalogue_release.as_deref(),
    );
    push_optional_header(
        &mut out,
        "source_catalogue_license",
        provenance.source_catalogue_license.as_deref(),
    );
    push_optional_header(
        &mut out,
        "source_catalogue_checksum",
        provenance.source_catalogue_checksum.as_deref(),
    );
    push_optional_header(
        &mut out,
        "magnitude_limit",
        provenance.magnitude_limit.as_deref(),
    );
    push_header(&mut out, "photometry_model", &provenance.photometry_model);
    push_header(&mut out, "band_definition", &provenance.band_definition);
    push_optional_header(&mut out, "smoothing_fwhm_deg", provenance.smoothing.as_deref());
    out.push_str("healpix_index,integrated_ph_cm2_ns_sr,b_s10,v_s10\n");
    for (index, value) in map.values().iter().enumerate() {
        out.push_str(&format!(
            "{},{:.17e},{:.17e},{:.17e}\n",
            index, value.integrated_ph_cm2_ns_sr, value.b_s10, value.v_s10
        ));
    }
    out
}

fn push_header(out: &mut String, key: &str, value: &str) {
    out.push_str(&format!("# {key}={}\n", sanitize_header_value(value)));
}

fn push_optional_header(out: &mut String, key: &str, value: Option<&str>) {
    if let Some(value) = value {
        push_header(out, key, value);
    }
}

fn sanitize_header_value(value: &str) -> String {
    value.replace('\r', " ").replace('\n', " ")
}
