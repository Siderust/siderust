//! CSV serialization for generated stellar surface-brightness maps.
//!
//! The CSV representation is intentionally simple and deterministic so it can
//! be consumed by downstream tools that build runtime data assets. Each row
//! stores the zero-based HEALPix pixel index and the stellar brightness values
//! associated with that pixel.

use crate::starlight::StellarSurfaceBrightnessMap;

/// Serialize a generated stellar surface-brightness map to CSV text.
///
/// Rows are emitted in the map grid ordering. The current schema is:
/// `healpix_index,integrated_ph_cm2_ns_sr,b_s10,v_s10`.
#[must_use]
pub fn to_csv(map: &StellarSurfaceBrightnessMap) -> String {
    let mut out = String::from("healpix_index,integrated_ph_cm2_ns_sr,b_s10,v_s10\n");
    for (index, value) in map.values().iter().enumerate() {
        let line = format!(
            "{index},{},{},{}\n",
            value.integrated_ph_cm2_ns_sr, value.b_s10, value.v_s10
        );
        out.push_str(&line);
    }
    out
}
