use crate::starlight::StellarSurfaceBrightnessMap;

#[must_use]
pub fn to_csv(map: &StellarSurfaceBrightnessMap) -> String {
    let mut out = String::new();
    for (index, value) in map.values().iter().enumerate() {
        let line = format!("{index},{},{},{}\n", value.integrated_ph_cm2_ns_sr, value.b_s10, value.v_s10);
        out.push_str(&line);
    }
    out
}
