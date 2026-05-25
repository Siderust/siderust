// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Generate Sun-Earth Lagrange Chebyshev source files.
//!
//! This development binary is enabled by the `lagrange-centers` feature. It
//! loads a selected ephemeris backend, fits L1-L5 over the requested Julian-Date
//! span, and writes Rust source fragments into the requested output directory.

#[cfg(feature = "lagrange-centers")]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    use siderust::ephemeris::lagrange::fit::FitConfig;
    use siderust::ephemeris::{RuntimeEphemeris, Vsop87Ephemeris};
    use siderust::qtty::Seconds;
    use std::env;
    use std::fs;
    use std::path::PathBuf;

    let args: Vec<String> = env::args().collect();
    let source = value_arg(&args, "--source").unwrap_or("vsop87");
    let bsp_path = value_arg(&args, "--bsp-path");
    let out =
        PathBuf::from(value_arg(&args, "--out").unwrap_or("src/embedded_data/lagrange/generated"));
    let from = parse_jd(value_arg(&args, "--from"), 2_451_545.0)?;
    let to = parse_jd(value_arg(&args, "--to"), 2_451_577.0)?;

    fs::create_dir_all(&out)?;
    let config = FitConfig {
        from,
        to,
        block: Seconds::new(32.0 * siderust::qtty::time::SECONDS_PER_DAY),
        validation_step: Seconds::new(6.0 * 3_600.0),
        ..FitConfig::default()
    };

    match source {
        "vsop87" => generate_for(&Vsop87Ephemeris, &out, config)?,
        "de440" | "de441" => {
            let path = bsp_path.ok_or("--bsp-path is required for de440/de441")?;
            let eph = RuntimeEphemeris::from_bsp(path)?;
            generate_for(&eph, &out, config)?;
        }
        other => return Err(format!("unsupported --source value: {other}").into()),
    }
    Ok(())
}

#[cfg(not(feature = "lagrange-centers"))]
fn main() {
    eprintln!("generate-lagrange-cheby requires --features lagrange-centers");
}

#[cfg(feature = "lagrange-centers")]
fn value_arg<'a>(args: &'a [String], name: &str) -> Option<&'a str> {
    args.windows(2)
        .find_map(|pair| (pair[0] == name).then_some(pair[1].as_str()))
}

#[cfg(feature = "lagrange-centers")]
fn parse_jd(
    value: Option<&str>,
    default: f64,
) -> Result<siderust::JulianDate, Box<dyn std::error::Error>> {
    let raw = value.map_or(Ok(default), str::parse::<f64>)?;
    siderust::time::try_jd_f64(raw).map_err(|err| err.into())
}

#[cfg(feature = "lagrange-centers")]
fn generate_for(
    ephemeris: &dyn siderust::ephemeris::DynEphemeris,
    out: &std::path::Path,
    config: siderust::ephemeris::lagrange::fit::FitConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    use siderust::ephemeris::lagrange::SunEarthLagrangePoint;
    use std::fs;

    let points = [
        (SunEarthLagrangePoint::L1, "l1"),
        (SunEarthLagrangePoint::L2, "l2"),
        (SunEarthLagrangePoint::L3, "l3"),
        (SunEarthLagrangePoint::L4, "l4"),
        (SunEarthLagrangePoint::L5, "l5"),
    ];
    let mut metadata = String::new();
    for (point, name) in points {
        let fitted =
            siderust::ephemeris::lagrange::fit::fit_sun_earth_lagrange(ephemeris, point, config)?;
        let const_name = format!("RECORDS_{}", point.label());
        let body = format!(
            "pub const {const_name}: &[f64] = &{:?};\npub const NCOEFF_{}: usize = {};\n",
            fitted.records,
            point.label(),
            fitted.ncoeff
        );
        fs::write(out.join(format!("{name}.rs")), body)?;
        metadata.push_str(&format!(
            "{} max_abs_m={} rms_m={} ncoeff={}\n",
            point.label(),
            fitted.stats.max_abs_error.value(),
            fitted.stats.rms_error.value(),
            fitted.ncoeff
        ));
    }
    fs::write(out.join("metadata.txt"), metadata)?;
    Ok(())
}