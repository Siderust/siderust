// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use std::{
    collections::BTreeMap,
    env,
    fs::{self, File},
    io::{BufRead, BufReader, Read},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result};
use regex::Regex;
use std::fmt::Write;

/// Compute a stable SHA-256 fingerprint of the ELP2000 input dataset
/// directory by hashing every `ELP*` file in sorted name order.
///
/// Returned as `(hex_digest, total_bytes)`. The hex digest is checked
/// against `SIDERUST_ELP2000_SHA256` when set, allowing reproducible
/// builds to refuse unverified coefficient data.
fn dataset_sha256(data_dir: &Path) -> Result<(String, u64)> {
    use sha2::{Digest, Sha256};

    let mut paths: Vec<PathBuf> = fs::read_dir(data_dir)
        .with_context(|| format!("read-dir {data_dir:?}"))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .map(|s| s.starts_with("ELP"))
                .unwrap_or(false)
        })
        .collect();
    paths.sort();

    let mut hasher = Sha256::new();
    let mut total: u64 = 0;
    for path in &paths {
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            hasher.update(name.as_bytes());
            hasher.update(b"\0");
        }
        let mut f = File::open(path).with_context(|| format!("open {path:?}"))?;
        let mut buf = [0u8; 64 * 1024];
        loop {
            let n = f.read(&mut buf).with_context(|| format!("read {path:?}"))?;
            if n == 0 {
                break;
            }
            total += n as u64;
            hasher.update(&buf[..n]);
        }
    }
    let digest = hasher.finalize();
    let mut s = String::with_capacity(64);
    for b in digest {
        let _ = write!(s, "{:02x}", b);
    }
    Ok((s, total))
}

fn check_pinned_sha256(actual: &str, env_var: &str) -> Result<()> {
    if let Ok(expected) = env::var(env_var) {
        let expected = expected.trim().to_ascii_lowercase();
        if !expected.is_empty() && expected != actual {
            anyhow::bail!(
                "ELP2000: SHA-256 mismatch.\n  expected: {expected}\n  actual:   {actual}\n\
                 Refusing to regenerate tables from unverified data."
            );
        }
    }
    Ok(())
}

/// ---------------------------  PARSER  ---------------------------

#[derive(Debug)]
struct Entry {
    ints: Vec<i64>,
    floats: Vec<String>, // We use String to avoid formating the values in generate_rust
}

fn file_format(key: &str) -> Option<(usize, usize)> {
    let idx: u32 = key.trim_start_matches("ELP").parse().ok()?;
    match idx {
        1..=3 => Some((4, 7)),
        10..=21 => Some((11, 3)),
        4..=9 | 22..=36 => Some((5, 3)),
        _ => None,
    }
}

fn parse_line(line: &str, n_ints: usize, n_floats: usize, token_re: &Regex) -> Result<Entry> {
    let tokens: Vec<&str> = token_re.find_iter(line).map(|m| m.as_str()).collect();
    if tokens.len() != n_ints + n_floats {
        anyhow::bail!(
            "Unexpected column count {} (wanted {})",
            tokens.len(),
            n_ints + n_floats
        );
    }
    let ints = tokens[..n_ints]
        .iter()
        .map(|t| t.parse::<i64>().map_err(|e| anyhow::anyhow!(e)))
        .collect::<Result<Vec<_>>>()?;
    let floats = tokens[n_ints..].iter().map(|t| t.to_string()).collect();
    Ok(Entry { ints, floats })
}

fn parse_file(path: &Path, n_ints: usize, n_floats: usize) -> Result<Vec<Entry>> {
    let file = File::open(path).with_context(|| format!("open {path:?}"))?;
    let reader = BufReader::new(file);
    let line_re = Regex::new(r"^\s*[-+]?\d").unwrap();
    let token_re = Regex::new(r"[-+]?\d+\.\d+|[-+]?\d+").unwrap();

    reader
        .lines()
        .map_while(Result::ok)
        .filter(|l| line_re.is_match(l))
        .map(|l| parse_line(&l, n_ints, n_floats, &token_re))
        .collect()
}

fn parse_all_elps(dir: &Path) -> Result<BTreeMap<String, Vec<Entry>>> {
    let mut map = BTreeMap::new();

    let mut paths: Vec<_> = fs::read_dir(dir)
        .with_context(|| format!("read-dir {dir:?}"))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .map(|s| s.starts_with("ELP"))
                .unwrap_or(false)
        })
        .collect();

    paths.sort();

    for path in paths {
        let key = path.file_stem().unwrap().to_string_lossy().to_uppercase();
        if let Some((n_ints, n_floats)) = file_format(&key) {
            let entries =
                parse_file(&path, n_ints, n_floats).with_context(|| format!("Parsing {key}"))?;
            map.insert(key, entries);
        } else {
            println!("cargo:warning=Skipping {key}: unknown format");
        }
    }
    Ok(map)
}

/// ----------------------  RUST CODE GENERATOR  -------------------

/// Provenance fields embedded as a comment header in the generated file.
#[derive(Debug, Clone, Default)]
struct ElpProvenance {
    source_url: String,
    sha256_hex: String,
    retrieved_at: String,
    byte_count: u64,
}

fn generate_rust(data: &BTreeMap<String, Vec<Entry>>, prov: Option<&ElpProvenance>) -> Result<String> {
    let mut out = String::new();

    writeln!(out, "// ───────────────────────────────────────────────────────────────────")?;
    writeln!(out, "// **AUTOGENERATED** by build.rs – DO NOT EDIT BY HAND")?;
    if let Some(p) = prov {
        writeln!(out, "// Source URL  : {}", p.source_url)?;
        writeln!(out, "// Source bytes: {}", p.byte_count)?;
        writeln!(out, "// Source SHA256: {}", p.sha256_hex)?;
        writeln!(out, "// Retrieved at : {}", p.retrieved_at)?;
        writeln!(out, "// Generator    : siderust scripts/elp2000/mod.rs")?;
    }
    writeln!(out, "// ───────────────────────────────────────────────────────────────────")?;
    writeln!(out)?;
    writeln!(out, "use crate::calculus::elp2000::elp_structs::*;\n")?;

    for (name, entries) in data {
        let idx: u32 = name.trim_start_matches("ELP").parse().unwrap();
        let typ = match idx {
            1..=3 => "MainProblem",
            4..=9 | 22..=36 => "EarthPert",
            10..=21 => "PlanetPert",
            _ => continue,
        };

        writeln!(out, "pub static {name}: &[{typ}] = &[")?;
        for e in entries {
            match typ {
                "MainProblem" => {
                    let ilu = e
                        .ints
                        .iter()
                        .map(|v| v.to_string())
                        .collect::<Vec<_>>()
                        .join(", ");
                    let a = &e.floats[0];
                    let b = e.floats[1..].join(", ");
                    writeln!(out, "    {typ} {{ ilu: [{ilu}], a: {a}, b: [{b}] }},")?;
                }
                "EarthPert" => {
                    let iz = format!("{}{}", e.ints[0], ".0"); // "1.0", "-2.0", …
                    let ilu = e.ints[1..5]
                        .iter()
                        .map(|v| v.to_string())
                        .collect::<Vec<_>>()
                        .join(", ");
                    let (o, a, p) = (&e.floats[0], &e.floats[1], &e.floats[2]);
                    writeln!(
                        out,
                        "    {typ} {{ iz: {iz}, ilu: [{ilu}], o: {o}, a: {a}, p: {p} }},"
                    )?;
                }
                "PlanetPert" => {
                    let ipla = e
                        .ints
                        .iter()
                        .map(|v| v.to_string())
                        .collect::<Vec<_>>()
                        .join(", ");
                    let (theta, o, p) = (&e.floats[0], &e.floats[1], &e.floats[2]);
                    writeln!(
                        out,
                        "    {typ} {{ ipla: [{ipla}], theta: {theta}, o: {o}, p: {p} }},"
                    )?;
                }
                _ => unreachable!(),
            }
        }
        writeln!(out, "];\n")?;
    }
    Ok(out)
}

/// ----------------------  DATASET HANDLING  ----------------------
/// Ensure the 36 `ELPn` files are present under `dir`.
///
fn ensure_dataset(dir: &Path) -> Result<()> {
    use reqwest::blocking::Client;

    fs::create_dir_all(dir)?;

    let base = "https://cdsarc.cds.unistra.fr/ftp/VI/79/";
    let client = Client::builder()
        .user_agent("ELP build script (rust)")
        .build()?;

    for n in 1..=36 {
        let name = format!("ELP{n}");
        let path = dir.join(&name);
        if path.exists() {
            continue;
        }
        let url = format!("{base}/{name}");
        println!("cargo:info=Downloading {url}");
        let bytes = client
            .get(&url)
            .send()
            .with_context(|| format!("GET {url}"))?
            .error_for_status()?
            .bytes()?;
        fs::write(&path, &bytes).with_context(|| format!("write {path:?}"))?;
    }
    Ok(())
}

/// ---------------------------  ENTRYPOINT  -----------------------
#[allow(dead_code)]
pub fn run(data_dir: &Path) -> Result<()> {
    ensure_dataset(data_dir)?;

    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed={}", data_dir.display());

    let parsed = parse_all_elps(data_dir)?;
    let code = generate_rust(&parsed, None)?;

    let out_dir = PathBuf::from(env::var("OUT_DIR")?);
    fs::write(out_dir.join("elp_data.rs"), code.as_bytes())?;
    println!("cargo:info=elp_data.rs generated");

    Ok(())
}

/// Like [`run`] but writes `elp_data.rs` to `gen_dir` instead of `OUT_DIR`.
///
/// Used by `build.rs` when `SIDERUST_REGEN=1` to overwrite the committed
/// table in `src/generated/`. Verifies the input dataset SHA-256 against
/// the `SIDERUST_ELP2000_SHA256` environment variable when set, and
/// embeds the SHA-256 + retrieval timestamp as a comment header in the
/// emitted `elp_data.rs`.
pub fn run_regen(data_dir: &Path, gen_dir: &Path) -> Result<()> {
    ensure_dataset(data_dir)?;

    let (sha, bytes) = dataset_sha256(data_dir)?;
    eprintln!(
        "ELP2000: dataset SHA-256 = {sha} ({bytes} bytes across all files)"
    );
    check_pinned_sha256(&sha, "SIDERUST_ELP2000_SHA256")?;

    let prov = ElpProvenance {
        source_url: "https://cdsarc.cds.unistra.fr/ftp/VI/79/".to_string(),
        sha256_hex: sha,
        retrieved_at: iso8601_now(),
        byte_count: bytes,
    };

    let parsed = parse_all_elps(data_dir)?;
    let code = generate_rust(&parsed, Some(&prov))?;
    fs::create_dir_all(gen_dir)?;
    fs::write(gen_dir.join("elp_data.rs"), code.as_bytes())?;
    Ok(())
}

fn iso8601_now() -> String {
    use std::time::SystemTime;
    let secs = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = (secs / 86_400) as i64;
    let rem = (secs % 86_400) as u32;
    let hour = (rem / 3_600) as u8;
    let minute = ((rem % 3_600) / 60) as u8;
    let second = (rem % 60) as u8;
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u8;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u8;
    let year = (y + if m <= 2 { 1 } else { 0 }) as i32;
    format!(
        "{year:04}-{m:02}-{d:02}T{hour:02}:{minute:02}:{second:02}Z"
    )
}
