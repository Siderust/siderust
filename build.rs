// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Build script: downloads JPL DE440/DE441 kernels on demand.
#[cfg(any(feature = "regen-data", feature = "de440", feature = "de441"))]
use std::env;
#[cfg(any(feature = "regen-data", feature = "de440", feature = "de441"))]
use std::path::PathBuf;

#[cfg(feature = "regen-data")]
#[path = "scripts/vsop87/mod.rs"]
mod vsop87_build;

#[cfg(feature = "regen-data")]
#[path = "scripts/elp2000/mod.rs"]
mod elp2000_build;

#[cfg(feature = "de440")]
#[path = "scripts/jpl/de440/mod.rs"]
mod de440_build;

#[cfg(feature = "de441")]
#[path = "scripts/jpl/de441/mod.rs"]
mod de441_build;

#[cfg(any(feature = "de440", feature = "de441"))]
#[path = "scripts/jpl/daf.rs"]
pub(crate) mod jpl_daf;

#[cfg(any(feature = "de440", feature = "de441"))]
#[path = "scripts/jpl/pipeline.rs"]
pub(crate) mod jpl_pipeline;

#[cfg(any(feature = "de440", feature = "de441"))]
#[path = "scripts/jpl/spk.rs"]
pub(crate) mod jpl_spk;

fn main() {
    println!("cargo:rerun-if-env-changed=SIDERUST_DATASETS_DIR");
    println!("cargo:rerun-if-env-changed=SIDERUST_REGEN");
    println!("cargo:rerun-if-env-changed=SIDERUST_JPL_STUB");
    println!("cargo:rustc-check-cfg=cfg(siderust_mock_de440)");
    println!("cargo:rustc-check-cfg=cfg(siderust_mock_de441)");
    println!("cargo:rustc-check-cfg=cfg(siderust_archive_present)");

    #[cfg(feature = "regen-data")]
    regen_tables();

    #[cfg(not(feature = "regen-data"))]
    eprintln!(
        "siderust build: using committed generated VSOP87/ELP2000 tables in src/embedded_data/. \
         Enable the `regen-data` feature and set SIDERUST_REGEN=1 to refresh them."
    );

    #[cfg(feature = "de440")]
    build_de440();

    #[cfg(feature = "de441")]
    build_de441();

    archive::emit();
}

mod archive {
    use std::env;
    use std::fs;
    use std::path::PathBuf;

    /// Always called from `build.rs::main`.
    ///
    /// * Writes an `archive_registry.rs` stub into `OUT_DIR` so the runtime
    ///   layer always compiles, even when the submodule is absent.
    /// * If `feature = "archive-data"` is enabled and `archive/MANIFEST.toml`
    ///   is present, parses the manifest and populates the registry with
    ///   concrete entries pulled from each family `manifest.toml`.
    /// * Emits `cargo:rustc-cfg=siderust_archive_present` when the submodule
    ///   is checked out, so downstream code can opportunistically rely on
    ///   archive data without requiring the feature.
    pub(super) fn emit() {
        let out_dir = PathBuf::from(env::var_os("OUT_DIR").expect("OUT_DIR must be set"));
        let manifest_dir = PathBuf::from(
            env::var_os("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR must be set"),
        );

        // The archive is no longer a git submodule. Consumers that want
        // archive-backed builds either:
        //   * point `SIDERUST_ARCHIVE_ROOT` at a local checkout of
        //     https://github.com/Siderust/archive, or
        //   * keep a sibling `../archive` directory, or
        //   * place an `archive/` directory next to this crate.
        // The fallback chain below tries each in order.
        println!("cargo:rerun-if-env-changed=SIDERUST_ARCHIVE_ROOT");
        let archive_root = env::var_os("SIDERUST_ARCHIVE_ROOT")
            .map(PathBuf::from)
            .filter(|p| p.join("MANIFEST.toml").is_file())
            .or_else(|| {
                let candidates = [
                    manifest_dir.join("archive"),
                    manifest_dir.join("../archive"),
                ];
                candidates
                    .into_iter()
                    .find(|p| p.join("MANIFEST.toml").is_file())
            })
            .unwrap_or_else(|| manifest_dir.join("archive"));
        let top_manifest = archive_root.join("MANIFEST.toml");

        println!("cargo:rerun-if-changed={}", top_manifest.display());

        let archive_present = top_manifest.is_file();
        if archive_present {
            println!("cargo:rustc-cfg=siderust_archive_present");
        }

        let registry_path = out_dir.join("archive_registry.rs");
        let registry = if cfg!(feature = "archive-data") && archive_present {
            match build_registry(&archive_root) {
                Ok(rendered) => rendered,
                Err(err) => {
                    println!(
                        "cargo:warning=archive-data: failed to parse archive manifests: {err}"
                    );
                    empty_registry()
                }
            }
        } else if cfg!(feature = "archive-data") && !archive_present {
            println!(
                "cargo:warning=archive-data feature is enabled but no archive checkout was found. \
                 Set SIDERUST_ARCHIVE_ROOT to a clone of https://github.com/Siderust/archive, \
                 or place an `archive/` directory next to this crate."
            );
            empty_registry()
        } else {
            empty_registry()
        };

        if let Err(err) = fs::write(&registry_path, registry) {
            panic!(
                "failed to write archive_registry.rs to {}: {err}",
                registry_path.display()
            );
        }

        emit_lagrange_byte_paths(&out_dir, &archive_root, archive_present);
    }

    /// Emits `OUT_DIR/lagrange_paths.rs` declaring `L1_BYTES..L5_BYTES`
    /// `&'static [u8]` constants. When the `lagrange-centers` feature is
    /// enabled and the archive checkout has the required `.sck` files, the
    /// constants are wired up with `include_bytes!`. Otherwise a build-time
    /// `compile_error!` is emitted with an actionable diagnostic so users
    /// can't accidentally enable the feature without supplying the data.
    fn emit_lagrange_byte_paths(
        out_dir: &std::path::Path,
        archive_root: &std::path::Path,
        archive_present: bool,
    ) {
        let path = out_dir.join("lagrange_paths.rs");
        let lagrange_dir = archive_root.join("lagrange").join("vsop87");
        let required = ["l1.sck", "l2.sck", "l3.sck", "l4.sck", "l5.sck"];
        let all_present = archive_present
            && required
                .iter()
                .all(|name| lagrange_dir.join(name).is_file());

        let body = if cfg!(feature = "lagrange-centers") {
            if all_present {
                let mut s = String::from(
                    "/// Auto-generated by build.rs. Embeds Lagrange SCK kernels \
                     from the archive checkout.\n",
                );
                for name in required {
                    let abs = lagrange_dir.join(name);
                    println!("cargo:rerun-if-changed={}", abs.display());
                    let upper = name.split('.').next().unwrap().to_uppercase();
                    s.push_str(&format!(
                        "static {upper}_BYTES: &[u8] = include_bytes!(\"{}\");\n",
                        abs.display()
                    ));
                }
                s
            } else {
                let dir_str = lagrange_dir.display().to_string().replace('\\', "\\\\");
                format!(
                    "compile_error!(\"lagrange-centers feature requires the Siderust \
                     Archive at {dir_str} (missing one of l1.sck..l5.sck). Clone \
                     https://github.com/Siderust/archive and either: (a) set \
                     SIDERUST_ARCHIVE_ROOT to the checkout, (b) place an `archive/` \
                     directory next to this crate, or (c) place a sibling \
                     `../archive/` directory.\");\n"
                )
            }
        } else {
            String::from(
                "/// Auto-generated by build.rs. `lagrange-centers` feature is disabled.\n",
            )
        };

        if let Err(err) = fs::write(&path, body) {
            panic!(
                "failed to write lagrange_paths.rs to {}: {err}",
                path.display()
            );
        }
    }

    fn empty_registry() -> String {
        let mut s = String::new();
        s.push_str("/// Auto-generated by build.rs (no archive entries).\n");
        s.push_str("pub const ARCHIVE_ENTRIES: &[ArchiveEntry] = &[];\n");
        s
    }

    #[cfg(feature = "archive-data")]
    fn build_registry(archive_root: &std::path::Path) -> Result<String, String> {
        use std::fmt::Write;

        let top_text = fs::read_to_string(archive_root.join("MANIFEST.toml"))
            .map_err(|e| format!("reading MANIFEST.toml: {e}"))?;
        let top: toml::Value =
            toml::from_str(&top_text).map_err(|e| format!("parsing MANIFEST.toml: {e}"))?;

        let families = top
            .get("family")
            .and_then(|v| v.as_array())
            .cloned()
            .unwrap_or_default();

        let mut entries: Vec<RegistryEntry> = Vec::new();
        for fam in families {
            let id = fam
                .get("id")
                .and_then(|v| v.as_str())
                .ok_or("family without id")?
                .to_string();
            let manifest_rel = fam
                .get("manifest")
                .and_then(|v| v.as_str())
                .ok_or_else(|| format!("family `{id}` without manifest path"))?
                .to_string();
            let family_dir = std::path::Path::new(&manifest_rel)
                .parent()
                .map(|p| p.to_string_lossy().into_owned())
                .unwrap_or_default();
            let kind = fam
                .get("kind")
                .and_then(|v| v.as_str())
                .unwrap_or("unknown")
                .to_string();

            let fam_manifest = archive_root.join(&manifest_rel);
            println!("cargo:rerun-if-changed={}", fam_manifest.display());
            if !fam_manifest.is_file() {
                entries.push(RegistryEntry {
                    family: id.clone(),
                    dataset_id: id.clone(),
                    kind: kind.clone(),
                    relative_path: family_dir,
                    manifest_path: manifest_rel,
                    valid_from_jd: f64::NAN,
                    valid_to_jd: f64::NAN,
                    checksum_sha256: String::new(),
                });
                continue;
            }

            let fam_text = fs::read_to_string(&fam_manifest)
                .map_err(|e| format!("reading {}: {e}", fam_manifest.display()))?;
            let fam_value: toml::Value = toml::from_str(&fam_text)
                .map_err(|e| format!("parsing {}: {e}", fam_manifest.display()))?;

            let dataset_id = fam_value
                .get("dataset_id")
                .and_then(|v| v.as_str())
                .unwrap_or(&id)
                .to_string();
            let valid_from = fam_value
                .get("valid_from_jd")
                .and_then(|v| v.as_float().or_else(|| v.as_integer().map(|i| i as f64)))
                .unwrap_or(f64::NAN);
            let valid_to = fam_value
                .get("valid_to_jd")
                .and_then(|v| v.as_float().or_else(|| v.as_integer().map(|i| i as f64)))
                .unwrap_or(f64::NAN);

            entries.push(RegistryEntry {
                family: id,
                dataset_id,
                kind,
                relative_path: family_dir,
                manifest_path: manifest_rel,
                valid_from_jd: valid_from,
                valid_to_jd: valid_to,
                checksum_sha256: String::new(),
            });
        }

        let mut s = String::new();
        s.push_str("/// Auto-generated by build.rs from archive/MANIFEST.toml\n");
        let _ = writeln!(s, "pub const ARCHIVE_ENTRIES: &[ArchiveEntry] = &[");
        for e in &entries {
            let _ = writeln!(
                s,
                "    ArchiveEntry {{ family: {family:?}, dataset_id: {ds:?}, kind: {kind:?}, \
                 relative_path: {path:?}, manifest_path: {mp:?}, valid_from_jd: {vf}_f64, \
                 valid_to_jd: {vt}_f64, checksum_sha256: {ck:?} }},",
                family = e.family,
                ds = e.dataset_id,
                kind = e.kind,
                path = e.relative_path,
                mp = e.manifest_path,
                vf = if e.valid_from_jd.is_nan() {
                    "f64::NAN".to_string()
                } else {
                    format!("{:?}", e.valid_from_jd)
                },
                vt = if e.valid_to_jd.is_nan() {
                    "f64::NAN".to_string()
                } else {
                    format!("{:?}", e.valid_to_jd)
                },
                ck = e.checksum_sha256,
            );
        }
        s.push_str("];\n");
        Ok(s)
    }

    #[cfg(not(feature = "archive-data"))]
    fn build_registry(_: &std::path::Path) -> Result<String, String> {
        Ok(empty_registry())
    }

    #[cfg(feature = "archive-data")]
    struct RegistryEntry {
        family: String,
        dataset_id: String,
        kind: String,
        relative_path: String,
        manifest_path: String,
        valid_from_jd: f64,
        valid_to_jd: f64,
        checksum_sha256: String,
    }
}

// ── Shared helpers ────────────────────────────────────────────────────────────

/// Returns the base directory for dataset storage.
/// Prefers `SIDERUST_DATASETS_DIR`; falls back to `OUT_DIR`.
#[cfg(any(feature = "regen-data", feature = "de440", feature = "de441"))]
fn datasets_base_dir() -> PathBuf {
    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR not set by Cargo"));
    env::var_os("SIDERUST_DATASETS_DIR")
        .map(PathBuf::from)
        .unwrap_or(out_dir)
}

// ── `regen-data` feature ──────────────────────────────────────────────────────

/// Regenerates the committed VSOP87 and ELP2000 tables from source data.
///
/// Only runs when `SIDERUST_REGEN=1` (or `true`/`yes`) is set. Requires the
/// `regen-data` build feature so that `reqwest` is compiled only when needed.
///
/// The preferred way to invoke this is via the helper script:
///   scripts/update_generated_tables.sh
#[cfg(feature = "regen-data")]
fn regen_tables() {
    let regen = env::var("SIDERUST_REGEN")
        .map(|v| {
            let v = v.trim().to_ascii_lowercase();
            v == "1" || v == "true" || v == "yes"
        })
        .unwrap_or(false);

    if !regen {
        eprintln!(
            "siderust build (regen-data): SIDERUST_REGEN not set, skipping table regeneration."
        );
        return;
    }

    let base = datasets_base_dir();
    let manifest_dir =
        PathBuf::from(env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR not set"));
    let gen_dir = manifest_dir.join("src/embedded_data");

    eprintln!("Regenerating VSOP87 data...");
    vsop87_build::run_regen(base.join("vsop87_dataset").as_path(), &gen_dir)
        .unwrap_or_else(|e| panic!("VSOP87 codegen failed: {e}"));
    eprintln!("VSOP87 regeneration complete");

    eprintln!("Regenerating ELP2000 data...");
    elp2000_build::run_regen(base.join("elp2000_dataset").as_path(), &gen_dir)
        .unwrap_or_else(|e| panic!("ELP2000 codegen failed: {e}"));
    eprintln!("ELP2000 regeneration complete");
}

// ── `de440` feature ───────────────────────────────────────────────────────────

/// Builds DE440 coefficient data from the JPL BSP file.
///
/// When `SIDERUST_JPL_STUB=all` (or `de440`), emits `siderust_mock_de440` so
/// the library falls back to `Vsop87Ephemeris` without downloading the BSP.
#[cfg(feature = "de440")]
fn build_de440() {
    if jpl_stub_active_for("de440") {
        println!("cargo:rustc-cfg=siderust_mock_de440");
    }
    let de440_dir = datasets_base_dir().join("de440_dataset");
    eprintln!("Building DE440 data...");
    de440_build::run(de440_dir.as_path()).unwrap_or_else(|e| panic!("DE440 codegen failed: {e}"));
    eprintln!("DE440 data generation complete");
}

// ── `de441` feature ───────────────────────────────────────────────────────────

/// Builds DE441 coefficient data from the JPL BSP file.
///
/// When `SIDERUST_JPL_STUB=all` (or `de441`), emits `siderust_mock_de441` so
/// the library falls back to `Vsop87Ephemeris` without downloading the BSP.
#[cfg(feature = "de441")]
fn build_de441() {
    if jpl_stub_active_for("de441") {
        println!("cargo:rustc-cfg=siderust_mock_de441");
    }
    let de441_dir = datasets_base_dir().join("de441_dataset");
    eprintln!("Building DE441 data...");
    de441_build::run(de441_dir.as_path()).unwrap_or_else(|e| panic!("DE441 codegen failed: {e}"));
    eprintln!("DE441 data generation complete");
}

// ── Shared JPL helpers ────────────────────────────────────────────────────────

#[cfg(any(feature = "de440", feature = "de441"))]
fn jpl_stub_active_for(dataset: &str) -> bool {
    let Ok(raw) = env::var("SIDERUST_JPL_STUB") else {
        return false;
    };
    let lower = raw.trim().to_ascii_lowercase();
    if lower.is_empty() {
        return false;
    }
    if matches!(lower.as_str(), "all" | "1" | "true" | "yes" | "on") {
        return true;
    }
    lower
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|s| !s.is_empty())
        .any(|tok| tok == dataset)
}
