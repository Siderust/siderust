//! Validate tag, Cargo, and changelog consistency before publication.

use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use serde::Deserialize;

use crate::changelog::{extract_version_section, ChangelogError};

/// Release metadata is not ready for publication.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidationError(pub String);

impl fmt::Display for ValidationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ValidationError {}

impl From<ChangelogError> for ValidationError {
    fn from(value: ChangelogError) -> Self {
        Self(value.0)
    }
}

impl From<std::io::Error> for ValidationError {
    fn from(value: std::io::Error) -> Self {
        Self(value.to_string())
    }
}

impl From<serde_json::Error> for ValidationError {
    fn from(value: serde_json::Error) -> Self {
        Self(value.to_string())
    }
}

/// Map `vX.Y.Z` to `X.Y.Z`.
pub fn parse_tag(tag: &str) -> Result<String, ValidationError> {
    let tag = tag.trim();
    let Some(version) = tag.strip_prefix('v') else {
        return Err(ValidationError(format!(
            "malformed release tag {tag:?}; expected vX.Y.Z"
        )));
    };
    if !is_release_version(version) {
        return Err(ValidationError(format!(
            "malformed release tag {tag:?}; expected vX.Y.Z (optionally with pre-release)"
        )));
    }
    Ok(version.to_owned())
}

fn is_release_version(value: &str) -> bool {
    let mut parts = value.splitn(2, ['-', '+']);
    let core = parts.next().unwrap_or("");
    let mut nums = core.split('.');
    let (Some(major), Some(minor), Some(patch)) = (nums.next(), nums.next(), nums.next()) else {
        return false;
    };
    if nums.next().is_some() {
        return false;
    }
    !major.is_empty()
        && major.chars().all(|c| c.is_ascii_digit())
        && !minor.is_empty()
        && minor.chars().all(|c| c.is_ascii_digit())
        && !patch.is_empty()
        && patch.chars().all(|c| c.is_ascii_digit())
}

#[derive(Debug, Deserialize)]
struct Metadata {
    packages: Vec<MetadataPackage>,
    workspace_members: Vec<String>,
}

#[derive(Debug, Deserialize)]
struct MetadataPackage {
    name: String,
    version: String,
    id: String,
    publish: Option<PublishField>,
    dependencies: Vec<MetadataDependency>,
    manifest_path: String,
}

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum PublishField {
    Flag(bool),
    Registries(Vec<String>),
}

#[derive(Debug, Deserialize)]
struct MetadataDependency {
    name: String,
    source: Option<String>,
    path: Option<String>,
    kind: Option<String>,
}

/// Name / version / publishability for one manifest.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackageMeta {
    pub name: String,
    pub version: String,
    pub publishable: bool,
}

fn is_publishable(publish: &Option<PublishField>) -> bool {
    match publish {
        None => true,
        Some(PublishField::Flag(flag)) => *flag,
        Some(PublishField::Registries(regs)) => !regs.is_empty(),
    }
}

fn cargo_metadata_for_manifest(manifest: &Path) -> Result<Metadata, ValidationError> {
    let output = Command::new("cargo")
        .args([
            "metadata",
            "--format-version",
            "1",
            "--no-deps",
            "--manifest-path",
        ])
        .arg(manifest)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()?;
    if !output.status.success() {
        return Err(ValidationError(format!(
            "cargo metadata failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )));
    }
    Ok(serde_json::from_slice(&output.stdout)?)
}

/// Read package identity from a Cargo manifest via `cargo metadata`.
pub fn package_meta(manifest: &Path) -> Result<PackageMeta, ValidationError> {
    let metadata = cargo_metadata_for_manifest(manifest)?;
    let manifest_canon = fs::canonicalize(manifest).unwrap_or_else(|_| manifest.to_path_buf());
    let package = metadata
        .packages
        .into_iter()
        .find(|package| {
            fs::canonicalize(&package.manifest_path)
                .map(|path| path == manifest_canon)
                .unwrap_or_else(|_| Path::new(&package.manifest_path) == manifest)
        })
        .ok_or_else(|| {
            ValidationError(format!(
                "no package found for manifest {}",
                manifest.display()
            ))
        })?;
    Ok(PackageMeta {
        name: package.name,
        version: package.version,
        publishable: is_publishable(&package.publish),
    })
}

fn root_crate_version(root: &Path) -> Result<String, ValidationError> {
    Ok(package_meta(&root.join("Cargo.toml"))?.version)
}

fn validate_tag_version(tag: Option<&str>, version: &str) -> Result<(), ValidationError> {
    let Some(tag) = tag else {
        return Ok(());
    };
    let tagged = parse_tag(tag)?;
    if tagged != version {
        return Err(ValidationError(format!(
            "tag {tag} does not match crate version {version}"
        )));
    }
    Ok(())
}

fn dependency_is_blocked(dep: &MetadataDependency) -> Option<String> {
    if dep.path.is_some() && dep.source.is_none() {
        return Some("path dependency without crates.io source".to_owned());
    }
    if let Some(source) = &dep.source {
        if source.starts_with("git+") {
            return Some(format!("git dependency ({source})"));
        }
    }
    None
}

fn validate_workspace(root: &Path, expected_version: &str) -> Result<Vec<String>, ValidationError> {
    let mut notes = Vec::new();
    let root_version = root_crate_version(root)?;
    if root_version != expected_version {
        return Err(ValidationError(format!(
            "root crate version {root_version} != expected {expected_version}"
        )));
    }

    let ffi_manifest = root.join("siderust-ffi").join("Cargo.toml");
    if ffi_manifest.exists() {
        let ffi = package_meta(&ffi_manifest)?;
        if ffi.version != expected_version {
            return Err(ValidationError(format!(
                "siderust-ffi version {} != expected {expected_version}",
                ffi.version
            )));
        }
        if ffi.publishable {
            notes.push("siderust-ffi is publishable".to_owned());
        } else {
            notes.push("siderust-ffi remains publish = false".to_owned());
        }
    }

    let root_toml = fs::read_to_string(root.join("Cargo.toml"))?;
    if root_toml.lines().any(|line| {
        let trimmed = line.trim();
        trimmed == "[patch.crates-io]" || trimmed.starts_with("[patch.")
    }) {
        return Err(ValidationError(
            "Cargo.toml [patch] is present; remove patches before publishing".to_owned(),
        ));
    }

    let metadata = cargo_metadata_for_manifest(&root.join("Cargo.toml"))?;
    let members: std::collections::BTreeSet<_> =
        metadata.workspace_members.iter().cloned().collect();
    for package in metadata.packages {
        if !members.contains(&package.id) {
            continue;
        }
        if matches!(package.name.as_str(), "siderust" | "siderust-ffi")
            && package.version != expected_version
        {
            return Err(ValidationError(format!(
                "{} metadata version {} != expected {expected_version}",
                package.name, package.version
            )));
        }
        if !is_publishable(&package.publish) || package.name != "siderust" {
            continue;
        }
        for dep in &package.dependencies {
            if matches!(dep.kind.as_deref(), Some("dev") | Some("build")) {
                continue;
            }
            if let Some(reason) = dependency_is_blocked(dep) {
                return Err(ValidationError(format!(
                    "publishable crate {} has blocked dependency {}: {reason}",
                    package.name, dep.name
                )));
            }
        }
    }
    notes.push(format!(
        "workspace packages consistent at {expected_version}"
    ));
    Ok(notes)
}

fn validate_changelog(root: &Path, version: &str) -> Result<String, ValidationError> {
    let text = fs::read_to_string(root.join("CHANGELOG.md"))?;
    let section = extract_version_section(&text, version)?;
    Ok(section.body)
}

/// Return true when crates.io already hosts `name/version`.
pub fn crates_io_has_version(name: &str, version: &str) -> Result<bool, ValidationError> {
    let url = format!("https://crates.io/api/v1/crates/{name}/{version}");
    let output = Command::new("curl")
        .args([
            "-sS",
            "-o",
            "/dev/null",
            "-w",
            "%{http_code}",
            "-A",
            "siderust-release-validate (https://github.com/Siderust/siderust)",
            &url,
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()?;
    if !output.status.success() {
        return Err(ValidationError(format!(
            "crates.io lookup for {name} {version} failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )));
    }
    let code = String::from_utf8_lossy(&output.stdout);
    match code.trim() {
        "200" => Ok(true),
        "404" => Ok(false),
        other => Err(ValidationError(format!(
            "crates.io lookup for {name} {version} failed: HTTP {other}"
        ))),
    }
}

fn resolve_version(
    root: &Path,
    tag: Option<&str>,
    version: Option<&str>,
) -> Result<String, ValidationError> {
    match (tag, version) {
        (Some(tag), Some(version)) => {
            let tagged = parse_tag(tag)?;
            if tagged != version {
                return Err(ValidationError(format!(
                    "--version {version} does not match --tag {tag}"
                )));
            }
            Ok(version.to_owned())
        }
        (Some(tag), None) => parse_tag(tag),
        (None, Some(version)) => {
            let version = version
                .strip_prefix('v')
                .or_else(|| version.strip_prefix('V'))
                .unwrap_or(version);
            Ok(version.to_owned())
        }
        (None, None) => root_crate_version(root),
    }
}

/// Run the full preflight / tag validation suite.
pub fn run_validation(
    root: &Path,
    tag: Option<&str>,
    version: Option<&str>,
    require_tag: bool,
    check_published: bool,
) -> Result<Vec<String>, ValidationError> {
    if require_tag && tag.is_none() {
        return Err(ValidationError(
            "--require-tag was set but no tag was provided".to_owned(),
        ));
    }
    let resolved = resolve_version(root, tag, version)?;
    validate_tag_version(tag, &resolved)?;
    let mut notes = validate_workspace(root, &resolved)?;
    let body = validate_changelog(root, &resolved)?;
    notes.push(format!(
        "CHANGELOG.md contains nonempty [{resolved}] ({} lines)",
        body.lines().count()
    ));
    if check_published {
        let present = crates_io_has_version("siderust", &resolved)?;
        notes.push(if present {
            format!("siderust {resolved} already on crates.io")
        } else {
            format!("siderust {resolved} not yet on crates.io")
        });
    }
    Ok(notes)
}

/// Resolve repository root from this crate's manifest directory.
pub fn default_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accepts_semver_tag() {
        assert_eq!(parse_tag("v0.15.0").unwrap(), "0.15.0");
        assert_eq!(parse_tag("v0.15.0-rc.1").unwrap(), "0.15.0-rc.1");
    }

    #[test]
    fn rejects_malformed_tags() {
        for tag in ["0.15.0", "v0.15", "vv0.15.0", "release-0.15.0"] {
            assert!(parse_tag(tag).is_err(), "{tag}");
        }
    }

    #[test]
    fn tag_must_match_version() {
        assert!(validate_tag_version(Some("v0.15.0"), "0.14.1").is_err());
        assert!(validate_tag_version(Some("v0.15.0"), "0.15.0").is_ok());
    }

    #[test]
    fn require_tag_enforced() {
        let root = default_root();
        let err = run_validation(&root, None, None, true, false).unwrap_err();
        assert!(err.0.contains("require-tag"));
    }

    #[test]
    fn repo_preflight_without_tag() {
        let root = default_root();
        let notes = run_validation(&root, None, None, false, false).expect("repo ready");
        assert!(notes.iter().any(|note| note.contains("CHANGELOG.md")));
    }

    #[test]
    fn mismatched_tag_and_version_rejected() {
        let root = default_root();
        let err = run_validation(&root, Some("v0.15.0"), Some("0.14.1"), false, false)
            .expect_err("mismatch");
        assert!(err.0.contains("does not match"));
    }
}
