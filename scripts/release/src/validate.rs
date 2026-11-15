//! Validate tag, Cargo, and changelog consistency before publication.

use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use serde::Deserialize;

use crate::changelog::{extract_version_section, normalize_version, ChangelogError};
use crate::crates_io_has_version;

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
    tag.strip_prefix('v')
        .and_then(|version| semver::Version::parse(version).ok())
        .map(|version| version.to_string())
        .ok_or_else(|| {
            ValidationError(format!(
                "malformed release tag {tag:?}; expected v followed by a SemVer version"
            ))
        })
}

#[derive(Debug, Deserialize)]
pub(crate) struct Metadata {
    pub(crate) packages: Vec<MetadataPackage>,
    pub(crate) workspace_members: Vec<String>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct MetadataPackage {
    pub(crate) name: String,
    pub(crate) version: String,
    pub(crate) id: String,
    publish: Option<PublishField>,
    dependencies: Vec<MetadataDependency>,
    pub(crate) manifest_path: PathBuf,
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

impl MetadataPackage {
    pub(crate) fn publishable(&self) -> bool {
        match &self.publish {
            None | Some(PublishField::Flag(true)) => true,
            Some(PublishField::Flag(false)) => false,
            Some(PublishField::Registries(registries)) => {
                registries.iter().any(|registry| registry == "crates-io")
            }
        }
    }
}

pub(crate) fn cargo_metadata_for_manifest(manifest: &Path) -> Result<Metadata, ValidationError> {
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
    let manifest = fs::canonicalize(manifest)?;
    let metadata = cargo_metadata_for_manifest(&manifest)?;
    let package = metadata
        .packages
        .into_iter()
        .find(|package| fs::canonicalize(&package.manifest_path).is_ok_and(|path| path == manifest))
        .ok_or_else(|| ValidationError(format!("no package found for {}", manifest.display())))?;
    let publishable = package.publishable();
    Ok(PackageMeta {
        name: package.name,
        version: package.version,
        publishable,
    })
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

fn validate_workspace(
    root: &Path,
    expected_version: &str,
    publishing: bool,
    metadata: &Metadata,
) -> Result<Vec<String>, ValidationError> {
    let mut notes = Vec::new();
    let root_toml = fs::read_to_string(root.join("Cargo.toml"))?;
    let manifest: toml::Table = toml::from_str(&root_toml)
        .map_err(|err| ValidationError(format!("invalid Cargo.toml: {err}")))?;
    let has_patch = manifest.contains_key("patch");
    if publishing && has_patch {
        return Err(ValidationError(
            "Cargo.toml [patch] is present; remove patches before publishing".to_owned(),
        ));
    }

    if has_patch {
        notes.push("Workspace checks use temporary Cargo patches".to_owned());
    }

    for package in &metadata.packages {
        if !metadata.workspace_members.contains(&package.id) {
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
        if package.name == "siderust-ffi" {
            notes.push(if package.publishable() {
                "siderust-ffi permits crates.io publication".to_owned()
            } else {
                "siderust-ffi is restricted from crates.io publication".to_owned()
            });
        }
        if !package.publishable() || package.name != "siderust" {
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

/// Resolve explicit CLI version parameters using the same rules for every command.
pub fn explicit_version(
    tag: Option<&str>,
    version: Option<&str>,
) -> Result<Option<String>, ValidationError> {
    let version = version
        .map(|value| {
            normalize_version(value)
                .map(|version| version.to_string())
                .ok_or_else(|| ValidationError(format!("invalid release version {value:?}")))
        })
        .transpose()?;
    if let Some(tag) = tag {
        let tagged = parse_tag(tag)?;
        if let Some(version) = version.filter(|version| version != &tagged) {
            return Err(ValidationError(format!(
                "--version {version} does not match --tag {tag}"
            )));
        }
        return Ok(Some(tagged));
    }
    Ok(version)
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
    let (resolved, notes, _) = preflight(root, tag, version, tag.is_some())?;
    let mut notes = notes;
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

pub(crate) fn preflight(
    root: &Path,
    tag: Option<&str>,
    version: Option<&str>,
    publishing: bool,
) -> Result<(String, Vec<String>, Metadata), ValidationError> {
    let explicit = explicit_version(tag, version)?;
    let manifest = fs::canonicalize(root.join("Cargo.toml"))?;
    let metadata = cargo_metadata_for_manifest(&manifest)?;
    let primary = metadata
        .packages
        .iter()
        .find(|package| {
            fs::canonicalize(&package.manifest_path).is_ok_and(|path| path == manifest)
                && package.name == "siderust"
                && metadata.workspace_members.contains(&package.id)
        })
        .ok_or_else(|| ValidationError("workspace root must be the siderust package".to_owned()))?;
    let resolved = explicit.unwrap_or_else(|| primary.version.clone());
    let mut notes = validate_workspace(root, &resolved, publishing, &metadata)?;
    let body = validate_changelog(root, &resolved)?;
    notes.push(format!(
        "CHANGELOG.md contains nonempty [{resolved}] ({} lines)",
        body.lines().count()
    ));
    Ok((resolved, notes, metadata))
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
        for tag in [
            "0.15.0",
            "v0.15",
            "vv0.15.0",
            "release-0.15.0",
            "v00.15.0",
            "v0.15.0-",
            "v0.15.0-rc.01",
            "v0.15.0+",
            "v0.15.0-rc..1",
            " v0.15.0",
        ] {
            assert!(parse_tag(tag).is_err(), "{tag}");
        }
    }

    #[test]
    fn tag_must_match_version() {
        let root = default_root();
        assert!(run_validation(&root, None, Some("0.0.0"), false, false).is_err());
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
    fn patch_policy_and_malformed_workspace() {
        let dir = tempfile::tempdir().unwrap();
        let root = dir.path();
        fs::create_dir(root.join("src")).unwrap();
        fs::write(root.join("src/lib.rs"), "").unwrap();
        fs::write(root.join("CHANGELOG.md"), "## [1.2.3]\n\n- notes").unwrap();
        let manifest = root.join("Cargo.toml");
        assert!(preflight(root, None, None, false).is_err());
        fs::write(
            &manifest,
            "[package]\nname = 'siderust'\nversion = '1.2.3'\n[patch.'crates-io']\n",
        )
        .unwrap();
        assert!(preflight(root, None, None, false).is_ok());
        assert!(preflight(root, Some("v1.2.3"), None, true)
            .unwrap_err()
            .0
            .contains("remove patches"));
        fs::write(manifest, "[workspace]\n").unwrap();
        assert!(preflight(root, None, None, false)
            .unwrap_err()
            .0
            .contains("workspace root"));
    }

    #[test]
    fn mismatched_tag_and_version_rejected() {
        let root = default_root();
        let err = run_validation(&root, Some("v0.15.0"), Some("0.14.1"), false, false)
            .expect_err("mismatch");
        assert!(err.0.contains("does not match"));
    }
}
