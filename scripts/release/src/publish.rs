//! crates.io lookup and the ordered workspace publication workflow.

use std::env;
use std::path::Path;
use std::process::Command;
use std::time::Duration;

use serde::Deserialize;

use crate::validate::{preflight, Metadata, MetadataPackage};
use crate::ValidationError;

/// Return whether the exact version is on crates.io, including yanked versions.
pub fn crates_io_has_version(name: &str, version: &str) -> Result<bool, ValidationError> {
    if name.is_empty()
        || !name
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_'))
        || semver::Version::parse(version).is_err()
    {
        return Err(ValidationError(
            "expected a crate name and a SemVer version".into(),
        ));
    }
    let agent: ureq::Agent = ureq::Agent::config_builder()
        .https_only(true)
        .timeout_global(Some(Duration::from_secs(20)))
        .max_redirects(0)
        .http_status_as_error(false)
        .build()
        .into();
    let url = format!("https://crates.io/api/v1/crates/{name}/{version}");
    let result = agent
        .get(&url)
        .header(
            "User-Agent",
            "siderust-release (https://github.com/Siderust/siderust)",
        )
        .header("Accept", "application/json")
        .call()
        .map_err(|error| error.to_string())
        .and_then(|mut response| {
            let status = response.status().as_u16();
            if !matches!(status, 200 | 404) {
                return Err(format!("unexpected HTTP {status}"));
            }
            let bytes = response
                .body_mut()
                .with_config()
                .limit(64 * 1024)
                .read_to_vec()
                .map_err(|error| format!("cannot read HTTP {status} response: {error}"))?;
            let body = serde_json::from_slice(&bytes)
                .map_err(|error| format!("invalid HTTP {status} response body: {error}"))?;
            version_status(status, body, name, version)
        });
    result.map_err(|message| {
        ValidationError(format!(
            "cannot verify crates.io status for {name} {version}: {message}"
        ))
    })
}

#[derive(Deserialize)]
#[serde(untagged)]
enum VersionResponse {
    Found { version: RegistryVersion },
    Error { errors: Vec<RegistryError> },
}

#[derive(Deserialize)]
struct RegistryVersion {
    #[serde(rename = "crate")]
    name: String,
    num: String,
}

#[derive(Deserialize)]
struct RegistryError {
    detail: String,
}

fn version_status(
    status: u16,
    body: VersionResponse,
    name: &str,
    version: &str,
) -> Result<bool, String> {
    match (status, body) {
        (200, VersionResponse::Found { version: found })
            if found.name == name && found.num == version =>
        {
            Ok(true)
        }
        (404, VersionResponse::Error { errors }) if errors.len() == 1 => {
            let detail = &errors[0].detail;
            if *detail == format!("crate `{name}` does not exist")
                || *detail == format!("crate `{name}` does not have a version `{version}`")
            {
                return Ok(false);
            }
            Err("unrecognized crates.io HTTP 404 response".into())
        }
        _ => Err(format!("unexpected crates.io HTTP {status} response")),
    }
}

fn publish_candidates(
    metadata: &Metadata,
    confirm_ffi: bool,
) -> Result<Vec<&MetadataPackage>, ValidationError> {
    let mut candidates = Vec::new();
    for name in ["siderust", "siderust-ffi"] {
        if name == "siderust-ffi" && !confirm_ffi {
            continue;
        }
        let package = metadata
            .packages
            .iter()
            .find(|package| {
                package.name == name && metadata.workspace_members.contains(&package.id)
            })
            .ok_or_else(|| ValidationError(format!("{name} is missing from the workspace")))?;
        if package.publishable() {
            candidates.push(package);
        } else if name == "siderust-ffi" {
            return Err(ValidationError(
                "--confirm-ffi requires siderust-ffi to permit crates.io publication in Cargo.toml"
                    .into(),
            ));
        } else {
            println!("skip {name}: manifest restricts crates.io publication");
        }
    }
    Ok(candidates)
}

fn publish_missing(
    candidates: &[&MetadataPackage],
    dry_run: bool,
    lookup: impl Fn(&str, &str) -> Result<bool, ValidationError>,
    mut cargo: impl FnMut(&MetadataPackage, bool) -> Result<(), ValidationError>,
) -> Result<(), ValidationError> {
    let mut pending = Vec::new();
    for package in candidates {
        if lookup(&package.name, &package.version)? {
            println!(
                "skip {} {}: already on crates.io",
                package.name, package.version
            );
        } else {
            pending.push(*package);
        }
    }
    for package in pending {
        cargo(package, dry_run)?;
    }
    Ok(())
}

fn cargo_publish_command(root: &Path, package: &MetadataPackage, dry_run: bool) -> Command {
    let mut command = Command::new("cargo");
    command
        .current_dir(root)
        .args(["publish", "--registry", "crates-io", "--manifest-path"])
        .arg(&package.manifest_path)
        .arg("--allow-dirty");
    if dry_run {
        command.arg("--dry-run").env_remove("CARGO_REGISTRY_TOKEN");
    }
    command
}

fn cargo_publish(
    root: &Path,
    package: &MetadataPackage,
    dry_run: bool,
) -> Result<(), ValidationError> {
    let name = &package.name;
    let version = &package.version;
    println!(
        "{} {name} {version}",
        if dry_run { "dry-run" } else { "publish" }
    );
    let mut command = cargo_publish_command(root, package, dry_run);
    let status = command.status().map_err(|err| {
        ValidationError(format!("could not run Cargo for {name} {version}: {err}"))
    })?;
    if !status.success() {
        return Err(ValidationError(format!(
            "cargo publish for {name} {version} failed ({status}); rerun to recheck crates.io"
        )));
    }
    Ok(())
}

/// Validate the release, check all selected versions, then publish missing packages.
pub fn publish(root: &Path, dry_run: bool, confirm_ffi: bool) -> Result<(), ValidationError> {
    let (_, _, metadata) = preflight(root, None, None, true)?;
    let candidates = publish_candidates(&metadata, confirm_ffi)?;
    if !dry_run && env::var_os("CARGO_REGISTRY_TOKEN").is_none_or(|token| token.is_empty()) {
        return Err(ValidationError(
            "CARGO_REGISTRY_TOKEN must be set for live publication".into(),
        ));
    }
    publish_missing(
        &candidates,
        dry_run,
        crates_io_has_version,
        |package, dry_run| cargo_publish(root, package, dry_run),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use std::cell::RefCell;

    fn workspace(ffi_publish: serde_json::Value) -> Metadata {
        let packages = [
            ("siderust", json!(null)),
            ("siderust-ffi", ffi_publish),
            ("siderust-release", json!([])),
            ("unreviewed", json!(null)),
        ]
        .map(|(name, publish)| {
            json!({
                "name": name, "id": name, "version": "1.2.3",
                "manifest_path": format!("{name}/Cargo.toml"),
                "publish": publish, "dependencies": []
            })
        });
        serde_json::from_value(json!({
            "packages": packages,
            "workspace_members": ["siderust", "siderust-ffi", "siderust-release", "unreviewed"]
        }))
        .unwrap()
    }

    #[test]
    fn only_confirmed_ffi_can_follow_primary() {
        for restriction in [json!([]), json!(false), json!(["internal"])] {
            let metadata = workspace(restriction);
            let candidates = publish_candidates(&metadata, false).unwrap();
            assert_eq!(candidates.len(), 1);
            assert_eq!(candidates[0].name, "siderust");
            assert!(publish_candidates(&metadata, true).is_err());
        }
        let metadata = workspace(json!(["crates-io"]));
        let candidates = publish_candidates(&metadata, true).unwrap();
        assert_eq!(candidates.len(), 2);
        assert_eq!(candidates[1].name, "siderust-ffi");
        assert!(serde_json::from_value::<Metadata>(json!({"packages": []})).is_err());
        let malformed = json!({
            "name": "broken", "version": "1.2.3", "id": "broken",
            "publish": "crates-io", "manifest_path": "Cargo.toml", "dependencies": []
        });
        assert!(serde_json::from_value::<MetadataPackage>(malformed).is_err());
    }

    #[test]
    fn verifies_status_and_response_identity() {
        let found = || VersionResponse::Found {
            version: RegistryVersion {
                name: "siderust".into(),
                num: "1.2.3".into(),
            },
        };
        let missing = || VersionResponse::Error {
            errors: vec![RegistryError {
                detail: "crate `siderust` does not have a version `1.2.3`".into(),
            }],
        };
        assert_eq!(version_status(200, found(), "siderust", "1.2.3"), Ok(true));
        assert_eq!(
            version_status(404, missing(), "siderust", "1.2.3"),
            Ok(false)
        );
        assert!(version_status(200, found(), "siderust", "1.2.4").is_err());
        assert!(version_status(404, missing(), "siderust", "1.2.4").is_err());
        for status in [301, 401, 403, 429, 500, 503] {
            assert!(version_status(status, missing(), "siderust", "1.2.3").is_err());
        }
    }

    #[test]
    fn cargo_dry_run_has_no_publication_token() {
        let metadata = workspace(json!(false));
        let candidates = publish_candidates(&metadata, false).unwrap();
        let dry_run = cargo_publish_command(Path::new("."), candidates[0], true);
        assert!(dry_run.get_args().any(|arg| arg == "--dry-run"));
        assert!(dry_run
            .get_envs()
            .any(|(key, value)| key == "CARGO_REGISTRY_TOKEN" && value.is_none()));
        let live = cargo_publish_command(Path::new("."), candidates[0], false);
        assert!(!live.get_args().any(|arg| arg == "--dry-run"));
    }

    #[test]
    fn dry_run_and_rerun_skip_existing_versions() {
        let metadata = workspace(json!(true));
        let candidates = publish_candidates(&metadata, true).unwrap();
        let published = RefCell::new(vec!["siderust".to_owned()]);
        for (dry_run, expected_calls) in [(true, 1), (false, 1), (false, 0)] {
            let mut calls = 0;
            publish_missing(
                &candidates,
                dry_run,
                |name, _| Ok(published.borrow().iter().any(|item| item == name)),
                |package, dry_run| {
                    assert_eq!(package.name, "siderust-ffi");
                    assert_eq!(calls, 0);
                    calls += 1;
                    if !dry_run {
                        published.borrow_mut().push(package.name.clone());
                    }
                    Ok(())
                },
            )
            .unwrap();
            assert_eq!(calls, expected_calls);
        }
        assert_eq!(published.into_inner(), ["siderust", "siderust-ffi"]);
    }

    #[test]
    fn failures_stop_before_any_further_publication() {
        let metadata = workspace(json!(true));
        let candidates = publish_candidates(&metadata, true).unwrap();
        let error = publish_missing(
            &candidates,
            false,
            |name, _| {
                if name == "siderust-ffi" {
                    Err(ValidationError("lookup timed out".into()))
                } else {
                    Ok(false)
                }
            },
            |_, _| panic!("must verify all versions before publishing"),
        )
        .unwrap_err();
        assert_eq!(error.0, "lookup timed out");

        let error = publish_missing(
            &candidates,
            false,
            |_, _| Ok(false),
            |package, _| {
                assert_eq!(package.name, "siderust");
                Err(ValidationError("Cargo failed".into()))
            },
        )
        .unwrap_err();
        assert_eq!(error.0, "Cargo failed");
    }
}
