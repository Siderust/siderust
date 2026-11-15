//! CLI for Siderust release validation, changelog extraction, and publication.
//!
//! ```text
//! cargo run -p siderust-release -- validate
//! cargo run -p siderust-release -- validate --tag v0.15.0 --require-tag
//! cargo run -p siderust-release -- extract-changelog --tag v0.15.0 -o notes.md
//! cargo run -p siderust-release -- package-meta Cargo.toml
//! cargo run -p siderust-release -- already-published siderust 0.14.1
//! cargo run -p siderust-release -- publish --dry-run
//! ```

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::slice::Iter;

use siderust_release::{
    crates_io_has_version, default_root, explicit_version, extract_version_section,
    format_release_notes, package_meta, publish, run_validation,
};

const USAGE: &str = "\
usage:
  siderust-release validate [--root DIR] [--tag vX.Y.Z] [--version X.Y.Z] [--require-tag] [--check-published]
  siderust-release extract-changelog [--root DIR] (--tag vX.Y.Z | --version X.Y.Z) [-o FILE]
  siderust-release package-meta <MANIFEST>
  siderust-release already-published <NAME> <VERSION>
  siderust-release publish [--root DIR] [--dry-run] [--confirm-ffi]

already-published exits 0 for present, 1 for absent, 2 for usage, 3 for unverifiable.
publish uploads missing versions unless --dry-run is set.
Live publication requires CARGO_REGISTRY_TOKEN. FFI requires explicit confirmation.";

fn help() -> ExitCode {
    println!("{USAGE}");
    ExitCode::SUCCESS
}

fn usage() -> ! {
    eprintln!("{USAGE}");
    std::process::exit(2);
}

fn take_value<'a>(args: &mut Iter<'a, String>, flag: &str) -> &'a str {
    match args.next() {
        Some(value) if !value.starts_with('-') => value.as_str(),
        _ => {
            eprintln!("error: {flag} requires a value");
            usage();
        }
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() {
        usage();
    }
    if matches!(args[0].as_str(), "-h" | "--help" | "help")
        || matches!(args.get(1).map(String::as_str), Some("-h" | "--help"))
    {
        return help();
    }
    let mut iter = args.iter();
    let cmd = iter.next().map(String::as_str).unwrap_or_default();
    match cmd {
        "validate" => cmd_validate(&mut iter),
        "extract-changelog" => cmd_extract(&mut iter),
        "package-meta" => cmd_package_meta(&mut iter),
        "already-published" => cmd_already_published(&mut iter),
        "publish" => cmd_publish(&mut iter),
        other => {
            eprintln!("error: unknown command {other}");
            usage();
        }
    }
}

fn cmd_validate(args: &mut Iter<'_, String>) -> ExitCode {
    let mut root = None;
    let mut tag = None;
    let mut version = None;
    let mut require_tag = false;
    let mut check_published = false;
    while let Some(arg) = args.next().map(String::as_str) {
        match arg {
            "--root" => root = Some(PathBuf::from(take_value(args, "--root"))),
            "--tag" => tag = Some(take_value(args, "--tag").to_owned()),
            "--version" => version = Some(take_value(args, "--version").to_owned()),
            "--require-tag" => require_tag = true,
            "--check-published" => check_published = true,
            "-h" | "--help" => return help(),
            other => {
                eprintln!("error: unknown argument {other}");
                usage();
            }
        }
    }
    let root = root.unwrap_or_else(default_root);
    match run_validation(
        &root,
        tag.as_deref(),
        version.as_deref(),
        require_tag,
        check_published,
    ) {
        Ok(notes) => {
            for note in notes {
                println!("ok: {note}");
            }
            println!("release validation passed");
            ExitCode::SUCCESS
        }
        Err(err) => {
            eprintln!("error: {err}");
            ExitCode::FAILURE
        }
    }
}

fn cmd_extract(args: &mut Iter<'_, String>) -> ExitCode {
    let mut root = None;
    let mut tag = None;
    let mut version = None;
    let mut output = None;
    while let Some(arg) = args.next().map(String::as_str) {
        match arg {
            "--root" => root = Some(PathBuf::from(take_value(args, "--root"))),
            "--tag" => tag = Some(take_value(args, "--tag").to_owned()),
            "--version" => version = Some(take_value(args, "--version").to_owned()),
            "-o" | "--output" => output = Some(PathBuf::from(take_value(args, "-o"))),
            "-h" | "--help" => return help(),
            other => {
                eprintln!("error: unknown argument {other}");
                usage();
            }
        }
    }
    let root = root.unwrap_or_else(default_root);
    let resolved = match explicit_version(tag.as_deref(), version.as_deref()) {
        Ok(Some(version)) => version,
        Ok(None) => {
            eprintln!("error: provide --tag and/or --version");
            usage();
        }
        Err(err) => {
            eprintln!("error: {err}");
            return ExitCode::FAILURE;
        }
    };
    let text = match fs::read_to_string(root.join("CHANGELOG.md")) {
        Ok(text) => text,
        Err(err) => {
            eprintln!("error: {err}");
            return ExitCode::FAILURE;
        }
    };
    let section = match extract_version_section(&text, &resolved) {
        Ok(section) => section,
        Err(err) => {
            eprintln!("error: {err}");
            return ExitCode::FAILURE;
        }
    };
    let notes = format_release_notes(&section, &resolved);
    if let Some(path) = output {
        if let Err(err) = fs::write(&path, &notes) {
            eprintln!("error: {err}");
            return ExitCode::FAILURE;
        }
    } else {
        print!("{notes}");
    }
    ExitCode::SUCCESS
}

fn cmd_package_meta(args: &mut Iter<'_, String>) -> ExitCode {
    let Some(manifest) = args.next() else {
        eprintln!("error: package-meta requires a manifest path");
        usage();
    };
    if args.next().is_some() {
        eprintln!("error: unexpected extra arguments");
        usage();
    }
    match package_meta(Path::new(manifest)) {
        Ok(meta) => {
            println!(
                "{} {} {}",
                meta.name,
                meta.version,
                if meta.publishable {
                    "publishable"
                } else {
                    "private"
                }
            );
            ExitCode::SUCCESS
        }
        Err(err) => {
            eprintln!("error: {err}");
            ExitCode::FAILURE
        }
    }
}

fn cmd_already_published(args: &mut Iter<'_, String>) -> ExitCode {
    let name = args.next().map(String::as_str);
    let version = args.next().map(String::as_str);
    if name.is_none() || version.is_none() || args.next().is_some() {
        eprintln!("error: already-published requires <NAME> <VERSION>");
        usage();
    }
    match crates_io_has_version(name.unwrap(), version.unwrap()) {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(err) => {
            eprintln!("error: {err}");
            ExitCode::from(3)
        }
    }
}

fn cmd_publish(args: &mut Iter<'_, String>) -> ExitCode {
    let mut root = None;
    let mut dry_run = false;
    let mut confirm_ffi = false;
    while let Some(arg) = args.next().map(String::as_str) {
        match arg {
            "--root" => root = Some(PathBuf::from(take_value(args, "--root"))),
            "--dry-run" => dry_run = true,
            "--confirm-ffi" => confirm_ffi = true,
            "-h" | "--help" => return help(),
            other => {
                eprintln!("error: unknown argument {other}");
                usage();
            }
        }
    }
    match publish(&root.unwrap_or_else(default_root), dry_run, confirm_ffi) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("error: {err}");
            ExitCode::FAILURE
        }
    }
}
