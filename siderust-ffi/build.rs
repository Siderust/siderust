use std::env;
use std::path::PathBuf;
use std::process::Command;

fn nightly_rustc_path() -> Option<String> {
    let output = Command::new("rustup")
        .args(["which", "--toolchain", "nightly", "rustc"])
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let path = String::from_utf8(output.stdout).ok()?;
    let path = path.trim();
    if path.is_empty() {
        None
    } else {
        Some(path.to_owned())
    }
}

fn main() {
    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let include_dir = PathBuf::from(&crate_dir).join("include");

    let config =
        cbindgen::Config::from_file("cbindgen.toml").expect("Unable to read cbindgen.toml");

    let previous_rustc = env::var_os("RUSTC");
    if let Some(nightly_rustc) = nightly_rustc_path() {
        env::set_var("RUSTC", nightly_rustc);
    } else {
        println!(
            "cargo:warning=nightly rustc not found; falling back to the checked-in siderust_ffi.h if cbindgen macro expansion fails"
        );
    }

    let bindings_result = cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_config(config)
        .generate();

    match previous_rustc {
        Some(value) => env::set_var("RUSTC", value),
        None => env::remove_var("RUSTC"),
    }

    let bindings = match bindings_result {
        Ok(bindings) => bindings,
        Err(err) => {
            let checked_in_header = include_dir.join("siderust_ffi.h");
            if checked_in_header.exists() {
                println!(
                    "cargo:warning=unable to regenerate siderust_ffi.h ({err}); copying the checked-in header instead"
                );
                std::fs::create_dir_all(&include_dir).expect("Unable to create include directory");
                std::fs::copy(&checked_in_header, out_dir.join("siderust_ffi.h"))
                    .expect("Unable to copy checked-in header to OUT_DIR");
                return;
            }
            panic!("Unable to generate C bindings: {err}");
        }
    };

    bindings.write_to_file(out_dir.join("siderust_ffi.h"));

    std::fs::create_dir_all(&include_dir).expect("Unable to create include directory");
    bindings.write_to_file(include_dir.join("siderust_ffi.h"));
}
