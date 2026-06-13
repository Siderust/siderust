//! Build script for `siderust-ffi`, generating the exported C header.

use std::env;
use std::path::{Path, PathBuf};

fn copy_checked_in_header(include_dir: &Path, out_dir: &Path) -> bool {
    let checked_in_header = include_dir.join("siderust_ffi.h");
    if !checked_in_header.exists() {
        return false;
    }

    std::fs::create_dir_all(include_dir).expect("Unable to create include directory");
    std::fs::copy(&checked_in_header, out_dir.join("siderust_ffi.h"))
        .expect("Unable to copy checked-in header to OUT_DIR");
    true
}

fn main() {
    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let include_dir = PathBuf::from(&crate_dir).join("include");

    let config =
        cbindgen::Config::from_file("cbindgen.toml").expect("Unable to read cbindgen.toml");

    let bindings_result = cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_config(config)
        .generate();

    let bindings: cbindgen::Bindings = match bindings_result {
        Ok(bindings) => bindings,
        Err(err) => {
            if copy_checked_in_header(&include_dir, &out_dir) {
                println!(
                    "cargo:warning=unable to regenerate siderust_ffi.h ({err}); copying the checked-in header instead"
                );
                return;
            }
            panic!("Unable to generate C bindings: {err}");
        }
    };

    bindings.write_to_file(out_dir.join("siderust_ffi.h"));

    std::fs::create_dir_all(&include_dir).expect("Unable to create include directory");
    bindings.write_to_file(include_dir.join("siderust_ffi.h"));
}
