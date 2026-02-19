# siderust-ffi

`siderust-ffi` is the C ABI layer for the Rust crate `siderust` (high-precision astronomy).
It exposes a flat, C-compatible API that is consumed by `siderust-cpp` and can also be used
directly from C or other languages that can call a C ABI.

The public C header `include/siderust_ffi.h` is generated from the Rust sources using
`cbindgen` (see `build.rs` and `cbindgen.toml`). Do not edit the header manually; edit the
Rust types/functions and regenerate instead.

## What This Crate Contains

- Rust `extern "C"` functions (FFI entrypoints) in `src/*.rs`.
- C-compatible types in `src/types.rs` and status codes in `src/error.rs`.
- Generated C header in `include/siderust_ffi.h` (via `cbindgen`).
- A shared library and/or static library artifact (see `Cargo.toml` `crate-type`).

## Dependencies

This crate is meant to be built from a workspace checkout where these path dependencies exist:

- `siderust` (Rust implementation) at `../siderust`
- `tempoch-ffi` (C ABI for time types used by the API) at `../siderust/tempoch/tempoch-ffi`
- `tempoch`, `affn`, `qtty` (Rust crates used internally)

If you cloned via git submodules, make sure they are initialized:

```bash
git submodule update --init --recursive
```

## Features

Cargo features mirror `siderust` options:

- `serde`: enables serialization-related APIs where applicable
- `de440`: enables DE440 ephemeris support
- `de441`: enables DE441 ephemeris support

Example:

```bash
cargo build --release --features de441
```

## Build

Building the crate generates the header and produces the library artifacts under `target/`:

```bash
cargo build --release
```

Outputs (platform dependent) are typically:

- Linux: `target/release/libsiderust_ffi.so`
- macOS: `target/release/libsiderust_ffi.dylib`
- Windows: `target/release/siderust_ffi.dll` (and an import library)

## Install / Use From C or C++

Minimum pieces to consume the ABI:

- Header: `include/siderust_ffi.h` (also requires `tempoch_ffi.h`)
- Library: the built `siderust_ffi` shared/static library

If you are using this as part of the `siderust-cpp` superproject, the top-level CMake build
invokes Cargo for you and installs both headers and libraries via `cmake --install`.

## Maintenance Notes (Avoiding Drift)

- Treat Rust as the source of truth for enums/structs/constants.
- `include/siderust_ffi.h` is generated; regenerate it by running `cargo build`.
- When changing exported signatures or the meaning/values of public enums, consider ABI
  compatibility and bump `siderust_ffi_version()` in `src/lib.rs` accordingly.

