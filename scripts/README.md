# Developer tooling

Run these commands from the repository root unless the tool's instructions
specify another directory.

| Tool | Purpose | Prerequisites |
| --- | --- | --- |
| [Release CLI](#releases) | Maintainer preflight, notes, and crates.io publication | Cargo |
| [CSPICE references](cspice/README.md) | Reproduce or verify independent lunar and surface references | Python, pinned SpiceyPy and NumPy, local NAIF kernels |
| [Compact and SGP4 references](validation/README.md) | Reproduce independent DE440s/ERFA and Vallado vectors | Python and the pinned reference libraries |
| [ESP32 harness](../benches/reports/compact_ephemeris.md#esp32) | Build the compact workload for hardware measurements | Espressif Xtensa toolchain |
| [Markdown link checker](check-doc-links.py) | Check relative links and heading fragments in the manuals and reports | Python 3 standard library |

The Markdown checker uses a small standard-library scanner for the
repository's inline links and GitHub heading fragments, so the docs job
requires no Markdown library installation. Run it with
`python3 scripts/check-doc-links.py`.

## Releases

```sh
cargo run --locked -p siderust-release -- validate
cargo run --locked -p siderust-release -- validate --tag v0.19.0 --require-tag
cargo run --locked -p siderust-release -- extract-changelog --tag v0.19.0 -o /tmp/siderust-notes.md
cargo run --locked -p siderust-release -- publish --dry-run
```

Local validation checks workspace versions, dependency sources, and the
matching changelog section. Tagged validation and live publication reject
Cargo patches.
Use `--check-published` to include a read-only crates.io lookup. The
`already-published NAME VERSION` subcommand exits with 0 for an existing
version, 1 for a confirmed absence, 2 for incorrect usage, and 3 when the
lookup cannot be verified.

The dry run queries crates.io and runs Cargo packaging checks for missing
versions. A live `publish` requires `CARGO_REGISTRY_TOKEN` and can upload the
primary crate. Both modes skip versions already on crates.io and stop on
unverifiable registry responses. The FFI crate requires `--confirm-ffi`
and a manifest that permits crates.io publication. Publish jobs in GitHub
Actions run only on release tag pushes.

## Adding tools

Use Rust for repository automation and keep tooling out of the library's
runtime dependencies. Set `publish = false` on tooling crates. Document the
inputs, external effects, and verification command beside each tool.
