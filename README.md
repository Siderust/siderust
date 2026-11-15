# Siderust

[![Crates.io](https://img.shields.io/crates/v/siderust.svg)](https://crates.io/crates/siderust)
[![Docs.rs](https://docs.rs/siderust/badge.svg)](https://docs.rs/siderust)
[![CI](https://github.com/Siderust/siderust/actions/workflows/github-ci.yml/badge.svg)](https://github.com/Siderust/siderust/actions/workflows/github-ci.yml)
[![License: AGPL-3.0-only](https://img.shields.io/badge/license-AGPL--3.0--only-blue)](LICENSE)

Siderust is a Rust library for positional astronomy and satellite mechanics.
It provides typed coordinate transformations, analytical and runtime
ephemerides, observing geometry and event searches, and orbit propagation
from Two-Line Elements (TLEs). Positions carry reference centers, frames,
and units in their types.

## Quickstart

Run `cargo add siderust`, then follow the
[Earth–Mars separation tutorial](doc/getting-started/quickstart.md).
The default build includes analytical planetary and lunar models and
requires no external ephemeris kernels.

## Documentation

- [Siderust manual](doc/README.md) for task guides and scientific conventions.
- [Rust API reference](https://docs.rs/siderust) and [runnable examples](examples/README.md).
- [Cargo features and datasets](doc/reference/features-and-data.md), including
  compact offline builds and runtime JPL kernels.
- [C ABI documentation](siderust-ffi/README.md).
- [Changelog](CHANGELOG.md) and [benchmark methodology](benches/README.md).

## Contributing

See the [developer tooling guide](scripts/README.md) for release checks,
reference generation, and the ESP32 harness.

See the [contributor architecture](doc/contributing/architecture.md) and
[documentation guide](doc/contributing/documentation.md). For numerical
discrepancies, include the epoch and scale, frame and center, units,
selected model, input data, and an independent reference if available.
[Report an issue](https://github.com/Siderust/siderust/issues).

## Citation and license

Cite the software using [CITATION.cff](CITATION.cff) and record the
software revision and the models and datasets used in the analysis.

Siderust is licensed under the [GNU Affero General Public License v3.0](LICENSE).
For commercial, proprietary, or differently licensed use, contact the
authors about dual licensing.

## Acknowledgments

Thanks to [Màrius Montón](https://github.com/mariusmm), whose Rust course at
the Universitat Autònoma de Barcelona helped start this project.
