# VSOP87 Scripts

This directory contains the build-time scripts for handling the VSOP87 planetary theory data used by Siderust.

## What is VSOP87?

VSOP87 (Variations Séculaires des Orbites Planétaires 1987) is a high-precision analytical model of planetary motions, developed by Bretagnon & Francou at the Bureau des Longitudes. It expresses the positions of the major planets as trigonometric series, with coefficients derived from astronomical observations and numerical integrations. The model is widely used in astronomy for ephemeris calculations, spacecraft navigation, and astronomical software.

## Purpose of These Scripts

The scripts in this directory automate the process of:

1. **Fetching** the official VSOP87 data files from the IMCCE FTP server.
2. **Parsing** and collecting the coefficients into an in-memory data structure.
3. **Generating** Rust source code from the parsed data, producing static arrays for each planet, coordinate, and time exponent.
4. **Writing** the generated code into Rust modules (e.g., `vsop87a.rs`, `vsop87e.rs`) in the build output directory.

This process is run automatically at build time by the [`build.rs`](../../build.rs) script. The generated Rust code is then included in the main crate, allowing Siderust to compute planetary positions efficiently and with high accuracy.

## Why is This Needed?

- **Automated Updates:** The scripts ensure that the latest VSOP87 data can be fetched and integrated without manual intervention.
- **Compile-Time Generation:** By generating Rust code at build time, the planetary data is embedded directly into the binary, eliminating runtime parsing and improving performance.
- **Reproducibility:** The scripts guarantee that the generated code is deterministic and consistent across builds, aiding in reproducibility and version control.

## Script Overview

- [`fetch.rs`](fetch.rs): Downloads the VSOP87 data files if not already present.
- [`collect.rs`](collect.rs): Parses the downloaded files and organizes the coefficients.
- [`codegen.rs`](codegen.rs): Converts the parsed data into Rust source code.
- [`io.rs`](io.rs): Handles writing the generated code to files.
- [`mod.rs`](mod.rs): Orchestrates the build pipeline and exposes the main entry point for the build script.

## Usage

You do not need to run these scripts manually. They are invoked automatically by Cargo during the build process via [`build.rs`](../../build.rs).

---
*For more details on the VSOP87 theory, see the [IMCCE documentation](https://www.imcce.fr/inpop/ephemerides/vsop87).*