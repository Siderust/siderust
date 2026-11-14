# Cargo features and scientific data

## Cargo features

[Cargo.toml](../../Cargo.toml) defines the complete feature graph. The
features below change the available models, platform requirements, or data.

| Feature | Effect |
| --- | --- |
| `std` (default) | Standard library and filesystem APIs. Without it the library uses `no_std` with `alloc`. |
| `serde` (default) | Serialization for supported types. Filesystem-backed catalogs also require `std`. |
| `full-ephemeris` (default) | Full VSOP87A/E and ELP2000-82B. Takes precedence in the default context when both ephemeris features are enabled. |
| `embedded-compact` | Checked, truncated VSOP87E and ELP2000-82B. Becomes the default ephemeris if `full-ephemeris` is disabled. |
| `full-nutation` (default) | Enables the full IAU 2000A and 2006A models and makes IAU 2006A the default. Without it, IAU 2000B is the default. |
| `lagrange-centers` | Archive-backed Sun–Earth Lagrange support. Implies `serde`. |
| `runtime-data` | Dataset download/cache helpers. Implies `std`. |
| `spice` | High-level kernel context and lunar observer. Implies `std`. |
| `atmosphere`, `photometry` | Atmospheric optics and photometric interfaces. |
| `pod`, `pod-parquet`, `pod-doris` | Precise orbit determination, optional Parquet output, and optional DORIS input. `pod` implies `std` and `serde`. |

Select at least one ephemeris feature when disabling defaults. Cargo
unifies features requested for the same dependency, so inspect
`cargo tree -e features -i siderust` if the final application enables more
than intended. See [the compact guide](../guides/embedded.md) for
deployment considerations. `bench-internals` exposes benchmark helpers
outside the stable public API.

## Data ownership and coverage

| Data | Source and loading | Coverage and provenance |
| --- | --- | --- |
| VSOP87, ELP2000 and compact coefficients | Bundled accessors from [Siderust Archive](https://github.com/Siderust/archive). The underlying series come from IMCCE. | Record the Archive version and model. Compact validity is [bounded in TDB](../guides/embedded.md#evaluate-a-supported-epoch). Full analytical series have model-dependent accuracy. |
| Nutation, gravity, atmosphere, Pluto, Lagrange SCK | Archive snapshots and feature-gated accessors. | Use the Archive family manifests for upstream sources, checksums, and domain. |
| Leap seconds, ΔT and EOP | [Tempoch](https://docs.rs/tempoch) owns time interpretation and active bundles. Archive supplies bundled time data. | The bundled fallback includes monthly ΔT and leap-second data. Daily IERS EOP requires suitable active data. Inspect the actual table date range and observed/predicted status. |
| JPL DE and supplemental SPK, PCK and FK files | Load at runtime. Acquire manually or use Archive fetch helpers where supported. | Coverage depends on the selected segment and all required chain edges. Record identities, hashes, and load order. |

With `runtime-data`, `siderust_archive::jpl::DatasetManager` can acquire and
cache supported JPL datasets. Its default cache is `~/.siderust/data/`.
Set `SIDERUST_DATA_DIR` to override it. Check the data manager’s cache path
when troubleshooting a missing kernel. Archive owns dataset regeneration and
validation tooling. Use its [repository and manifests](https://github.com/Siderust/archive)
for maintenance.

For body-center offsets, check the kernel summaries and overlap. Examples
of supplemental SPKs referenced by this project are `mar099.bsp` (Mars
499), `jup365.bsp` (Jupiter 599), `sat441l.bsp` (Saturn 699),
`ura184.bsp` (Uranus 799), and `nep097.bsp` (Neptune 899). Use
[the runtime guide](../guides/ephemeris-backends.md#match-the-physical-target)
to choose target and origin chains.
