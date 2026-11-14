# Benchmarking Siderust

Siderust uses Criterion to measure optimized host workloads. This catalog lists
the Cargo targets declared in [Cargo.toml](../Cargo.toml). For timing boundaries,
validation evidence, and open coverage gaps, see the
[report library](reports/README.md).

```sh
cargo bench --locked -p siderust --bench converters
cargo bench --locked -p siderust --bench moon_altitude --features bench-internals -- --list
cargo bench --locked -p siderust --bench moon_altitude --features bench-internals -- \
  --exact 'public_api/moon/30d/above_threshold/horizon'
cargo bench --locked -p siderust --benches --features bench-internals,embedded-compact --no-run
```

A Cargo benchmark **target** is an executable, selected with `--bench TARGET`.
Inside it, Criterion registers **groups** and **cases**. For example, the executable
`moon_altitude` registers the group `public_api/moon/30d` and case
`above_threshold/horizon`. Criterion joins them with `/`. Pass a regex after
Cargo's `--` to filter cases, or use `--exact` for a literal ID. Use `--list` with
the features and dataset for the intended run. Harness setup, including opening
a supplied BSP, can run during listing or filtering. With `--test`, Criterion
executes each selected case once for a smoke check. Use `-- --help` on an
individual target to see its options.

## Catalog

Features below are Cargo requirements. `std` and `full-ephemeris` are defaults,
along with `serde` and `full-nutation`. Keep the defaults for the commands shown.
Plain `cargo bench` skips the three targets gated by `bench-internals`, omits
the compact cases, and skips the BSP-dependent groups without a kernel.
A compact-only library build selects a different default ephemeris.
That configuration requires its own baselines. The cross-body suite requires
the full tables because its single-evaluation epoch is outside compact coverage.

| Target | Measured operation and scope | Required features | Data | Execution and report |
| --- | --- | --- | --- | --- |
| [`altitude_comparison`](01_altitude_comparison.rs) | Public Sun, Moon and Sirius single-altitude calls and 7-, 30-, 365-day searches | `full-ephemeris` | Bundled | Host, [events](reports/coverage.md#observation-events) |
| [`converters`](02_converters.rs) | Public Cartesian ↔ spherical ICRS position conversions | None | None | Host, [coordinates](reports/coverage.md#coordinates-and-time) |
| [`de441`](03_de441.rs) | Public runtime DE4xx positions and Earth velocity, plus J2000 sidereal cases | `std` | Local BSP for all cases | Dataset host, [SPK](reports/coverage.md#spk-and-satellite-propagation) |
| [`elp2000`](04_elp2000.rs) | Public geocentric Moon position at four epochs | `full-ephemeris` | Bundled | Host, [ephemerides](reports/coverage.md#ephemerides) |
| [`ephemeris_comparison`](05_ephemeris_comparison.rs) | Public full ephemeris queries, with optional compact position pairs | `full-ephemeris`, add `embedded-compact` for pairs | Bundled | Host, [ephemerides](reports/coverage.md#ephemerides) |
| [`moon_altitude`](06_moon_altitude.rs) | Public horizon and range searches, with internal scan baseline | `bench-internals` | Bundled | Host, [Moon](reports/moon_altitude_perf_report.md) |
| [`sgp4`](07_sgp4.rs) | Public SGP4/SDP4 initialization, single and repeated TEME propagation | None | Checked-in Vallado TLEs | Host, [satellites](reports/coverage.md#spk-and-satellite-propagation) |
| [`sidereal`](08_sidereal.rs) | Public GMST and GAST with internal or precomputed nutation | None | Bundled time/nutation data | Host, [time](reports/coverage.md#coordinates-and-time) |
| [`solar_altitude`](09_solar_altitude.rs) | Public twilight, horizon and range searches, plus internal engine and twilight batch comparisons | `bench-internals` | Bundled | Host, [Sun](reports/solar_altitude_perf_report.md) |
| [`spk_perf`](10_spk_perf.rs) | Public synthetic SPK query paths, in-memory BSP parsing and real kernel queries | `std` | Synthetic built in, local BSP for optional cases | Host or dataset host, [SPK](reports/coverage.md#spk-and-satellite-propagation) |
| [`star_altitude`](11_star_altitude.rs) | Public Sirius altitude and event searches, with internal scan baseline | `bench-internals` | Bundled time/nutation data | Host, [events](reports/coverage.md#observation-events) |
| [`vsop87`](12_vsop87.rs) | Public VSOP87A/E trait calls for eight planets | `full-ephemeris` | Bundled | Host, [ephemerides](reports/coverage.md#ephemerides) |

All targets can be checked by ordinary CI. Dataset-free cases can also be
smoke-tested there. Use a controlled host for timing. The
[ESP32 harness](../tools/compact-esp32) is separate from Criterion and requires
an Xtensa toolchain and board for device measurements.

## Kernel-backed runs

Set `SIDERUST_BSP_PATH` to an existing compatible DE4xx BSP covering J2000 TT
with the Sun, Earth-chain and Moon segments. The filename `de441` is a target
name, not a requirement to use DE441 specifically. For example, with a local DE440s file:

```sh
SIDERUST_BSP_PATH=/path/to/de440s.bsp cargo bench --locked -p siderust --bench de441
SIDERUST_BSP_PATH=/path/to/de440s.bsp cargo bench --locked -p siderust --bench spk_perf -- --test
```

Without that variable, `de441` registers no cases and `spk_perf` runs its
synthetic group. A supplied file that fails to load or cover a selected query
fails the run. The harnesses never fetch a kernel. Record its checksum and
provenance before collecting measurements.

Criterion's artifacts normally go to `target/criterion/`. The locked Criterion
version has HTML reports disabled by default. Terminal estimates and saved
measurements are available without HTML. See the
[measurement and reproduction workflow](reports/methodology.md) and the
[Linux profiling guide](reports/profiling.md).
