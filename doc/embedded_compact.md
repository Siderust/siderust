# Compact offline ephemerides

Enable the profile with `default-features = false, features = ["embedded-compact"]`.
With defaults disabled, select either this profile or `full-ephemeris`.
Add `full-nutation` to select the full IAU 2006A default. The ordinary
default still includes full VSOP87A/E, ELP2000-82B, and IAU 2006A.
Both ephemeris features can coexist. In that case the ordinary default
context selects the full backend and callers can choose
`AstroContext::<CompactEphemeris, _>` explicitly.

Use the checked `Ephemeris::try_*` methods for input dates. The input is TT
and Tempoch converts it to TDB before the Archive evaluator receives a
scalar Julian date. Both models cover `[2458849.5, 2470172.5)` TDB, from
2020-01-01 up to 2051-01-01. The existing infallible trait methods panic
outside coverage. Tempoch rejects non-finite values when constructing typed
epochs. The Sun and planets are barycentric in AU, and the Moon is
geocentric in km. Each model uses its own native-to-canonical rotation.
DE440s comparisons for Mars through Neptune use the system-barycenter
reference segments. The compact theory supplies no separate satellite
correction for those bodies.

The compact Earth velocity uses a 675-second centered difference with a
second-order one-sided stencil at the coverage edges. The position and
observation path allocates no heap and performs no I/O. The ESP32 harness
links a null allocator and reserves zero heap.

[Example 25](../examples/25_embedded_compact.rs) shows the migration for
the lunar-altitude, phase, and bright-limb workload in issue #140. Compute
the checked Moon vector and the geocentric Sun (negative heliocentric Earth),
convert them once to topocentric true-of-date coordinates, then reuse them
with `moon_phase_from_topocentric`, `moon_bright_limb_from_topocentric` and
`equatorial_to_horizontal_true_of_date_with_ctx`.

The example includes WGS84 observer parallax, IAU 2006 precession, IAU 2000B
nutation, and Earth rotation. Light-time, aberration, gravitational deflection,
and refraction are outside this geometric example. The independent regression
also tests annual aberration with the compact Earth velocity. For high-accuracy
pointing, supply an EOP provider and verify the epoch with
`AstroContext::try_eop_at_tt`. The Archive snapshot bundles monthly ΔT from
MJD 41714 through 63871 (observed through 61314) and no daily IERS EOP.
The example explicitly chooses `NullEop` (UT1 = UTC and zero polar motion).

## Independent comparison

The [generator](../tools/validation/compact_reference.py) uses DE440s and
ERFA. Regenerate the fixture with Python and `numpy pyerfa jplephem`.
The SHA-256 of `de440s.bsp` is
`c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2`.

```sh
python tools/validation/compact_reference.py /path/to/de440s.bsp > tests/fixtures/compact_reference.csv
cargo test --locked -p siderust --test compact_ephemeris --no-default-features --features embedded-compact -- --nocapture
```

Across 486 body/site/epoch cases, including both ends of the supported
interval and targets within 0.051° of the horizon:

| Maximum discrepancy | Measured |
| --- | ---: |
| Geocentric angular separation | 3.89″ |
| Topocentric angular separation | 4.16″ |
| Raw azimuth difference | 23.57″ |
| Altitude difference | 3.99″ |
| Topocentric separation with annual aberration | 4.16″ |
| Lunar radial distance | 0.91 km |
| Earth barycentric velocity | 0.182 m/s |
| First-order aberration bound from velocity difference | 0.000126″ |
| TT–TDB offset difference against ERFA | 8.02 µs |
| IAU 2000B vs ERFA IAU 2006A longitude / obliquity | 1.73 / 0.65 mas |
| Changing synthetic EOP to UT1=UTC and zero polar motion | 6.07″ |

The geometric reference uses shared synthetic DUT1 = 0.4 seconds, xp = 0.15″
and yp = −0.25″, including for future dates. This checks transformations
independently of the availability of measured EOP. The comparison led to
fixing the observer-parallax Earth rotation sign and the CIO rotation in
the shared transformation pipeline.

## Binary measurements

Linux reproduction [loftyinclination/moon](https://github.com/loftyinclination/moon)
at `1e7c973808695074edf3326e1291d3c60487489e`, x86-64, Rust nightly
`3ef0c22cd`. Use the reporter's `cargo +nightly build --release -Z
build-std=std,panic_abort` with its existing `opt-level = "z"`, LTO,
one codegen unit and stripping. For the migration, point the dependency
at this checkout, enable `embedded-compact`, add the Archive and Tempoch
patches to the reproduction's root manifest and substitute the lunar workflow
from example 25. The reporter's location passed 51.507222 as longitude
and −0.1275 as latitude. That order was retained in the measurement.

| Linux configuration | .text bytes | .rodata bytes | PT_LOAD file bytes | PT_LOAD memory bytes |
| --- | ---: | ---: | ---: | ---: |
| Registry 0.18.1 | 71,981 | 4,631,066 | 4,712,192 | 4,714,760 |
| 0.18.1, ordinary inlining throughout full lunar evaluator | 47,712 | 4,631,466 | 4,688,712 | 4,690,176 |
| 0.18.1, no inlining throughout full lunar evaluator | 47,908 | 4,631,466 | 4,689,984 | 4,694,280 |
| Compact, IAU 2006A | 29,329 | 162,624 | 199,296 | 200,952 |
| Compact, IAU 2000B | 28,373 | 62,856 | 98,576 | 102,648 |

Use `size -A` for the sections and sum the file and memory sizes from
`readelf -Wl` for the loadable segments. The compact executable retains
248 bytes of .bss and 8 bytes of .data. These measurements exclude loader,
shared-library, and process-stack overhead. The compact workload took
39.8 µs per evaluation versus 13.9 ms for the original on an Intel Core
Ultra 5 325, using the median of three 2,000-evaluation trials with
successive one-minute offsets and output outside the timed loop. The
original evaluates multiple full lunar positions while the migration
reuses its vectors. Ordinary and never inlining took 10.5 ms and 8.8 ms
respectively in those trials. Changing inlining primarily reduced code
size, while truncating the coefficients removed the dominant read-only data.

## ESP32

Use the [ESP32 harness](../tools/compact-esp32) and the Xtensa toolchain
from espup (tested with `rustc 1.97.0-nightly 8ea53bcd7`). Source the espup
environment so Cargo can find `xtensa-esp32-elf-gcc`.

```sh
cd tools/compact-esp32
cargo +esp build --locked --release
elf=target/xtensa-esp32-none-elf/release/siderust-compact-esp32
xtensa-esp32-elf-size -A "$elf"
readelf -Wl "$elf"
esptool --chip esp32 elf2image --flash-mode dio --flash-size 4MB --flash-freq 40m -o compact.bin "$elf"
espflash flash --monitor --chip esp32 "$elf"
```

The all-body harness produced a 189,888-byte flash image with esptool 5.6.0,
68,800 bytes .rodata, 58,726 bytes flash .text, 4,512 bytes RAM .rwtext,
2,008 bytes .data, and 344 bytes .bss. The linker reserves 194,256 bytes
of the DRAM range for the stack. This reservation does not measure peak
stack use. The UART reports seven 128-epoch trials at 80 MHz when run on
hardware. Measure stack high-water on the actual board with the application's
memory partition and concurrent tasks. ESP32 device latency and peak RAM
remain unmeasured.

## Dependency pins and publication follow-up

The original locked graph used Archive 0.1.5 through both Siderust and
Tempoch Core 0.11.2. Full VSOP, ELP, and nutation were unconditional.
The new graph uses one Archive 0.2.0 at
`5ec0343d4623b84ffbe5f990d15f0d48986b425d` and the Tempoch 0.11.3 family at
`1a28a91de777b21caf3499b4d5f950545381f0f1`. Archive's blocking package-version
mismatch and lack of independently selectable IAU 2000B tables were fixed
on the pinned branch. [Tempoch #27](https://github.com/Siderust/tempoch/pull/27)
updates the Archive time-data dependency and coordinates the FFI crate.

Check out Archive and Tempoch next to the Siderust repository. Cargo's
temporary root patches and the ESP32 harness patches point at those
directories. CI populates the same layout from exact commits.

After the maintainers publish the upstream versions, remove the patches
from both manifests and the sibling checkouts from both projects' CI.
Restore the Siderust publishing workflow's package dry run. Regenerate
both lockfiles, check the Archive dependency with
`cargo tree -i siderust-archive`, repeat tests and size measurements, and
verify the packages before publishing Siderust.
