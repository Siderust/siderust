# Compact ephemeris footprint and timing

Measurements from the compact integration at Siderust commit `c15ccf1`.
See the [compact profile](../../doc/guides/embedded.md) for current setup.

## Linux

Linux reproduction [loftyinclination/moon](https://github.com/loftyinclination/moon)
at `1e7c973808695074edf3326e1291d3c60487489e`, x86-64, Rust nightly
`3ef0c22cd`. Use the reporter's `cargo +nightly build --release -Z
build-std=std,panic_abort` with its existing `opt-level = "z"`, LTO,
one codegen unit and stripping. The compact measurements used the pinned
dependency graph at the integration revision, enabled `embedded-compact`,
and substituted the lunar workflow from example 25. The reporter's location passed 51.507222 as longitude
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

Use the [ESP32 harness](../../tools/compact-esp32) and the Xtensa toolchain
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
