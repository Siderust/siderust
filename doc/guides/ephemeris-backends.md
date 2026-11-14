# Choose an ephemeris backend

| Need | Backend | Setup |
| --- | --- | --- |
| Bundled analytical positions | `Vsop87Ephemeris`, using VSOP87A/E and ELP2000-82B | Default Cargo features. No kernels. |
| A small offline build with a bounded observing interval | `CompactEphemeris` | Enable [`embedded-compact`](embedded.md). Check the epoch before computing dependent quantities. |
| Numerical JPL Sun, Earth, and Moon states | `RuntimeEphemeris` through `DynEphemeris` | Load a compatible BSP file or its bytes. Use `DynAstroContext` where accepted. |
| Planet-center chains or kernel-defined lunar orientation | `SpkKernelSet` for compatible SPK chains, or `SpiceContext` for [lunar observation](spice-kernels.md) | Supply the required target, center, and orientation kernels. |

`AstroContext` selects the static `Ephemeris` backend through a type
parameter. If both analytical features are enabled, the full backend is the
default. Loading a kernel never changes that default. Direct theory calls
such as `Mars::vsop87a` select that theory independently of context.

## Load a JPL kernel

From the checkout, run the [runtime example](../../examples/12_runtime_ephemeris.rs)
with a local modern DE BSP that covers J2000:

```sh
cargo run --example 12_runtime_ephemeris -- /path/to/de440.bsp
```

`RuntimeEphemeris::from_bsp` requires `std`. `from_bytes` accepts bytes
supplied by the application. Enable `runtime-data` to inspect the cache
with `DatasetManager`. To download DE440 when it is missing, run:

```sh
cargo run --features runtime-data --example 12_runtime_ephemeris -- --download
```

See [features and data](../reference/features-and-data.md) for provenance
and cache configuration.

Check parsing and evaluation separately. A valid file can lack coverage at
the query epoch. Prefer the checked `DynEphemeris::try_*` methods. The
runtime front-end needs Sun/SSB, Earth–Moon barycenter (EMB)/SSB and Moon/EMB
chains, where SSB is the Solar System barycenter. It uses the explicit Earth/EMB chain when present, otherwise it
derives the Earth offset from Moon/EMB with a fixed mass ratio. This
interface converts TT to ET using a built-in periodic approximation and
exposes Sun, Earth, and Moon states. Use `SpkKernelSet` for other supported
target/center chains. See [kernel validation and precedence](../reference/guarantees.md#runtime-kernel-selection).

## Match the physical target

For geocentric comparisons, subtract the same definition of Earth's origin
from each target. A planetary system barycenter and the planet's center of
mass can be distinct SPK targets (for example, Mars system `4` and Mars
`499`). Use a compatible satellite SPK if the planetary kernel only has the
system barycenter and the calculation requires the body center. Check
coverage of both kernels over the requested dates. The compact Mars through
Neptune comparisons use system-barycenter references, and the compact
series supplies no separate satellite offsets.

Static Mercury-through-Neptune origin shifts request the selected backend's
planet capability. The lunar origin combines geocentric Moon and
barycentric Earth states from that same backend. Pluto uses a coarse catalog
orbit plus the Sun in the full default backend, while a custom backend must
supply its own Pluto capability. A missing capability returns an error on
checked paths. For comparisons across backends, also match
[axes, time, and correction settings](../concepts/accuracy-and-models.md).
