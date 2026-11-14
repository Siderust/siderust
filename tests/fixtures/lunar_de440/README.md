# Independent CSPICE validation of lunar geometry

The [lunar integration tests](../../spice_lunar.rs) use 282 independently
generated observations over 32 days. The sites include the equator,
mid-latitudes, both poles, a representative Commercial Lunar Payload
Services site, and both sides of the ±180-degree boundary. Sun horizon
crossings are sampled one second before and after the transition. Tests
cover PA, ME, and a synthetic 30/50/20-degree FK frame, as well as rotation
inverses and determinants, interpolation boundaries, temporal continuity,
precedence, and invalid data. Azimuth comparisons use circular differences
at nonsingular points.

Only 16 KiB of original SPK/BPC coefficient excerpts are committed, alongside
the FK, LSK, references and hashes. CSPICE `SPKSUB` and `PCKW02` preserve the
original coefficients. An extra BPC record preserves the original kernel's
selection at the excerpt's final endpoint. The generator queries the original
full kernels through SpiceyPy 8.3.0 and CSPICE N0067. CSPICE `pxform`,
`spkpos(...,"NONE",...)`, and `latrec` supply matrices, positions, and sites.
Independent NumPy ENU dot products give angles, ranges, and center visibility.

Observed maximum residuals on Linux x86-64 (both excerpts and full SPK/BPC):

| Quantity | Maximum residual | Regression bound |
| --- | --- | --- |
| Rotation matrix element | 2.23e-16 | 5e-14 |
| Fixed target/site-relative component | 2.98e-8 km (0.030 mm) | 2e-7 km |
| Slant range | 5.96e-8 km (0.060 mm) | 2e-7 km |
| Altitude/circular azimuth | 5.69e-14 degree | 1e-10 degree |

The separate time conversion test covers epochs in 1972, 1999, 2016,
2017 and 2026, including the 2016 leap second. The test allows residuals
below 1 microsecond for TT and 50 microseconds for Tempoch's TDB versus
CSPICE DELTET. Geometry fixtures use the same ET in both implementations.

The bounds account for floating-point evaluation and platform rounding.
Angle reduction follows NAIF PCKE02's signed D_MOD subtraction. These
residuals describe numerical agreement when both implementations use
the same model and kernels. Full-BPC coverage-edge rotation
checks have a 5e-12 element bound, including epochs centuries from J2000.
A Horizons comparison requires matching corrections, epochs, sites, and
frame models.

To regenerate the references, obtain the [original kernels](../../../doc/guides/spice-kernels.md#kernel-set) and run these commands from the repository root. The
[provenance manifest](provenance.json) records the hashes.

```sh
uv venv /tmp/lunar-reference-env
uv pip install --python /tmp/lunar-reference-env/bin/python spiceypy==8.3.0 numpy==2.5.4
KERNEL_DIR=/tmp/lunar-kernels /tmp/lunar-reference-env/bin/python \
  tools/cspice_refs/cspice_lunar_de440_references.py
cargo test --features spice --test spice_lunar
SIDERUST_LUNAR_KERNEL_DIR=/tmp/lunar-kernels \
  cargo test --features spice --test spice_lunar cspice_external -- --ignored --nocapture
```

The optional external test spans approximately 1850–2150 and tests the original
BPC's coverage endpoints. No large kernels or network access are needed for
ordinary CI.
