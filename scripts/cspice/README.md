# Independent CSPICE references

The generators use the N0067 CSPICE toolkit bundled with SpiceyPy 8.3.0
for frames, time conversion, positions, and kernel excerpts. Independent
NumPy ENU calculations supply surface angles and ranges. SpiceyPy provides
the existing validated CSPICE calls without a separate native toolchain
setup. [Rust CSPICE bindings](https://docs.rs/rust-spice/1.0.1/spice/)
are an alternative that would need native toolkit setup and parity checks
for the written kernels and numerical references.

Download the [five NAIF kernels](../../doc/guides/spice-kernels.md#kernel-set)
into a local directory. Both scripts verify the files they use against
[the checked-in SHA-256 manifest](../../tests/fixtures/lunar_de440/provenance.json)
before loading them. The commands below use the Python 3.14 environment
used to verify byte-for-byte reproduction.

```sh
uv venv --python 3.14 /tmp/lunar-reference-env
uv pip install --python /tmp/lunar-reference-env/bin/python -r scripts/cspice/requirements.txt
kernel_dir=/path/to/naif-kernels
/tmp/lunar-reference-env/bin/python scripts/cspice/cspice_lunar_de440_references.py \
  --kernel-dir "$kernel_dir" --verify
```

The generator queries the original DE440s SPK and lunar DE440 orientation
kernels using geometric positions (aberration correction `NONE`). It
writes excerpts with the original coefficients, ephemeris-time samples
and a separate UTC/TT/TDB comparison. See the
[fixture notes](../../tests/fixtures/lunar_de440/README.md) for frames,
sampling, and numerical tolerances. Replace `--verify` with
`--output-dir /tmp/lunar-de440` to inspect a regenerated set, or omit it to
update the checked-in fixture directory after complete generation.

For the surface tests, use the two text kernels. The output prints the
IAU Moon orientation at J2000 TT, independent lunar surface geometry,
and Mars planetographic longitude conversions for comparison with the
[stored test literals](../../tests/body_surface_horizontal.rs).

```sh
CSPICE_PCK="$kernel_dir/pck00011.tpc" CSPICE_LSK="$kernel_dir/naif0012.tls" \
  /tmp/lunar-reference-env/bin/python scripts/cspice/cspice_surface_references.py

cargo test --locked --features spice --test spice_lunar --test body_surface_horizontal
SIDERUST_LUNAR_KERNEL_DIR="$kernel_dir" \
  cargo test --locked --features spice --test spice_lunar cspice_external -- --ignored
```
