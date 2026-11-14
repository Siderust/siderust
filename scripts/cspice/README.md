# Independent CSPICE references

The [generators and pinned environment in siderust-py](https://github.com/Siderust/siderust-py/tree/fb8069a2d5710d450d40f1abcbbc3df3960a8f76/scripts/validation/cspice)
reproduce the lunar DE440 artifacts and the surface test literals using CSPICE
N0067 and independent NumPy ENU calculations. The
[lunar fixture notes](../../tests/fixtures/lunar_de440/README.md) describe
the frames, sampling and tolerances. The generators read the fixture
provenance and verify the hashes of the
[five NAIF kernels](../../doc/guides/spice-kernels.md#kernel-set).

Run the Rust regression tests with the checked-in artifacts:

```sh
cargo test --locked --features spice --test spice_lunar --test body_surface_horizontal
SIDERUST_LUNAR_KERNEL_DIR=/path/to/naif-kernels \
  cargo test --locked --features spice --test spice_lunar cspice_external -- --ignored
```
