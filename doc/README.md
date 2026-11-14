# Siderust manual

Siderust provides typed coordinates, time conversions, ephemerides, and
observation geometry for astronomy and satellite mechanics. This manual explains
how to run a calculation, choose the physical models and data it needs, and
interpret the result.

- **Start here.** [Compute the Earth–Mars separation](getting-started/quickstart.md)
  with the default build and a fixed epoch.
- **Do a calculation.** [Transform coordinates](guides/coordinate-transformations.md),
  [choose an ephemeris](guides/ephemeris-backends.md), or
  [observe from a lunar site](guides/lunar-observations.md).
  For specialized deployments, see [compact offline builds](guides/embedded.md)
  and [lunar geometry with SPICE kernels](guides/spice-kernels.md).
- **Understand the result.** Read about [coordinates and reference frames](concepts/coordinates-and-frames.md),
  [time and Earth orientation](concepts/time-and-earth-orientation.md), and
  [model accuracy and reproducibility](concepts/accuracy-and-models.md).
- **Look up a contract.** See [features and data](reference/features-and-data.md),
  [coordinate and evaluation guarantees](reference/guarantees.md), and
  [Rust API documentation](https://docs.rs/siderust).
- **Contribute.** Start with the [crate boundaries and transform invariants](contributing/architecture.md)
  and the [documentation maintenance guide](contributing/documentation.md).

The [example gallery](../examples/README.md) covers event searches, target
tracking, orbit propagation, and other workflows. The [C ABI documentation](../siderust-ffi/README.md)
describes the bindings shipped in this checkout.
