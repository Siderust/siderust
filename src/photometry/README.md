# `photometry`

## Scientific scope

`siderust::photometry` owns astronomical photometric concepts: photometric
systems, spectral response curves (passbands), throughput unit, synthetic
photometry, and color indices. Generic sampled-spectrum infrastructure lives
in `optica::spectrum`.

## Technical scope

What you will find here:

- `passbands/`: curated built-in filter datasets such as Bessell 1990 UBVRI,
  with lazily-initialised, statically-cached `SampledSpectrum<Nanometer, Throughput>`
  constants and `johnson_b()` / `johnson_v()` convenience accessors.
- `mod.rs`: thin re-export of [`Throughput`].

Use `optica::spectrum::{SampledSpectrum, Interpolation, SpectrumError}` and
`optica::spectrum::loaders` for generic spectrum handling.
