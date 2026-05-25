# `spectra`

## Scientific scope

`siderust::spectra` now contains only astronomical spectral datasets that are
specific to Siderust, chiefly photometric passbands and their throughput unit.
Generic sampled-spectrum infrastructure lives in `optica::spectrum`.

## Technical scope

What remains here:

- `passbands/`: curated built-in filter datasets such as Bessell 1990.
- `mod.rs`: thin astronomical wrapper re-exporting [`Throughput`].

Use `optica::spectrum::{SampledSpectrum, Interpolation, SpectrumError}` and
`optica::spectrum::loaders` for generic spectrum handling.
