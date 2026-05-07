# `data`

## Scientific scope

Large modern astronomy datasets such as JPL ephemerides and SPICE kernels are
too large and too frequently revised to be treated like small compile-time
constants. This module exists to manage those external reference products in a
reproducible way: a caller can identify a dataset, fetch it into a local cache,
verify its integrity, and then hand the resulting bytes to the scientific code
that understands the format.

The scientific content of the data is not defined here; this module is the
logistics layer that keeps access to authoritative upstream products explicit,
auditable, and byte-stable.

## Technical scope

What you will find here:

- `registry.rs`: the pinned list of supported dataset identifiers and metadata.
- `manager.rs`, `download.rs`, `cache.rs`: dataset acquisition and cache
  management when the `runtime-data` feature is enabled.
- `daf.rs`: Double Precision Array File parsing helpers.
- `spk.rs`: SPK-specific parsing primitives.
- `DataError`: the common error type for IO, integrity, and parsing failures.

Without `runtime-data`, the registry and parser-facing types remain available,
so dependent code can reason about datasets without pulling in network or cache
machinery.

## References

- NASA NAIF. SPK Required Reading.
  <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html>
- Acton, C. H., Bachman, N., Semenov, B., and Wright, E. (2018). "A look
  towards the future in the handling of space science mission geometry."
  Planetary and Space Science 150, 9-12.
- Astropy Project documentation. Downloadable data management.
  <https://docs.astropy.org/en/latest/utils/data.html>
