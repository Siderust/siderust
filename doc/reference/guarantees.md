# Coordinate and evaluation contracts

## Type and runtime checks

- Position operations require compatible center, frame, and unit types.
  Positions and displacements have different affine roles. Direction types
  carry a frame and omit distance. Callers must label imported coordinates
  correctly, and most coordinate types do not carry the observation epoch.
- Operations requiring a common parameterized origin also compare the
  parameters. Use the checked operation where available for independent
  input sites or origins. Operator and other convenience paths may panic
  on a parameter mismatch.
- Use checked ephemeris and center-shift entry points to handle unavailable
  coverage or capabilities. Look up each method's error and panic contract
  in [Rustdoc](https://docs.rs/siderust).
- Coordinate Serde data contains components and any center parameters.
  Preserve the frame, unit, and epoch in the surrounding format and
  deserialize into the intended type.
- Siderust forbids unsafe Rust in the main crate. The separate C ABI has
  its own boundary and ownership rules. See the [C ABI documentation](../../siderust-ffi/README.md).

Frame rotations act on positions and independent vectors. To transform a
coordinate velocity between rotating axes, include the frame-rate term
`Rdot · r`. Changes between moving origins also require the relative origin
velocity. A block-diagonal rotation of a state supplies neither term.

## Runtime kernel selection

The ICRF-oriented loaders (`RuntimeEphemeris`, `SpkKernelSet`, and the
legacy body-chain parser) validate frame ID 1 and SPK types 2/3 across the
input, including unrelated segments. A producer-documented ICRF realization
is a caller requirement. Parsed constructors without frame metadata assume
ICRF. Use the indexed or stacked loaders when a kernel has multiple segments
for a fundamental target. The legacy `parse_bsp` representation cannot
store those stacks.

Within declared coverage, later segments and later loaded kernels take
priority. If the selected data or center chain fails, the query fails.
Outside a later segment's coverage, an earlier covering segment remains
eligible. See [NAIF SPK data precedence](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Data%20Precedence).

The `SpiceContext` path has its own validation. It supports J2000 state
chains for lunar observation, type 2 binary PCK orientation, and fixed
FK rotations in MATRIX and three-axis ANGLES form (radians, degrees, or
arcseconds). Binary PCK supports both little- and big-endian IEEE DAF.
Type 3 and 20 binary PCK segments can be indexed but fail when selected. General dynamic
frames, CK/SCLK composition, DSK, EK, meta-kernel expansion, and
kernel-pool redefinition semantics are outside this interface. The
lunar observer exposes position geometry. Derivatives and state-frame
transforms are outside the lunar observer interface. See [lunar kernel setup](../guides/spice-kernels.md).
