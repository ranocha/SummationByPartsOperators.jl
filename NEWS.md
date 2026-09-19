# Changelog

SummationByPartsOperators.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.


## Changes in the v0.5 lifecycle

#### Changed

- The minimum Julia version was updated to 1.10 in version 0.5.91.
- `mul!` with vectors of composite element types such as `ForwardDiff.Dual`s,
  `Complex` numbers, and `StaticVector`s is significantly faster when using
  the default `FastMode()`.

#### Fixed

- `mul!` does not throw an error anymore for element types that cannot be
  handled by LoopVectorization.jl, e.g., `SVector{2, BigFloat}`,
  `MVector{2, Float64}`, and user-defined scalar types.
- `mul!` does not throw an error anymore for scaling factors that are no
  native numbers, e.g., `BigFloat`s.

#### Deprecated

- The (keyword) argument `parallel::Union{Val{:serial}, Val{:threads}}`
  is deprecated in favor of `mode` with possible values
  `FastMode()` (default), `SafeMode()`, and `ThreadedMode()`
- The non-exported struct `SumOfDerivativeOperators` is deprecated in favor of
  `LinearlyCombinedDerivativeOperators`.


## Breaking changes from v0.4.x to v0.5

- Switch from British English to American English consistently, e.g.,
  `semidiscretise` → `semidiscretize`
- `add_transpose_derivative_left!` and `add_transpose_derivative_right!`
  were replaced by the more general functions
  `mul_transpose_derivative_left!` and `mul_transpose_derivative_right!`,
  which use the same interface as `mul!`
- The number of nodes passed to `periodic_central_derivative_operator`, and
  `periodic_derivative_operator` changed from the number of visualization nodes
  to the number of compute nodes (= number of visualization nodes minus one),
  in accordance with `fourier_derivative_operator`
