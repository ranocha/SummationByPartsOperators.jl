# Changelog

SummationByPartsOperators.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.


## Changes in the v0.5 lifecycle

#### Added

- The boundary optimized operators of `MattssonAlmquistVanDerWeide2018Minimal`
  and `MattssonAlmquistVanDerWeide2018Accurate` are now also available with
  interior accuracy orders 10 and 12 (in addition to 4, 6, and 8).

#### Changed

- `MattssonAlmquistVanDerWeide2018Minimal` and
  `MattssonAlmquistVanDerWeide2018Accurate` can now be constructed on smaller
  grids. Previously, the equispaced grid points following the non-uniform ones
  were listed explicitly, which required more nodes than the boundary closures
  actually need. The grid points themselves are unchanged up to round-off.
- The minimum Julia version was updated to 1.10 in version 0.5.91.
- `mul!` with vectors of composite element types such as `ForwardDiff.Dual`s,
  `Complex` numbers, and `StaticVector`s is significantly faster when using
  the default `FastMode()`.

#### Fixed

- The bandwidths reported for the fourth-order accurate variable-coefficient
  second-derivative operator of `Mattsson2012()` were too small. Thus,
  converting such an operator to a `BandedMatrix` silently dropped some
  coefficients of the boundary closure (unless the variable coefficient was
  constant).
- `mul!` does not throw an error anymore for element types that cannot be
  handled by LoopVectorization.jl, e.g., `SVector{2, BigFloat}`,
  `MVector{2, Float64}`, and user-defined scalar types.
- `mul!` does not throw an error anymore for scaling factors that are no
  native numbers, e.g., `BigFloat`s.
- Some coefficients of `Mattsson2012`, `Mattsson2014`, `MattssonNordström2004`,
  `MattssonSvärdNordström2004`, `MattssonSvärdShoeybi2008`, and
  `SharanBradyLivescu2022` were rounded to `Float64` before being converted to
  the requested element type. This did not change anything for `Float64` but
  led to inexact coefficients for other element types such as `BigFloat` or
  `Rational`. In particular, `SharanBradyLivescu2022` with accuracy order 6
  did not work at all for exact element types.
- The third-derivative operator of `Mattsson2014` with accuracy order 2 was
  missing the fourth row of its boundary closure and did thus not satisfy the
  SBP property. The missing row `(1//16, -5//8, 17//16, 0, -1, 1//2)` is
  determined by the antisymmetric matrix `R` given in the reference. Since the
  boundary closure got wider, this operator requires at least 9 nodes now
  (instead of 7).

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
