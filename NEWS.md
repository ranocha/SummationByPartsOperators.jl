# Changelog

SummationByPartsOperators.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.


## Changes in the v0.5 lifecycle

#### Added

- The boundary-optimized SBP operators for second derivatives with variable
  coefficients of `StiernströmAlmquistMattsson2023` are now available with
  interior accuracy orders 4, 6, 8, 10, and 12. They are defined on the
  non-uniform grid of `MattssonAlmquistVanDerWeide2018Accurate` and are fully
  compatible with the corresponding first-derivative operators, which are
  available via the same source of coefficients.
- The boundary optimized operators of `MattssonAlmquistVanDerWeide2018Minimal`
  and `MattssonAlmquistVanDerWeide2018Accurate` are now also available with
  interior accuracy orders 10 and 12 (in addition to 4, 6, and 8).
- `integrate` now also works for the variable coefficient derivative operators
  created by `var_coef_derivative_operator`, for periodic dissipation
  operators, and for the lazy rational operators, operator quotients, and
  (super) spectral viscosity operators. Previously, these threw a `MethodError`
  or an error about a missing field `Δx`. The wrapping operators use the
  quadrature rule of the derivative operator they are built from.

#### Changed

- `MattssonAlmquistVanDerWeide2018Minimal` and
  `MattssonAlmquistVanDerWeide2018Accurate` can now be constructed on smaller
  grids. Previously, the equispaced grid points following the non-uniform ones
  were listed explicitly, which required more nodes than the boundary closures
  actually need. The grid points themselves are unchanged up to round-off.
- The coefficients of `MattssonAlmquistVanDerWeide2018Minimal` and
  `MattssonAlmquistVanDerWeide2018Accurate` are now stored as the diagonal norm
  `H` and the antisymmetric `Q` printed in the paper, read as exact rational
  numbers, instead of the truncated decimals of `Q[i, j] / H[i, i]`. Thus,
  these operators now satisfy the SBP property exactly for exact element types
  such as `Rational{BigInt}`. In `Float64`, all but eight of the 224 boundary
  coefficients of the accuracy orders 4, 6, and 8 available before are
  bit-identical to the previous ones; the remaining eight change by at most
  101 ulp since the paper prints the interior stencil coefficients they are
  coupled to (`1/12`, `1/60`, `4/105`, `1/280`) as truncated decimals. For
  element types with more precision than `Float64`, e.g., `BigFloat`, the
  coefficients are more accurate than before since they are no longer rounded
  to `Float64` first.
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
