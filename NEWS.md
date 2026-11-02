# Changelog

SummationByPartsOperators.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.


## Changes in the v0.5 lifecycle

#### Added

- The first-derivative operators of `MattssonAlmquistCarpenter2014Optimal()`
  are available with accuracy order 8 now. Two entries of `Q` listed in
  Appendix B.4 of the reference are erroneous and have been corrected, see
  [#14](https://github.com/ranocha/SummationByPartsOperators.jl/issues/14).
- The boundary optimized operators of `MattssonAlmquistVanDerWeide2018Minimal`
  and `MattssonAlmquistVanDerWeide2018Accurate` are now also available with
  interior accuracy orders 10 and 12 (in addition to 4, 6, and 8).

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

- `D + combi` for a derivative operator `D` and a linear combination `combi`
  of derivative operators paired the coefficients with the wrong operators,
  e.g., `D3 + (D1 - D2)` computed `D1 + D2 - D3`.
- `ldiv!` and `\` for rational Fourier derivative operators used a wrong
  wavenumber for the last Fourier mode.
- `rat + λ * I`, `rat - λ * I`, and `λ * I - rat` ignored `λ` for rational
  periodic finite difference operators `rat`, i.e., they were only correct for
  `λ == 1`. This also affected `D ± λ * I` for periodic derivative operators.
- `fourier_derivative_matrix` returned wrong results for an odd number of
  nodes and did not initialize the diagonal, which could result in garbage
  values or an `UndefRefError` for non-bits types such as `BigFloat`.
- `mul_transpose_derivative_right!` ignored `β`. The methods of
  `mul_transpose_derivative_left!` and `mul_transpose_derivative_right!` for
  `Val(0)` computed `(α + β) * u` at the boundary node instead of `α + β * u`.
- `MattssonAlmquistCarpenter2014Optimal()` with accuracy order 2 and seven
  nodes resulted in a grid spacing of zero.
- `evaluate_coefficients` for Legendre operators returned plotting nodes in
  `[-1, 1]` instead of the physical interval.
- The Jacobians of coupled operators were computed from cell bounds that are
  shrunk slightly for floating point types, which resulted in a systematic
  relative error of the mass matrix and the derivative of a few ulp of the
  cell bounds divided by the cell length, e.g., `1.5e-3` for `Float32` cells of
  length `0.01` at `x = 100`.
- `integrate(func, u, cD)` for coupled operators `cD` does not accept vectors
  `u` of the wrong length anymore.
- `mul!` for finite difference operators throws an `ArgumentError` for arrays
  with offset axes instead of silently computing wrong results.
  `integrate(func, u, D)` for finite difference operators and
  `mul_transpose_derivative_left!` support arrays with offset axes now.
- `I - D` for Fourier derivative operators `D` resulted in a
  `StackOverflowError`.
- `issymmetric` for rational Fourier derivative operators returned `true` for
  some non-symmetric operators, e.g., `(I + D) / (I + 2 * D)`.
- `show` and `summary` of upwind operators with different accuracy orders
  displayed the accuracy order of the minus operator instead of the plus
  operator.
- `mass_matrix_boundary(D, dim)` threw an `ArgumentError` for a
  `MultidimensionalMatrixDerivativeOperator` whose `boundary_indices` contain a
  node several times, e.g., a corner shared by two boundary faces, although such
  operators are explicitly allowed. Now, the contributions of all entries of a
  node are summed, so the result does not depend on the order of the entries
  and is consistent with `integrate_boundary`. For a `TensorProductOperator`,
  the result is unchanged; its special method, which selected the entries of the
  corners based on their position in `boundary_indices`, has been removed.
- `tensor_product_operator_2D(D_x, D_y)` threw a `MethodError` if only one of
  the 1D operators has a specific source of coefficients, e.g., for a finite
  difference operator combined with a `LegendreDerivativeOperator`. Now,
  `SourceOfCoefficientsCombination` accepts any object returned by
  `source_of_coefficients`.
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
- The kernels of the Fourier spectral viscosity operators
  `MadayTadmor1989()`, `TadmorWaagan2012Convergent()`, and `Tadmor1993()`
  were shifted by one wavenumber relative to the references, e.g.,
  `Tadmor1993()` used `1 - ((m + 1) / (k + 1))^(2s)` instead of
  `1 - (m / k)^(2s)` for the wavenumber `k` and the cutoff wavenumber
  `m = cutoff - 1`. Now, they match the references. The keyword argument
  `cutoff` of `dissipation_operator` for Fourier derivative operators is
  documented now. Moreover, the kernels of these operators and of
  `TadmorWaagan2012Standard()` are computed in the requested element type
  instead of `Float64`.

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
