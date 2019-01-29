# SummationByPartsOperators

[![Build Status](https://travis-ci.org/ranocha/SummationByPartsOperators.jl.svg?branch=master)](https://travis-ci.org/ranocha/SummationByPartsOperators.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/i1saoodeqrepiodl?svg=true)](https://ci.appveyor.com/project/ranocha/SummationByPartsOperators-jl)
[![Coverage Status](https://coveralls.io/repos/github/ranocha/SummationByPartsOperators.jl/badge.svg?branch=master)](https://coveralls.io/github/ranocha/SummationByPartsOperators.jl?branch=master)
[![codecov.io](http://codecov.io/github/ranocha/SummationByPartsOperators.jl/coverage.svg?branch=master)](http://codecov.io/github/ranocha/SummationByPartsOperators.jl?branch=master)

A library of classical summation-by-parts (SBP) operators used in finite difference
methods to get provably stable semidiscretisations, paying special attention to
boundary conditions.


## Basic Operators

The following derivative operators are implemented as "lazy operators", i.e. no matrix is formed explicitly.


### Periodic Domains

- `periodic_derivative_operator(derivative_order, accuracy_order, xmin, xmax, N)`

  These are classical central finite difference operators using `N` nodes on the
  interval `[xmin, xmax]`.

- `periodic_derivative_operator(Holoborodko2008(), derivative_order, accuracy_order, xmin, xmax, N)`

  These are central finite difference operators using `N` nodes on the
  interval `[xmin, xmax]` and the coefficients of [Pavel Holoborodko](http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/).

- `fourier_derivative_operator(xmin, xmax, N)`

  Fourier derivative operators are implemented using the fast Fourier transform of [FFTW.jl](https://github.com/JuliaMath/FFTW.jl).


### Finite/Nonperiodic Domains

- `derivative_operator(source_of_coefficients, derivative_order, accuracy_order, xmin, xmax, N)`

  Finite difference SBP operators for first and second derivatives can be obained by using `MattssonNordström2004()` as `source_of_coefficients`.
  Other sources of coefficients are implemented as well. To obtain a full list for all operators, use `subtypes(SourceOfCoefficients)`.

- `legendre_derivative_operator(xmin, xmax, N)`

  Use Lobatto Legendre polynomial collocation schemes on `N`, i.e.
  polynomials of degree `N-1`, implemented via [PolynomialBases.jl](https://github.com/ranocha/PolynomialBases.jl).


### Dissipation Operators

Additionally, some artificial dissipation/viscosity operators are implemented. The most basic usage is `Di = dissipation_operator(D)`,
where `D` can be a (periodic, Fourier, Legendre, SBP FD) derivative
operator. Use `?dissipation_operator` for more details.


### Conversion to Other Forms

Sometimes, it can be convenient to obtain an explicit (sparse, banded) matrix form of the operators. Therefore, some conversion functions are supplied, e.g.
```julia
julia> using SummationByPartsOperators

julia> D = derivative_operator(MattssonNordström2004(), 1, 2, 0., 1., 5)
SBP 1st derivative operator of order 2 {T=Float64, Parallel=Val{:serial}}
on a grid in [0.0, 1.0] using 5 nodes
and coefficients given in
  Mattsson, Nordström (2004)
  Summation by parts operators for finite difference approximations of second
    derivaties.
  Journal of Computational Physics 199, pp.503-540.


julia> Matrix(D)
5×5 Array{Float64,2}:
 -4.0   4.0   0.0   0.0  0.0
 -2.0   0.0   2.0   0.0  0.0
  0.0  -2.0   0.0   2.0  0.0
  0.0   0.0  -2.0   0.0  2.0
  0.0   0.0   0.0  -4.0  4.0

julia> using SparseArrays

julia> sparse(D)
5×5 SparseMatrixCSC{Float64,Int64} with 10 stored entries:
  [1, 1]  =  -4.0
  [2, 1]  =  -2.0
  [1, 2]  =  4.0
  [3, 2]  =  -2.0
  [2, 3]  =  2.0
  [4, 3]  =  -2.0
  [3, 4]  =  2.0
  [5, 4]  =  -4.0
  [4, 5]  =  2.0
  [5, 5]  =  4.0

julia> using BandedMatrices

julia> BandedMatrix(D)
5×5 BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:
 -4.0   4.0    ⋅     ⋅    ⋅
 -2.0   0.0   2.0    ⋅    ⋅
   ⋅   -2.0   0.0   2.0   ⋅
   ⋅     ⋅   -2.0   0.0  2.0
   ⋅     ⋅     ⋅   -4.0  4.0
```

## Documentation

Examples can e found in the directory [`notebooks`](https://github.com/ranocha/SummationByPartsOperators.jl/tree/master/notebooks). In particular, examples of complete discretisations of
[the linear advection equation](https://github.com/ranocha/SummationByPartsOperators.jl/blob/master/notebooks/Advection_equation.ipynb),
[the heat equation](https://github.com/ranocha/SummationByPartsOperators.jl/blob/master/notebooks/Heat_equation.ipynb),
and the [wave equation](https://github.com/ranocha/SummationByPartsOperators.jl/blob/master/notebooks/Wave_equation.ipynb) are supplied.
Further examples are supplied as [tests](https://github.com/ranocha/SummationByPartsOperators.jl/tree/master/test).
