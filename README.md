# SummationByPartsOperators.jl: A Julia library of provably stable semidiscretization techniques with mimetic properties

[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ranocha.github.io/SummationByPartsOperators.jl/stable)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ranocha.github.io/SummationByPartsOperators.jl/dev)
[![Build Status](https://github.com/ranocha/SummationByPartsOperators.jl/workflows/CI/badge.svg)](https://github.com/ranocha/SummationByPartsOperators.jl/actions)
[![Codecov](http://codecov.io/github/ranocha/SummationByPartsOperators.jl/coverage.svg?branch=master)](http://codecov.io/github/ranocha/SummationByPartsOperators.jl?branch=master)
[![Coveralls](https://coveralls.io/repos/github/ranocha/SummationByPartsOperators.jl/badge.svg?branch=master)](https://coveralls.io/github/ranocha/SummationByPartsOperators.jl?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4773575.svg)](https://doi.org/10.5281/zenodo.4773575)
[![GitHub commits since tagged version](https://img.shields.io/github/commits-since/ranocha/SummationByPartsOperators.jl/v0.4.5.svg?style=social&logo=github)](https://github.com/ranocha/SummationByPartsOperators.jl)
<!-- [![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/S/SummationByPartsOperators.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html) -->

[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
is a Julia library of summation-by-parts (SBP) operators, which are discrete
derivative operators developed to get provably stable semidiscretizations,
paying special attention to boundary conditions. Discretizations included in this
framework are finite difference, Fourier pseudospectral, continuous Galerkin,
and discontinuous Galerkin methods. The main aim of
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
is to be useful for researchers and students to learn the basic concepts by
providing a unified framework of all of these seemingly different discretizations.
At the same time, the implementation is optimized to achieve good performance
without sacrificing flexibility.


## Installation

[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
is a registered Julia package. Thus, you can install it from the Julia REPL via
```julia
julia> using Pkg; Pkg.add("SummationByPartsOperators")
```


## Basic example





## Brief overview

The following derivative operators are implemented as "lazy"/matrix-free
operators, i.e. no large (size of the computational grid) matrix is formed
explicitly.


### Periodic domains

- `periodic_derivative_operator(derivative_order, accuracy_order, xmin, xmax, N)`

  These are classical central finite difference operators using `N` nodes on the
  interval `[xmin, xmax]`.

- `periodic_derivative_operator(Holoborodko2008(), derivative_order, accuracy_order, xmin, xmax, N)`

  These are central finite difference operators using `N` nodes on the
  interval `[xmin, xmax]` and the coefficients of [Pavel Holoborodko](http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/).

- `fourier_derivative_operator(xmin, xmax, N)`

  Fourier derivative operators are implemented using the fast Fourier transform of [FFTW.jl](https://github.com/JuliaMath/FFTW.jl).


### Finite (nonperiodic) domains

- `derivative_operator(source_of_coefficients, derivative_order, accuracy_order, xmin, xmax, N)`

  Finite difference SBP operators for first and second derivatives can be obained by using `MattssonNordström2004()` as `source_of_coefficients`.
  Other sources of coefficients are implemented as well. To obtain a full list for all operators, use `subtypes(SourceOfCoefficients)`.

- `legendre_derivative_operator(xmin, xmax, N)`

  Use Lobatto Legendre polynomial collocation schemes on `N`, i.e.
  polynomials of degree `N-1`, implemented via [PolynomialBases.jl](https://github.com/ranocha/PolynomialBases.jl).


### Dissipation operators

Additionally, some artificial dissipation/viscosity operators are implemented.
The most basic usage is `Di = dissipation_operator(D)`,
where `D` can be a (periodic, Fourier, Legendre, SBP FD) derivative
operator. Use `?dissipation_operator` for more details.


### Conversion to other forms

Sometimes, it can be convenient to obtain an explicit (sparse, banded) matrix form
of the operators. Therefore, some conversion functions are supplied, e.g.
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

Examples can be found in the directory
[`notebooks`](https://github.com/ranocha/SummationByPartsOperators.jl/tree/master/notebooks).
In particular, examples of complete discretisations of
[the linear advection equation](https://github.com/ranocha/SummationByPartsOperators.jl/blob/master/notebooks/Advection_equation.ipynb),
[the heat equation](https://github.com/ranocha/SummationByPartsOperators.jl/blob/master/notebooks/Heat_equation.ipynb),
and the [wave equation](https://github.com/ranocha/SummationByPartsOperators.jl/blob/master/notebooks/Wave_equation.ipynb) are supplied.
Further examples are supplied as
[tests](https://github.com/ranocha/SummationByPartsOperators.jl/tree/master/test).

The latest documentation is available
[online](https://ranocha.github.io/SummationByPartsOperators.jl/stable)
and under [`docs/src`](docs/src).


## Referencing

If you use
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
for your research, please cite it using the bibtex entry
```bibtex
@misc{ranocha2021sbp,
  title={{SummationByPartsOperators.jl}: {A} {J}ulia library of provably stable
         semidiscretization techniques with mimetic properties},
  author={Ranocha, Hendrik},
  year={2021},
  howpublished={\url{https://github.com/ranocha/SummationByPartsOperators.jl},
  doi={10.5281/zenodo.4773575}
}
```


## License and contributing

This project is licensed under the MIT license (see [LICENSE.md](LICENSE.md)).
Since it is an open-source project, we are very happy to accept contributions
from the community. Please refer to [CONTRIBUTING.md](CONTRIBUTING.md) for more
details.
