# Introduction

Summation-by-parts (SBP) operators are discrete derivative operators designed to
enable (semi-) discrete stability proofs mimicking the energy method from the
continuous level. To do so, SBP operators mimic integration-by-parts discretely.
Here, we will briefly explain the basic concepts. If you want to learn more about
this subject, the classical review articles of [^SvärdNordström2014] and
[^FernándezHickenZingg2014] are good starting points. More recent references and
applications of SBP operators from many classes implemented in
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
are given by [^RanochaMitsotakisKetcheson2021].

Since SBP operators are designed to mimic integration-by-parts, they need a notion
of derivatives and integrals. Here, derivatives are interpreted as linear operators `D`
(derivative matrices) and integrals are interpreted as discrete inner products,
represented by the associated mass/norm matrices `M`. Thus, the discrete derivative
of a grid function `u` is `D * u` and the discrete inner product of two grid functions
`u` and `v` is `dot(u, M, v)`, where `M = mass_matrix(D)`. Here, we have already
introduced some basic interfaces provided by
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl):
- Derivative operators act as linear operators implementing `*` (and `mul!` for
  more efficient in-place updates avoiding allocations).
- The mass matrix associated to an SBP derivative operator can be retrieved via
  [`mass_matrix`](@ref).


## Periodic domains

Periodic (central) SBP operators mimic the properties of differential operators
on periodic domains. Hence, they are

- skew-symmetric if they approximate odd derivatives
- symmetric and semi-definite if they approximate even derivatives;
  second-derivative operators are negative semi-definite,
  fourth-derivative operators are positive semi-definite etc.

Classical central finite difference operators on periodic domains are periodic
SBP operators. They can be constructed via [`periodic_derivative_operator`](@ref).
Similarly, Fourier collocation methods can be interpreted as periodic SBP operators,
which can be constructed via [`fourier_derivative_operator`](@ref).

```jldoctest
julia> using SummationByPartsOperators, LinearAlgebra

julia> D = periodic_derivative_operator(derivative_order=1, accuracy_order=2,
                                        xmin=0.0, xmax=2.0, N=20)
Periodic first-derivative operator of order 2 on a grid in [0.0, 2.0] using 20 nodes,
stencils with 1 nodes to the left, 1 nodes to the right, and coefficients of Fornberg (1998)
  Calculation of Weights in Finite Difference Formulas.
  SIAM Rev. 40.3, pp. 685-691.

julia> M = mass_matrix(D)
UniformScaling{Float64}
0.1*I

julia> M * Matrix(D) + Matrix(D)' * M |> norm
0.0

julia> D = fourier_derivative_operator(xmin=0.0, xmax=2.0, N=20)
Periodic 1st derivative Fourier operator {T=Float64}
on a grid in [0.0, 2.0] using 20 nodes and 11 modes

julia> M = mass_matrix(D)
UniformScaling{Float64}
0.1*I

julia> norm(M * Matrix(D) + Matrix(D)' * M) < 10 * eps(eltype(D))
true
```

As you have seen above, conversion methods to other common types such as `Matrix`,
`sparse` from the standard library SparseArrays, and `BandedMatrix` from
[BandedMatrices.jl](https://github.com/JuliaMatrices/BandedMatrices.jl) are
available.


## Non-periodic domains

On non-periodic domains, additional boundary terms appear. Thus, the basic
symmetry properties of SBP operators are the same as the ones of periodic SBP
operators modulo boundary terms. Note that the correct handling of boundary terms
is the basic reason of the success of SBP operators. In particular for hyperbolic
problems, other boundary treatments that might appear senseful can result in
catastrophic failure.

### First-derivative operators

First-derivative SBP operators need to mimic
```math
  \int_{x_\mathrm{min}}^{x_\mathrm{max}} u(x) \bigl( \partial_x v(x) \bigr) \mathrm{d}x
+ \int_{x_\mathrm{min}}^{x_\mathrm{max}} \bigl( \partial_x u(x) \bigr) v(x) \mathrm{d}x
= u(x_\mathrm{max}) v(x_\mathrm{max}) - u(x_\mathrm{min}) v(x_\mathrm{min}).
```
Thus, a discrete evaluation at the boundary of the domain is necessary. For
SBP operators with a grid including the boundary nodes, this can be achieved
by simply picking the first/last nodal coefficient of a grid function `u`.
If boundary nodes are not included, some interpolation is necessary in general.
Nevertheless, getting a boundary value is a linear functional that is often
represented in the literature using (transposed) vectors `tL, tR`. Then,
an SBP operator has to satisfy `M * D + D' * M == tR * tR' - tL * tL'`.
The boundary operators are represented matrix-free via
[`derivative_left`](@ref) and [`derivative_right`](@ref) for zeroth-order
derivatives.

```@repl
using SummationByPartsOperators, LinearAlgebra

D = derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=2,
                        xmin=0//1, xmax=1//1, N=9)
tL = zeros(eltype(D), size(D, 1)); tL[1] = 1; tL'
tR = zeros(eltype(D), size(D, 1)); tR[end] = 1; tR'
M = mass_matrix(D)

M * Matrix(D) + Matrix(D)' * M == tR * tR' - tL * tL'
u = randn(size(grid(D))); derivative_left(D, u, Val(0)) == u[begin]
u = randn(size(grid(D))); derivative_right(D, u, Val(0)) == u[end]
```

Here, we have introduced some additional features. Firstly, exact rational
coefficients are provided, based on the type of `xmin` and `xmax` (if available).
Secondly, a [`source_of_coefficients`](@ref) has to be provided when constructing
the SBP operator. You can list them using
```@example
using InteractiveUtils, SummationByPartsOperators
subtypes(SourceOfCoefficients)
```
Here and in the following, the order of accuracy of (finite difference) SBP
operators refers to the local order of accuracy in the interior, cf.
[`accuracy_order`](@ref).

A special case of first-derivative SBP operators are polynomial derivative operators
on Lobatto-Legendre nodes, implemented in [`legendre_derivative_operator`](@ref).

### Second-derivative operators

To mimic integration-by-parts of second derivatives,
```math
  \int_{x_\mathrm{min}}^{x_\mathrm{max}} u(x) \bigl( \partial_x^2 v(x) \bigr) \mathrm{d}x
= - \int_{x_\mathrm{min}}^{x_\mathrm{max}} \bigl( \partial_x u(x) \bigr) \bigl( \partial_x v(x) \bigr) \mathrm{d}x
  + u(x_\mathrm{max}) \bigl( \partial_x v(x_\mathrm{max}) \bigr)
  - \bigl( \partial_x u(x_\mathrm{min})) v(x_\mathrm{min}),
```
the evaluation of the first derivative at the boundaries is necessary. These
linear functionals are available as [`derivative_left`](@ref) and
[`derivative_right`](@ref). In the literature, they are often called `dL` and
`dR`. Then, a second-derivative SBP operator has to be of the form
`M * D == -A + tR * dR' - tL * dL'`, where `A` is symmetric and positive
semidefinite.

```@repl
using SummationByPartsOperators, LinearAlgebra

D = derivative_operator(MattssonNordström2004(), derivative_order=2, accuracy_order=2,
                        xmin=0//1, xmax=1//1, N=9)

M = mass_matrix(D)
tL = derivative_left(D, Val(0)); tL'
tR = derivative_right(D, Val(0)); tR'
dL = derivative_left(D, Val(1)); dL'
dR = derivative_right(D, Val(1)); dR'

A = -M * Matrix(D) + tR * dR' - tL * dL'
isposdef(A)
```
Usually, there is no need to form `dL, dR` explicitly. Instead, you can use the
matrix-free variants [`derivative_left`](@ref) and [`derivative_right`](@ref).
Some procedures imposing boundary conditions weakly require adding the transposed
boundary derivatives to a grid function, which can be achieved by
[`mul_transpose_derivative_left!`](@ref) and [`mul_transpose_derivative_right!`](@ref).
You can find applications of these operators in the source code of
[`WaveEquationNonperiodicSemidiscretization`](@ref).

A special case of second-derivative SBP operators are polynomial derivative operators
on Lobatto-Legendre nodes, implemented in [`legendre_second_derivative_operator`](@ref).


## [Upwind operators](@id intro-upwind-operators)

Upwind SBP operators were introduced by [`Mattsson2017`](@ref). They combine
two derivative operators `Dp` (`:plus`) and `Dm` (`:minus`) such that
`M * Dp + Dm' * M == tR * tR' - tL * tL'` and `M * (Dp - Dm)` is negative
semidefinite.

```@repl
using SummationByPartsOperators, LinearAlgebra

Dp = derivative_operator(Mattsson2017(:plus), derivative_order=1, accuracy_order=2,
                         xmin=0//1, xmax=1//1, N=9)
Matrix(Dp)
Dm = derivative_operator(Mattsson2017(:minus), derivative_order=1, accuracy_order=2,
                         xmin=0//1, xmax=1//1, N=9)
Matrix(Dm)

M = mass_matrix(Dp)
M * Matrix(Dp) + Matrix(Dm)' * M
minimum(eigvals(-M * (Matrix(Dp) - Matrix(Dm)))) # > 0 up to floating point tolerances
```

You can also set up fully periodic upwind operators by setting the argument
`left_offset` of [`periodic_derivative_operator`](@ref) appropriately. For example,

```@repl
using SummationByPartsOperators, LinearAlgebra

Dp = periodic_derivative_operator(derivative_order=1, accuracy_order=2, left_offset=0,
                                  xmin=0//1, xmax=1//1, N=8)
Matrix(Dp)
Dm = periodic_derivative_operator(derivative_order=1, accuracy_order=2, left_offset=-2,
                                  xmin=0//1, xmax=1//1, N=8)
Matrix(Dm)

M = mass_matrix(Dp)
M * Matrix(Dp) + Matrix(Dm)' * M |> iszero
minimum(eigvals(-M * (Matrix(Dp) - Matrix(Dm)))) # > 0 up to floating point tolerances
```

Note that we used `N=8` here, i.e., one node less than for the non-periodic
example. This is necessary since the additional node at the right boundary is
identified with the left boundary node for periodic operators.

To create all upwind operators for a single setup, you can use
[`upwind_operators`](@ref).

```@repl
using SummationByPartsOperators

D = upwind_operators(Mattsson2017, derivative_order=1, accuracy_order=2,
                     xmin=0, xmax=1//1, N=9)
Matrix(D.plus)
Matrix(D.minus)
```

You can also couple upwind operators continuously across elements using
[`couple_continuously`](@ref) to obtain global upwind operators, see
[below](@ref intro-CGSEM) and Theorem 2.4 of [^RanochaMitsotakisKetcheson2021].

Similarly, you can couple classical and upwind operators discontinuously across
elements using [`couple_discontinuously`](@ref) to obtain global upwind operators,
see [below](@ref intro-DGSEM) and Theorem 2.2 of [^RanochaMitsotakisKetcheson2021].


## [Continuous Galerkin methods](@id intro-CGSEM)

SBP operators can be coupled to obtain (nodal) continuous Galerkin (CG) methods.
If the underlying SBP operators are [`LegendreDerivativeOperator`](@ref)s,
these are CG spectral element methods (CGSEM). However, a continuous coupling
of arbitrary SBP operators is supported.

```@repl
using SummationByPartsOperators, LinearAlgebra

D = couple_continuously(
        legendre_derivative_operator(xmin=-1.0, xmax=1.0, N=3),
        UniformMesh1D(xmin=0.0, xmax=1.0, Nx=3))
Matrix(D)
mass_matrix(D)
```


## [Discontinuous Galerkin methods](@id intro-DGSEM)

SBP operators can also be coupled as in discontinuous Galerkin (DG) methods.
Using a central numerical flux results in central SBP operators; upwind fluxes
yield upwind SBP operators. If [`LegendreDerivativeOperator`](@ref)s are used,
the discontinuous coupling yields DG spectral element methods (DGSEM).

```@repl
using SummationByPartsOperators, LinearAlgebra

D = couple_discontinuously(
        legendre_derivative_operator(xmin=-1.0, xmax=1.0, N=3),
        UniformPeriodicMesh1D(xmin=0.0, xmax=1.0, Nx=3),
        Val(:central))

M = mass_matrix(D);
M * Matrix(D) + Matrix(D)' * M |> iszero
```

Right now, only uniform meshes [`UniformMesh1D`](@ref) and [`UniformPeriodicMesh1D`](@ref)
are implemented.

You can also specify a different coupling than `Val(:central)` to obtain
[upwind operators](@ref intro-upwind-operators).

```@repl
using SummationByPartsOperators, LinearAlgebra

Dp = couple_discontinuously(
        legendre_derivative_operator(xmin=-1.0, xmax=1.0, N=3),
        UniformPeriodicMesh1D(xmin=0.0, xmax=1.0, Nx=3),
        Val(:plus))

Matrix(Dp)

Dm = couple_discontinuously(
        legendre_derivative_operator(xmin=-1.0, xmax=1.0, N=3),
        UniformPeriodicMesh1D(xmin=0.0, xmax=1.0, Nx=3),
        Val(:minus))

Matrix(Dm)
```


## Basic interfaces and additional features

To actually compute and plot the discrete grid functions, a few additional ingredients
are necessary.
- The discrete coefficients of a function on the [`grid`](@ref) of an SBP
  operator can usually be computed as `x = grid(D); u = u_function.(x)`,
  at least for nodal bases. In general, [`compute_coefficients`](@ref)
  (or the in-place version [`compute_coefficients!`](@ref))
  can also be used for this task.
- To get a grid and discrete values suitable for plotting, you can use
  [`evaluate_coefficients`](@ref) (or the in-place version
  [`evaluate_coefficients!`](@ref)).
  The plot nodes returned from [`evaluate_coefficients`](@ref) can be different
  from the nodes of the [`grid`](@ref) associated to an SBP operator.
- To implement boundary procedures, the weights of the mass matrix at the boundary
  are often needed. These can be obtained without forming `M = mass_matrix(D)`
  explicitly via [`left_boundary_weight`](@ref) and [`right_boundary_weight`](@ref).
- Instead of forming a mass matrix explicitly, discrete integrals can be evaluated
  efficiently using [`integrate`](@ref).
- Dissipation operators based on the same discrete inner product as SBP derivative
  operators can be obtained via [`dissipation_operator`](@ref).


## Next steps

If you are familiar with SBP operators in general, this introduction might already
be enough for you to apply
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
to your problems. Otherwise, you might want to have a look at the references,
the tutorials coming next,
or some ready-to-use semidiscretizations of the following partial differential
equations (PDEs). These are shipped with this package and you are encouraged to
look at their source code to learn more about it.

- Linear scalar advection with variable coefficient:
  [`VariableLinearAdvectionNonperiodicSemidiscretization`](@ref), [`VariableLinearAdvectionPeriodicSemidiscretization`](@ref)
- Burgers' equation (inviscid):
  [`BurgersPeriodicSemidiscretization`](@ref), [`BurgersNonperiodicSemidiscretization`](@ref)
- Scalar conservation law with cubic flux:
  [`CubicPeriodicSemidiscretization`](@ref), [`CubicNonperiodicSemidiscretization`](@ref)
- A scalar conservation law with quartic, non-convex flux:
  [`QuarticNonconvexPeriodicSemidiscretization`](@ref)
- The second-order wave equation:
  [`WaveEquationNonperiodicSemidiscretization`](@ref)

Some additional examples are included as [Jupyter](https://jupyter.org) notebooks
in the directory [`notebooks`](https://github.com/ranocha/SummationByPartsOperators.jl/tree/main/notebooks).
Even more examples and research articles making use of
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
are listed in the section [Applications](@ref).
If you want to know even more, you can have a look at the
[test](https://github.com/ranocha/SummationByPartsOperators.jl/tree/main/test).


## References

[^SvärdNordström2014]:
    Svärd, Nordström (2014).
    Review of summation-by-parts schemes for initial–boundary-value problems.
    [DOI: 10.1016/j.jcp.2014.02.031](https://doi.org/10.1016/j.jcp.2014.02.031)

[^FernándezHickenZingg2014]:
    Fernández, Hicken, Zingg (2014).
    Review of summation-by-parts operators with simultaneous approximation terms
    for the numerical solution of partial differential equations.
    [DOI: 10.1016/j.compfluid.2014.02.016](https://doi.org/10.1016/j.compfluid.2014.02.016)

[^RanochaMitsotakisKetcheson2021]:
    Ranocha, Mitsotakis, Ketcheson (2021).
    A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations.
    [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)

