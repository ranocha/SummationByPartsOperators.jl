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
                                        xmin=0.0, xmax=2.0, N=21)
Periodic 1st derivative operator of order 2 {T=Float64, Parallel=Val{:serial}}
on a grid in [0.0, 2.0] using 21 nodes,
stencils with 1 nodes to the left, 1 nodes to the right, and coefficients from
  Fornberg (1998)
  Calculation of Weights in Finite Difference Formulas.
  SIAM Rev. 40.3, pp. 685-691.


julia> M = mass_matrix(D)
UniformScaling{Float64}
0.1*I

julia> M * Matrix(D) + Matrix(D)' * M |> norm
0.0

julia> D = fourier_derivative_operator(xmin=0.0, xmax=2.0, N=20)
Periodic 1st derivative Fourier operator {T=Float64}
on a grid in [0.0, 2.0] using 20 nodes and 11 modes.


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
= u(x_\mathrm{max}) v(x_\mathrm{max}) - u(x_\mathrm{min}) - v(x_\mathrm{min}).
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

```jldoctest
julia> using SummationByPartsOperators, LinearAlgebra

julia> D = derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=2,
                               xmin=0//1, xmax=1//1, N=9)
SBP 1st derivative operator of order 2 {T=Rational{Int64}, Parallel=Val{:serial}}
on a grid in [0//1, 1//1] using 9 nodes
and coefficients given in
  Mattsson, Nordström (2004)
  Summation by parts operators for finite difference approximations of second
    derivatives.
  Journal of Computational Physics 199, pp. 503-540.


julia> tL = zeros(eltype(D), size(D, 1)); tL[1] = 1; tL'
1×9 adjoint(::Vector{Rational{Int64}}) with eltype Rational{Int64}:
 1//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1

julia> tR = zeros(eltype(D), size(D, 1)); tR[end] = 1; tR'
1×9 adjoint(::Vector{Rational{Int64}}) with eltype Rational{Int64}:
 0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1  1//1

julia> M = mass_matrix(D)
9×9 Diagonal{Rational{Int64}, Vector{Rational{Int64}}}:
 1//16   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅
  ⋅     1//8   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅
  ⋅      ⋅    1//8   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅
  ⋅      ⋅     ⋅    1//8   ⋅     ⋅     ⋅     ⋅     ⋅
  ⋅      ⋅     ⋅     ⋅    1//8   ⋅     ⋅     ⋅     ⋅
  ⋅      ⋅     ⋅     ⋅     ⋅    1//8   ⋅     ⋅     ⋅
  ⋅      ⋅     ⋅     ⋅     ⋅     ⋅    1//8   ⋅     ⋅
  ⋅      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    1//8   ⋅
  ⋅      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    1//16

julia> M * Matrix(D) + Matrix(D)' * M == tR * tR' - tL * tL'
true

julia> u = randn(size(grid(D))); derivative_left(D, u, Val(0)) == u[begin]
true

julia> u = randn(size(grid(D))); derivative_right(D, u, Val(0)) == u[end]
true
```

Here, we have introduced two additional features. Firstly, exact rational
coefficients are provided, based on the type of `xmin` and `xmax` (if available).
Secondly, a [`source_of_coefficients`](@ref) has to be provided when constructing
the SBP operator. You can list them using
```@example
subtypes(SourceOfCoefficients)
```

### Second-derivative operators

To mimic integration-by-parts of second derivatives, the evaluation of the first
derivative at the boundaries is necessary. These linear functionals are available
as [`derivative_left`](@ref) and [`derivative_right`](@ref). For example,
```jldoctest
julia> using SummationByPartsOperators, LinearAlgebra

julia> D = derivative_operator(MattssonNordström2004(), derivative_order=2, accuracy_order=4,
                               xmin=0//1, xmax=1//1, N=9)
SBP 2nd derivative operator of order 4 {T=Rational{Int64}, Parallel=Val{:serial}}
on a grid in [0//1, 1//1] using 9 nodes
and coefficients given in
  Mattsson, Nordström (2004)
  Summation by parts operators for finite difference approximations of second
    derivatives.
  Journal of Computational Physics 199, pp. 503-540.


julia> derivative_left(D, grid(D), Val(1))
1//1
```
This is correct since `grid(D)` represents the identity function on the domain,
whose derivative is unity everywhere. Some procedures imposing boundary conditions
weakly require adding the transposed boundary derivatives to a grid function,
which can be achieved by [`add_transpose_derivative_left!`](@ref) and
[`add_transpose_derivative_right!`](@ref).
You can find applications of these operators in the source code of
[`WaveEquationNonperiodicSemidiscretisation`](@ref).


## Upwind operators

TODO


## Basic interfaces

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

