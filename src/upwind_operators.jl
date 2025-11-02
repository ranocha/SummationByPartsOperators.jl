
"""
    UpwindOperators
    UpwindOperators(D_minus, D_central, D_plus)

A struct bundling the individual operators available for non-periodic
upwind SBP operators. The individual operators are available as
`D.minus`, `D.plus` (and optionally `D.central`, if provided), where
`D::UpwindOperators`.

The combined struct behaves as much as possible as an operator itself as long
as no ambiguities arise. For example, upwind operators need to use the same
grid and mass matrix, so [`mass_matrix`](@ref), [`grid`](@ref),
[`xmin`](@ref), [`xmax`](@ref) etc. are available but `mul!` is not.

It is recommended to construct an instance of `UpwindOperators` using
[`upwind_operators`](@ref). An instance can also be constructed manually
by passing the operators in the order `D_minus, D_central, D_plus`.

See also [`upwind_operators`](@ref), [`PeriodicUpwindOperators`](@ref)
"""
@auto_hash_equals struct UpwindOperators{T,
                                         Minus <: AbstractNonperiodicDerivativeOperator{T},
                                         Central <:
                                         AbstractNonperiodicDerivativeOperator{T},
                                         Plus <: AbstractNonperiodicDerivativeOperator{T}} <:
                         AbstractNonperiodicDerivativeOperator{T}
    minus::Minus
    central::Central
    plus::Plus

    function UpwindOperators(minus::Minus, central::Central,
                             plus::Plus) where {T,
                                                Minus <:
                                                AbstractNonperiodicDerivativeOperator{T},
                                                Central <:
                                                AbstractNonperiodicDerivativeOperator{T},
                                                Plus <:
                                                AbstractNonperiodicDerivativeOperator{T}}
        @argcheck grid(minus) == grid(central) == grid(plus)
        @argcheck size(minus) == size(central) == size(plus)
        @argcheck derivative_order(minus) == derivative_order(central) ==
                  derivative_order(plus)
        # Note: We do not check more expensive properties such as
        #       mass_matrix(minus) == mass_matrix(central) == mass_matrix(plus)
        #       central == (minus + plus) / 2
        #       M * (plus - minus) is negative semidefinite
        new{T, Minus, Central, Plus}(minus, central, plus)
    end
end

"""
    PeriodicUpwindOperators
    PeriodicUpwindOperators(D_minus, D_central, D_plus)

A struct bundling the individual operators available for periodic
upwind SBP operators. The individual operators are available as
`D.minus`, `D.plus` (and optionally `D.central`, if provided), where
`D::PeriodicUpwindOperators`.

The combined struct behaves as much as possible as an operator itself as long
as no ambiguities arise. For example, upwind operators need to use the same
grid and mass matrix, so [`mass_matrix`](@ref), [`grid`](@ref),
[`xmin`](@ref), [`xmax`](@ref) etc. are available but `mul!` is not.

It is recommended to construct an instance of `PeriodicUpwindOperators` using
[`upwind_operators`](@ref). An instance can also be constructed manually
by passing the operators in the order `D_minus, D_central, D_plus`.

See also [`upwind_operators`](@ref), [`UpwindOperators`](@ref)
"""
@auto_hash_equals struct PeriodicUpwindOperators{T,
                                                 Minus <:
                                                 AbstractPeriodicDerivativeOperator{T},
                                                 Central <:
                                                 AbstractPeriodicDerivativeOperator{T},
                                                 Plus <:
                                                 AbstractPeriodicDerivativeOperator{T}} <:
                         AbstractPeriodicDerivativeOperator{T}
    minus::Minus
    central::Central
    plus::Plus

    function PeriodicUpwindOperators(minus::Minus, central::Central,
                                     plus::Plus) where {T,
                                                        Minus <:
                                                        AbstractPeriodicDerivativeOperator{T},
                                                        Central <:
                                                        AbstractPeriodicDerivativeOperator{T},
                                                        Plus <:
                                                        AbstractPeriodicDerivativeOperator{T}
                                                        }
        @argcheck grid(minus) == grid(central) == grid(plus)
        @argcheck size(minus) == size(central) == size(plus)
        @argcheck derivative_order(minus) == derivative_order(central) ==
                  derivative_order(plus)
        # Note: We do not check more expensive properties such as
        #       mass_matrix(minus) == mass_matrix(central) == mass_matrix(plus)
        #       central == (minus + plus) / 2
        #       M * (plus - minus) is negative semidefinite
        new{T, Minus, Central, Plus}(minus, central, plus)
    end
end

function Base.show(io::IO, D::Union{UpwindOperators, PeriodicUpwindOperators})
    if derivative_order(D) == 1
        print(io, "Upwind SBP first-derivative operators")
    elseif derivative_order(D) == 2
        print(io, "Upwind SBP second-derivative operators")
    elseif derivative_order(D) == 3
        print(io, "Upwind SBP third-derivative operators")
    else
        print(io, "Upwind SBP ", derivative_order(D),
              "-derivative operators")
    end
    acc = accuracy_order(D.minus), accuracy_order(D.central), accuracy_order(D.minus)
    if all(==(first(acc)), acc)
        print(io, " of order ", first(acc))
    else
        print(io, " of orders ", acc)
    end
    if get(io, :compact, false) == false
        print(io, " on a grid in [", first(grid(D)), ", ", last(grid(D)),
              "] using ", length(grid(D)), " nodes \n")
        print(io, "and coefficients")
    end
    # Assumes that the same source is used for all coefficients
    print(io, " of ", nameof(typeof(source_of_coefficients(D.minus))))
end

function Base.summary(io::IO, D::Union{UpwindOperators, PeriodicUpwindOperators})
    acc = accuracy_order(D.minus), accuracy_order(D.central), accuracy_order(D.minus)
    if all(==(first(acc)), acc)
        acc_string = string(first(acc))
    else
        acc_string = string(acc)
    end
    print(io, nameof(typeof(D)), "(derivative:", derivative_order(D),
          ", accuracy:", acc_string, ")")
end

function derivative_order(D::Union{UpwindOperators, PeriodicUpwindOperators})
    derivative_order(D.minus)
end

grid(D::Union{UpwindOperators, PeriodicUpwindOperators}) = grid(D.minus)
xmin(D::Union{UpwindOperators, PeriodicUpwindOperators}) = xmin(D.minus)
xmax(D::Union{UpwindOperators, PeriodicUpwindOperators}) = xmax(D.minus)

mass_matrix(D::Union{UpwindOperators, PeriodicUpwindOperators}) = mass_matrix(D.minus)
left_boundary_weight(D::UpwindOperators) = left_boundary_weight(D.minus)
right_boundary_weight(D::UpwindOperators) = right_boundary_weight(D.minus)
function integrate(func::Func, u::AbstractVector,
                   D::Union{UpwindOperators, PeriodicUpwindOperators}) where {Func}
    integrate(func, u, D.minus)
end

function scale_by_mass_matrix!(u::AbstractVector,
                               D::Union{UpwindOperators, PeriodicUpwindOperators})
    scale_by_mass_matrix!(u, D.minus)
end

function scale_by_inverse_mass_matrix!(u::AbstractVector,
                                       D::Union{UpwindOperators, PeriodicUpwindOperators})
    scale_by_inverse_mass_matrix!(u, D.minus)
end

"""
    upwind_operators(source_type, args...; derivative_order = 1, kwargs...)

Create [`UpwindOperators`](@ref) from the given source type.
The positional arguments `args...` and keyword arguments `kwargs...` are passed
directly to [`derivative_operator`](@ref).

## Examples

```jldoctest
julia> D = upwind_operators(Mattsson2017, accuracy_order = 2,
                            xmin = 0//1, xmax = 9//1, N = 10)
Upwind SBP first-derivative operators of order 2 on a grid in [0//1, 9//1] using 10 nodes
and coefficients of Mattsson2017

julia> D.minus
SBP first-derivative operator of order 2 on a grid in [0//1, 9//1] using 10 nodes
and coefficients of Mattsson (2017)
  Diagonal-norm upwind SBP operators.
  Journal of Computational Physics 335, pp. 283-310.
  (upwind coefficients minus)

julia> D.plus
SBP first-derivative operator of order 2 on a grid in [0//1, 9//1] using 10 nodes
and coefficients of Mattsson (2017)
  Diagonal-norm upwind SBP operators.
  Journal of Computational Physics 335, pp. 283-310.
  (upwind coefficients plus)
```

```julia
julia> Matrix(D.central)
10×10 Matrix{Rational{Int64}}:
 -2//1   3//1  -1//1   0//1   0//1   0//1   0//1   0//1   0//1   0//1
 -3//5   0//1   4//5  -1//5   0//1   0//1   0//1   0//1   0//1   0//1
  1//4  -1//1   0//1   1//1  -1//4   0//1   0//1   0//1   0//1   0//1
  0//1   1//4  -1//1   0//1   1//1  -1//4   0//1   0//1   0//1   0//1
  0//1   0//1   1//4  -1//1   0//1   1//1  -1//4   0//1   0//1   0//1
  0//1   0//1   0//1   1//4  -1//1   0//1   1//1  -1//4   0//1   0//1
  0//1   0//1   0//1   0//1   1//4  -1//1   0//1   1//1  -1//4   0//1
  0//1   0//1   0//1   0//1   0//1   1//4  -1//1   0//1   1//1  -1//4
  0//1   0//1   0//1   0//1   0//1   0//1   1//5  -4//5   0//1   3//5
  0//1   0//1   0//1   0//1   0//1   0//1   0//1   1//1  -3//1   2//1
```
"""
function upwind_operators(source_type, args...; derivative_order = 1, kwargs...)
    Dm = derivative_operator(source_type(:minus), args...; derivative_order, kwargs...)
    Dc = derivative_operator(source_type(:central), args...; derivative_order, kwargs...)
    Dp = derivative_operator(source_type(:plus), args...; derivative_order, kwargs...)
    return UpwindOperators(Dm, Dc, Dp)
end

"""
    upwind_operators(periodic_derivative_operator;
                     derivative_order = 1, accuracy_order,
                     xmin, xmax, N,
                     mode = FastMode()))

Create [`PeriodicUpwindOperators`](@ref) from operators constructed by
[`periodic_derivative_operator`](@ref). The keyword arguments are passed
directly to [`periodic_derivative_operator`](@ref).

## Examples

```jldoctest
julia> D = upwind_operators(periodic_derivative_operator, accuracy_order = 2,
                            xmin = 0//1, xmax = 8//1, N = 8)
Upwind SBP first-derivative operators of order 2 on a grid in [0//1, 7//1] using 8 nodes
and coefficients of Fornberg1998

julia> D.minus
Periodic first-derivative operator of order 2 on a grid in [0//1, 8//1] using 8 nodes,
stencils with 2 nodes to the left, 0 nodes to the right, and coefficients of Fornberg (1998)
  Calculation of Weights in Finite Difference Formulas.
  SIAM Rev. 40.3, pp. 685-691.

julia> D.plus
Periodic first-derivative operator of order 2 on a grid in [0//1, 8//1] using 8 nodes,
stencils with 0 nodes to the left, 2 nodes to the right, and coefficients of Fornberg (1998)
  Calculation of Weights in Finite Difference Formulas.
  SIAM Rev. 40.3, pp. 685-691.
```

```julia
julia> Matrix(D.central)
8×8 Matrix{Rational{Int64}}:
  0//1   1//1  -1//4   0//1   0//1   0//1   1//4  -1//1
 -1//1   0//1   1//1  -1//4   0//1   0//1   0//1   1//4
  1//4  -1//1   0//1   1//1  -1//4   0//1   0//1   0//1
  0//1   1//4  -1//1   0//1   1//1  -1//4   0//1   0//1
  0//1   0//1   1//4  -1//1   0//1   1//1  -1//4   0//1
  0//1   0//1   0//1   1//4  -1//1   0//1   1//1  -1//4
 -1//4   0//1   0//1   0//1   1//4  -1//1   0//1   1//1
  1//1  -1//4   0//1   0//1   0//1   1//4  -1//1   0//1
```
"""
function upwind_operators(::typeof(periodic_derivative_operator);
                          derivative_order = 1, accuracy_order,
                          xmin, xmax, N,
                          mode = FastMode())
    Dm = periodic_derivative_operator(; derivative_order, accuracy_order,
                                      left_offset = -(accuracy_order ÷ 2 + 1),
                                      xmin, xmax, N, mode)
    Dp = periodic_derivative_operator(; derivative_order, accuracy_order,
                                      left_offset = -(accuracy_order - 1) ÷ 2,
                                      xmin, xmax, N, mode)
    # central coefficients obtained by averaging
    upper_coef_central = widening_plus(Dm.coefficients.upper_coef,
                                       Dp.coefficients.upper_coef) / 2
    central_coef_central = (Dm.coefficients.central_coef + Dp.coefficients.central_coef) / 2
    lower_coef_central = widening_plus(Dm.coefficients.lower_coef,
                                       Dp.coefficients.lower_coef) / 2
    coef_central = PeriodicDerivativeCoefficients(lower_coef_central, central_coef_central,
                                                  upper_coef_central,
                                                  mode, derivative_order, accuracy_order,
                                                  source_of_coefficients(Dm))
    Dc = PeriodicDerivativeOperator(coef_central, Dm.grid_evaluate)
    return PeriodicUpwindOperators(Dm, Dc, Dp)
end

"""
    upwind_operators(couple_discontinuously, D, mesh)

Create [`UpwindOperators`](@ref) (for a non-periodic `mesh`) or
[`PeriodicUpwindOperators`](@ref) (for a periodic `mesh`) from
an existing first-derivative SBP operator `D`. The individual operators
are constructed using [`couple_discontinuously`](@ref). If `D` is a
[`LegendreDerivativeOperator`](@ref) (constructed using
[`legendre_derivative_operator`](@ref)), the resulting upwind operators
are local discontinuous Galerkin (LDG) operators.

## References

- Ranocha, Mitsotakis, Ketcheson (2021).
  A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations.
  [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)

## Examples

```jldoctest
julia> D = upwind_operators(couple_discontinuously,
                            legendre_derivative_operator(xmin = -1.0, xmax = 1.0, N = 3),
                            UniformMesh1D(xmin = 0.0, xmax = 1.0, Nx = 4))
Upwind SBP first-derivative operators of order 2 on a grid in [5.0e-324, 0.9999999999999999] using 12 nodes
and coefficients of Module
```

```julia
julia> sparse(D.minus)
12×12 SparseArrays.SparseMatrixCSC{Float64, Int64} with 35 stored entries:
 -12.0   16.0   -4.0    ⋅      ⋅      ⋅     ⋅      ⋅      ⋅     ⋅      ⋅     ⋅
  -4.0     ⋅     4.0    ⋅      ⋅      ⋅     ⋅      ⋅      ⋅     ⋅      ⋅     ⋅
   4.0  -16.0   12.0    ⋅      ⋅      ⋅     ⋅      ⋅      ⋅     ⋅      ⋅     ⋅
    ⋅      ⋅   -24.0  12.0   16.0   -4.0    ⋅      ⋅      ⋅     ⋅      ⋅     ⋅
    ⋅      ⋅      ⋅   -4.0     ⋅     4.0    ⋅      ⋅      ⋅     ⋅      ⋅     ⋅
    ⋅      ⋅      ⋅    4.0  -16.0   12.0    ⋅      ⋅      ⋅     ⋅      ⋅     ⋅
    ⋅      ⋅      ⋅     ⋅      ⋅   -24.0  12.0   16.0   -4.0    ⋅      ⋅     ⋅
    ⋅      ⋅      ⋅     ⋅      ⋅      ⋅   -4.0     ⋅     4.0    ⋅      ⋅     ⋅
    ⋅      ⋅      ⋅     ⋅      ⋅      ⋅    4.0  -16.0   12.0    ⋅      ⋅     ⋅
    ⋅      ⋅      ⋅     ⋅      ⋅      ⋅     ⋅      ⋅   -24.0  12.0   16.0  -4.0
    ⋅      ⋅      ⋅     ⋅      ⋅      ⋅     ⋅      ⋅      ⋅   -4.0     ⋅    4.0
    ⋅      ⋅      ⋅     ⋅      ⋅      ⋅     ⋅      ⋅      ⋅    4.0  -16.0  12.0

julia> sparse(D.plus)
12×12 SparseArrays.SparseMatrixCSC{Float64, Int64} with 35 stored entries:
 -12.0   16.0   -4.0     ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅     ⋅
  -4.0     ⋅     4.0     ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅     ⋅
   4.0  -16.0  -12.0   24.0     ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅     ⋅
    ⋅      ⋅      ⋅   -12.0   16.0   -4.0     ⋅      ⋅      ⋅      ⋅      ⋅     ⋅
    ⋅      ⋅      ⋅    -4.0     ⋅     4.0     ⋅      ⋅      ⋅      ⋅      ⋅     ⋅
    ⋅      ⋅      ⋅     4.0  -16.0  -12.0   24.0     ⋅      ⋅      ⋅      ⋅     ⋅
    ⋅      ⋅      ⋅      ⋅      ⋅      ⋅   -12.0   16.0   -4.0     ⋅      ⋅     ⋅
    ⋅      ⋅      ⋅      ⋅      ⋅      ⋅    -4.0     ⋅     4.0     ⋅      ⋅     ⋅
    ⋅      ⋅      ⋅      ⋅      ⋅      ⋅     4.0  -16.0  -12.0   24.0     ⋅     ⋅
    ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅   -12.0   16.0  -4.0
    ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅    -4.0     ⋅    4.0
    ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅      ⋅     4.0  -16.0  12.0
```
"""
function upwind_operators(::typeof(couple_discontinuously),
                          D::AbstractNonperiodicDerivativeOperator,
                          mesh::AbstractMesh1D)
    @assert derivative_order(D) == 1
    Dp = couple_discontinuously(D, mesh, Val{:plus}())
    Dc = couple_discontinuously(D, mesh, Val{:central}())
    Dm = couple_discontinuously(D, mesh, Val{:minus}())
    if isperiodic(mesh)
        return PeriodicUpwindOperators(Dm, Dc, Dp)
    else
        return UpwindOperators(Dm, Dc, Dp)
    end
end
