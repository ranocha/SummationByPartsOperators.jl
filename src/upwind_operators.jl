
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
"""
@auto_hash_equals struct UpwindOperators{T, Minus   <: AbstractNonperiodicDerivativeOperator{T},
                                            Central <: AbstractNonperiodicDerivativeOperator{T},
                                            Plus    <: AbstractNonperiodicDerivativeOperator{T}} <: AbstractNonperiodicDerivativeOperator{T}
  minus::Minus
  central::Central
  plus::Plus

  function UpwindOperators(minus::Minus, central::Central, plus::Plus) where {T,
                                                                              Minus   <: AbstractNonperiodicDerivativeOperator{T},
                                                                              Central <: AbstractNonperiodicDerivativeOperator{T},
                                                                              Plus    <: AbstractNonperiodicDerivativeOperator{T}}
    @argcheck grid(minus) == grid(central) == grid(plus)
    @argcheck size(minus) == size(central) == size(plus)
    @argcheck derivative_order(minus) == derivative_order(central) == derivative_order(plus)
    # Note: We do not check more expensive properties such as
    #       mass_matrix(minus) == mass_matrix(central) == mass_matrix(plus)
    #       central == (minus + plus) / 2
    #       M * (plus - minus) is negative semidefinite
    new{T, Minus, Central, Plus}(minus, central, plus)
  end
end

function Base.show(io::IO, D::UpwindOperators)
  if derivative_order(D) == 1
      print(io, "Upwind SBP first-derivative operators")
  elseif  derivative_order(D) == 2
      print(io, "Upwind SBP second-derivative operators")
  elseif  derivative_order(D) == 3
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
  # TODO: Assumes that the same source is used for all coefficients
  print(io, " of ", nameof(typeof(source_of_coefficients(D.minus))))
end

function Base.summary(io::IO, D::UpwindOperators)
  acc = accuracy_order(D.minus), accuracy_order(D.central), accuracy_order(D.minus)
  if all(==(first(acc)), acc)
    acc_string = string(first(acc))
  else
    acc_string = string(acc)
  end
  print(io, nameof(typeof(D)), "(derivative:", derivative_order(D),
            ", accuracy:", acc_string, ")")
end

derivative_order(D::UpwindOperators) = derivative_order(D.minus)

grid(D::UpwindOperators) = grid(D.minus)
xmin(D::UpwindOperators) = xmin(D.minus)
xmax(D::UpwindOperators) = xmax(D.minus)

mass_matrix(D::UpwindOperators) = mass_matrix(D.minus)
left_boundary_weight(D::UpwindOperators) = left_boundary_weight(D.minus)
right_boundary_weight(D::UpwindOperators) = right_boundary_weight(D.minus)
function integrate(func, u::AbstractVector, D::UpwindOperators)
  integrate(func, u, D.minus)
end


"""
    upwind_operators(source_type, args...; kwargs...)

Create [`UpwindOperators`](@ref) from the given source type.
The positional arguments `args...` and keyword arguments `kwargs...` are passed
directly to [`derivative_operator`](@ref).

## Examples

```jldoctest
julia> D = upwind_operators(Mattsson2017, derivative_order=1, accuracy_order=2,
                            xmin=0//1, xmax=9//1, N=10)
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
function upwind_operators(source_type, args...; kwargs...)
  Dm = derivative_operator(source_type(:minus),   args...; kwargs...)
  Dc = derivative_operator(source_type(:central), args...; kwargs...)
  Dp = derivative_operator(source_type(:plus),    args...; kwargs...)
  return UpwindOperators(Dm, Dc, Dp)
end

