
"""
    VarCoefDerivativeOperator

A derivative operator with variable coefficients
on a nonperiodic finite difference grid.
See [`var_coef_derivative_operator`](@ref).
"""
@auto_hash_equals struct VarCoefDerivativeOperator{T,
                                                   Coefficients <:
                                                   VarCoefDerivativeCoefficients,
                                                   Grid} <:
                         AbstractVariableCoefficientNonperiodicDerivativeOperator{T}
    coefficients::Coefficients
    grid::Grid
    Δx::T
    factor::T
    b::Vector{T}

    function VarCoefDerivativeOperator(coefficients::Coefficients,
                                       grid::Grid,
                                       b::Vector{T}) where {T,
                                                            Coefficients <:
                                                            VarCoefDerivativeCoefficients,
                                                            Grid}
        @argcheck checkbounds(Bool, grid, coefficients.coefficient_cache) DimensionMismatch
        @argcheck length(grid) == length(b)

        Δx = step(grid)
        factor = inv(Δx^derivative_order(coefficients))
        new{T, Coefficients, Grid}(coefficients, grid, Δx, factor, b)
    end
end

function Base.show(io::IO, D::VarCoefDerivativeOperator{T}) where {T}
    print(io, "SBP variable coefficient ")
    if derivative_order(D) == 1
        print(io, "1st")
    elseif derivative_order(D) == 2
        print(io, "2nd")
    elseif derivative_order(D) == 3
        print(io, "3rd")
    else
        print(io, derivative_order(D), "th")
    end
    print(io, "derivative operator with order of accuracy ")
    print(io, accuracy_order(D), " {T=", T, ", ExecutionMode=", typeof(D.coefficients.mode),
          "} \n")
    print(io, "on a grid in [", first(grid(D)), ", ", last(grid(D)),
          "] using ", length(grid(D)), " nodes \n")
    print(io, "and coefficients given in \n")
    print(io, source_of_coefficients(D))
end

# Compute `α*D*u + β*dest` and store the result in `dest`.
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::VarCoefDerivativeOperator,
                                       u::AbstractVector, α, β)
    @boundscheck begin
        @argcheck size(D, 2)==length(u) DimensionMismatch
        @argcheck size(D, 1)==length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, D.b, D.factor * α, β)
end

# Compute `α*D*u` and store the result in `dest`.
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::VarCoefDerivativeOperator,
                                       u::AbstractVector, α)
    @boundscheck begin
        @argcheck size(D, 2)==length(u) DimensionMismatch
        @argcheck size(D, 1)==length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, D.b, D.factor * α)
end

"""
    var_coef_derivative_operator(source_of_coefficients, derivative_order, accuracy_order,
                                 xmin, xmax, N, bfunc,
                                 mode = FastMode())

Create a `VarCoefDerivativeOperator` approximating a `derivative_order`-th
derivative with variable coefficients `bfunc` on a grid between `xmin` and
`xmax` with `N` grid points up to order of accuracy `accuracy_order` with
coefficients given by `source_of_coefficients`.
The evaluation of the derivative can be parallelized using threads by choosing
`mode = ThreadedMode()`.

You can modify the variable coefficient provided by `bfunc` also after creating
the operator by modifying the field `b` of the returned operator. This allows
a manual evaluation of the variable coefficients on the [`grid`](@ref) of the
operator.

# Examples

```jldoctest
julia> D0 = var_coef_derivative_operator(Mattsson2012(), 2, 2, -1.0, 1.0, 11, zero);

julia> Matrix(D0)
11×11 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> D0.b .= grid(D0).^2;

julia> Matrix(D0)
11×11 Matrix{Float64}:
 34.0  -59.0   25.0   0.0   0.0   0.0   0.0   0.0    0.0    0.0   0.0
 20.5  -33.0   12.5   0.0   0.0   0.0   0.0   0.0    0.0    0.0   0.0
  0.0   12.5  -19.0   6.5   0.0   0.0   0.0   0.0    0.0    0.0   0.0
  0.0    0.0    6.5  -9.0   2.5   0.0   0.0   0.0    0.0    0.0   0.0
  0.0    0.0    0.0   2.5  -3.0   0.5   0.0   0.0    0.0    0.0   0.0
  0.0    0.0    0.0   0.0   0.5  -1.0   0.5   0.0    0.0    0.0   0.0
  0.0    0.0    0.0   0.0   0.0   0.5  -3.0   2.5    0.0    0.0   0.0
  0.0    0.0    0.0   0.0   0.0   0.0   2.5  -9.0    6.5    0.0   0.0
  0.0    0.0    0.0   0.0   0.0   0.0   0.0   6.5  -19.0   12.5   0.0
  0.0    0.0    0.0   0.0   0.0   0.0   0.0   0.0   12.5  -33.0  20.5
  0.0    0.0    0.0   0.0   0.0   0.0   0.0   0.0   25.0  -59.0  34.0

julia> D2 = var_coef_derivative_operator(Mattsson2012(), 2, 2, -1.0, 1.0, 11, abs2);

julia> Matrix(D2) ≈ Matrix(D0)
true
```
"""
function var_coef_derivative_operator(source_of_coefficients, derivative_order,
                                      accuracy_order,
                                      xmin, xmax, N, bfunc, mode = FastMode())
    if mode === Val(:serial) || mode === Val(:threads)
        # TODO: deprecated in v0.5
        Base.depwarn("Providing the argument `parallel` is deprecated." *
                     "Use `mode` instead.", :var_coef_derivative_operator)
        mode = _parallel_to_mode(mode)
    end
    grid = construct_grid(source_of_coefficients, accuracy_order, xmin, xmax, N)
    coefficients = var_coef_derivative_coefficients(source_of_coefficients,
                                                    derivative_order, accuracy_order, grid,
                                                    mode)
    VarCoefDerivativeOperator(coefficients, grid, bfunc.(grid))
end

"""
    mass_matrix(D::Union{DerivativeOperator,VarCoefDerivativeOperator})

Create the diagonal mass matrix for the SBP derivative operator `D`.
"""
function mass_matrix(D::Union{DerivativeOperator, VarCoefDerivativeOperator})
    m = fill(one(eltype(D)), length(grid(D)))
    @unpack left_weights, right_weights = D.coefficients

    m[1:length(left_weights)] = left_weights
    m[end:-1:(end - length(right_weights) + 1)] = right_weights
    Diagonal(D.Δx * m)
end

function scale_by_mass_matrix!(u::AbstractVector,
                               D::Union{DerivativeOperator, VarCoefDerivativeOperator},
                               factor = true)
    Base.require_one_based_indexing(u)
    N = size(D, 1)
    @boundscheck begin
        length(u) == N ||
            throw(DimensionMismatch("sizes of input vector and operator do not match"))
    end
    @unpack Δx = D
    @unpack left_weights, right_weights = D.coefficients

    @simd for i in eachindex(left_weights)
        @inbounds u[i] = factor * u[i] * (Δx * left_weights[i])
    end
    @simd for i in (length(left_weights) + 1):(N - length(right_weights))
        @inbounds u[i] = factor * u[i] * Δx
    end
    @simd for i in eachindex(right_weights)
        @inbounds u[end - i + 1] = factor * u[end - i + 1] * (Δx * right_weights[i])
    end

    return u
end

function scale_by_inverse_mass_matrix!(u::AbstractVector,
                                       D::Union{DerivativeOperator,
                                                VarCoefDerivativeOperator}, factor = true)
    Base.require_one_based_indexing(u)
    N = size(D, 1)
    @boundscheck begin
        length(u) == N ||
            throw(DimensionMismatch("sizes of input vector and operator do not match"))
    end
    @unpack Δx = D
    @unpack left_weights, right_weights = D.coefficients

    @simd for i in eachindex(left_weights)
        @inbounds u[i] = factor * u[i] / (Δx * left_weights[i])
    end
    @simd for i in (length(left_weights) + 1):(N - length(right_weights))
        @inbounds u[i] = factor * u[i] / Δx
    end
    @simd for i in eachindex(right_weights)
        @inbounds u[end - i + 1] = factor * u[end - i + 1] / (Δx * right_weights[i])
    end

    u
end

function get_weight(D::Union{DerivativeOperator, VarCoefDerivativeOperator}, i::Int)
    @unpack Δx = D
    @unpack left_weights, right_weights = D.coefficients
    N, _ = size(D)
    @boundscheck begin
        @argcheck 1 <= i <= N
    end
    @inbounds if i <= length(left_weights)
        ω = Δx * left_weights[i]
    elseif i > N - length(right_weights)
        ω = Δx * right_weights[N - i + 1]
    else
        ω = Δx
    end
    ω
end
