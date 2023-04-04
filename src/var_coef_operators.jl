
"""
    VarCoefDerivativeOperator

A dissipation operator on a nonperiodic finite difference grid.
"""
@auto_hash_equals struct VarCoefDerivativeOperator{T,CoefficientCache,ExecutionMode,SourceOfCoefficients,Grid} <: AbstractVariableCoefficientNonperiodicDerivativeOperator{T}
    coefficients::VarCoefDerivativeCoefficients{T,CoefficientCache,ExecutionMode,SourceOfCoefficients}
    grid::Grid
    Δx::T
    factor::T
    b::Vector{T}

    function VarCoefDerivativeOperator(coefficients::VarCoefDerivativeCoefficients{T,CoefficientCache,ExecutionMode,SourceOfCoefficients},
                                       grid::Grid, b::Vector{T}) where {T,CoefficientCache,ExecutionMode,SourceOfCoefficients,Grid}
        @argcheck checkbounds(Bool, grid, coefficients.coefficient_cache) DimensionMismatch
        @argcheck length(grid) == length(b)

        Δx = step(grid)
        factor = inv(Δx^derivative_order(coefficients))
        new{T,CoefficientCache,ExecutionMode,SourceOfCoefficients,Grid}(
            coefficients, grid, Δx, factor, b)
    end
end


function Base.show(io::IO, D::VarCoefDerivativeOperator{T}) where {T}
    print(io, "SBP variable coefficient ")
    if  derivative_order(D) == 1
        print(io, "1st")
    elseif  derivative_order(D) == 2
        print(io, "2nd")
    elseif  derivative_order(D) == 3
        print(io, "3rd")
    else
        print(io, derivative_order(D), "th")
    end
    print(io, "derivative operator with order of accuracy ")
    print(io, accuracy_order(D), " {T=", T, ", ExecutionMode=", typeof(D.coefficients.mode), "} \n")
    print(io, "on a grid in [", first(grid(D)), ", ", last(grid(D)),
                "] using ", length(grid(D)), " nodes \n")
    print(io, "and coefficients given in \n")
    print(io, source_of_coefficients(D))
end



# Compute `α*D*u + β*dest` and store the result in `dest`.
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::VarCoefDerivativeOperator,
                                        u::AbstractVector, α, β)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, D.b, D.factor*α, β)
end

# Compute `α*D*u` and store the result in `dest`.
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::VarCoefDerivativeOperator,
                                        u::AbstractVector, α)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, D.b, D.factor*α)
end


"""
    var_coef_derivative_operator(source_of_coefficients, derivative_order, accuracy_order,
                                 xmin, xmax, N, left_weights, right_weights, bfunc,
                                 mode=FastMode())

Create a `VarCoefDerivativeOperator` approximating a `derivative_order`-th
derivative with variable coefficients `bfunc` on a grid between `xmin` and
`xmax` with `N` grid points up to order of accuracy `accuracy_order` with
coefficients given by `source_of_coefficients`.
The evaluation of the derivative can be parallelized using threads by choosing
`mode=ThreadedMode()`.
"""
function var_coef_derivative_operator(source_of_coefficients, derivative_order, accuracy_order,
                                        xmin, xmax, N, bfunc, mode=FastMode())
    if mode === Val(:serial) || mode === Val(:threads)
        # TODO: deprecated in v0.5
        Base.depwarn("Providing the argument `parallel` is deprecated." *
                     "Use `mode` instead.", :var_coef_derivative_operator)
        mode = _parallel_to_mode(mode)
    end
    grid = construct_grid(source_of_coefficients, accuracy_order, xmin, xmax, N)
    coefficients = var_coef_derivative_coefficients(source_of_coefficients,
        derivative_order, accuracy_order, grid, mode)
    VarCoefDerivativeOperator(coefficients, grid, bfunc.(grid))
end



"""
    mass_matrix(D::Union{DerivativeOperator,VarCoefDerivativeOperator})

Create the diagonal mass matrix for the SBP derivative operator `D`.
"""
function mass_matrix(D::Union{DerivativeOperator,VarCoefDerivativeOperator})
    m = fill(one(eltype(D)), length(grid(D)))
    @unpack left_weights, right_weights = D.coefficients

    m[1:length(left_weights)] = left_weights
    m[end:-1:end-length(right_weights)+1] = right_weights
    Diagonal(D.Δx * m)
end

function scale_by_mass_matrix!(u::AbstractVector, D::Union{DerivativeOperator,VarCoefDerivativeOperator}, factor=true)
    @unpack Δx = D
    @unpack left_weights, right_weights = D.coefficients
    N, _ = size(D)
    @boundscheck begin
        @argcheck N == length(u)
    end

    @simd for i in eachindex(left_weights)
        @inbounds u[i] = factor * u[i] * (Δx * left_weights[i])
    end
    @simd for i in length(left_weights)+1:N-length(right_weights)
        @inbounds u[i] = factor * u[i] * Δx
    end
    @simd for i in eachindex(right_weights)
        @inbounds u[end-i+1] = factor * u[end-i+1] * (Δx * right_weights[i])
    end

    u
end

function scale_by_inverse_mass_matrix!(u::AbstractVector, D::Union{DerivativeOperator,VarCoefDerivativeOperator}, factor=true)
    @unpack Δx = D
    @unpack left_weights, right_weights = D.coefficients
    N, _ = size(D)
    @boundscheck begin
        @argcheck N == length(u)
    end

    @simd for i in eachindex(left_weights)
        @inbounds u[i] = factor * u[i] / (Δx * left_weights[i])
    end
    @simd for i in length(left_weights)+1:N-length(right_weights)
        @inbounds u[i] = factor * u[i] / Δx
    end
    @simd for i in eachindex(right_weights)
        @inbounds u[end-i+1] = factor * u[end-i+1] / (Δx * right_weights[i])
    end

    u
end

function get_weight(D::Union{DerivativeOperator,VarCoefDerivativeOperator}, i::Int)
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
