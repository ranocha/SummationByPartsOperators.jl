
"""
    VarCoefDerivativeOperator

A dissipation operator on a nonperiodic finite difference grid.
"""
struct VarCoefDerivativeOperator{T,CoefficientCache,Parallel,SourceOfCoefficients,Grid} <: AbstractVariableCoefficientNonperiodicDerivativeOperator{T}
    coefficients::VarCoefDerivativeCoefficients{T,CoefficientCache,Parallel,SourceOfCoefficients}
    grid::Grid
    Δx::T
    factor::T
    b::Vector{T}

    function VarCoefDerivativeOperator(coefficients::VarCoefDerivativeCoefficients{T,CoefficientCache,Parallel,SourceOfCoefficients},
                                       grid::Grid, b::Vector{T}) where {T,CoefficientCache,Parallel,SourceOfCoefficients,Grid}
        @argcheck checkbounds(Bool, grid, coefficients.coefficient_cache) DimensionMismatch
        @argcheck length(grid) == length(b)

        Δx = step(grid)
        factor = inv(Δx^derivative_order(coefficients))
        new{T,CoefficientCache,Parallel,SourceOfCoefficients,Grid}(
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
    print(io, accuracy_order(D), " {T=", T, ", Parallel=", typeof(D.coefficients.parallel), "} \n")
    print(io, "on a grid in [", first(grid(D)), ", ", last(grid(D)),
                "] using ", length(grid(D)), " nodes \n")
    print(io, "and coefficients given in \n")
    print(io, source_of_coeffcients(D))
end


"""
    mul!(dest::AbstractVector, D::VarCoefDerivativeOperator, u::AbstractVector, α, β)

Compute `α*D*u + β*dest` and store the result in `dest`.
"""
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::VarCoefDerivativeOperator,
                                        u::AbstractVector, α, β)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, D.b, D.factor*α, β)
end

"""
    mul!(dest::AbstractVector, D::VarCoefDerivativeOperator, u::AbstractVector, α)

Compute `α*D*u` and store the result in `dest`.
"""
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::VarCoefDerivativeOperator,
                                        u::AbstractVector, α)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, D.b, D.factor*α)
end


"""
    var_coef_derivative_operator(source_of_coefficients, derivative_order, accuracy_order, xmin, xmax, N, left_weights, right_weights, bfunc, parallel=Val{:serial}())

Create a `VarCoefDerivativeOperator` approximating a `derivative_order`-th
derivative with variable coefficients `bfunc` on a grid between `xmin` and
`xmax` with `N` grid points up to order of accuracy `accuracy_order` with
coefficients given by `source_of_coefficients`.
The evaluation of the derivative can be parallised using threads by chosing
`parallel=Val{:threads}())`.
"""
function var_coef_derivative_operator(source_of_coefficients, derivative_order, accuracy_order,
                                        xmin, xmax, N, bfunc, parallel=Val{:serial}())
    grid = construct_grid(source_of_coefficients, accuracy_order, xmin, xmax, N)
    coefficients = var_coef_derivative_coefficients(source_of_coefficients, derivative_order, accuracy_order, grid, parallel)
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

    @inbounds @simd for i in eachindex(left_weights)
        u[i] = factor * u[i] * (Δx * left_weights[i])
    end
    @inbounds @simd for i in length(left_weights)+1:N-length(right_weights)
        u[i] = factor * u[i] * Δx
    end
    @inbounds @simd for i in eachindex(right_weights)
        u[end-i+1] = factor * u[end-i+1] * (Δx * right_weights[i])
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

    @inbounds @simd for i in eachindex(left_weights)
        u[i] = factor * u[i] / (Δx * left_weights[i])
    end
    @inbounds @simd for i in length(left_weights)+1:N-length(right_weights)
        u[i] = factor * u[i] / Δx
    end
    @inbounds @simd for i in eachindex(right_weights)
        u[end-i+1] = factor * u[end-i+1] / (Δx * right_weights[i])
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
