
abstract type AbstractCoefficientCache{T} end

function Base.checkbounds(::Type{Bool}, u::AbstractVector, cache::AbstractCoefficientCache)
    length(u) > length(cache.inv_left_weights)+length(cache.inv_right_weights)
end

"""
    VarCoefDerivativeCoefficients

The coefficients of a variable coefficient derivative operator on a nonperiodic grid.
"""
struct VarCoefDerivativeCoefficients{T,CoefficientCache<:AbstractCoefficientCache{T},
                                     LeftWidth,RightWidth,Parallel,SourceOfCoefficients} <: AbstractDerivativeCoefficients{T}
    # coefficients defining the operator and its action
    coefficient_cache::CoefficientCache
    left_weights::SVector{LeftWidth, T}
    right_weights::SVector{RightWidth, T}
    parallel::Parallel
    # corresponding orders etc.
    derivative_order::Int
    accuracy_order  ::Int
    source_of_coeffcients::SourceOfCoefficients

    function VarCoefDerivativeCoefficients(coefficient_cache::CoefficientCache,
                                           left_weights::SVector{LeftWidth, T},
                                           right_weights::SVector{RightWidth, T},
                                           parallel::Parallel,
                                           derivative_order::Int, accuracy_order::Int,
                                           source_of_coeffcients::SourceOfCoefficients) where {T,CoefficientCache<:AbstractCoefficientCache{T},
                                                                                                LeftWidth,RightWidth,Parallel,SourceOfCoefficients}
        new{T,CoefficientCache,LeftWidth,RightWidth,Parallel,SourceOfCoefficients}(
            coefficient_cache, left_weights, right_weights, parallel, derivative_order, accuracy_order, source_of_coeffcients)
    end
end


function Base.show(io::IO, coefficients::VarCoefDerivativeCoefficients)
    print(io, "Coefficients of the variable coefficient ")
    if  derivative_order(coefficients) == 1
        print(io, "1st")
    elseif  derivative_order(coefficients) == 2
        print(io, "2nd")
    elseif  derivative_order(coefficients) == 3
        print(io, "3rd")
    else
        print(io, derivative_order(coefficients), "th")
    end
    print(io, "derivative operator with order of accuracy ",
            accuracy_order(coefficients), " given in \n")
    print(io, source_of_coeffcients(coefficients))
end


"""
    mul!(dest::AbstractVector, coefficients::VarCoefDerivativeCoefficients, u::AbstractVector, b::AbstractVector, α, β)

Compute `α*D*u + β*dest` using the coefficients `b` and store the result in `dest`.
"""
function mul!(dest::AbstractVector, coefficients::VarCoefDerivativeCoefficients,
                u::AbstractVector, b::AbstractVector, α, β)
    @unpack coefficient_cache, parallel = coefficients

    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(dest) == length(b) DimensionMismatch
        @argcheck checkbounds(Bool, dest, coefficient_cache) DimensionMismatch
    end

    convolve_boundary_coefficients!(dest, coefficient_cache, u, b, α, β)
    convolve_interior_coefficients!(dest, coefficient_cache, u, b, α, β, parallel)
end

"""
    mul!(dest::AbstractVector, coefficients::VarCoefDerivativeCoefficients, u::AbstractVector, b::AbstractVector, α)

Compute `α*D*u` using the coefficients `b` and store the result in `dest`.
"""
function mul!(dest::AbstractVector, coefficients::VarCoefDerivativeCoefficients,
                u::AbstractVector, b::AbstractVector, α)
    @unpack coefficient_cache, parallel = coefficients

    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(dest) == length(b) DimensionMismatch
        @argcheck checkbounds(Bool, dest, coefficient_cache) DimensionMismatch
    end

    convolve_boundary_coefficients!(dest, coefficient_cache, u, b, α)
    convolve_interior_coefficients!(dest, coefficient_cache, u, b, α, parallel)
end

@inline function left_length(cache::AbstractCoefficientCache)
    length(cache.inv_left_weights)
end

@inline function right_length(cache::AbstractCoefficientCache)
    length(cache.inv_right_weights)
end


function convolve_interior_coefficients!(dest::AbstractVector, cache::AbstractCoefficientCache,
                                         u::AbstractVector, b::AbstractVector, α, parallel)
    for i in (left_length(cache)+1):(length(dest)-right_length(cache)) @inbounds begin
        retval = convolve_interior_coefficients_loopbody(i, cache, u, b)
        dest[i] = α*retval
    end end
end

function convolve_interior_coefficients!(dest::AbstractVector, cache::AbstractCoefficientCache,
                                         u::AbstractVector, b::AbstractVector, α, ::Val{:threads})
    Threads.@threads for i in (left_length(cache)+1):(length(dest)-right_length(cache)) @inbounds begin
        retval = convolve_interior_coefficients_loopbody(i, cache, u, b)
        dest[i] = α*retval
    end end
end


function convolve_interior_coefficients!(dest::AbstractVector, cache::AbstractCoefficientCache,
                                         u::AbstractVector, b::AbstractVector, α, β, parallel)
    for i in (left_length(cache)+1):(length(dest)-right_length(cache)) @inbounds begin
        retval = convolve_interior_coefficients_loopbody(i, cache, u, b)
        dest[i] = α*retval + β*dest[i]
    end end
end

function convolve_interior_coefficients!(dest::AbstractVector, cache::AbstractCoefficientCache,
                                         u::AbstractVector, b::AbstractVector, α, β, ::Val{:threads})
    Threads.@threads for i in (left_length(cache)+1):(length(dest)-right_length(cache)) @inbounds begin
        retval = convolve_interior_coefficients_loopbody(i, cache, u, b)
        dest[i] = α*retval + β*dest[i]
    end end
end



abstract type AbstractVariableCoefficientDerivativeOperator{T} <: AbstractDerivativeOperator{T} end

@inline source_of_coeffcients(D::AbstractVariableCoefficientDerivativeOperator) = source_of_coeffcients(D.coefficients)

@inline function lower_bandwidth(D::AbstractVariableCoefficientDerivativeOperator)
    lower_bandwidth(D.coefficients.coefficient_cache)
end

@inline function upper_bandwidth(D::AbstractVariableCoefficientDerivativeOperator)
    upper_bandwidth((D.coefficients.coefficient_cache))
end



"""
    DissipationOperator

A dissipation operator on a nonperiodic finite difference grid.
"""
struct DissipationOperator{T,CoefficientCache,Parallel,SourceOfCoefficients,Grid} <: AbstractVariableCoefficientDerivativeOperator{T}
    factor::T
    coefficients::VarCoefDerivativeCoefficients{T,CoefficientCache,Parallel,SourceOfCoefficients}
    grid::Grid
    b::Vector{T}

    function DissipationOperator(factor::T, 
                                 coefficients::VarCoefDerivativeCoefficients{T,CoefficientCache,Parallel,SourceOfCoefficients},
                                 grid::Grid, b::Vector{T}) where {T,CoefficientCache,Parallel,SourceOfCoefficients,Grid}
        @argcheck checkbounds(Bool, grid, coefficients.coefficient_cache) DimensionMismatch
        @argcheck length(grid) == length(b)
        factor < 0 && @warn("Negative dissipation strength shouldn't be used.")

        new{T,CoefficientCache,Parallel,SourceOfCoefficients,Grid}(
            factor,coefficients, grid, b)
    end
end


function Base.show(io::IO, D::DissipationOperator{T}) where {T}
    if  derivative_order(D) == 2
        print(io, "SBP 2nd derivative dissipation operator of order ")
    else
        print(io, "SBP ", derivative_order(D), "th derivative dissipation operator of order ")
    end
    print(io, accuracy_order(D), " {T=", T, ", Parallel=", typeof(D.coefficients.parallel), "} \n")
    print(io, "on a grid in [", first(grid(D)), ", ", last(grid(D)),
                "] using ", length(grid(D)), " nodes \n")
    print(io, "and coefficients given in \n")
    print(io, source_of_coeffcients(D))
end


"""
    mul!(dest::AbstractVector, D::DissipationOperator, u::AbstractVector, α, β)

Compute `α*D*u + β*dest` and store the result in `dest`.
"""
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::DissipationOperator, 
                                       u::AbstractVector, α, β)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, D.b, D.factor*α, β)
end

"""
    mul!(dest::AbstractVector, D::DissipationOperator, u::AbstractVector, α)

Compute `α*D*u` and store the result in `dest`.
"""
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::DissipationOperator, 
                                       u::AbstractVector, α)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, D.b, D.factor*α)
end


"""
    dissipation_operator(source_of_coefficients, order, xmin, xmax, N, 
                         left_weights, right_weights, parallel=Val{:serial}())

Create a negative semidefinite `DissipationOperator` using undivided differences
approximating a weighted `order`-th derivative on a grid between `xmin` and 
`xmax` with `N` grid points up to order of accuracy 2 with coefficients given
by `source_of_coefficients`. 
The norm matrix is given by `left_weights` and `right_weights`.
The evaluation of the derivative can be parallised using threads by chosing
`parallel=Val{:threads}())`.
"""
function dissipation_operator(source_of_coefficients, order, xmin, xmax, N, 
                              left_weights, right_weights, 
                              strength=one(xmin+xmax),
                              parallel=Val{:serial}())
    grid = construct_grid(source_of_coefficients, order, xmin, xmax, N)
    coefficients, b = dissipation_coefficients(source_of_coefficients, order, grid, left_weights, right_weights, parallel)
    DissipationOperator(strength, coefficients, grid, b)
end

"""
    dissipation_operator(source_of_coefficients, D::DerivativeOperator{T};
                         strength=one(T),
                         order::Int=accuracy_order(D), 
                         parallel=D.coefficients.parallel)

Create a negative semidefinite `DissipationOperator` using undivided differences
approximating a weighted `order`-th derivative adapted to the derivative
operator `D` with coefficients given in `source_of_coefficients`.
The evaluation of the derivative can be parallised using threads by chosing
`parallel=Val{:threads}())`.
"""
function dissipation_operator(source_of_coefficients, D::DerivativeOperator{T};
                              strength=one(T),
                              order::Int=accuracy_order(D),
                              parallel=D.coefficients.parallel) where {T}
    x = grid(D)
    dissipation_operator(source_of_coefficients, order, first(x), last(x), length(x), D.coefficients.left_weights, D.coefficients.right_weights, strength, parallel)
end

"""
    dissipation_operator(D::DerivativeOperator; kwargs...)

Create a negative semidefinite `DissipationOperator` using undivided differences
approximating a weighted `order`-th derivative adapted to the derivative
operator `D`.
The evaluation of the derivative can be parallised using threads by chosing
`parallel=Val{:threads}())`.
"""
function dissipation_operator(D::DerivativeOperator; kwargs...)
    source = MattssonSvärdNordström2004()
    dissipation_operator(source, D; kwargs...)
end
