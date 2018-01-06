
abstract type AbstractCoefficientCache{T} end

function Base.checkbounds(::Type{Bool}, u::AbstractVector, cache::AbstractCoefficientCache)
    length(u) > length(cache.inv_left_weights)+length(cache.inv_right_weights)
end

"""
    DissipationCoefficients

The coefficients of a dissipation operator on a nonperiodic grid.
"""
struct DissipationCoefficients{T,CoefficientCache<:AbstractCoefficientCache{T},
                               Parallel,SourceOfCoefficients} <: AbstractDerivativeCoefficients{T}
    # coefficients defining the operator and its action
    coefficient_cache::CoefficientCache
    parallel::Parallel
    # corresponding orders etc.
    derivative_order::Int
    accuracy_order  ::Int
    source_of_coeffcients::SourceOfCoefficients

    function DissipationCoefficients(coefficient_cache::CoefficientCache,
                                    parallel::Parallel, derivative_order::Int, accuracy_order::Int,
                                    source_of_coeffcients::SourceOfCoefficients) where {T,CoefficientCache<:AbstractCoefficientCache{T},
                                                                                        Parallel,SourceOfCoefficients}
        new{T,CoefficientCache,Parallel,SourceOfCoefficients}(
            coefficient_cache, parallel, derivative_order, accuracy_order, source_of_coeffcients)
    end
end


function Base.show(io::IO, coefficients::DissipationCoefficients)
    if  derivative_order(coefficients) == 2
        print(io, "Coefficients of the SBP 2nd derivative dissipation operator of order ",
                    accuracy_order(coefficients), " given in \n")
    else
        print(io, "Coefficients of the SBP ", derivative_order(coefficients),
                    "th derivative dissipation operator of order ", accuracy_order(coefficients),
                    " given in \n")
    end
    print(io, source_of_coeffcients(coefficients))
end


"""
    mul!(dest::AbstractVector, coefficients::DissipationCoefficients, u::AbstractVector, b::AbstractVector, α, β)

Compute `α*D*u + β*dest` using the coefficients `b` and store the result in `dest`.
"""
function mul!(dest::AbstractVector, coefficients::DissipationCoefficients, u::AbstractVector, b::AbstractVector, α, β)
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
    mul!(dest::AbstractVector, coefficients::DissipationCoefficients, u::AbstractVector, b::AbstractVector, α)

Compute `α*D*u` using the coefficients `b` and store the result in `dest`.
"""
function mul!(dest::AbstractVector, coefficients::DissipationCoefficients, u::AbstractVector, b::AbstractVector, α)
    @unpack coefficient_cache, parallel = coefficients

    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(dest) == length(b) DimensionMismatch
        @argcheck checkbounds(Bool, dest, coefficient_cache) DimensionMismatch
    end

    convolve_boundary_coefficients!(dest, coefficient_cache, u, b, α)
    convolve_interior_coefficients!(dest, coefficient_cache, u, b, α, parallel)
end


function convolve_interior_coefficients!(dest::AbstractVector, cache::AbstractCoefficientCache, u::AbstractVector, b::AbstractVector, α, parallel)
    for i in (length(cache.inv_left_weights)+1):(length(dest)-length(cache.inv_right_weights)) @inbounds begin
        retval = convolve_interior_coefficients_loopbody(i, cache, u, b)
        dest[i] = α*retval
    end end
end

function convolve_interior_coefficients!(dest::AbstractVector, cache::AbstractCoefficientCache, u::AbstractVector, b::AbstractVector, α, parallel::Val{:threads})
    Threads.@threads for i in (length(cache.inv_left_weights)+1):(length(dest)-length(cache.inv_right_weights)) @inbounds begin
        retval = convolve_interior_coefficients_loopbody(i, cache, u, b)
        dest[i] = α*retval
    end end
end


function convolve_interior_coefficients!(dest::AbstractVector, cache::AbstractCoefficientCache, u::AbstractVector, b::AbstractVector, α, β, parallel)
    for i in (length(cache.inv_left_weights)+1):(length(dest)-length(cache.inv_right_weights)) @inbounds begin
        retval = convolve_interior_coefficients_loopbody(i, cache, u, b)
        dest[i] = α*retval + β*dest[i]
    end end
end

function convolve_interior_coefficients!(dest::AbstractVector, cache::AbstractCoefficientCache, u::AbstractVector, b::AbstractVector, α, β, parallel::Val{:threads})
    Threads.@threads for i in (length(cache.inv_left_weights)+1):(length(dest)-length(cache.inv_right_weights)) @inbounds begin
        retval = convolve_interior_coefficients_loopbody(i, cache, u, b)
        dest[i] = α*retval + β*dest[i]
    end end
end


"""
    DissipationOperator

A dissipation operator on a nonperiodic finite difference grid.
"""
struct DissipationOperator{T,CoefficientCache,Parallel,SourceOfCoefficients,Grid} <: AbstractDerivativeOperator{T}
    coefficients::DissipationCoefficients{T,CoefficientCache,Parallel,SourceOfCoefficients}
    grid::Grid
    b::Vector{T}

    function DissipationOperator(coefficients::DissipationCoefficients{T,CoefficientCache,Parallel,SourceOfCoefficients},
                                grid::Grid, b::Vector{T}) where {T,CoefficientCache,Parallel,SourceOfCoefficients,Grid}
        @argcheck checkbounds(Bool, grid, coefficients.coefficient_cache) DimensionMismatch
        @argcheck length(grid) == length(b)

        new{T,CoefficientCache,Parallel,SourceOfCoefficients,Grid}(
            coefficients, grid, b)
    end
end


@inline source_of_coeffcients(D::DissipationOperator) = source_of_coeffcients(D.coefficients)


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
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::DissipationOperator, u::AbstractVector, α, β)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, D.b, α, β)
end

"""
    mul!(dest::AbstractVector, D::DissipationOperator, u::AbstractVector, α)

Compute `α*D*u` and store the result in `dest`.
"""
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::DissipationOperator, u::AbstractVector, α)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, D.b, α)
end


"""
    dissipation_operator(source_of_coefficients, order, xmin, xmax, N, left_weights, right_weights, parallel=Val{:serial}())

Create a `DissipationOperator` approximating a weighted `order`-th derivative on
a grid between `xmin` and `xmax` with `N` grid points up to order of accuracy 2
with coefficients given by `source_of_coefficients`. The norm matrix is given by
`left_weights` and `right_weights`.
The evaluation of the derivative can be parallised using threads by chosing
`parallel=Val{:threads}())`.
"""
function dissipation_operator(source_of_coefficients, order, xmin, xmax, N, left_weights, right_weights, parallel=Val{:serial}())
    grid = construct_grid(source_of_coefficients, order, xmin, xmax, N)
    coefficients, b = dissipation_coefficients(source_of_coefficients, order, grid, left_weights, right_weights, parallel)
    DissipationOperator(coefficients, grid, b)
end

"""
    dissipation_operator(source_of_coefficients, D::DerivativeOperator, order::Int=accuracy_order(D), parallel=D.coefficients.parallel)

Create a `DissipationOperator` approximating a weighted `order`-th derivative
adapted to the derivative operator `D`.
The evaluation of the derivative can be parallised using threads by chosing
`parallel=Val{:threads}())`.
"""
function dissipation_operator(source_of_coefficients, D::DerivativeOperator, order::Int=accuracy_order(D), parallel=D.coefficients.parallel)
    dissipation_operator(source_of_coefficients, order, first(grid(D)), last(grid(D)), length(grid(D)), D.coefficients.left_weights, D.coefficients.right_weights, parallel)
end