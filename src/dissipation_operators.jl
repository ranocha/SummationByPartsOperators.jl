
abstract type AbstractCoefficientCache{T} end

function Base.checkbounds(::Type{Bool}, u::AbstractVector, cache::AbstractCoefficientCache)
    length(u) > length(cache.inv_left_weights)+length(cache.inv_right_weights)
end

"""
    VarCoefDerivativeCoefficients

The coefficients of a variable coefficient derivative operator on a nonperiodic grid.
"""
struct VarCoefDerivativeCoefficients{T,CoefficientCache<:AbstractCoefficientCache{T},
                                     LeftWidth,RightWidth,ExecutionMode,SourceOfCoefficients} <: AbstractDerivativeCoefficients{T}
    # coefficients defining the operator and its action
    coefficient_cache::CoefficientCache
    left_weights::SVector{LeftWidth, T}
    right_weights::SVector{RightWidth, T}
    mode::ExecutionMode
    # corresponding orders etc.
    derivative_order::Int
    accuracy_order  ::Int
    source_of_coefficients::SourceOfCoefficients

    function VarCoefDerivativeCoefficients(coefficient_cache::CoefficientCache,
                                           left_weights::SVector{LeftWidth, T},
                                           right_weights::SVector{RightWidth, T},
                                           mode::ExecutionMode,
                                           derivative_order::Int, accuracy_order::Int,
                                           source_of_coefficients::SourceOfCoefficients) where {T,CoefficientCache<:AbstractCoefficientCache{T},
                                                                                                LeftWidth,RightWidth,ExecutionMode,SourceOfCoefficients}
        new{T,CoefficientCache,LeftWidth,RightWidth,ExecutionMode,SourceOfCoefficients}(
            coefficient_cache, left_weights, right_weights, mode, derivative_order, accuracy_order, source_of_coefficients)
    end
end


function Base.show(io::IO, coefficients::VarCoefDerivativeCoefficients)
    if derivative_order(coefficients) == 1
        print(io, "Coefficients of the variable-coefficient first-derivative operator")
    elseif  derivative_order(coefficients) == 2
        print(io, "Coefficients of the variable-coefficient second-derivative operator")
    elseif  derivative_order(coefficients) == 3
        print(io, "Coefficients of the variable-coefficient third-derivative operator")
    else
        print(io, "Coefficients of the variable-coefficient ", derivative_order(coefficients),
              "-derivative operator")
    end
    print(io, " of order ", accuracy_order(coefficients), " of ")
    print(io, source_of_coefficients(coefficients))
end


# Compute `α*D*u + β*dest` using the coefficients `b` and store the result in `dest`.
function mul!(dest::AbstractVector, coefficients::VarCoefDerivativeCoefficients,
                u::AbstractVector, b::AbstractVector, α, β)
    @unpack coefficient_cache, mode = coefficients

    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(dest) == length(b) DimensionMismatch
        @argcheck checkbounds(Bool, dest, coefficient_cache) DimensionMismatch
    end

    convolve_boundary_coefficients!(dest, coefficient_cache, u, b, α, β)
    convolve_interior_coefficients!(dest, coefficient_cache, u, b, α, β, mode)
end

# Compute `α*D*u` using the coefficients `b` and store the result in `dest`.
function mul!(dest::AbstractVector, coefficients::VarCoefDerivativeCoefficients,
                u::AbstractVector, b::AbstractVector, α)
    @unpack coefficient_cache, mode = coefficients

    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(dest) == length(b) DimensionMismatch
        @argcheck checkbounds(Bool, dest, coefficient_cache) DimensionMismatch
    end

    convolve_boundary_coefficients!(dest, coefficient_cache, u, b, α)
    convolve_interior_coefficients!(dest, coefficient_cache, u, b, α, mode)
end

@inline function left_length(cache::AbstractCoefficientCache)
    length(cache.inv_left_weights)
end

@inline function right_length(cache::AbstractCoefficientCache)
    length(cache.inv_right_weights)
end


function convolve_interior_coefficients!(dest::AbstractVector, cache::AbstractCoefficientCache,
                                         u::AbstractVector, b::AbstractVector, α, mode)
    for i in (left_length(cache)+1):(length(dest)-right_length(cache)) @inbounds begin
        retval = convolve_interior_coefficients_loopbody(i, cache, u, b)
        dest[i] = α*retval
    end end
end

function convolve_interior_coefficients!(dest::AbstractVector, cache::AbstractCoefficientCache,
                                         u::AbstractVector, b::AbstractVector, α, mode::ThreadedMode)
    Threads.@threads for i in (left_length(cache)+1):(length(dest)-right_length(cache)) @inbounds begin
        retval = convolve_interior_coefficients_loopbody(i, cache, u, b)
        dest[i] = α*retval
    end end
end


function convolve_interior_coefficients!(dest::AbstractVector, cache::AbstractCoefficientCache,
                                         u::AbstractVector, b::AbstractVector, α, β, mode)
    for i in (left_length(cache)+1):(length(dest)-right_length(cache)) @inbounds begin
        retval = convolve_interior_coefficients_loopbody(i, cache, u, b)
        dest[i] = α*retval + β*dest[i]
    end end
end

function convolve_interior_coefficients!(dest::AbstractVector, cache::AbstractCoefficientCache,
                                         u::AbstractVector, b::AbstractVector, α, β, mode::ThreadedMode)
    Threads.@threads for i in (left_length(cache)+1):(length(dest)-right_length(cache)) @inbounds begin
        retval = convolve_interior_coefficients_loopbody(i, cache, u, b)
        dest[i] = α*retval + β*dest[i]
    end end
end



abstract type AbstractVariableCoefficientNonperiodicDerivativeOperator{T} <: AbstractNonperiodicDerivativeOperator{T} end
abstract type AbstractVariableCoefficientPeriodicDerivativeOperator{T} <: AbstractPeriodicDerivativeOperator{T} end


@inline source_of_coefficients(D::Union{AbstractVariableCoefficientNonperiodicDerivativeOperator,AbstractVariableCoefficientNonperiodicDerivativeOperator}) = source_of_coefficients(D.coefficients)

@inline function lower_bandwidth(D::AbstractVariableCoefficientNonperiodicDerivativeOperator)
    lower_bandwidth(D.coefficients.coefficient_cache)
end

@inline function upper_bandwidth(D::AbstractVariableCoefficientNonperiodicDerivativeOperator)
    upper_bandwidth((D.coefficients.coefficient_cache))
end



"""
    DissipationOperator

A dissipation operator on a nonperiodic finite difference grid.
See [`dissipation_operator`](@ref).
"""
struct DissipationOperator{T,Coefficients<:VarCoefDerivativeCoefficients{T},Grid} <: AbstractVariableCoefficientNonperiodicDerivativeOperator{T}
    factor::T
    coefficients::Coefficients
    grid::Grid
    b::Vector{T}

    function DissipationOperator(factor::T,
                                 coefficients::Coefficients,
                                 grid::Grid, b::Vector{T}) where {T,Coefficients<:VarCoefDerivativeCoefficients{T},Grid}
        @argcheck checkbounds(Bool, grid, coefficients.coefficient_cache) DimensionMismatch
        @argcheck length(grid) == length(b)
        factor < 0 && @warn("Negative dissipation strength shouldn't be used.")

        new{T,Coefficients,Grid}(factor,coefficients, grid, b)
    end
end


function Base.show(io::IO, D::DissipationOperator{T}) where {T}
    if  derivative_order(D) == 2
        print(io, "SBP second-derivative dissipation operator")
    else
        print(io, "SBP ", derivative_order(D),
              "-derivative dissipation operator")
    end
    print(io, " of order ", accuracy_order(D))
    if get(io, :compact, false) == false
        print(io, " on a grid in [", first(grid(D)), ", ", last(grid(D)),
                "] using ", length(grid(D)), " nodes \n")
        print(io, "and coefficients")
    end
    print(io, " of ", source_of_coefficients(D))
end



# Compute `α*D*u + β*dest` and store the result in `dest`.
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::DissipationOperator,
                                       u::AbstractVector, α, β)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, D.b, D.factor*α, β)
end

# Compute `α*D*u` and store the result in `dest`.
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
                         left_weights, right_weights, mode=FastMode())

Create a negative semidefinite `DissipationOperator` using undivided differences
approximating a weighted `order`-th derivative on a grid between `xmin` and
`xmax` with `N` grid points up to order of accuracy 2 with coefficients given
by `source_of_coefficients`.
The norm matrix is given by `left_weights` and `right_weights`.
The evaluation of the derivative can be parallized using threads by chosing
`mode=ThreadedMode()`.
"""
function dissipation_operator(source_of_coefficients, order, xmin, xmax, N,
                              left_weights, right_weights,
                              strength=one(xmin+xmax),
                              mode=FastMode(),
                              parallel=nothing)
    if parallel !== nothing
        # TODO: deprecated in v0.5
        Base.depwarn("Providing the keyword argument `parallel` is deprecated." *
                     "Use `mode` instead.", :dissipation_operator)
        mode = _parallel_to_mode(parallel)
    end
    grid = construct_grid(source_of_coefficients, order, xmin, xmax, N)
    coefficients, b = dissipation_coefficients(source_of_coefficients, order, grid, left_weights, right_weights, mode)
    DissipationOperator(strength, coefficients, grid, b)
end

"""
    dissipation_operator([source_of_coefficients=MattssonSvärdNordström2004()],
                         D::DerivativeOperator{T};
                         strength=one(T),
                         order::Int=accuracy_order(D),
                         mode=D.coefficients.mode)

Create a negative semidefinite `DissipationOperator` using undivided differences
approximating a weighted `order`-th derivative adapted to the derivative
operator `D` with coefficients given in `source_of_coefficients`.
The evaluation of the derivative can be parallized using threads by chosing
`mode=ThreadedMode()`.
"""
function dissipation_operator(source_of_coefficients, D::DerivativeOperator{T};
                              strength=one(T),
                              order::Int=accuracy_order(D),
                              mode=D.coefficients.mode,
                              parallel=nothing) where {T}
    if parallel !== nothing
        # TODO: deprecated in v0.5
        Base.depwarn("Providing the keyword argument `parallel` is deprecated." *
                     "Use `mode` instead.", :dissipation_operator)
        mode = _parallel_to_mode(parallel)
    end
    x = grid(D)
    dissipation_operator(source_of_coefficients, order, first(x), last(x), length(x), D.coefficients.left_weights, D.coefficients.right_weights, strength, mode)
end

function dissipation_operator(D::DerivativeOperator; kwargs...)
    source = MattssonSvärdNordström2004()
    dissipation_operator(source, D; kwargs...)
end
