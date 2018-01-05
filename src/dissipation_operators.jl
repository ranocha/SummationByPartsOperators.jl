
"""
    DissipationCoefficients

The coefficients of a dissipation operator on a nonperiodic grid.
"""
struct DissipationCoefficients{T,LeftBoundary,RightBoundary,LowerCoef,CentralCoef<:DerivativeCoefficientRow{T},
                               UpperCoef,Parallel,SourceOfCoefficients} <: AbstractDerivativeCoefficients{T}
    # coefficients defining the operator and its action
    left_boundary::LeftBoundary
    right_boundary::RightBoundary
    lower_coef::LowerCoef
    central_coef::CentralCoef
    upper_coef::UpperCoef
    parallel::Parallel
    # corresponding orders etc.
    derivative_order::Int
    accuracy_order  ::Int
    source_of_coeffcients::SourceOfCoefficients

    function DissipationCoefficients(left_boundary::LeftBoundary, right_boundary::RightBoundary,
                                    lower_coef::LowerCoef, central_coef::CentralCoef, upper_coef::UpperCoef,
                                    parallel::Parallel, derivative_order::Int, accuracy_order::Int,
                                    source_of_coeffcients::SourceOfCoefficients) where {T,LeftBoundary,RightBoundary,LowerCoef,CentralCoef<:DerivativeCoefficientRow{T},
                                                                                        UpperCoef,Parallel,SourceOfCoefficients}
        new{T,LeftBoundary,RightBoundary,LowerCoef,CentralCoef,UpperCoef,Parallel,SourceOfCoefficients}(
            left_boundary, right_boundary, lower_coef, central_coef, upper_coef,
            parallel, derivative_order, accuracy_order, source_of_coeffcients)
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
    @unpack left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel = coefficients

    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(dest) == length(b) DimensionMismatch
        @argcheck length(u) > length(lower_coef)+length(upper_coef) DimensionMismatch
        @argcheck length(u) > length(left_boundary) + length(right_boundary) DimensionMismatch
    end

    convolve_variable_boundary_coefficients!(dest, left_boundary, right_boundary, u, b, α, β)
    convolve_variable_interior_coefficients!(dest, lower_coef, central_coef, upper_coef, u, b, α, β, length(left_boundary), length(right_boundary), parallel)
end

"""
    mul!(dest::AbstractVector, coefficients::DissipationCoefficients, u::AbstractVector, b::AbstractVector, α)

Compute `α*D*u` using the coefficients `b` and store the result in `dest`.
"""
function mul!(dest::AbstractVector, coefficients::DissipationCoefficients, u::AbstractVector, b::AbstractVector, α)
    @unpack left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel = coefficients

    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(dest) == length(b) DimensionMismatch
        @argcheck length(u) > length(lower_coef)+length(upper_coef) DimensionMismatch
        @argcheck length(u) > length(left_boundary) + length(right_boundary) DimensionMismatch
    end

    convolve_variable_boundary_coefficients!(dest, left_boundary, right_boundary, u, b, α)
    convolve_variable_interior_coefficients!(dest, lower_coef, central_coef, upper_coef, u, b, α, length(left_boundary), length(right_boundary), parallel)
end


@inline function convolve_row(coef_row::DerivativeCoefficientRow{T,Offset,Length}, i::Int, u, b) where {T,Offset,Length}
    @inbounds begin
        tmp = coef_row.coef[1]*b[i+Offset]
        for j in 2:length(coef_row)
            tmp += coef_row.coef[j]*b[i+j-1+Offset]
        end
        retval = tmp*u[i]
    end
    retval
end
#=
@inline @unroll function convolve_row(coef_row::DerivativeCoefficientRow, i::Int, u, b)
    Offset = offset(coef_row)

    @inbounds begin
        tmp = zero(promote_type(eltype(coef_row.coef), eltype(b)))
        @unroll for j in 1:length(coef_row)
            tmp += coef_row.coef[j]*b[i+j-1+Offset]
        end
        retval = tmp*u[i]
    end
    retval
end
=#

@inline @unroll function convolve_variable_boundary_coefficients!(dest::AbstractVector, left_boundary, right_boundary, u::AbstractVector, b::AbstractVector, α)
    T = eltype(dest)

    @unroll for i in 1:length(left_boundary) @inbounds begin
        tmp = convolve_left_boundary(T, left_boundary[i], u, b)
        dest[i] = α*tmp
    end end

    @unroll for i in 1:length(right_boundary) @inbounds begin
        tmp = convolve_right_boundary(T, right_boundary[i], u, b)
        dest[end+1-i] = α*tmp
    end end
end

@inline @unroll function convolve_variable_boundary_coefficients!(dest::AbstractVector, left_boundary, right_boundary, u::AbstractVector, b::AbstractVector, α, β)
    T = eltype(dest)

    @unroll for i in 1:length(left_boundary) @inbounds begin
        tmp = convolve_left_boundary(T, left_boundary[i], u, b)
        dest[i] = α*tmp + β*dest[i]
    end end

    @unroll for i in 1:length(right_boundary) @inbounds begin
        tmp = convolve_right_boundary(T, right_boundary[i], u, b)
        dest[end+1-i] = α*tmp + β*dest[end+1-i]
    end end
end


@inline @unroll function convolve_left_boundary(T, boundary_coef, u, b)
    tmp = zero(T)
    @unroll for j in 1:length(boundary_coef)
        tmp += convolve_row(boundary_coef[j], j, u, b)
    end
    tmp
end

@inline @unroll function convolve_right_boundary(T, boundary_coef, u, b)
    tmp = zero(T)
    L = length(u)+1
    @unroll for j in 1:length(boundary_coef)
        tmp += convolve_row(boundary_coef[j], L-j, u, b)
    end
    tmp
end


@inline @unroll function convolve_variable_interior_coefficients!(dest::AbstractVector, lower_coef, central_coef, upper_coef, u::AbstractVector, b::AbstractVector, α, left_boundary_width, right_boundary_width, parallel)
    T = eltype(dest)
    for i in (left_boundary_width+1):(length(dest)-right_boundary_width) @inbounds begin
        tmp = zero(T)
        @unroll for j in 1:length(lower_coef)
            tmp += convolve_row(lower_coef[j], i-j, u, b)
        end
        tmp += convolve_row(central_coef, i, u, b)
        @unroll for j in 1:length(upper_coef)
            tmp += convolve_row(upper_coef[j], i+j, u, b)
        end
        dest[i] = α*tmp
    end end
end

@inline function convolve_variable_interior_coefficients!(dest::AbstractVector, lower_coef, central_coef, upper_coef, u::AbstractVector, b::AbstractVector, α, left_boundary_width, right_boundary_width, parallel::Val{:threads})
    Threads.@threads for i in (left_boundary_width+1):(length(dest)-right_boundary_width)
        convolve_variable_interior_coefficients_threaded_loop!(dest, i, lower_coef, central_coef, upper_coef, u, b, α)
    end
end

@noinline @unroll function convolve_variable_interior_coefficients_threaded_loop!(dest::AbstractVector, i, lower_coef, central_coef, upper_coef, u::AbstractVector, b::AbstractVector, α)
    T = eltype(dest)
    @inbounds begin
        tmp = zero(T)
        @unroll for j in 1:length(lower_coef)
            tmp += convolve_row(lower_coef[j], i-j, u, b)
        end
        tmp += convolve_row(central_coef, i, u, b)
        @unroll for j in 1:length(upper_coef)
            tmp += convolve_row(upper_coef[j], i+j, u, b)
        end
        dest[i] = α*tmp
    end
end

@inline @unroll function convolve_variable_interior_coefficients!(dest::AbstractVector, lower_coef, central_coef, upper_coef, u::AbstractVector, b::AbstractVector, α, β, left_boundary_width, right_boundary_width, parallel)
    T = eltype(dest)
    for i in (left_boundary_width+1):(length(dest)-right_boundary_width) @inbounds begin
        tmp = zero(T)
        @unroll for j in 1:length(lower_coef)
            tmp += convolve_row(lower_coef[j], i-j, u, b)
        end
        tmp += convolve_row(central_coef, i, u, b)
        @unroll for j in 1:length(upper_coef)
            tmp += convolve_row(upper_coef[j], i+j, u, b)
        end
        dest[i] = α*tmp + β*dest[i]
    end end
end

@inline function convolve_variable_interior_coefficients!(dest::AbstractVector, lower_coef, central_coef, upper_coef, u::AbstractVector, b::AbstractVector, α, β, left_boundary_width, right_boundary_width, parallel::Val{:threads})
    Threads.@threads for i in (left_boundary_width+1):(length(dest)-right_boundary_width)
        convolve_variable_interior_coefficients_threaded_loop!(dest, i, lower_coef, central_coef, upper_coef, u, b, α, β)
    end
end

@noinline @unroll function convolve_variable_interior_coefficients_threaded_loop!(dest::AbstractVector, i, lower_coef, central_coef, upper_coef, u::AbstractVector, b::AbstractVector, α, β)
    T = eltype(dest)
    @inbounds begin
        tmp = zero(T)
        @unroll for j in 1:length(lower_coef)
            tmp += convolve_row(lower_coef[j], i-j, u, b)
        end
        tmp += convolve_row(central_coef, i, u, b)
        @unroll for j in 1:length(upper_coef)
            tmp += convolve_row(upper_coef[j], i+j, u, b)
        end
        dest[i] = α*tmp + β*dest[i]
    end
end



"""
    DissipationOperator

A dissipation operator on a nonperiodic finite difference grid.
"""
struct DissipationOperator{T,LeftBoundary,RightBoundary,LowerCoef,CentralCoef,
                           UpperCoef,Parallel,SourceOfCoefficients,Grid} <: AbstractDerivativeOperator{T}
    coefficients::DissipationCoefficients{T,LeftBoundary,RightBoundary,LowerCoef,CentralCoef,
                                          UpperCoef,Parallel,SourceOfCoefficients}
    grid::Grid
    b::Vector{T}

    function DissipationOperator(coefficients::DissipationCoefficients{T,LeftBoundary,RightBoundary,LowerCoef,CentralCoef,
                                                                       UpperCoef,Parallel,SourceOfCoefficients},
                                grid::Grid, b::Vector{T}) where {T,LeftBoundary,RightBoundary,LowerCoef,CentralCoef,
                                                                 UpperCoef,Parallel,SourceOfCoefficients,Grid}
        @argcheck length(grid) > length(coefficients.lower_coef)+length(coefficients.upper_coef) DimensionMismatch
        @argcheck length(grid) > length(coefficients.left_boundary) + length(coefficients.right_boundary) DimensionMismatch
        @argcheck length(grid) == length(b)

        new{T,LeftBoundary,RightBoundary,LowerCoef,CentralCoef,UpperCoef,Parallel,SourceOfCoefficients,Grid}(
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
    dissipation_operator(source_of_coefficients, D::DerivativeOperator, order::Int=accuracy_order(D), parallel=Val{:serial}())

Create a `DissipationOperator` approximating a weighted `order`-th derivative
adapted to the derivative operator `D`.
The evaluation of the derivative can be parallised using threads by chosing
`parallel=Val{:threads}())`.
"""
function dissipation_operator(source_of_coefficients, D::DerivativeOperator, order::Int=accuracy_order(D), parallel=D.coefficients.parallel)
    dissipation_operator(source_of_coefficients, order, first(grid(D)), last(grid(D)), length(grid(D)), D.coefficients.left_weights, D.coefficients.right_weights, parallel)
end
