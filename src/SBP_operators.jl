
"""
    DerivativeCoefficients{T,BoundaryWidth,LowerOffset,UpperOffset,Parallel}

The coefficients of a derivative operator on a periodic grid.
"""
struct DerivativeCoefficients{T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients} <: AbstractDerivativeCoefficients{T}
    # coefficients defining the operator and its action
    left_boundary::LeftBoundary
    right_boundary::RightBoundary
    lower_coef::SVector{LowerOffset, T}
    central_coef::T
    upper_coef::SVector{UpperOffset, T}
    parallel::Parallel
    # corresponding orders etc.
    derivative_order::Int
    accuracy_order  ::Int
    symmetric       ::Bool
    source_of_coeffcients::SourceOfCoefficients

    function DerivativeCoefficients(left_boundary::LeftBoundary, right_boundary::RightBoundary, lower_coef::SVector{LowerOffset, T}, central_coef::T, upper_coef::SVector{UpperOffset, T}, parallel::Parallel, derivative_order::Int, accuracy_order::Int, source_of_coeffcients::SourceOfCoefficients) where {T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients}
        symmetric = LowerOffset == UpperOffset
        if symmetric
            @inbounds for i in Base.OneTo(LowerOffset)
                symmetric = symmetric && lower_coef[i] == upper_coef[i]
            end
        end
        #TODO: check boundary coefficients
        new{T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients}(left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel, derivative_order, accuracy_order, symmetric, source_of_coeffcients)
    end
end


@inline source_of_coeffcients(coefficients::DerivativeCoefficients) = coefficients.source_of_coeffcients


function Base.show(io::IO, coefficients::DerivativeCoefficients)
    if derivative_order(coefficients) == 1
        print(io, "Coefficients of the SBP 1st derivative operator of order ",
                    accuracy_order(coefficients), " given in \n")
    elseif  derivative_order(coefficients) == 2
        print(io, "Coefficients of the SBP 2nd derivative operator of order ",
                    accuracy_order(coefficients), " given in \n")
    elseif  derivative_order(coefficients) == 3
        print(io, "Coefficients of the SBP 3rd derivative operator of order ",
                    accuracy_order(coefficients), " given in \n")
    else
        print(io, "Coefficients of the SBP ", derivative_order(coefficients),
                    "th derivative operator of order ", accuracy_order(coefficients),
                    " given in \n")
    end
    print(io, source_of_coeffcients(coefficients))
end


"""
    mul!(dest::AbstractVector, coefficients::DerivativeCoefficients, u::AbstractVector, α, β)

Compute `α*D*u + β*dest` and store the result in `dest`.
"""
function mul!(dest::AbstractVector, coefficients::DerivativeCoefficients{T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset}, u::AbstractVector, α, β) where {T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset}
    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(u) > LowerOffset+UpperOffset DimensionMismatch
        @argcheck length(u) > length(coefficients.left_boundary) + length(coefficients.right_boundary) DimensionMismatch
    end

    @unpack left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel = coefficients
    convolve_boundary_coefficients!(dest, left_boundary, right_boundary, u, α, β)
    convolve_interior_coefficients!(dest, lower_coef, central_coef, upper_coef, u, α, β, length(left_boundary), length(right_boundary), parallel)
end

"""
    mul!(dest::AbstractVector, coefficients::DerivativeCoefficients, u::AbstractVector, α)

Compute `α*D*u` and store the result in `dest`.
"""
function mul!(dest::AbstractVector, coefficients::DerivativeCoefficients{T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset}, u::AbstractVector, α) where {T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset}
    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(u) > LowerOffset+UpperOffset DimensionMismatch
        @argcheck length(u) > length(coefficients.left_boundary) + length(coefficients.right_boundary) DimensionMismatch
    end

    @unpack left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel = coefficients
    convolve_boundary_coefficients!(dest, left_boundary, right_boundary, u, α)
    convolve_interior_coefficients!(dest, lower_coef, central_coef, upper_coef, u, α, length(left_boundary), length(right_boundary), parallel)
end


"""
    DerivativeCoefficientRow{T,Start,Length}

A struct representing a row in the boundary block of an SBP derivative operator
with scalar type `T`.
"""
struct DerivativeCoefficientRow{T,Start,Length}
    coef::SVector{Length, T}
end

function -(coef_row::DerivativeCoefficientRow{T,Start,Length}) where {T,Start,Length}
    DerivativeCoefficientRow{T,Start,Length}(-coef_row.coef)
end

@inline function convolve_left_row(coef_row::DerivativeCoefficientRow{T,Start,Length}, u) where {T,Start,Length}
    @inbounds tmp = coef_row.coef[1]*u[Start]
    @inbounds for i in 2:Length
        tmp += coef_row.coef[i]*u[Start+i-1]
    end
    tmp
end

@inline function convolve_right_row(coef_row::DerivativeCoefficientRow{T,Start,Length}, u) where {T,Start,Length}
    @inbounds tmp = coef_row.coef[1]*u[end-Start+1]
    @inbounds for i in 2:Length
        tmp += coef_row.coef[i]*u[end-(Start+i-2)]
    end
    tmp
end


@unroll function convolve_boundary_coefficients!(dest, left_boundary, right_boundary, u, α, β)
    @unroll for i in 1:length(left_boundary)
        tmp = convolve_left_row(left_boundary[i], u)
        @inbounds dest[i] = β*dest[i] + α*tmp
    end
    @unroll for i in 1:length(right_boundary)
        tmp = convolve_right_row(right_boundary[i], u)
        @inbounds dest[end-i+1] = β*dest[end-i+1] + α*tmp
    end
end

@unroll function convolve_boundary_coefficients!(dest, left_boundary, right_boundary, u, α)
    @unroll for i in 1:length(left_boundary)
        tmp = convolve_left_row(left_boundary[i], u)
        @inbounds dest[i] = α*tmp
    end
    @unroll for i in 1:length(right_boundary)
        tmp = convolve_right_row(right_boundary[i], u)
        @inbounds dest[end-i+1] = α*tmp
    end
end


@generated function convolve_interior_coefficients!(dest::AbstractVector, lower_coef::SVector{LowerOffset}, central_coef, upper_coef::SVector{UpperOffset}, u::AbstractVector, α, β, left_boundary_width, right_boundary_width, parallel) where {LowerOffset, UpperOffset}
    ex = :( lower_coef[$LowerOffset]*u[i-$LowerOffset] )
    for j in LowerOffset-1:-1:1
        ex = :( $ex + lower_coef[$j]*u[i-$j] )
    end
    ex = :( $ex + central_coef*u[i] )
    for j in 1:UpperOffset
        ex = :( $ex + upper_coef[$j]*u[i+$j] )
    end

    quote
        @inbounds for i in (left_boundary_width+1):(length(dest)-right_boundary_width)
            dest[i] = β*dest[i] + α*$ex
        end
    end
end

function convolve_interior_coefficients!(dest::AbstractVector{T}, lower_coef::SVector{LowerOffset}, central_coef, upper_coef::SVector{UpperOffset}, u::AbstractVector, α, β, left_boundary_width, right_boundary_width, parallel::Val{:threads}) where {T, LowerOffset, UpperOffset}
    Threads.@threads for i in (left_boundary_width+1):(length(dest)-right_boundary_width) @inbounds begin
        tmp = zero(T)
        for j in Base.OneTo(LowerOffset)
            tmp += lower_coef[j]*u[i-j]
        end
        tmp += central_coef*u[i]
        for j in Base.OneTo(UpperOffset)
            tmp += upper_coef[j]*u[i+j]
        end
        dest[i] = β*dest[i] + α*tmp
    end end
end

@generated function convolve_interior_coefficients!(dest::AbstractVector, lower_coef::SVector{LowerOffset}, central_coef, upper_coef::SVector{UpperOffset}, u::AbstractVector, α, left_boundary_width, right_boundary_width, parallel) where {LowerOffset, UpperOffset}
    ex = :( lower_coef[$LowerOffset]*u[i-$LowerOffset] )
    for j in LowerOffset-1:-1:1
        ex = :( $ex + lower_coef[$j]*u[i-$j] )
    end
    ex = :( $ex + central_coef*u[i] )
    for j in 1:UpperOffset
        ex = :( $ex + upper_coef[$j]*u[i+$j] )
    end

    quote
        @inbounds for i in (left_boundary_width+1):(length(dest)-right_boundary_width)
            dest[i] = α*$ex
        end
    end
end

function convolve_interior_coefficients!(dest::AbstractVector{T}, lower_coef::SVector{LowerOffset}, central_coef, upper_coef::SVector{UpperOffset}, u::AbstractVector, α, left_boundary_width, right_boundary_width, parallel::Val{:threads}) where {T, LowerOffset, UpperOffset}
    Threads.@threads for i in (left_boundary_width+1):(length(dest)-right_boundary_width) @inbounds begin
        tmp = zero(T)
        for j in Base.OneTo(LowerOffset)
            tmp += lower_coef[j]*u[i-j]
        end
        tmp += central_coef*u[i]
        for j in Base.OneTo(UpperOffset)
            tmp += upper_coef[j]*u[i+j]
        end
        dest[i] = α*tmp
    end end
end



"""
    DerivativeOperator{T,StencilWidth,Parallel}

A derivative operator on a finite difference grid.
"""
struct DerivativeOperator{T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients,Grid} <: AbstractDerivativeOperator{T}
    coefficients::DerivativeCoefficients{T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients}
    grid::Grid
    factor::T

    function DerivativeOperator(coefficients::DerivativeCoefficients{T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients}, grid::Grid) where {T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients,Grid}
        @argcheck length(grid) > LowerOffset+UpperOffset DimensionMismatch
        @argcheck length(grid) > length(coefficients.left_boundary) + length(coefficients.right_boundary) DimensionMismatch

        Δx = step(grid)
        factor = inv(Δx^derivative_order(coefficients))
        new{T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients,Grid}(coefficients, grid, factor)
    end
end


@inline source_of_coeffcients(D::DerivativeOperator) = source_of_coeffcients(D.coefficients)


function Base.show(io::IO, D::DerivativeOperator{T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel}) where {T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel}
    if derivative_order(D) == 1
        print(io, "SBP 1st derivative operator of order ")
    elseif  derivative_order(D) == 2
        print(io, "SBP 2nd derivative operator of order ")
    elseif  derivative_order(D) == 3
        print(io, "SBP 3rd derivative operator of order ")
    else
        print(io, "SBP ", derivative_order(D), "th derivative operator of order ")
    end
    print(io, accuracy_order(D), " {T=", T, ", Parallel=", Parallel, "} \n")
    print(io, "on a grid in [", first(grid(D)), ", ", last(grid(D)),
                "] using ", length(grid(D)), " nodes \n")
    print(io, "and coefficients given in \n")
    print(io, source_of_coeffcients(D))
end


Base.@propagate_inbounds function Base.A_mul_B!(dest, D::DerivativeOperator, u)
    mul!(dest, D, u, one(eltype(dest)))
end

@noinline function *(D::DerivativeOperator, u)
    T = promote_type(eltype(D), eltype(u))
    dest = similar(u, T); fill!(dest, zero(T))
    @inbounds A_mul_B!(dest, D, u)
    dest
end

"""
    mul!(dest::AbstractVector, D::DerivativeOperator, u::AbstractVector, α, β)

Compute `α*D*u + β*dest` and store the result in `dest`.
"""
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::DerivativeOperator, u::AbstractVector, α, β)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, α*D.factor, β)
end

"""
    mul!(dest::AbstractVector, D::erivativeOperator, u::AbstractVector, α)

Compute `α*D*u` and store the result in `dest`.
"""
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::DerivativeOperator, u::AbstractVector, α)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, α*D.factor)
end


"""
    derivative_operator(source_of_coefficients, derivative_order, accuracy_order,
                        xmin, xmax, N, parallel=Val{:serial}())

Create a `DerivativeOperator` approximating the `derivative_order`-th derivative
on a grid between `xmin` and `xmax` with `N` grid points up to order of accuracy
`accuracy_order`. with coefficients given by `source_of_coefficients`.
The evaluation of the derivative can be parallised using threads by chosing
`parallel=Val{:threads}())`.
"""
function derivative_operator(source_of_coefficients, derivative_order, accuracy_order, xmin, xmax, N, parallel=Val{:serial}())
    grid = construct_grid(source_of_coefficients, accuracy_order, xmin, xmax, N)
    if derivative_order == 1
        coefficients = first_derivative_coefficients(source_of_coefficients, accuracy_order, eltype(grid), parallel)
    elseif derivative_order == 2
        coefficients = second_derivative_coefficients(source_of_coefficients, accuracy_order, eltype(grid), parallel)
    elseif derivative_order == 3
        coefficients = third_derivative_coefficients(source_of_coefficients, accuracy_order, eltype(grid), parallel)
    elseif derivative_order == 4
        coefficients = fourth_derivative_coefficients(source_of_coefficients, accuracy_order, eltype(grid), parallel)
    else
        throw(ArgumentError("Derivative order $derivative_order not implemented."))
    end
    DerivativeOperator(coefficients, grid)
end


@inline construct_grid(source_of_coefficients, accuracy_order, xmin, xmax, N) = linspace(xmin, xmax, N)
