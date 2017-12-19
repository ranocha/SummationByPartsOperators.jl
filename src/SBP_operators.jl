
"""
    DerivativeCoefficients{T,BoundaryWidth,LowerOffset,UpperOffset,Parallel}

The coefficients of a derivative operator on a periodic grid.
"""
struct DerivativeCoefficients{T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel} <: AbstractDerivativeCoefficients{T}
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

    function DerivativeCoefficients(left_boundary::LeftBoundary, right_boundary::RightBoundary, lower_coef::SVector{LowerOffset, T}, central_coef::T, upper_coef::SVector{UpperOffset, T}, parallel::Parallel, derivative_order::Int, accuracy_order::Int) where {T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel}
        symmetric = LowerOffset == UpperOffset
        if symmetric
            @inbounds for i in Base.OneTo(LowerOffset)
                symmetric = symmetric && lower_coef[i] == upper_coef[i]
            end
        end
        #TODO: check boundary coefficients
        new{T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel}(left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel, derivative_order, accuracy_order, symmetric)
    end
end


"""
    mul!(dest::AbstractVector, coefficients::DerivativeCoefficients, u::AbstractVector, α, β)

Compute `α*D*u + β*dest` and store the result in `dest`.
"""
function mul!(dest::AbstractVector, coefficients::DerivativeCoefficients{T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel}, u::AbstractVector, α, β) where {T,LeftBoundary,RightBoundary,LowerOffset,UpperOffset,Parallel}
    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(u) > LowerOffset+UpperOffset DimensionMismatch
        @argcheck length(u) > length(coefficients.left_boundary) + length(coefficients.right_boundary) DimensionMismatch
    end

    @unpack left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel = coefficients
    convolve_boundary_coefficients!(dest, left_boundary, right_boundary, u, α, β)
    convolve_interior_coefficients!(dest, lower_coef, central_coef, upper_coef, u, α, β, length(left_boundary), length(right_boundary), parallel)
end


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


# Coefficients
struct MattssonSvärdShoeybi2008 end

function first_derivative_coefficients(::MattssonSvärdShoeybi2008, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            DerivativeCoefficientRow{T,1,2}(SVector(-one(T), one(T))),
        )
        right_boundary = (
            -left_boundary[1],
        )
        upper_coef = SVector(T(1//2))
        central_coef = zero(T)
        lower_coef = -upper_coef

        DerivativeCoefficients(left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel, 1, order)
    elseif order == 4
        left_boundary = (
            DerivativeCoefficientRow{T,1,4}(SVector(T(-24//17), T(59//34), T(-4//17), T(-3//34))),
            DerivativeCoefficientRow{T,1,3}(SVector(T(-1//2), T(0), T(1//2))),
            DerivativeCoefficientRow{T,1,5}(SVector(T(4//43), T(-59//86), T(0), T(59//86), T(-4//43))),
            DerivativeCoefficientRow{T,1,6}(SVector(T(3//98), T(0), T(-59//98), T(0), T(32//49), T(-4//49))),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(2//3), T(-1//12))
        central_coef = zero(T)
        lower_coef = -upper_coef

        DerivativeCoefficients(left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel, 1, order)
    #=elseif order == 6
        x1 = 342523 // 518400
        left_boundary = (
            DerivativeCoefficientRow{T,1,6}(SVector(T(-21600//13649),
                                                    T(104009//54596),
                                                    T(30443//81894),
                                                    T(-33311//27298),
                                                    T(16863//27298),
                                                    T(-15025//163788))),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(3//4), T(-3//20), T(1//60))
        central_coef = zero(T)
        lower_coef = -upper_coef

        DerivativeCoefficients(left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel, 1, order)
    =#
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end
