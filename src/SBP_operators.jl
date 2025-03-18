
"""
    DerivativeCoefficients

The coefficients of a derivative operator on a nonperiodic grid.
"""
@auto_hash_equals struct DerivativeCoefficients{
    T,
    LeftBoundary,
    RightBoundary,
    LeftBoundaryDerivatives,
    RightBoundaryDerivatives,
    LowerOffset,
    UpperOffset,
    LeftWidth,
    RightWidth,
    ExecutionMode,
    SourceOfCoefficients,
} <: AbstractDerivativeCoefficients{T}
    # coefficients defining the operator and its action
    left_boundary::LeftBoundary
    right_boundary::RightBoundary
    left_boundary_derivatives::LeftBoundaryDerivatives
    right_boundary_derivatives::RightBoundaryDerivatives
    lower_coef::SVector{LowerOffset,T}
    central_coef::T
    upper_coef::SVector{UpperOffset,T}
    left_weights::SVector{LeftWidth,T}
    right_weights::SVector{RightWidth,T}
    mode::ExecutionMode
    # corresponding orders etc.
    derivative_order::Int
    accuracy_order::Int
    symmetric::Bool
    source_of_coefficients::SourceOfCoefficients

    function DerivativeCoefficients(
        left_boundary::LeftBoundary,
        right_boundary::RightBoundary,
        left_boundary_derivatives::LeftBoundaryDerivatives,
        right_boundary_derivatives::RightBoundaryDerivatives,
        lower_coef::SVector{LowerOffset,T},
        central_coef::T,
        upper_coef::SVector{UpperOffset,T},
        left_weights::SVector{LeftWidth,T},
        right_weights::SVector{RightWidth,T},
        mode::ExecutionMode,
        derivative_order::Int,
        accuracy_order::Int,
        source_of_coefficients::SourceOfCoefficients,
    ) where {
        T,
        LeftBoundary,
        RightBoundary,
        LeftBoundaryDerivatives,
        RightBoundaryDerivatives,
        LowerOffset,
        UpperOffset,
        LeftWidth,
        RightWidth,
        ExecutionMode,
        SourceOfCoefficients,
    }
        @argcheck length(left_boundary) == LeftWidth
        @argcheck length(right_boundary) == RightWidth
        @argcheck length(left_boundary_derivatives) == length(right_boundary_derivatives)
        symmetric = LowerOffset == UpperOffset
        if symmetric
            @inbounds for i in Base.OneTo(LowerOffset)
                symmetric = symmetric && lower_coef[i] == upper_coef[i]
            end
        end
        symmetric = symmetric && length(left_boundary[end]) == length(left_boundary)
        #TODO: check boundary coefficients
        if derivative_order - 1 != length(left_boundary_derivatives)
            @warn(
                "Derivative coefficients of degree $derivative_order should provide $(derivative_order-1) boundary derivatives."
            )
        end
        new{
            T,
            LeftBoundary,
            RightBoundary,
            LeftBoundaryDerivatives,
            RightBoundaryDerivatives,
            LowerOffset,
            UpperOffset,
            LeftWidth,
            RightWidth,
            ExecutionMode,
            SourceOfCoefficients,
        }(
            left_boundary,
            right_boundary,
            left_boundary_derivatives,
            right_boundary_derivatives,
            lower_coef,
            central_coef,
            upper_coef,
            left_weights,
            right_weights,
            mode,
            derivative_order,
            accuracy_order,
            symmetric,
            source_of_coefficients,
        )
    end
end


@inline source_of_coefficients(coefficients::AbstractDerivativeCoefficients) =
    coefficients.source_of_coefficients


function Base.show(io::IO, coefficients::DerivativeCoefficients)
    if derivative_order(coefficients) == 1
        print(io, "Coefficients of the first-derivative operator")
    elseif derivative_order(coefficients) == 2
        print(io, "Coefficients of the second-derivative operator")
    elseif derivative_order(coefficients) == 3
        print(io, "Coefficients of the third-derivative operator")
    else
        print(
            io,
            "Coefficients of the ",
            derivative_order(coefficients),
            "-derivative operator",
        )
    end
    print(io, " of order ", accuracy_order(coefficients), " of ")
    print(io, source_of_coefficients(coefficients))
end


# Compute `α*D*u + β*dest` and store the result in `dest`.
function mul!(
    dest::AbstractVector,
    coefficients::DerivativeCoefficients,
    u::AbstractVector,
    α,
    β,
)
    @unpack left_boundary, right_boundary, lower_coef, central_coef, upper_coef, mode =
        coefficients

    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(u) > length(lower_coef) + length(upper_coef) DimensionMismatch
        @argcheck length(u) > length(left_boundary) + length(right_boundary) DimensionMismatch
    end

    convolve_boundary_coefficients!(dest, left_boundary, right_boundary, u, α, β, mode)
    convolve_interior_coefficients!(
        dest,
        lower_coef,
        central_coef,
        upper_coef,
        u,
        α,
        β,
        static_length(left_boundary),
        static_length(right_boundary),
        mode,
    )
end

# Compute `α*D*u` and store the result in `dest`.
function mul!(
    dest::AbstractVector,
    coefficients::DerivativeCoefficients,
    u::AbstractVector,
    α,
)
    @unpack left_boundary, right_boundary, lower_coef, central_coef, upper_coef, mode =
        coefficients

    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(u) > length(lower_coef) + length(upper_coef) DimensionMismatch
        @argcheck length(u) >
                  length(coefficients.left_boundary) + length(coefficients.right_boundary) DimensionMismatch
    end

    convolve_boundary_coefficients!(dest, left_boundary, right_boundary, u, α, mode)
    convolve_interior_coefficients!(
        dest,
        lower_coef,
        central_coef,
        upper_coef,
        u,
        α,
        static_length(left_boundary),
        static_length(right_boundary),
        mode,
    )
end


"""
    DerivativeCoefficientRow{T,Start,Length}

A struct representing a row in the boundary block of an SBP derivative operator
with scalar type `T`.
"""
@auto_hash_equals struct DerivativeCoefficientRow{T,Start,Length}
    coef::SVector{Length,T}
end

function Base.length(::DerivativeCoefficientRow{T,Start,Length}) where {T,Start,Length}
    Length
end

function offset(::DerivativeCoefficientRow{T,Start,Length}) where {T,Start,Length}
    Start
end

function Base.:-(coef_row::DerivativeCoefficientRow{T,Start,Length}) where {T,Start,Length}
    DerivativeCoefficientRow{T,Start,Length}(-coef_row.coef)
end

function widening_plus(
    row1::SVector{Length1,T},
    row2::SVector{Length2,T},
) where {T,Length1,Length2}
    Length = max(Length1, Length2)
    row = SVector{Length,T}(ntuple(i -> zero(T), Length))
    for i = 1:Length1
        row = Base.setindex(row, row[i] + row1[i], i)
    end
    for i = 1:Length2
        row = Base.setindex(row, row[i] + row2[i], i)
    end
    row
end

function Base.:+(
    coef_row1::DerivativeCoefficientRow{T,Start1,Length1},
    coef_row2::DerivativeCoefficientRow{T,Start2,Length2},
) where {T,Start1,Length1,Start2,Length2}
    End = max(Start1 + Length1, Start2 + Length2)
    Start = min(Start1, Start2)
    Length = End - Start
    row = SVector{Length,T}(ntuple(i -> zero(T), Length))
    for i = 1:Length1
        j = i + Start1 - Start
        row = Base.setindex(row, row[j] + coef_row1.coef[i], j)
    end
    for i = 1:Length2
        j = i + Start2 - Start
        row = Base.setindex(row, row[j] + coef_row2.coef[i], j)
    end
    DerivativeCoefficientRow{T,Start,Length}(row)
end

function Base.:/(
    coef_row::DerivativeCoefficientRow{T,Start,Length},
    α::Number,
) where {T,Start,Length}
    DerivativeCoefficientRow{T,Start,Length}(coef_row.coef / α)
end

@inline function convolve_left_row(
    coef_row::DerivativeCoefficientRow{T,Start,Length},
    u,
    ::SafeMode,
) where {T,Start,Length}
    @inbounds tmp = coef_row.coef[1] * u[Start]
    @inbounds for i = 2:Length
        tmp += coef_row.coef[i] * u[Start+i-1]
    end
    tmp
end
@inline function convolve_left_row(
    coef_row::DerivativeCoefficientRow{T,Start,Length},
    u,
    ::Union{FastMode,ThreadedMode},
) where {T,Start,Length}
    @inbounds tmp = coef_row.coef[1] * u[Start]
    @inbounds for i = 2:Length
        @muladd tmp = tmp + coef_row.coef[i] * u[Start+i-1]
    end
    tmp
end

@inline function convolve_right_row(
    coef_row::DerivativeCoefficientRow{T,Start,Length},
    u,
    ::SafeMode,
) where {T,Start,Length}
    @inbounds tmp = coef_row.coef[1] * u[end-Start+1]
    @inbounds for i = 2:Length
        tmp += coef_row.coef[i] * u[end-(Start+i-2)]
    end
    tmp
end
@inline function convolve_right_row(
    coef_row::DerivativeCoefficientRow{T,Start,Length},
    u,
    ::Union{FastMode,ThreadedMode},
) where {T,Start,Length}
    @inbounds tmp = coef_row.coef[1] * u[end-Start+1]
    @inbounds for i = 2:Length
        @muladd tmp = tmp + coef_row.coef[i] * u[end-(Start+i-2)]
    end
    tmp
end


@unroll function convolve_boundary_coefficients!(
    dest::AbstractVector,
    left_boundary,
    right_boundary,
    u::AbstractVector,
    α,
    β,
    mode::SafeMode,
)
    @unroll for i = 1:length(left_boundary)
        tmp = convolve_left_row(left_boundary[i], u, mode)
        @inbounds dest[i] = β * dest[i] + α * tmp
    end
    @unroll for i = 1:length(right_boundary)
        tmp = convolve_right_row(right_boundary[i], u, mode)
        @inbounds dest[end-i+1] = β * dest[end-i+1] + α * tmp
    end
end
@unroll function convolve_boundary_coefficients!(
    dest::AbstractVector,
    left_boundary,
    right_boundary,
    u::AbstractVector,
    α,
    β,
    mode::Union{FastMode,ThreadedMode},
)
    @unroll for i = 1:length(left_boundary)
        tmp = convolve_left_row(left_boundary[i], u, mode)
        @inbounds @muladd dest[i] = β * dest[i] + α * tmp
    end
    @unroll for i = 1:length(right_boundary)
        tmp = convolve_right_row(right_boundary[i], u, mode)
        @inbounds @muladd dest[end-i+1] = β * dest[end-i+1] + α * tmp
    end
end

@unroll function convolve_boundary_coefficients!(
    dest::AbstractVector,
    left_boundary,
    right_boundary,
    u::AbstractVector,
    α,
    mode::SafeMode,
)
    @unroll for i = 1:length(left_boundary)
        tmp = convolve_left_row(left_boundary[i], u, mode)
        @inbounds dest[i] = α * tmp
    end
    @unroll for i = 1:length(right_boundary)
        tmp = convolve_right_row(right_boundary[i], u, mode)
        @inbounds dest[end-i+1] = α * tmp
    end
end
@unroll function convolve_boundary_coefficients!(
    dest::AbstractVector,
    left_boundary,
    right_boundary,
    u::AbstractVector,
    α,
    mode::Union{FastMode,ThreadedMode},
)
    @unroll for i = 1:length(left_boundary)
        tmp = convolve_left_row(left_boundary[i], u, mode)
        @inbounds @muladd dest[i] = α * tmp
    end
    @unroll for i = 1:length(right_boundary)
        tmp = convolve_right_row(right_boundary[i], u, mode)
        @inbounds @muladd dest[end-i+1] = α * tmp
    end
end


@generated function convolve_interior_coefficients!(
    dest::AbstractVector,
    lower_coef::SVector{LowerOffset},
    central_coef,
    upper_coef::SVector{UpperOffset},
    u::AbstractVector,
    α,
    β,
    ::StaticInt{left_boundary_width},
    ::StaticInt{right_boundary_width},
    mode,
) where {LowerOffset,UpperOffset,left_boundary_width,right_boundary_width}
    if LowerOffset > 0
        ex = :(lower_coef[$LowerOffset] * u[i-$LowerOffset])
        for j = LowerOffset-1:-1:1
            ex = :($ex + lower_coef[$j] * u[i-$j])
        end
        ex = :($ex + central_coef * u[i])
    else
        ex = :(central_coef * u[i])
    end
    for j = 1:UpperOffset
        ex = :($ex + upper_coef[$j] * u[i+$j])
    end

    if mode <: ThreadedMode
        quote
            Base.@_inline_meta
            @tturbo warn_check_args = false for i =
                                                (left_boundary_width+1):(length(
                dest,
            )-right_boundary_width)
                dest[i] = β * dest[i] + α * $ex
            end
        end
    elseif mode <: SafeMode
        quote
            Base.@_inline_meta
            @inbounds for i = (left_boundary_width+1):(length(dest)-right_boundary_width)
                dest[i] = β * dest[i] + α * $ex
            end
        end
    else
        quote
            Base.@_inline_meta
            @turbo warn_check_args = false for i =
                                               (left_boundary_width+1):(length(
                dest,
            )-right_boundary_width)
                dest[i] = β * dest[i] + α * $ex
            end
        end
    end
end

@generated function convolve_interior_coefficients!(
    dest::AbstractVector,
    lower_coef::SVector{LowerOffset},
    central_coef,
    upper_coef::SVector{UpperOffset},
    u::AbstractVector,
    α,
    ::StaticInt{left_boundary_width},
    ::StaticInt{right_boundary_width},
    mode,
) where {LowerOffset,UpperOffset,left_boundary_width,right_boundary_width}
    if LowerOffset > 0
        ex = :(lower_coef[$LowerOffset] * u[i-$LowerOffset])
        for j = LowerOffset-1:-1:1
            ex = :($ex + lower_coef[$j] * u[i-$j])
        end
        ex = :($ex + central_coef * u[i])
    else
        ex = :(central_coef * u[i])
    end
    for j = 1:UpperOffset
        ex = :($ex + upper_coef[$j] * u[i+$j])
    end

    if mode <: ThreadedMode
        quote
            Base.@_inline_meta
            @tturbo warn_check_args = false for i =
                                                (firstindex(dest)+$left_boundary_width):(lastindex(
                dest,
            )-$right_boundary_width)
                dest[i] = α * $ex
            end
        end
    elseif mode <: SafeMode
        quote
            Base.@_inline_meta
            @simd ivdep for i =
                            (firstindex(dest)+$left_boundary_width):(lastindex(
                dest,
            )-$right_boundary_width)
                @inbounds dest[i] = α * $ex
            end
        end
    else
        quote
            Base.@_inline_meta
            @turbo warn_check_args = false for i =
                                               (firstindex(dest)+$left_boundary_width):(lastindex(
                dest,
            )-$right_boundary_width)
                dest[i] = α * $ex
            end
        end
    end
end


@static if VERSION >= v"1.6" # `reinterpret(reshape, ...)` was introduced in Julia v1.6
    # Specialized for vectors of `StaticVector`s
    @generated function convolve_interior_coefficients!(
        _dest::AbstractVector{<:StaticVector{N,T1}},
        lower_coef::SVector{LowerOffset},
        central_coef,
        upper_coef::SVector{UpperOffset},
        _u::AbstractVector{<:StaticVector{N,T2}},
        α,
        β,
        ::StaticInt{left_boundary_width},
        ::StaticInt{right_boundary_width},
        mode,
    ) where {LowerOffset,UpperOffset,N,T1,T2,left_boundary_width,right_boundary_width}
        if LowerOffset > 0
            ex = :(lower_coef[$LowerOffset] * u[v, i-$LowerOffset])
            for j = LowerOffset-1:-1:1
                ex = :($ex + lower_coef[$j] * u[v, i-$j])
            end
            ex = :($ex + central_coef * u[v, i])
        else
            ex = :(central_coef * u[v, i])
        end
        for j = 1:UpperOffset
            ex = :($ex + upper_coef[$j] * u[v, i+$j])
        end

        if N == 1
            ex_dest = :(reshape(reinterpret($T1, _dest), 1, length(_dest)))
            ex_u = :(reshape(reinterpret($T2, _u), 1, length(_u)))
        else
            ex_dest = :(reinterpret(reshape, $T1, _dest))
            ex_u = :(reinterpret(reshape, $T2, _u))
        end

        if mode <: ThreadedMode
            quote
                Base.@_inline_meta
                dest = $ex_dest
                u = $ex_u
                fi = firstindex(dest, 2) + left_boundary_width
                la = lastindex(dest, 2) - right_boundary_width
                @tturbo warn_check_args = false for i = fi:la
                    for v in LoopVectorization.indices((dest, u), (1, 1))
                        dest[v, i] = β * dest[v, i] + α * $ex
                    end
                end
            end
        elseif mode <: SafeMode
            quote
                Base.@_inline_meta
                dest = $ex_dest
                u = $ex_u
                fi = firstindex(dest, 2) + left_boundary_width
                la = lastindex(dest, 2) - right_boundary_width
                @inbounds for i = fi:la
                    for v in LoopVectorization.indices((dest, u), (1, 1))
                        dest[v, i] = β * dest[v, i] + α * $ex
                    end
                end
            end
        else
            quote
                Base.@_inline_meta
                dest = $ex_dest
                u = $ex_u
                fi = firstindex(dest, 2) + left_boundary_width
                la = lastindex(dest, 2) - right_boundary_width
                @turbo warn_check_args = false for i = fi:la
                    for v in LoopVectorization.indices((dest, u), (1, 1))
                        dest[v, i] = β * dest[v, i] + α * $ex
                    end
                end
            end
        end
    end

    @generated function convolve_interior_coefficients!(
        _dest::AbstractVector{<:StaticVector{N,T1}},
        lower_coef::SVector{LowerOffset},
        central_coef,
        upper_coef::SVector{UpperOffset},
        _u::AbstractVector{<:StaticVector{N,T2}},
        α,
        ::StaticInt{left_boundary_width},
        ::StaticInt{right_boundary_width},
        mode,
    ) where {LowerOffset,UpperOffset,N,T1,T2,left_boundary_width,right_boundary_width}
        if LowerOffset > 0
            ex = :(lower_coef[$LowerOffset] * u[v, i-$LowerOffset])
            for j = LowerOffset-1:-1:1
                ex = :($ex + lower_coef[$j] * u[v, i-$j])
            end
            ex = :($ex + central_coef * u[v, i])
        else
            ex = :(central_coef * u[v, i])
        end
        for j = 1:UpperOffset
            ex = :($ex + upper_coef[$j] * u[v, i+$j])
        end

        if N == 1
            ex_dest = :(reshape(reinterpret($T1, _dest), 1, length(_dest)))
            ex_u = :(reshape(reinterpret($T2, _u), 1, length(_u)))
        else
            ex_dest = :(reinterpret(reshape, $T1, _dest))
            ex_u = :(reinterpret(reshape, $T2, _u))
        end

        if mode <: ThreadedMode
            quote
                Base.@_inline_meta
                dest = $ex_dest
                u = $ex_u
                fi = firstindex(dest, 2) + left_boundary_width
                la = lastindex(dest, 2) - right_boundary_width
                @tturbo warn_check_args = false for i = fi:la
                    for v in LoopVectorization.indices((dest, u), (1, 1))
                        dest[v, i] = α * $ex
                    end
                end
            end
        elseif mode <: SafeMode
            quote
                Base.@_inline_meta
                dest = $ex_dest
                u = $ex_u
                fi = firstindex(dest, 2) + left_boundary_width
                la = lastindex(dest, 2) - right_boundary_width
                @inbounds for i = fi:la
                    for v in LoopVectorization.indices((dest, u), (1, 1))
                        dest[v, i] = α * $ex
                    end
                end
            end
        else
            quote
                Base.@_inline_meta
                dest = $ex_dest
                u = $ex_u
                fi = firstindex(dest, 2) + left_boundary_width
                la = lastindex(dest, 2) - right_boundary_width
                @turbo warn_check_args = false for i = fi:la
                    for v in LoopVectorization.indices((dest, u), (1, 1))
                        dest[v, i] = α * $ex
                    end
                end
            end
        end
    end
end


"""
    DerivativeOperator

A derivative operator on a nonperiodic finite difference grid.
See [`derivative_operator`](@ref).
"""
@auto_hash_equals struct DerivativeOperator{
    T,
    LeftBoundary,
    RightBoundary,
    LeftBoundaryDerivatives,
    RightBoundaryDerivatives,
    LowerOffset,
    UpperOffset,
    LeftWidth,
    RightWidth,
    ExecutionMode,
    SourceOfCoefficients,
    Grid,
} <: AbstractNonperiodicDerivativeOperator{T}
    coefficients::DerivativeCoefficients{
        T,
        LeftBoundary,
        RightBoundary,
        LeftBoundaryDerivatives,
        RightBoundaryDerivatives,
        LowerOffset,
        UpperOffset,
        LeftWidth,
        RightWidth,
        ExecutionMode,
        SourceOfCoefficients,
    }
    grid::Grid
    Δx::T
    factor::T

    function DerivativeOperator(
        coefficients::DerivativeCoefficients{
            T,
            LeftBoundary,
            RightBoundary,
            LeftBoundaryDerivatives,
            RightBoundaryDerivatives,
            LowerOffset,
            UpperOffset,
            LeftWidth,
            RightWidth,
            ExecutionMode,
            SourceOfCoefficients,
        },
        grid::Grid,
    ) where {
        T,
        LeftBoundary,
        RightBoundary,
        LeftBoundaryDerivatives,
        RightBoundaryDerivatives,
        LowerOffset,
        UpperOffset,
        LeftWidth,
        RightWidth,
        ExecutionMode,
        SourceOfCoefficients,
        Grid,
    }
        @argcheck length(grid) > LowerOffset + UpperOffset DimensionMismatch
        @argcheck length(grid) >
                  length(coefficients.left_boundary) + length(coefficients.right_boundary) DimensionMismatch

        Δx = step(grid)
        factor = inv(Δx^derivative_order(coefficients))
        new{
            T,
            LeftBoundary,
            RightBoundary,
            LeftBoundaryDerivatives,
            RightBoundaryDerivatives,
            LowerOffset,
            UpperOffset,
            LeftWidth,
            RightWidth,
            ExecutionMode,
            SourceOfCoefficients,
            Grid,
        }(
            coefficients,
            grid,
            Δx,
            factor,
        )
    end
end


@inline source_of_coefficients(D::DerivativeOperator) =
    source_of_coefficients(D.coefficients)


function Base.show(io::IO, D::DerivativeOperator)
    if derivative_order(D) == 1
        print(io, "SBP first-derivative operator")
    elseif derivative_order(D) == 2
        print(io, "SBP second-derivative operator")
    elseif derivative_order(D) == 3
        print(io, "SBP third-derivative operator")
    else
        print(io, "SBP ", derivative_order(D), "-derivative operator")
    end
    print(io, " of order ", accuracy_order(D))
    if get(io, :compact, false) == false
        print(
            io,
            " on a grid in [",
            first(grid(D)),
            ", ",
            last(grid(D)),
            "] using ",
            length(grid(D)),
            " nodes \n",
        )
        print(io, "and coefficients")
    end
    print(io, " of ", source_of_coefficients(D))
end


"""
    derivative_left(D::DerivativeOperator, u, der_order::Val{N})

Compute the `N`-th derivative of the function given by the coefficients `u` at
the left boundary of the grid.
"""
@inline function derivative_left(D::DerivativeOperator, u, der_order::Val{N}) where {N}
    convolve_left_row(D.coefficients.left_boundary_derivatives[N], u, D.coefficients.mode) /
    D.Δx^N
end

@inline function derivative_left(
    D::AbstractNonperiodicDerivativeOperator,
    u,
    der_order::Val{0},
)
    u[begin]
end
@inline function derivative_left(D::DerivativeOperator, u, der_order::Val{0})
    u[begin]
end

"""
    derivative_right(D::DerivativeOperator, u, der_order::Val{N})

Compute the `N`-th derivative of the function given by the coefficients `u` at
the right boundary of the grid.
"""
@inline function derivative_right(D::DerivativeOperator, u, der_order::Val{N}) where {N}
    convolve_right_row(
        D.coefficients.right_boundary_derivatives[N],
        u,
        D.coefficients.mode,
    ) / D.Δx^N
end

@inline function derivative_right(
    D::AbstractNonperiodicDerivativeOperator,
    u,
    der_order::Val{0},
)
    u[end]
end
@inline function derivative_right(D::DerivativeOperator, u, der_order::Val{0})
    u[end]
end


"""
    mul!(du, D::DerivativeOperator, u, α=true, β=false)

Efficient in-place version of `du = α * D * u + β * du`. Note that `du` must not
be aliased with `u`.
"""
function mul! end


"""
    mul_transpose_derivative_left!(u, D::DerivativeOperator, der_order::Val{N}, α=true, β=false)

Set the grid function `u` to `α` times the transposed `N`-th derivative functional
applied to `u` plus `β` times `u` in the domain of the `N`-th derivative functional
at the left boundary of the grid.
Thus, the coefficients `α, β` have the same meaning as in [`mul!`](@ref).
"""
@inline function mul_transpose_derivative_left!(
    u::AbstractVector,
    D::DerivativeOperator,
    der_order::Val{N},
    α = true,
    β = false,
) where {N}
    factor = α / D.Δx^N
    coef = D.coefficients.left_boundary_derivatives[N].coef
    for i in eachindex(coef)
        u[i] = factor * coef[i] + β * u[i]
    end
end

@inline function mul_transpose_derivative_left!(
    u::AbstractVector,
    D::DerivativeOperator,
    der_order::Val{0},
    α = true,
    β = false,
)
    @inbounds u[begin] = α * u[begin] + β * u[begin]
    return nothing
end

"""
    mul_transpose_derivative_right!(u, D::DerivativeOperator, der_order::Val{N}, α=true, β=false)

Set the grid function `u` to `α` times the transposed `N`-th derivative functional
applied to `u` plus `β` times `u` in the domain of the `N`-th derivative functional
at the right boundary of the grid.
Thus, the coefficients `α, β` have the same meaning as in [`mul!`](@ref).
"""
@inline function mul_transpose_derivative_right!(
    u::AbstractVector,
    D::DerivativeOperator,
    der_order::Val{N},
    α = true,
    β = false,
) where {N}
    factor = α / D.Δx^N
    coef = D.coefficients.right_boundary_derivatives[N].coef
    for i in eachindex(coef)
        u[end-i+1] += factor * coef[i]
    end
end

@inline function mul_transpose_derivative_right!(
    u::AbstractVector,
    D::DerivativeOperator,
    der_order::Val{0},
    α = true,
    β = false,
)
    @inbounds u[end] = α * u[end] + β * u[end]
    return nothing
end


@auto_hash_equals struct OneHot <: AbstractVector{Bool}
    i::Int
    length::Int
end

Base.getindex(v::OneHot, i::Integer) = i == v.i
Base.size(v::OneHot) = (v.length,)

"""
    derivative_left(D::AbstractNonperiodicDerivativeOperator, der_order)

Get a representation of the linear functional evaluation the `N`th derivative at
the left boundary node as (dense) vector.
"""
function derivative_left(D::AbstractNonperiodicDerivativeOperator, der_order)
    [derivative_left(D, OneHot(i, size(D, 2)), der_order) for i = 1:size(D, 2)]
end

"""
    derivative_right(D::AbstractNonperiodicDerivativeOperator, der_order)

Get a representation of the linear functional evaluation the `N`th derivative at
the right boundary node as (dense) vector.
"""
function derivative_right(D::AbstractNonperiodicDerivativeOperator, der_order)
    [derivative_right(D, OneHot(i, size(D, 2)), der_order) for i = 1:size(D, 2)]
end



# Compute `α*D*u + β*dest` and store the result in `dest`.
Base.@propagate_inbounds function mul!(
    dest::AbstractVector,
    D::DerivativeOperator,
    u::AbstractVector,
    α,
    β,
)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, α * D.factor, β)
end

# Compute `α*D*u` and store the result in `dest`.
Base.@propagate_inbounds function mul!(
    dest::AbstractVector,
    D::DerivativeOperator,
    u::AbstractVector,
    α,
)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, α * D.factor)
end


"""
    integrate(func, u, D::DerivativeOperator)

Map the function `func` to the coefficients `u` and integrate with respect to
the quadrature rule associated with the SBP derivative operator `D`.
"""
function integrate(func, u::AbstractVector, D::DerivativeOperator)
    @boundscheck begin
        length(u) == length(grid(D)) ||
            throw(DimensionMismatch("sizes of input vector and operator do not match"))
    end
    @unpack Δx = D
    @unpack left_weights, right_weights = D.coefficients

    @inbounds res =
        sum(func, view(u, 1+length(left_weights):length(u)-length(right_weights)))
    @inbounds for i in Base.OneTo(length(left_weights))
        res += left_weights[i] * func(u[i])
    end
    @inbounds for i in Base.OneTo(length(right_weights))
        res += right_weights[i] * func(u[end-i+1])
    end

    Δx * res
end

function scale_by_mass_matrix!(u::AbstractVector, D::DerivativeOperator)
    Base.require_one_based_indexing(u)
    @boundscheck begin
        length(u) == size(D, 2) ||
            throw(DimensionMismatch("sizes of input vector and operator do not match"))
    end
    @unpack Δx = D
    @unpack left_weights, right_weights = D.coefficients

    @inbounds for i = (1+length(left_weights)):(length(u)-length(right_weights))
        u[i] *= Δx
    end
    @inbounds for i in Base.OneTo(length(left_weights))
        u[i] *= left_weights[i] * Δx
    end
    @inbounds for i in Base.OneTo(length(right_weights))
        u[end-i+1] *= right_weights[i] * Δx
    end

    return u
end

function scale_by_inverse_mass_matrix!(u::AbstractVector, D::DerivativeOperator)
    Base.require_one_based_indexing(u)
    @boundscheck begin
        length(u) == size(D, 2) ||
            throw(DimensionMismatch("sizes of input vector and operator do not match"))
    end
    @unpack Δx = D
    @unpack left_weights, right_weights = D.coefficients

    @inbounds for i = (1+length(left_weights)):(length(u)-length(right_weights))
        u[i] /= Δx
    end
    @inbounds for i in Base.OneTo(length(left_weights))
        u[i] /= left_weights[i] * Δx
    end
    @inbounds for i in Base.OneTo(length(right_weights))
        u[end-i+1] /= right_weights[i] * Δx
    end

    return u
end


"""
    derivative_operator(source_of_coefficients,
                        derivative_order, accuracy_order,
                        xmin, xmax, N, mode=FastMode())
    derivative_operator(source_of_coefficients;
                        derivative_order, accuracy_order,
                        xmin, xmax, N, mode=FastMode())

Create a [`DerivativeOperator`](@ref) approximating the `derivative_order`-th derivative
on a grid between `xmin` and `xmax` with `N` grid points up to order of accuracy
`accuracy_order`. with coefficients given by `source_of_coefficients`.
The evaluation of the derivative can be parallelized using threads by choosing
`mode=ThreadedMode()`.
"""
function derivative_operator(
    source_of_coefficients,
    derivative_order,
    accuracy_order,
    xmin,
    xmax,
    N,
    mode = FastMode(),
)
    if mode === Val(:serial) || mode === Val(:threads)
        # TODO: deprecated in v0.5
        Base.depwarn(
            "Providing the argument `parallel` is deprecated." * "Use `mode` instead.",
            :derivative_operator,
        )
        mode = _parallel_to_mode(mode)
    end
    grid = construct_grid(source_of_coefficients, accuracy_order, xmin, xmax, N)
    if derivative_order == 1
        coefficients = first_derivative_coefficients(
            source_of_coefficients,
            accuracy_order,
            eltype(grid),
            mode,
        )
    elseif derivative_order == 2
        coefficients = second_derivative_coefficients(
            source_of_coefficients,
            accuracy_order,
            eltype(grid),
            mode,
        )
    elseif derivative_order == 3
        coefficients = third_derivative_coefficients(
            source_of_coefficients,
            accuracy_order,
            eltype(grid),
            mode,
        )
    elseif derivative_order == 4
        coefficients = fourth_derivative_coefficients(
            source_of_coefficients,
            accuracy_order,
            eltype(grid),
            mode,
        )
    else
        throw(ArgumentError("Derivative order $derivative_order not implemented."))
    end
    DerivativeOperator(coefficients, grid)
end

function derivative_operator(
    source_of_coefficients;
    derivative_order,
    accuracy_order,
    xmin,
    xmax,
    N,
    mode = FastMode(),
    parallel = nothing,
)
    if parallel !== nothing
        # TODO: deprecated in v0.5
        Base.depwarn(
            "Providing the keyword argument `parallel` is deprecated." *
            "Use `mode` instead.",
            :derivative_operator,
        )
        mode = _parallel_to_mode(parallel)
    end
    derivative_operator(
        source_of_coefficients,
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
        mode,
    )
end


@inline construct_grid(source_of_coefficients, accuracy_order, xmin, xmax, N) =
    range(xmin, stop = xmax, length = N)


function lower_bandwidth(D::DerivativeOperator)
    @unpack left_boundary, right_boundary, lower_coef = D.coefficients

    l = length(lower_coef)
    for (i, coef_row) in enumerate(left_boundary)
        l = max(l, i - offset(coef_row))
    end
    for (i, coef_row) in enumerate(right_boundary)
        l = max(l, 1 + length(coef_row) - offset(coef_row) - i)
    end
    l
end

function upper_bandwidth(D::DerivativeOperator)
    @unpack left_boundary, right_boundary, upper_coef = D.coefficients

    u = length(upper_coef)
    for (i, coef_row) in enumerate(left_boundary)
        u = max(u, 1 + length(coef_row) - offset(coef_row) - i)
    end
    for (i, coef_row) in enumerate(right_boundary)
        u = max(u, i - offset(coef_row))
    end
    u
end

function left_boundary_weight(D::DerivativeOperator)
    @inbounds retval = D.Δx * D.coefficients.left_weights[1]
    retval
end

function right_boundary_weight(D::DerivativeOperator)
    @inbounds retval = D.Δx * D.coefficients.right_weights[1]
    retval
end
