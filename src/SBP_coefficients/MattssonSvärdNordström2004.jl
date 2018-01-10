
"""
    MattssonSvärdNordström2004

Coefficients of the SBP operators given in
  Mattsson, Svärd, Nordström (2004)
  Stable and Accurate Artificial Dissipation.
  Journal of Scientific Computing 21.1, pp. 57-79.
"""
struct MattssonSvärdNordström2004 <: SourceOfCoefficients end

function Base.show(io::IO, ::MattssonSvärdNordström2004)
    print(io,
        "  Mattsson, Svärd, Nordström (2004) \n",
        "  Stable and Accurate Artificial Dissipation. \n",
        "  Journal of Scientific Computing 21.1, pp. 57-79. \n")
end


function dissipation_coefficients(source::MattssonSvärdNordström2004, order::Int, grid, left_weights, right_weights, parallel=Val{:serial}())
    T = promote_type(eltype(grid), eltype(left_weights), eltype(right_weights))
    if order == 2
        inv_left_weights = one(T) ./ left_weights
        inv_right_weights = one(T) ./ right_weights
        coefficient_cache = MattssonSvärdNordström2004Cache2(inv_left_weights, inv_right_weights)
        b = ones(grid)
        b[1] = T(0)
    elseif order == 4
        inv_left_weights = one(T) ./ left_weights
        inv_right_weights = one(T) ./ right_weights
        coefficient_cache = MattssonSvärdNordström2004Cache4(inv_left_weights, inv_right_weights)
        b = ones(grid)
        b[1] = b[end] = T(0)
    elseif order == 6
        inv_left_weights = one(T) ./ left_weights
        inv_right_weights = one(T) ./ right_weights
        coefficient_cache = MattssonSvärdNordström2004Cache6(inv_left_weights, inv_right_weights)
        b = ones(grid)
        b[1] = b[2] = b[end] = T(0)
    elseif order == 8
        inv_left_weights = one(T) ./ left_weights
        inv_right_weights = one(T) ./ right_weights
        coefficient_cache = MattssonSvärdNordström2004Cache8(inv_left_weights, inv_right_weights)
        b = ones(grid)
        b[1] = b[2] = b[end] = b[end-1] = T(0)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end


    DissipationCoefficients(coefficient_cache, parallel, order, 2, source), b
end



struct MattssonSvärdNordström2004Cache2{T,LeftWidth,RightWidth} <: AbstractCoefficientCache{T}
    inv_left_weights::SVector{LeftWidth,T}
    inv_right_weights::SVector{RightWidth,T}

    function MattssonSvärdNordström2004Cache2(inv_left_weights::SVector{LeftWidth,T}, inv_right_weights::SVector{RightWidth,T}) where {T,LeftWidth,RightWidth}
        boundary_length = 2

        if length(inv_left_weights) < boundary_length
            inv_left_weights = SVector(inv_left_weights..., ntuple(j->one(T), boundary_length-length(inv_left_weights))...)
        end
        if length(inv_right_weights) < boundary_length
            inv_right_weights = SVector(inv_right_weights..., ntuple(j->one(T), boundary_length-length(inv_right_weights))...)
        end

        new{T,length(inv_left_weights),length(inv_right_weights)}(inv_left_weights, inv_right_weights)
    end
end

lower_bandwidth(cache::MattssonSvärdNordström2004Cache2) = 1
upper_bandwidth(cache::MattssonSvärdNordström2004Cache2) = 1

function convolve_boundary_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache2, u::AbstractVector, b::AbstractVector, α)
    @unpack inv_left_weights, inv_right_weights = cache

    @inbounds begin
        dest[1] = α * inv_left_weights[1] * (
                    (b[1] + b[2]) * u[1]
                    + (b[1] - b[2]) * u[2]
                )
        dest[2] = α * inv_left_weights[2] * (
                    (b[1] - b[2]) * u[1]
                    + (b[1] + b[2] + b[3]) * u[2]
                    - b[3] * u[3]
                )
        for i in 3:length(inv_left_weights)
            dest[i] = α * inv_left_weights[i] * (
                        -b[i] * u[i-1]
                        + (b[i] + b[i+1]) * u[i]
                        - b[i+1] * u[i+1]
                    )
        end

        dest[end] = α * inv_right_weights[1] * (
                        b[end] * u[end]
                        - b[end] * u[end-1]
                    )
        dest[end-1] = α * inv_right_weights[2] * (
                        -b[end-1] * u[end-2]
                        + (b[end-1] + b[end]) * u[end-1]
                        - b[end] * u[end]
                    )
        for i in 2:length(inv_right_weights)-1
            dest[end-i] = α * inv_right_weights[i+1] * (
                        -b[end-i] * u[end-(i-1)]
                        + (b[end-i] + b[end-(i+1)]) * u[end-i]
                        - b[end-(i+1)] * u[end-(i+1)]
                    )
        end
    end
end

function convolve_boundary_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache2, u::AbstractVector, b::AbstractVector, α, β)
    @unpack inv_left_weights, inv_right_weights = cache

    @inbounds begin
        dest[1] = α * inv_left_weights[1] * (
                    (b[1] + b[2]) * u[1]
                    + (b[1] - b[2]) * u[2]
                ) + β*dest[1]
        dest[2] = α * inv_left_weights[2] * (
                    (b[1] - b[2]) * u[1]
                    + (b[1] + b[2] + b[3]) * u[2]
                    - b[3] * u[3]
                ) + β*dest[2]
        for i in 3:length(inv_left_weights)
            dest[i] = α * inv_left_weights[i] * (
                        -b[i] * u[i-1]
                        + (b[i] + b[i+1]) * u[i]
                        - b[i+1] * u[i+1]
                    ) + β*dest[i]
        end

        dest[end] = α * inv_right_weights[1] * (
                        b[end] * u[end]
                        - b[end] * u[end-1]
                    ) + β*dest[end]
        dest[end-1] = α * inv_right_weights[2] * (
                        -b[end-1] * u[end-2]
                        + (b[end-1] + b[end]) * u[end-1]
                        - b[end] * u[end]
                    ) + β*dest[end-1]
        for i in 2:length(inv_right_weights)-1
            dest[end-i] = α * inv_right_weights[i+1] * (
                        -b[end-i] * u[end-(i-1)]
                        + (b[end-i] + b[end-(i+1)]) * u[end-i]
                        - b[end-(i+1)] * u[end-(i+1)]
                    ) + β*dest[end-i]
        end
    end
end

@inline function convolve_interior_coefficients_loopbody(i, cache::MattssonSvärdNordström2004Cache2, u, b)
    @inbounds begin
        b_i = b[i]
        b_ip1 = b[i+1]

        retval = -b_i * u[i-1] + (b_i + b_ip1) * u[i] - b_ip1 * u[i+1]
    end

    retval
end



struct MattssonSvärdNordström2004Cache4{T,LeftWidth,RightWidth} <: AbstractCoefficientCache{T}
    inv_left_weights::SVector{LeftWidth,T}
    inv_right_weights::SVector{RightWidth,T}

    function MattssonSvärdNordström2004Cache4(inv_left_weights::SVector{LeftWidth,T}, inv_right_weights::SVector{RightWidth,T}) where {T,LeftWidth,RightWidth}
        boundary_length = 3

        if length(inv_left_weights) < boundary_length
            inv_left_weights = SVector(inv_left_weights..., ntuple(j->one(T), boundary_length-length(inv_left_weights))...)
        end
        if length(inv_right_weights) < boundary_length
            inv_right_weights = SVector(inv_right_weights..., ntuple(j->one(T), boundary_length-length(inv_right_weights))...)
        end

        new{T,length(inv_left_weights),length(inv_right_weights)}(inv_left_weights, inv_right_weights)
    end
end

lower_bandwidth(cache::MattssonSvärdNordström2004Cache4) = 2
upper_bandwidth(cache::MattssonSvärdNordström2004Cache4) = 2

function convolve_boundary_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache4, u::AbstractVector, b::AbstractVector, α)
    @unpack inv_left_weights, inv_right_weights = cache

    @inbounds begin
        dest[1] = α * inv_left_weights[1] * (
                    (b[1] + b[2]) * u[1]
                    - 2*(b[1] + b[2]) * u[2]
                    + (b[1] + b[2]) * u[3]
                )
        dest[2] = α * inv_left_weights[2] * (
                    -2*(b[1] + b[2]) * u[1]
                    + (4*b[1] + 4*b[2] + b[3]) * u[2]
                    -2*(b[1] + b[2] + b[3]) * u[3]
                    + b[3] * u[4]
                )
        dest[3] = α * inv_left_weights[3] * (
                    (b[1] + b[2]) * u[1]
                    -2*(b[1] + b[2] + b[3]) * u[2]
                    + (b[1] + b[2] + 4*b[3] + b[4]) * u[3]
                    -2*(b[3] + b[4]) * u[4]
                    + b[4] * u[5]
                )
        for i in 4:length(inv_left_weights)
            dest[i] = α * inv_left_weights[i] * (
                        b[i-1] * u[i-2]
                        - 2*(b[i-1] + b[i]) * u[i-1]
                        + (b[i-1] + 4*b[i] + b[i+1]) * u[i]
                        - 2*(b[i] + b[i+1]) * u[i+1]
                        + b[i+1] * u[i+2]
                    )
        end

        dest[end] = α * inv_right_weights[1] * (
                        (b[end] + b[end-1]) * u[end]
                        -2*(b[end] + b[end-1]) * u[end-1]
                        + (b[end] + b[end-1]) * u[end-2]
                    )
        dest[end-1] = α * inv_right_weights[2] * (
                        -2*(b[end] + b[end-1]) * u[end]
                        + (b[end-2] + 4*b[end-1] + 4*b[end]) * u[end-1]
                        -2*(b[end-2] + b[end-1] + b[end]) * u[end-2]
                        + b[end-2] * u[end-3]
                    )
        dest[end-2] = α * inv_right_weights[3] * (
                        (b[end-1] + b[end]) * u[end]
                        -2*(b[end-2] + b[end-1] + b[end]) * u[end-1]
                        + (b[end-3] + 4*b[end-2] + b[end-1] + b[end]) * u[end-2]
                        -2*(b[end-3] + b[end-2]) * u[end-3]
                        + b[end-3] * u[end-4]
                    )
        for i in 3:length(inv_right_weights)-1
            dest[end-i] = α * inv_right_weights[i+1] * (
                        b[end-(i-1)] * u[end-(i-2)]
                        - 2*(b[end-(i-1)] + b[i]) * u[end-(i-1)]
                        + (b[end-(i-1)] + 4*b[end-i] + b[end-(i+1)]) * u[end-i]
                        - 2*(b[end-i] + b[end-(i+1)]) * u[end-(i+1)]
                        + b[end-(i+1)] * u[end-(i+2)]
                    )
        end
    end
end

function convolve_boundary_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache4, u::AbstractVector, b::AbstractVector, α, β)
    @unpack inv_left_weights, inv_right_weights = cache

    @inbounds begin
        dest[1] = α * inv_left_weights[1] * (
                    (b[1] + b[2]) * u[1]
                    - 2*(b[1] + b[2]) * u[2]
                    + (b[1] + b[2]) * u[3]
                ) + β*dest[1]
        dest[2] = α * inv_left_weights[2] * (
                    -2*(b[1] + b[2]) * u[1]
                    + (4*b[1] + 4*b[2] + b[3]) * u[2]
                    -2*(b[1] + b[2] + b[3]) * u[3]
                    + b[3] * u[4]
                ) + β*dest[2]
        dest[3] = α * inv_left_weights[3] * (
                    (b[1] + b[2]) * u[1]
                    -2*(b[1] + b[2] + b[3]) * u[2]
                    + (b[1] + b[2] + 4*b[3] + b[4]) * u[3]
                    -2*(b[3] + b[4]) * u[4]
                    + b[4] * u[5]
                ) + β*dest[3]
        for i in 4:length(inv_left_weights)
            dest[i] = α * inv_left_weights[i] * (
                        b[i-1] * u[i-2]
                        - 2*(b[i-1] + b[i]) * u[i-1]
                        + (b[i-1] + 4*b[i] + b[i+1]) * u[i]
                        - 2*(b[i] + b[i+1]) * u[i+1]
                        + b[i+1] * u[i+2]
                    ) + β*dest[i]
        end

        dest[end] = α * inv_right_weights[1] * (
                        (b[end] + b[end-1]) * u[end]
                        -2*(b[end] + b[end-1]) * u[end-1]
                        + (b[end] + b[end-1]) * u[end-2]
                    ) + β*dest[end]
        dest[end-1] = α * inv_right_weights[2] * (
                        -2*(b[end] + b[end-1]) * u[end]
                        + (b[end-2] + 4*b[end-1] + 4*b[end]) * u[end-1]
                        -2*(b[end-2] + b[end-1] + b[end]) * u[end-2]
                        + b[end-2] * u[end-3]
                    ) + β*dest[end-1]
        dest[end-2] = α * inv_right_weights[3] * (
                        (b[end-1] + b[end]) * u[end]
                        -2*(b[end-2] + b[end-1] + b[end]) * u[end-1]
                        + (b[end-3] + 4*b[end-2] + b[end-1] + b[end]) * u[end-2]
                        -2*(b[end-3] + b[end-2]) * u[end-3]
                        + b[end-3] * u[end-4]
                    ) + β*dest[end-2]
        for i in 3:length(inv_right_weights)-1
            dest[end-i] = α * inv_right_weights[i+1] * (
                        b[end-(i-1)] * u[end-(i-2)]
                        - 2*(b[end-(i-1)] + b[i]) * u[end-(i-1)]
                        + (b[end-(i-1)] + 4*b[end-i] + b[end-(i+1)]) * u[end-i]
                        - 2*(b[end-i] + b[end-(i+1)]) * u[end-(i+1)]
                        + b[end-(i+1)] * u[end-(i+2)]
                    ) + β*dest[end-i]
        end
    end
end

@inline function convolve_interior_coefficients_loopbody(i, cache::MattssonSvärdNordström2004Cache4, u, b)
    @inbounds begin
        b_im1 = b[i-1]
        b_i   = b[i]
        b_ip1 = b[i+1]

        retval = (
                    b_im1 * u[i-2]
                    - 2*(b_im1 + b_i) * u[i-1]
                    + (b_im1 + 4b_i + b_ip1) * u[i]
                    - 2*(b_i + b_ip1) * u[i+1]
                    + b_ip1 * u[i+2]
                )
    end

    retval
end



struct MattssonSvärdNordström2004Cache6{T,LeftWidth,RightWidth} <: AbstractCoefficientCache{T}
    inv_left_weights::SVector{LeftWidth,T}
    inv_right_weights::SVector{RightWidth,T}

    function MattssonSvärdNordström2004Cache6(inv_left_weights::SVector{LeftWidth,T}, inv_right_weights::SVector{RightWidth,T}) where {T,LeftWidth,RightWidth}
        boundary_length = 4

        if length(inv_left_weights) < boundary_length
            inv_left_weights = SVector(inv_left_weights..., ntuple(j->one(T), boundary_length-length(inv_left_weights))...)
        end
        if length(inv_right_weights) < boundary_length
            inv_right_weights = SVector(inv_right_weights..., ntuple(j->one(T), boundary_length-length(inv_right_weights))...)
        end

        new{T,length(inv_left_weights),length(inv_right_weights)}(inv_left_weights, inv_right_weights)
    end
end

lower_bandwidth(cache::MattssonSvärdNordström2004Cache6) = 3
upper_bandwidth(cache::MattssonSvärdNordström2004Cache6) = 3

function convolve_boundary_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache6, u::AbstractVector, b::AbstractVector, α)
    @unpack inv_left_weights, inv_right_weights = cache

    @inbounds begin
        dest[1] = α * inv_left_weights[1] * (
                    (b[1] + b[2] + b[3]) * u[1]
                    - 3*(b[1] + b[2] + b[3]) * u[2]
                    + 3*(b[1] + b[2] + b[3]) * u[3]
                    - (b[1] + b[2] + b[3]) * u[4]
                )
        dest[2] = α * inv_left_weights[2] * (
                    - 3*(b[1] + b[2] + b[3]) * u[1]
                    + (9*b[1] + 9*b[2] + 9*b[3] + b[4]) * u[2]
                    - 3*(3*b[1] + 3*b[2] + 3*b[3] + b[4]) * u[3]
                    + 3*(b[1] + b[2] + b[3] + b[4]) * u[4]
                    - b[4] * u[5]
                )
        dest[3] = α * inv_left_weights[3] * (
                    3*(b[1] + b[2] + b[3]) * u[1]
                    - 3*(3*b[1] + 3*b[2] + 3*b[3] + b[4]) * u[2]
                    + (9*b[1] + 9*b[2] + 9*b[3] + 9*b[4] + b[5]) * u[3]
                    - 3*(b[1] + b[2] + b[3] + 3*b[4] + b[5]) * u[4]
                    + 3*(b[4] + b[5]) * u[5]
                    - b[5] * u[6]
                )
        dest[4] = α * inv_left_weights[4] * (
                    - (b[1] + b[2] + b[3]) * u[1]
                    + 3*(b[1] + b[2] + b[3] + b[4]) * u[2]
                    - 3*(b[1] + b[2] + b[3] + 3*b[4] + b[5]) * u[3]
                    + (b[1] + b[2] + b[3] + 9*b[4] + 9*b[5] + b[6]) * u[4]
                    - 3*(b[4] + 3*b[5] + b[6]) * u[5]
                    + 3*(b[5] + b[6]) * u[6]
                    - b[6] * u[7]
                )
        for i in 5:length(inv_left_weights)
            dest[i] = α * inv_left_weights[i] * (
                        - b[i-1] * u[i-3]
                        + 3*(b[i-1] + b[i]) * u[i-2]
                        - 3*(b[i-1] + 3*b[i] + b[i+1]) * u[i-1]
                        + (b[i-1] + 9*b[i] + 9*b[i+1] + b[i+2]) * u[i]
                        - 3*(b[i] + 3*b[i+1] + b[i+2]) * u[i+1]
                        + 3*(b[i+1] + b[i+2]) * u[i+2]
                        - b[i+2] * u[i+3]
                    )
        end

        dest[end] = α * inv_right_weights[1] * (
                        (b[end] + b[end-1]) * u[end]
                        - 3*(b[end] + b[end-1]) * u[end-1]
                        + 3*(b[end] + b[end-1]) * u[end-2]
                        - (b[end] + b[end-1]) * u[end-3]
                    )
        dest[end-1] = α * inv_right_weights[2] * (
                        - 3*(b[end] + b[end-1]) * u[end]
                        + (b[end-2] + 9*b[end-1] + 9*b[end]) * u[end-1]
                        - 3*(b[end-2] + 3*b[end-1] + 3*b[end]) * u[end-2]
                        + 3*(b[end-2] + b[end-1] + b[end]) * u[end-3]
                        - b[end-2] * u[end-4]
                    )
        dest[end-2] = α * inv_right_weights[3] * (
                        3*(b[end] + b[end-1]) * u[end]
                        - 3*(b[end-2] + 3*b[end-1] + 3*b[end]) * u[end-1]
                        + (b[end-3] + 9*b[end-2] + 9*b[end-1] + 9*b[end]) * u[end-2]
                        - 3*(b[end-3] + 3*b[end-2] + b[end-1] + b[end]) * u[end-3]
                        + 3*(b[end-3] + b[end-2]) * u[end-4]
                        - b[end-3] * u[end-5]
                    )
        dest[end-3] = α * inv_right_weights[4] * (
                        - (b[end] + b[end-1]) * u[end]
                        + 3*(b[end-2] + b[end-1] + b[end]) * u[end-1]
                        - 3*(b[end-3] + 3*b[end-2] + b[end-1] + b[end]) * u[end-2]
                        + (b[end-4] + 9*b[end-3] + 9*b[end-2] + b[end-1] + b[end]) * u[end-3]
                        - 3*(b[end-4] + 3*b[end-3] + b[end-2]) * u[end-4]
                        + 3*(b[end-4] + b[end-3]) * u[end-5]
                        - b[end-4] * u[end-6]
                    )
        for i in 4:length(inv_right_weights)-1
            dest[end-i] = α * inv_right_weights[i+1] * (
                        - b[end-(i+1)] * u[end-(i+3)]
                        + 3*(b[end-(i+1)] + b[end-i]) * u[end-(i+2)]
                        - 3*(b[end-(i+1)] + 3*b[end-i] + b[end-(i-1)]) * u[end-(i+1)]
                        + (b[end-(i+1)] + 9*b[end-i] + 9*b[end-(i-1)] + b[end-(i-2)]) * u[end-i]
                        - 3*(b[end-i] + 3*b[end-(i-1)] + b[end-(i-2)]) * u[end-(i-1)]
                        + 3*(b[end-(i-1)] + b[end-(i-2)]) * u[end-(i-2)]
                        - b[end-(i-2)] * u[end-(i-3)]
                    )
        end
    end
end

function convolve_boundary_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache6, u::AbstractVector, b::AbstractVector, α, β)
    @unpack inv_left_weights, inv_right_weights = cache

    @inbounds begin
        dest[1] = α * inv_left_weights[1] * (
                    (b[1] + b[2] + b[3]) * u[1]
                    - 3*(b[1] + b[2] + b[3]) * u[2]
                    + 3*(b[1] + b[2] + b[3]) * u[3]
                    - (b[1] + b[2] + b[3]) * u[4]
                ) + β*dest[1]
        dest[2] = α * inv_left_weights[2] * (
                    - 3*(b[1] + b[2] + b[3]) * u[1]
                    + (9*b[1] + 9*b[2] + 9*b[3] + b[4]) * u[2]
                    - 3*(3*b[1] + 3*b[2] + 3*b[3] + b[4]) * u[3]
                    + 3*(b[1] + b[2] + b[3] + b[4]) * u[4]
                    - b[4] * u[5]
                ) + β*dest[2]
        dest[3] = α * inv_left_weights[3] * (
                    3*(b[1] + b[2] + b[3]) * u[1]
                    - 3*(3*b[1] + 3*b[2] + 3*b[3] + b[4]) * u[2]
                    + (9*b[1] + 9*b[2] + 9*b[3] + 9*b[4] + b[5]) * u[3]
                    - 3*(b[1] + b[2] + b[3] + 3*b[4] + b[5]) * u[4]
                    + 3*(b[4] + b[5]) * u[5]
                    - b[5] * u[6]
                ) + β*dest[3]
        dest[4] = α * inv_left_weights[4] * (
                    - (b[1] + b[2] + b[3]) * u[1]
                    + 3*(b[1] + b[2] + b[3] + b[4]) * u[2]
                    - 3*(b[1] + b[2] + b[3] + 3*b[4] + b[5]) * u[3]
                    + (b[1] + b[2] + b[3] + 9*b[4] + 9*b[5] + b[6]) * u[4]
                    - 3*(b[4] + 3*b[5] + b[6]) * u[5]
                    + 3*(b[5] + b[6]) * u[6]
                    - b[6] * u[7]
                ) + β*dest[4]
        for i in 5:length(inv_left_weights)
            dest[i] = α * inv_left_weights[i] * (
                        - b[i-1] * u[i-3]
                        + 3*(b[i-1] + b[i]) * u[i-2]
                        - 3*(b[i-1] + 3*b[i] + b[i+1]) * u[i-1]
                        + (b[i-1] + 9*b[i] + 9*b[i+1] + b[i+2]) * u[i]
                        - 3*(b[i] + 3*b[i+1] + b[i+2]) * u[i+1]
                        + 3*(b[i+1] + b[i+2]) * u[i+2]
                        - b[i+2] * u[i+3]
                    ) + β*dest[i]
        end

        dest[end] = α * inv_right_weights[1] * (
                        (b[end] + b[end-1]) * u[end]
                        - 3*(b[end] + b[end-1]) * u[end-1]
                        + 3*(b[end] + b[end-1]) * u[end-2]
                        - (b[end] + b[end-1]) * u[end-3]
                    ) + β*dest[end]
        dest[end-1] = α * inv_right_weights[2] * (
                        - 3*(b[end] + b[end-1]) * u[end]
                        + (b[end-2] + 9*b[end-1] + 9*b[end]) * u[end-1]
                        - 3*(b[end-2] + 3*b[end-1] + 3*b[end]) * u[end-2]
                        + 3*(b[end-2] + b[end-1] + b[end]) * u[end-3]
                        - b[end-2] * u[end-4]
                    ) + β*dest[end-1]
        dest[end-2] = α * inv_right_weights[3] * (
                        3*(b[end] + b[end-1]) * u[end]
                        - 3*(b[end-2] + 3*b[end-1] + 3*b[end]) * u[end-1]
                        + (b[end-3] + 9*b[end-2] + 9*b[end-1] + 9*b[end]) * u[end-2]
                        - 3*(b[end-3] + 3*b[end-2] + b[end-1] + b[end]) * u[end-3]
                        + 3*(b[end-3] + b[end-2]) * u[end-4]
                        - b[end-3] * u[end-5]
                    ) + β*dest[end-2]
        dest[end-3] = α * inv_right_weights[4] * (
                        - (b[end] + b[end-1]) * u[end]
                        + 3*(b[end-2] + b[end-1] + b[end]) * u[end-1]
                        - 3*(b[end-3] + 3*b[end-2] + b[end-1] + b[end]) * u[end-2]
                        + (b[end-4] + 9*b[end-3] + 9*b[end-2] + b[end-1] + b[end]) * u[end-3]
                        - 3*(b[end-4] + 3*b[end-3] + b[end-2]) * u[end-4]
                        + 3*(b[end-4] + b[end-3]) * u[end-5]
                        - b[end-4] * u[end-6]
                    ) + β*dest[end-3]
        for i in 4:length(inv_right_weights)-1
            dest[end-i] = α * inv_right_weights[i+1] * (
                        - b[end-(i+1)] * u[end-(i+3)]
                        + 3*(b[end-(i+1)] + b[end-i]) * u[end-(i+2)]
                        - 3*(b[end-(i+1)] + 3*b[end-i] + b[end-(i-1)]) * u[end-(i+1)]
                        + (b[end-(i+1)] + 9*b[end-i] + 9*b[end-(i-1)] + b[end-(i-2)]) * u[end-i]
                        - 3*(b[end-i] + 3*b[end-(i-1)] + b[end-(i-2)]) * u[end-(i-1)]
                        + 3*(b[end-(i-1)] + b[end-(i-2)]) * u[end-(i-2)]
                        - b[end-(i-2)] * u[end-(i-3)]
                    ) + β*dest[end-i]
        end
    end
end

@inline function convolve_interior_coefficients_loopbody(i, cache::MattssonSvärdNordström2004Cache6, u, b)
    @inbounds begin
        b_im1 = b[i-1]
        b_i   = b[i]
        b_ip1 = b[i+1]
        b_ip2 = b[i+2]

        retval = (
                    -b_im1 * u[i-3]
                    + 3*(b_im1 + b_i) * u[i-2]
                    - 3*(b_im1 + 3*b_i + b_ip1) * u[i-1]
                    + (b_im1 + 9*b_i + 9*b_ip1 + b_ip2) * u[i]
                    - 3*(b_i + 3*b_ip1 + b_ip2) * u[i+1]
                    + 3*(b_ip1 + b_ip2) * u[i+2]
                    - b_ip2 * u[i+3]
                )
    end

    retval
end



struct MattssonSvärdNordström2004Cache8{T,LeftWidth,RightWidth} <: AbstractCoefficientCache{T}
    inv_left_weights::SVector{LeftWidth,T}
    inv_right_weights::SVector{RightWidth,T}

    function MattssonSvärdNordström2004Cache8(inv_left_weights::SVector{LeftWidth,T}, inv_right_weights::SVector{RightWidth,T}) where {T,LeftWidth,RightWidth}
        boundary_length = 5

        if length(inv_left_weights) < boundary_length
            inv_left_weights = SVector(inv_left_weights..., ntuple(j->one(T), boundary_length-length(inv_left_weights))...)
        end
        if length(inv_right_weights) < boundary_length
            inv_right_weights = SVector(inv_right_weights..., ntuple(j->one(T), boundary_length-length(inv_right_weights))...)
        end

        new{T,length(inv_left_weights),length(inv_right_weights)}(inv_left_weights, inv_right_weights)
    end
end

lower_bandwidth(cache::MattssonSvärdNordström2004Cache8) = 4
upper_bandwidth(cache::MattssonSvärdNordström2004Cache8) = 4

function convolve_boundary_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache8, u::AbstractVector, b::AbstractVector, α)
    @unpack inv_left_weights, inv_right_weights = cache

    @inbounds begin
        dest[1] = α * inv_left_weights[1] * (
                    (b[1] + b[2] + b[3]) * u[1]
                    - 4*(b[1] + b[2] + b[3]) * u[2]
                    + 6*(b[1] + b[2] + b[3]) * u[3]
                    - 4*(b[1] + b[2] + b[3]) * u[4]
                    + (b[1] + b[2] + b[3]) * u[5]
                )
        dest[2] = α * inv_left_weights[2] * (
                    - 4*(b[1] + b[2] + b[3]) * u[1]
                    + (16*b[1] + 16*b[2] + 16*b[3] + b[4]) * u[2]
                    - 4*(6*b[1] + 6*b[2] + 6*b[3] + b[4]) * u[3]
                    + (16*b[1] + 16*b[2] + 16*b[3] + 6*b[4]) * u[4]
                    - 4*(b[1] + b[2] + b[3] + b[4]) * u[5]
                    + b[4] * u[6]
                )
        dest[3] = α * inv_left_weights[3] * (
                    6*(b[1] + b[2] + b[3]) * u[1]
                    - 4*(6*b[1] + 6*b[2] + 6*b[3] + b[4]) * u[2]
                    + (36*b[1] + 36*b[2] + 36*b[3] + 16*b[4] + b[5]) * u[3]
                    - 4*(6*b[1] + 6*b[2] + 6*b[3] + 6*b[4] + b[5]) * u[4]
                    + (6*b[1] + 6*b[2] + 6*b[3] + 16*b[4] + 6*b[5]) * u[5]
                    - 4*(b[4] + b[5]) * u[6]
                    + b[5] * u[7]
                )
        dest[4] = α * inv_left_weights[4] * (
                    - 4*(b[1] + b[2] + b[3]) * u[1]
                    + (16*b[1] + 16*b[2] + 16*b[3] + 6*b[4]) * u[2]
                    - 4*(6*b[1] + 6*b[2] + 6*b[3] + 6*b[4] + b[5]) * u[3]
                    + (16*b[1] + 16*b[2] + 16*b[3] + 36*b[4] + 16*b[5] + b[6]) * u[4]
                    - 4*(b[1] + b[2] + b[3] + 6*b[4] + 6*b[5] + b[6]) * u[5]
                    + (6*b[4] + 16*b[5] + 6*b[6]) * u[6]
                    - 4*(b[5] + b[6]) * u[7]
                    + b[6] * u[8]
                )
        dest[5] = α * inv_left_weights[5] * (
                    (b[1] + b[2] + b[3]) * u[1]
                    - 4*(b[1] + b[2] + b[3] + b[4]) * u[2]
                    + (6*b[1] + 6*b[2] + 6*b[3] + 16*b[4] + 6*b[5]) * u[3]
                    - 4*(b[1] + b[2] + b[3] + 6*b[4] + 6*b[5] + b[6]) * u[4]
                    + (b[1] + b[2] + b[3] + 16*b[4] + 36*b[5] + 16*b[6] + b[7]) * u[5]
                    - 4*(b[4] + 6*b[5] + 6*b[6] + b[7]) * u[6]
                    + (6*b[5] + 16*b[6] + 6*b[7]) * u[7]
                    - 4*(b[6] + b[7]) * u[8]
                    + b[7] * u[9]
                )
        for i in 6:length(inv_left_weights)
            dest[i] = α * inv_left_weights[i] * (
                        b[i-2] * u[i-4]
                        - 4*(b[i-2] + b[i-1]) * u[i-3]
                        + (6*b[i-2] + 16*b[i-1] + 6*b[i]) * u[i-2]
                        - 4*(b[i-2] + 6*b[i-1] + 6*b[i] + b[i+1]) * u[i-1]
                        + (b[i-2] + 16*b[i-1] + 36*b[i] + 16*b[i+1] + b[i+2]) * u[i]
                        - 4*(b[i-1] + 6*b[i] + 6*b[i+1] + b[i+2]) * u[i+1]
                        + (6*b[i] + 16*b[i+1] + 6*b[i+2]) * u[i+2]
                        - 4*(b[i+1] + b[i+2]) * u[i+3]
                        + b[i+2] * u[i+4]
                    )
        end

        dest[end] = α * inv_right_weights[1] * (
                        (b[end-2] + b[end-1] + b[end]) * u[end]
                        - 4*(b[end-2] + b[end-1] + b[end]) * u[end-1]
                        + 6*(b[end-2] + b[end-1] + b[end]) * u[end-2]
                        - 4*(b[end-2] + b[end-1] + b[end]) * u[end-3]
                        + (b[end-2] + b[end-1] + b[end]) * u[end-4]
                    )
        dest[end-1] = α * inv_right_weights[2] * (
                        - 4*(b[end-2] + b[end-1] + b[end]) * u[end]
                        + (b[end-3] + 16*b[end-2] + 16*b[end-1] + 16*b[end]) * u[end-1]
                        - 4*(b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-2]
                        + (6*b[end-3] + 16*b[end-2] + 16*b[end-1] + 16*b[end]) * u[end-3]
                        - 4*(b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-4]
                        + b[end-3] * u[end-5]
                    )
        dest[end-2] = α * inv_right_weights[3] * (
                        6*(b[end-2] + b[end-1] + b[end]) * u[end]
                        - 4*(b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-1]
                        + (b[end-4] + 16*b[end-3] + 36*b[end-2] + 36*b[end-1] + 36*b[end]) * u[end-2]
                        - 4*(b[end-4] + 6*b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-3]
                        + (6*b[end-4] + 16*b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-4]
                        - 4*(b[end-4] + b[end-3]) * u[end-5]
                        + b[end-4] * u[end-6]
                    )
        dest[end-3] = α * inv_right_weights[4] * (
                        - 4*(b[end-2] + b[end-1] + b[end]) * u[end]
                        + (6*b[end-3] + 16*b[end-2] + 16*b[end-1] + 16*b[end]) * u[end-1]
                        - 4*(b[end-4] + 6*b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-2]
                        + (b[end-5] + 16*b[end-4] + 36*b[end-3] + 16*b[end-2] + 16*b[end-1] + 16*b[end]) * u[end-3]
                        - 4*(b[end-5] + 6*b[end-4] + 6*b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-4]
                        + (6*b[end-5] + 16*b[end-4] + 6*b[end-3]) * u[end-5]
                        - 4*(b[end-5] + b[end-4]) * u[end-6]
                        + b[end-5] * u[end-7]
                    )
        dest[end-4] = α * inv_right_weights[5] * (
                        (b[end-2] + b[end-1] + b[end]) * u[end]
                        - 4*(b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-1]
                        + (6*b[end-4] + 16*b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-2]
                        - 4*(b[end-5] + 6*b[end-4] + 6*b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-3]
                        + (b[end-6] + 16*b[end-5] + 36*b[end-4] + 16*b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-4]
                        - 4*(b[end-6] + 6*b[end-5] + 6*b[end-4] + b[end-3]) * u[end-5]
                        + (6*b[end-6] + 16*b[end-5] + 6*b[end-4]) * u[end-6]
                        - 4*(b[end-6] + b[end-5]) * u[end-7]
                        + b[end-6] * u[end-8]
                    )
        for i in 5:length(inv_right_weights)-1
            dest[end-i] = α * inv_right_weights[i+1] * (
                        b[end-(i+2)] * u[end-(i+4)]
                        - 4*(b[end-(i+2)] + b[end-(i+1)]) * u[end-(i+3)]
                        + (6*b[end-(i+2)] + 16*b[end-(i+1)] + 6*b[end-i]) * u[end-(i+2)]
                        - 4*(b[end-(i+2)] + 6*b[end-(i+1)] + 6*b[end-i] + b[end-(i-1)]) * u[end-(i+1)]
                        + (b[end-(i+2)] + 16*b[end-(i+1)] + 36*b[end-i] + 16*b[end-(i-1)] + b[end-(i-2)]) * u[end-i]
                        - 4*(b[end-(i+1)] + 6*b[end-i] + 6*b[end-(i-1)] + b[end-(i-2)]) * u[end-(i-1)]
                        + (6*b[end-i] + 16*b[end-(i-1)] + 6*b[end-(i-2)]) * u[end-(i-2)]
                        - 4*(b[end-(i-1)] + b[end-(i-2)]) * u[end-(i-3)]
                        + b[end-(i-2)] * u[end-(i-4)]
                    )
        end
    end
end

function convolve_boundary_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache8, u::AbstractVector, b::AbstractVector, α, β)
    @unpack inv_left_weights, inv_right_weights = cache

    @inbounds begin
        dest[1] = α * inv_left_weights[1] * (
                    (b[1] + b[2] + b[3]) * u[1]
                    - 4*(b[1] + b[2] + b[3]) * u[2]
                    + 6*(b[1] + b[2] + b[3]) * u[3]
                    - 4*(b[1] + b[2] + b[3]) * u[4]
                    + (b[1] + b[2] + b[3]) * u[5]
                ) + β*dest[1]
        dest[2] = α * inv_left_weights[2] * (
                    - 4*(b[1] + b[2] + b[3]) * u[1]
                    + (16*b[1] + 16*b[2] + 16*b[3] + b[4]) * u[2]
                    - 4*(6*b[1] + 6*b[2] + 6*b[3] + b[4]) * u[3]
                    + (16*b[1] + 16*b[2] + 16*b[3] + 6*b[4]) * u[4]
                    - 4*(b[1] + b[2] + b[3] + b[4]) * u[5]
                    + b[4] * u[6]
                ) + β*dest[2]
        dest[3] = α * inv_left_weights[3] * (
                    6*(b[1] + b[2] + b[3]) * u[1]
                    - 4*(6*b[1] + 6*b[2] + 6*b[3] + b[4]) * u[2]
                    + (36*b[1] + 36*b[2] + 36*b[3] + 16*b[4] + b[5]) * u[3]
                    - 4*(6*b[1] + 6*b[2] + 6*b[3] + 6*b[4] + b[5]) * u[4]
                    + (6*b[1] + 6*b[2] + 6*b[3] + 16*b[4] + 6*b[5]) * u[5]
                    - 4*(b[4] + b[5]) * u[6]
                    + b[5] * u[7]
                ) + β*dest[3]
        dest[4] = α * inv_left_weights[4] * (
                    - 4*(b[1] + b[2] + b[3]) * u[1]
                    + (16*b[1] + 16*b[2] + 16*b[3] + 6*b[4]) * u[2]
                    - 4*(6*b[1] + 6*b[2] + 6*b[3] + 6*b[4] + b[5]) * u[3]
                    + (16*b[1] + 16*b[2] + 16*b[3] + 36*b[4] + 16*b[5] + b[6]) * u[4]
                    - 4*(b[1] + b[2] + b[3] + 6*b[4] + 6*b[5] + b[6]) * u[5]
                    + (6*b[4] + 16*b[5] + 6*b[6]) * u[6]
                    - 4*(b[5] + b[6]) * u[7]
                    + b[6] * u[8]
                ) + β*dest[4]
        dest[5] = α * inv_left_weights[5] * (
                    (b[1] + b[2] + b[3]) * u[1]
                    - 4*(b[1] + b[2] + b[3] + b[4]) * u[2]
                    + (6*b[1] + 6*b[2] + 6*b[3] + 16*b[4] + 6*b[5]) * u[3]
                    - 4*(b[1] + b[2] + b[3] + 6*b[4] + 6*b[5] + b[6]) * u[4]
                    + (b[1] + b[2] + b[3] + 16*b[4] + 36*b[5] + 16*b[6] + b[7]) * u[5]
                    - 4*(b[4] + 6*b[5] + 6*b[6] + b[7]) * u[6]
                    + (6*b[5] + 16*b[6] + 6*b[7]) * u[7]
                    - 4*(b[6] + b[7]) * u[8]
                    + b[7] * u[9]
                ) + β*dest[5]
        for i in 6:length(inv_left_weights)
            dest[i] = α * inv_left_weights[i] * (
                        b[i-2] * u[i-4]
                        - 4*(b[i-2] + b[i-1]) * u[i-3]
                        + (6*b[i-2] + 16*b[i-1] + 6*b[i]) * u[i-2]
                        - 4*(b[i-2] + 6*b[i-1] + 6*b[i] + b[i+1]) * u[i-1]
                        + (b[i-2] + 16*b[i-1] + 36*b[i] + 16*b[i+1] + b[i+2]) * u[i]
                        - 4*(b[i-1] + 6*b[i] + 6*b[i+1] + b[i+2]) * u[i+1]
                        + (6*b[i] + 16*b[i+1] + 6*b[i+2]) * u[i+2]
                        - 4*(b[i+1] + b[i+2]) * u[i+3]
                        + b[i+2] * u[i+4]
                    ) + β*dest[i]
        end

        dest[end] = α * inv_right_weights[1] * (
                        (b[end-2] + b[end-1] + b[end]) * u[end]
                        - 4*(b[end-2] + b[end-1] + b[end]) * u[end-1]
                        + 6*(b[end-2] + b[end-1] + b[end]) * u[end-2]
                        - 4*(b[end-2] + b[end-1] + b[end]) * u[end-3]
                        + (b[end-2] + b[end-1] + b[end]) * u[end-4]
                    ) + β*dest[end]
        dest[end-1] = α * inv_right_weights[2] * (
                        - 4*(b[end-2] + b[end-1] + b[end]) * u[end]
                        + (b[end-3] + 16*b[end-2] + 16*b[end-1] + 16*b[end]) * u[end-1]
                        - 4*(b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-2]
                        + (6*b[end-3] + 16*b[end-2] + 16*b[end-1] + 16*b[end]) * u[end-3]
                        - 4*(b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-4]
                        + b[end-3] * u[end-5]
                    ) + β*dest[end-1]
        dest[end-2] = α * inv_right_weights[3] * (
                        6*(b[end-2] + b[end-1] + b[end]) * u[end]
                        - 4*(b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-1]
                        + (b[end-4] + 16*b[end-3] + 36*b[end-2] + 36*b[end-1] + 36*b[end]) * u[end-2]
                        - 4*(b[end-4] + 6*b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-3]
                        + (6*b[end-4] + 16*b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-4]
                        - 4*(b[end-4] + b[end-3]) * u[end-5]
                        + b[end-4] * u[end-6]
                    ) + β*dest[end-2]
        dest[end-3] = α * inv_right_weights[4] * (
                        - 4*(b[end-2] + b[end-1] + b[end]) * u[end]
                        + (6*b[end-3] + 16*b[end-2] + 16*b[end-1] + 16*b[end]) * u[end-1]
                        - 4*(b[end-4] + 6*b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-2]
                        + (b[end-5] + 16*b[end-4] + 36*b[end-3] + 16*b[end-2] + 16*b[end-1] + 16*b[end]) * u[end-3]
                        - 4*(b[end-5] + 6*b[end-4] + 6*b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-4]
                        + (6*b[end-5] + 16*b[end-4] + 6*b[end-3]) * u[end-5]
                        - 4*(b[end-5] + b[end-4]) * u[end-6]
                        + b[end-5] * u[end-7]
                    ) + β*dest[end-3]
        dest[end-4] = α * inv_right_weights[5] * (
                        (b[end-2] + b[end-1] + b[end]) * u[end]
                        - 4*(b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-1]
                        + (6*b[end-4] + 16*b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-2]
                        - 4*(b[end-5] + 6*b[end-4] + 6*b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-3]
                        + (b[end-6] + 16*b[end-5] + 36*b[end-4] + 16*b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-4]
                        - 4*(b[end-6] + 6*b[end-5] + 6*b[end-4] + b[end-3]) * u[end-5]
                        + (6*b[end-6] + 16*b[end-5] + 6*b[end-4]) * u[end-6]
                        - 4*(b[end-6] + b[end-5]) * u[end-7]
                        + b[end-6] * u[end-8]
                    ) + β*dest[end-4]
        for i in 5:length(inv_right_weights)-1
            dest[end-i] = α * inv_right_weights[i+1] * (
                        b[end-(i+2)] * u[end-(i+4)]
                        - 4*(b[end-(i+2)] + b[end-(i+1)]) * u[end-(i+3)]
                        + (6*b[end-(i+2)] + 16*b[end-(i+1)] + 6*b[end-i]) * u[end-(i+2)]
                        - 4*(b[end-(i+2)] + 6*b[end-(i+1)] + 6*b[end-i] + b[end-(i-1)]) * u[end-(i+1)]
                        + (b[end-(i+2)] + 16*b[end-(i+1)] + 36*b[end-i] + 16*b[end-(i-1)] + b[end-(i-2)]) * u[end-i]
                        - 4*(b[end-(i+1)] + 6*b[end-i] + 6*b[end-(i-1)] + b[end-(i-2)]) * u[end-(i-1)]
                        + (6*b[end-i] + 16*b[end-(i-1)] + 6*b[end-(i-2)]) * u[end-(i-2)]
                        - 4*(b[end-(i-1)] + b[end-(i-2)]) * u[end-(i-3)]
                        + b[end-(i-2)] * u[end-(i-4)]
                    ) + β*dest[end-i]
        end
    end
end

@inline function convolve_interior_coefficients_loopbody(i, cache::MattssonSvärdNordström2004Cache8, u, b)
    @inbounds begin
        b_im2 = b[i-2]
        b_im1 = b[i-1]
        b_i   = b[i]
        b_ip1 = b[i+1]
        b_ip2 = b[i+2]

        retval = (
                    b_im2 * u[i-4]
                    - 4*(b_im2 + b_im1) * u[i-3]
                    + (6*b_im2 + 16*b_im1 + 6*b_i) * u[i-2]
                    - 4*(b_im2 + 6*b_im1 + 6*b_i + b_ip1) * u[i-1]
                    + (b_im2 + 16*b_im1 + 36*b_i + 16*b_ip1 + b_ip2) * u[i]
                    - 4*(b_im1 + 6*b_i + 6*b_ip1 + b_ip2) * u[i+1]
                    + (6*b_i + 16*b_ip1 + 6*b_ip2) * u[i+2]
                    - 4*(b_ip1 + b_ip2) * u[i+3]
                    + b_ip2 * u[i+4]
                )
    end

    retval
end


function first_derivative_coefficients(source::MattssonSvärdNordström2004, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            DerivativeCoefficientRow{T,1,2}(SVector(-one(T), one(T))),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(1//2))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector(T(1//2))
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 1, order, source)
    #=elseif order == 4
        left_boundary = (
            # q1
            DerivativeCoefficientRow{T,1,4}(SVector(T(-24//17),
                                                    T(59//34),
                                                    T(-4//17),
                                                    T(-3//34) )),
            # q2
            DerivativeCoefficientRow{T,1,3}(SVector(T(-1//2),
                                                    T(0),
                                                    T(1//2))),
            # q3
            DerivativeCoefficientRow{T,1,5}(SVector(T(4//43),
                                                    T(-59//86),
                                                    T(0),
                                                    T(59//86),
                                                    T(-4//43))),
            # q4
            DerivativeCoefficientRow{T,1,6}(SVector(T(3//98),
                                                    T(0),
                                                    T(-59//98),
                                                    T(0),
                                                    T(32//49),
                                                    T(-4//49))),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(2//3), T(-1//12))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(17//48),
                                T(59//48),
                                T(43//48),
                                T(49//48) )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 1, order, source)
    elseif order == 6
        left_boundary = (
            # q1
            DerivativeCoefficientRow{T,1,6}(SVector(T(-21600//13649),
                                                    T(104009//54596),
                                                    T(30443//81894),
                                                    T(-33311//27298),
                                                    T(16863//27298),
                                                    T(-15025//163788) )),
            # q2
            DerivativeCoefficientRow{T,1,6}(SVector(T(-104009//240260),
                                                    T(0),
                                                    T(-311//72078),
                                                    T(20229//24026),
                                                    T(-24337//48052),
                                                    T(36661//360390) )),
            # q3
            DerivativeCoefficientRow{T,1,6}(SVector(T(-30443//162660),
                                                    T(311//32532),
                                                    T(0),
                                                    T(-11155//16266),
                                                    T(41287//32532),
                                                    T(-21999//54220) )),
            # q4
            DerivativeCoefficientRow{T,1,7}(SVector(T(33311//107180),
                                                    T(-20229//21436),
                                                    T(485//1398),
                                                    T(0),
                                                    T(4147//21436),
                                                    T(25427//321540),
                                                    T(72//5359) )),
            # q5
            DerivativeCoefficientRow{T,1,8}(SVector(T(-16863//78770),
                                                    T(24337//31508),
                                                    T(-41287//47262),
                                                    T(-4147//15754),
                                                    T(0),
                                                    T(342523//472620),
                                                    T(-1296//7877),
                                                    T(144//7877) )),
            # q6
            DerivativeCoefficientRow{T,1,9}(SVector(T(15025//525612),
                                                    T(-36661//262806),
                                                    T(21999//87602),
                                                    T(-25427//262806),
                                                    T(-342523//525612),
                                                    T(0),
                                                    T(32400//43801),
                                                    T(-6480//43801),
                                                    T(720//43801) ))
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(3//4), T(-3//20), T(1//60))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(13649//43200),
                                T(12013//8640),
                                T(2711//4320),
                                T(5359//4320),
                                T(7877//8640),
                                T(43801//43200) )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 1, order, source)
    elseif order == 8
        x1 =
        left_boundary = (
            # q1
            DerivativeCoefficientRow{T,1,8}(SVector(T(-2540160//1498139),
                                                    T(515174618457408334671//81327545493318772580),
                                                    T(-319653445960068136431//16265509098663754516),
                                                    T(1756838164798071830975//48796527295991263548),
                                                    T(-132855587210457907869//4066377274665938629),
                                                    T(205031990442925032261//16265509098663754516),
                                                    T(-147365687//1707878460),
                                                    T(-13621478277275954493//16265509098663754516) )),
            # q2
            DerivativeCoefficientRow{T,1,8}(SVector(T(-515174618457408334671//420776658856749814780),
                                                    T(0),
                                                    T(335158710375082986831//30055475632624986770),
                                                    T(-81614753349680015895//3005547563262498677),
                                                    T(996306452293020120805//36066570759149984124),
                                                    T(-338011221795904141851//30055475632624986770),
                                                    T(1109310//21038833),
                                                    T(503120655918727839389//631164988285124722170) )),
            # q3
            DerivativeCoefficientRow{T,1,8}(SVector(T(106551148653356045477//4733507987139058764),
                                                    T(-111719570125027662277//1690538566835378130),
                                                    T(0),
                                                    T(66611314247587471205//338107713367075626),
                                                    T(-177186717090009550283//676215426734151252),
                                                    T(13307564278128858831//112702571122358542),
                                                    T(-354462//1972295),
                                                    T(-7342453837193297117//788917997856509794) )),
            # q4
            DerivativeCoefficientRow{T,1,8}(SVector(T(-1756838164798071830975//297539771191584439476),
                                                    T(81614753349680015895//3542140133233148089),
                                                    T(-199833942742762413615//7084280266466296178),
                                                    T(0),
                                                    T(334165099234925485449//14168560532932592356),
                                                    T(-97726720448760690009//7084280266466296178),
                                                    T(-2934266//74384943),
                                                    T(65502786250599341781//49589961865264073246) )),
            # q5
            DerivativeCoefficientRow{T,1,9}(SVector(T(132855587210457907869//5691012984537517679),
                                                    T(-996306452293020120805//9756022259207173164),
                                                    T(531560151270028650849//3252007419735724388),
                                                    T(-334165099234925485449//3252007419735724388),
                                                    T(0),
                                                    T(68552420813742056079//3252007419735724388),
                                                    T(-2343561//22764052),
                                                    T(-197397151320811747355//68292155814450212148),
                                                    T(-2592//299527) )),
            # q6
            DerivativeCoefficientRow{T,1,10}(SVector(T(-9763428116329763441//3358016562304840404),
                                                     T(112670407265301380617//8395041405762101010),
                                                     T(-13307564278128858831//559669427050806734),
                                                     T(32575573482920230003//1679008281152420202),
                                                     T(-22850806937914018693//3358016562304840404),
                                                     T(0),
                                                     T(5346432//9794215),
                                                     T(1388463154122240//14728142817126493),
                                                     T(3072//103097),
                                                     T(-288//103097) )),
             # q7
             DerivativeCoefficientRow{T,1,11}(SVector(T(21052241//763903740),
                                                      T(-1109310//12731729),
                                                      T(3190158//63658645),
                                                      T(2934266//38195187),
                                                      T(2343561//50926916),
                                                      T(-48117888//63658645),
                                                      T(0),
                                                      T(52309152//63658645),
                                                      T(-145152//670091),
                                                      T(27648//670091),
                                                      T(-2592//670091) )),
              # q8
              DerivativeCoefficientRow{T,1,12}(SVector(T(13621478277275954493//55672594705880416916),
                                                       T(-503120655918727839389//417544460294103126870),
                                                       T(66082084534739674053//27836297352940208458),
                                                       T(-65502786250599341781//27836297352940208458),
                                                       T(197397151320811747355//167017784117641250748),
                                                       T(-87473178709701120//732534140866847591),
                                                       T(-366164064//487135205),
                                                       T(0),
                                                       T(4064256//5127739),
                                                       T(-1016064//5127739),
                                                       T(193536//5127739),
                                                       T(-18144//5127739) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(4//5), T(-1//5), T(4//105), T(-1//280))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(1498139//5080320),
                                T(1107307//725760),
                                T(20761//80640),
                                T(1304999//725760),
                                T(299527//725760),
                                T(103097//80640),
                                T(670091//725760),
                                T(5127739//5080320) )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 1, order, source)
    =#
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end


function second_derivative_coefficients(source::MattssonSvärdNordström2004, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,3}(SVector(T(1),
                                                    T(-2),
                                                    T(1) )),
        )
        right_boundary = left_boundary
        upper_coef = SVector(T(1))
        central_coef = T(-2)
        lower_coef = upper_coef
        left_weights = SVector(T(1//2))
        right_weights = left_weights
        left_boundary_derivatives = (
            DerivativeCoefficientRow{T,1,3}(SVector(T(-3//2),
                                                    T(2),
                                                    T(-1//2) )),
        )
        right_boundary_derivatives = (-left_boundary_derivatives[1], )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 2, order, source)
    #=elseif order == 4
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,4}(SVector(T(2),
                                                    T(-5),
                                                    T(4),
                                                    T(-1) )),
            # d2
            DerivativeCoefficientRow{T,1,3}(SVector(T(1),
                                                    T(-2),
                                                    T(1) )),
            # d3
            DerivativeCoefficientRow{T,1,5}(SVector(T(-4//43),
                                                    T(59//43),
                                                    T(-110//43),
                                                    T(59//43),
                                                    T(-4//43) )),
            # d4
            DerivativeCoefficientRow{T,1,6}(SVector(T(-1//49),
                                                    T(0),
                                                    T(59//49),
                                                    T(-118//49),
                                                    T(64//49),
                                                    T(-4//49) ))
        )
        right_boundary = left_boundary
        upper_coef = SVector(T(4//3), T(-1//12))
        central_coef = T(-5//2)
        lower_coef = upper_coef
        left_weights = SVector(T(17//48), T(59//48), T(43//48), T(49//48))
        right_weights = left_weights
        left_boundary_derivatives = (
            DerivativeCoefficientRow{T,1,4}(SVector(T(-11/6),
                                                    T(3),
                                                    T(-3//2),
                                                    T(1//3) )),
        )
        right_boundary_derivatives = (-left_boundary_derivatives[1], )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 2, order, source)
    elseif order == 6
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,6}(SVector(T(114170//40947),
                                                    T(-438107//54596),
                                                    T(336409//40947),
                                                    T(-276997//81894),
                                                    T(3747//13649),
                                                    T(21035//163788) )),
            # d2
            DerivativeCoefficientRow{T,1,6}(SVector(T(6173//5860),
                                                    T(-2066//879),
                                                    T(3283//1758),
                                                    T(-303//293),
                                                    T(2111//3516),
                                                    T(-601//4395) )),
            # d3
            DerivativeCoefficientRow{T,1,6}(SVector(T(-52391//81330),
                                                    T(134603//32532),
                                                    T(-21982//2711),
                                                    T(112915//16266),
                                                    T(-46969//16266),
                                                    T(30409//54220) )),
            # d4
            DerivativeCoefficientRow{T,1,7}(SVector(T(68603//321540),
                                                    T(-12423//10718),
                                                    T(112915//32154),
                                                    T(-75934//16077),
                                                    T(53369//21436),
                                                    T(-54899//160770),
                                                    T(48//5359) )),
            # d5
            DerivativeCoefficientRow{T,1,8}(SVector(T(-7053//39385),
                                                    T(86551//94524),
                                                    T(-46969//23631),
                                                    T(53369//15754),
                                                    T(-87904//23631),
                                                    T(820271//472620),
                                                    T(-1296//7877),
                                                    T(96//7877) )),
            # d6
            DerivativeCoefficientRow{T,1,9}(SVector(T(21035//525612),
                                                    T(-24641//131403),
                                                    T(30409//87602),
                                                    T(-54899//131403),
                                                    T(820271//525612),
                                                    T(-117600//43801),
                                                    T(64800//43801),
                                                    T(-6480//43801),
                                                    T(480//43801) )),
        )
        right_boundary = left_boundary
        upper_coef = SVector(T(3//2), T(-3//20), T(1//90))
        central_coef = T(-49//18)
        lower_coef = upper_coef
        left_weights = SVector( T(13649//43200),
                                T(12013//8640),
                                T(2711//4320),
                                T(5359//4320),
                                T(7877//8640),
                                T(43801//43200) )
        right_weights = left_weights
        left_boundary_derivatives = (
            DerivativeCoefficientRow{T,1,5}(SVector(T(-25//12),
                                                    T(4),
                                                    T(-3),
                                                    T(4//3),
                                                    T(-1//4) )),
        )
        right_boundary_derivatives = (-left_boundary_derivatives[1], )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 2, order, source)
    elseif order == 8
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,7}(SVector(T(4870382994799//1358976868290),
                                                    T(-893640087518//75498714905),
                                                    T(926594825119//60398971924),
                                                    T(-1315109406200//135897686829),
                                                    T(39126983272//15099742981),
                                                    T(12344491342//75498714905),
                                                    T(-451560522577//2717953736580) )),
            # d2
            DerivativeCoefficientRow{T,1,8}(SVector(T(333806012194//390619153855),
                                                    T(-154646272029//111605472530),
                                                    T(1168338040//33481641759),
                                                    T(82699112501//133926567036),
                                                    T(-171562838//11160547253),
                                                    T(-28244698346//167408208795),
                                                    T(11904122576//167408208795),
                                                    T(-2598164715//312495323084) )),
            # d3
            DerivativeCoefficientRow{T,1,8}(SVector(T(7838984095//52731029988),
                                                    T(1168338040//5649753213),
                                                    T(-88747895//144865467),
                                                    T(423587231//627750357),
                                                    T(-43205598281//22599012852),
                                                    T(4876378562//1883251071),
                                                    T(-5124426509//3766502142),
                                                    T(10496900965//39548272491) )),
            # d4
            DerivativeCoefficientRow{T,1,8}(SVector(T(-94978241528//828644350023),
                                                    T(82699112501//157837019052),
                                                    T(1270761693//13153084921),
                                                    T(-167389605005//118377764289),
                                                    T(48242560214//39459254763),
                                                    T(-31673996013//52612339684),
                                                    T(43556319241//118377764289),
                                                    T(-44430275135//552429566682) )),
            # d5
            DerivativeCoefficientRow{T,1,9}(SVector(T(1455067816//21132528431),
                                                    T(-171562838//3018932633),
                                                    T(-43205598281//36227191596),
                                                    T(48242560214//9056797899),
                                                    T(-52276055645//6037865266),
                                                    T(57521587238//9056797899),
                                                    T(-80321706377//36227191596),
                                                    T(8078087158//21132528431),
                                                    T(-1296//299527) )),
            # d6
            DerivativeCoefficientRow{T,1,10}(SVector(T(10881504334//327321118845),
                                                     T(-28244698346//140280479505),
                                                     T(4876378562//9352031967),
                                                     T(-10557998671//12469375956),
                                                     T(57521587238//28056095901),
                                                     T(-278531401019//93520319670),
                                                     T(73790130002//46760159835),
                                                     T(-137529995233//785570685228),
                                                     T(2048//103097),
                                                     T(-144//103097) )),
             # d7
             DerivativeCoefficientRow{T,1,11}(SVector(T(-135555328849//8509847458140),
                                                      T(11904122576//101307707835),
                                                      T(-5124426509//13507694378),
                                                      T(43556319241//60784624701),
                                                      T(-80321706377//81046166268),
                                                      T(73790130002//33769235945),
                                                      T(-950494905688//303923123505),
                                                      T(239073018673//141830790969),
                                                      T(-145152//670091),
                                                      T(18432//670091),
                                                      T(-1296//670091) )),
              # d8
              DerivativeCoefficientRow{T,1,12}(SVector(T(0),
                                                       T(-2598164715//206729925524),
                                                       T(10496900965//155047444143),
                                                       T(-44430275135//310094888286),
                                                       T(425162482//2720130599),
                                                       T(-137529995233//620189776572),
                                                       T(239073018673//155047444143),
                                                       T(-144648000000//51682481381),
                                                       T(8128512//5127739),
                                                       T(-1016064//5127739),
                                                       T(129024//5127739),
                                                       T(-9072//5127739) )),

        )
        right_boundary = left_boundary
        upper_coef = SVector(T(8//5), T(-1//5), T(8//315), T(-1//560))
        central_coef = T(-205//72)
        lower_coef = upper_coef
        left_weights = SVector( T(1498139//5080320),
                                T(1107307//725760),
                                T(20761//80640),
                                T(1304999//725760),
                                T(299527//725760),
                                T(103097//80640),
                                T(670091//725760),
                                T(5127739//5080320) )
        right_weights = left_weights
        left_boundary_derivatives = (
            DerivativeCoefficientRow{T,1,7}(SVector(T(-4723//2100),
                                                    T(839//175),
                                                    T(-157//35),
                                                    T(278//105),
                                                    T(-103//140),
                                                    T(-1//175),
                                                    T(6//175) )),
        )
        right_boundary_derivatives = (-left_boundary_derivatives[1], )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 2, order, source)
    =#
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end