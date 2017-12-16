
"""
    PeriodicDerivativeOperator{T,StencilWidth,Parallel}

A derivative operator on a periodic grid.
"""
struct PeriodicDerivativeOperator{T,LowerOffset,UpperOffset,Parallel} <: AbstractDerivativeOperator{T}
    lower_coef::SVector{LowerOffset, T}
    central_coef::T
    upper_coef::SVector{UpperOffset, T}
    parallel::Parallel
end


Base.A_mul_B!(dest, D::PeriodicDerivativeOperator, u) = mul!(dest, D, u, one(eltype(dest)), zero(eltype(dest)))

"""
    mul!(dest::AbstractVector, D::PeriodicDerivativeOperator, u::AbstractVector, α, β)

Compute `α*D*u + β*dest` and store the result in `dest`.
"""
function mul!(dest::AbstractVector, D::PeriodicDerivativeOperator{T,LowerOffset,UpperOffset}, u::AbstractVector, α, β) where {T,LowerOffset,UpperOffset}
    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(u) > LowerOffset+UpperOffset DimensionMismatch
    end

    @unpack lower_coef, central_coef, upper_coef, parallel = D
    convolve_boundary_coefficients!(dest, lower_coef, central_coef, upper_coef, u, α, β)
    convolve_interior_coefficients!(dest, lower_coef, central_coef, upper_coef, u, α, β, parallel)
end


@generated function convolve_boundary_coefficients!(dest::AbstractVector, lower_coef::SVector{LowerOffset}, central_coef, upper_coef::SVector{UpperOffset}, u::AbstractVector, α, β) where {LowerOffset, UpperOffset}
    ex_lower = :( nothing )
    for i in 1:LowerOffset
        ex = :( lower_coef[$LowerOffset]*u[end-$(LowerOffset-i)] )
        for j in LowerOffset-1:-1:1
            if i-j < 1
                ex = :( $ex + lower_coef[$j]*u[end-$(j-i)] )
            else
                ex = :( $(ex) + lower_coef[$j]*u[$(i-j)] )
            end
        end
        ex = :( $ex + central_coef*u[$i] )
        for j in 1:UpperOffset
            ex = :( $ex + upper_coef[$j]*u[$(i+j)] )
        end
        ex_lower = quote
            $ex_lower
            @inbounds dest[$i] = β*dest[$i] + α*$ex
        end
    end

    ex_upper = :( nothing )
    for i in (UpperOffset-1):-1:0
        ex = :( lower_coef[$LowerOffset]*u[end-$(i+LowerOffset)] )
        for j in LowerOffset-1:-1:1
            ex = :( $ex + lower_coef[$j]*u[end-$(j+i)] )
        end
        ex = :( $ex + central_coef*u[end-$i] )
        for j in 1:UpperOffset
            if i-j < 0
                ex = :( $ex + upper_coef[$j]*u[$(j-i)] )
            else
                ex = :( $ex + upper_coef[$j]*u[end-$(i-j)] )
            end
        end
        ex_upper = quote
            $ex_upper
            @inbounds dest[end-$i] = β*dest[end-$i] + α*$ex
        end
    end

    quote
        $ex_lower
        $ex_upper
    end
end


@generated function convolve_interior_coefficients!(dest::AbstractVector, lower_coef::SVector{LowerOffset}, central_coef, upper_coef::SVector{UpperOffset}, u::AbstractVector, α, β, parallel) where {LowerOffset, UpperOffset}
    ex = :( lower_coef[$LowerOffset]*u[i-$LowerOffset] )
    for j in LowerOffset-1:-1:1
        ex = :( $ex + lower_coef[$j]*u[i-$j] )
    end
    ex = :( $ex + central_coef*u[i] )
    for j in 1:UpperOffset
        ex = :( $ex + upper_coef[$j]*u[i+$j] )
    end

    quote
        @inbounds for i in $(LowerOffset+1):(length(dest)-$UpperOffset)
            dest[i] = β*dest[i] + α*$ex
        end
    end
end

function convolve_interior_coefficients!(dest::AbstractVector{T}, lower_coef::SVector{LowerOffset}, central_coef, upper_coef::SVector{UpperOffset}, u::AbstractVector, α, β, parallel::Val{:threads}) where {T, LowerOffset, UpperOffset}
    Threads.@threads for i in (LowerOffset+1):(length(dest)-UpperOffset) @inbounds begin
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
