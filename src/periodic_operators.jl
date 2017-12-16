
"""
    PeriodicDerivativeCoefficients{T,StencilWidth,Parallel}

The coefficients of a derivative operator on a periodic grid.
"""
struct PeriodicDerivativeCoefficients{T,LowerOffset,UpperOffset,Parallel} <: AbstractDerivativeCoefficients{T}
    # coefficients defining the operator and its action
    lower_coef::SVector{LowerOffset, T}
    central_coef::T
    upper_coef::SVector{UpperOffset, T}
    parallel::Parallel
    # corresponding orders etc.
    derivative_order::Int
    accuracy_order  ::Int
end

derivative_order(coefficients::AbstractDerivativeCoefficients) = coefficients.derivative_order
accuracy_order(coefficients::AbstractDerivativeCoefficients) = coefficients.accuracy_order
Base.eltype(coefficients::AbstractDerivativeCoefficients{T}) where {T} = T


function Base.A_mul_B!(dest, coefficients::PeriodicDerivativeCoefficients, u)
    mul!(dest, coefficients, u, one(eltype(dest)), zero(eltype(dest)))
end

"""
    mul!(dest::AbstractVector, coefficients::PeriodicDerivativeCoefficients, u::AbstractVector, α, β)

Compute `α*D*u + β*dest` and store the result in `dest`.
"""
function mul!(dest::AbstractVector, coefficients::PeriodicDerivativeCoefficients{T,LowerOffset,UpperOffset}, u::AbstractVector, α, β) where {T,LowerOffset,UpperOffset}
    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(u) > LowerOffset+UpperOffset DimensionMismatch
    end

    @unpack lower_coef, central_coef, upper_coef, parallel = coefficients
    convolve_periodic_boundary_coefficients!(dest, lower_coef, central_coef, upper_coef, u, α, β)
    convolve_interior_coefficients!(dest, lower_coef, central_coef, upper_coef, u, α, β, parallel)
end


@generated function convolve_periodic_boundary_coefficients!(dest::AbstractVector, lower_coef::SVector{LowerOffset}, central_coef, upper_coef::SVector{UpperOffset}, u::AbstractVector, α, β) where {LowerOffset, UpperOffset}
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


"""
    periodic_central_derivative_coefficients(derivative_order, accuracy_order, T=Float64, parallel=Val{:serial}())

Create the `PeriodicDerivativeCoefficients` approximating the `derivative_order`-th
derivative with an order of accuracy `accuracy_order` and scalar type `T`.
The evaluation of the derivative can be parallised using threads by chosing
`parallel=Val{:threads}())`.
"""
function periodic_central_derivative_coefficients(derivative_order, accuracy_order, T=Float64, parallel=Val{:serial}())
    @argcheck derivative_order > 0
    @argcheck accuracy_order > 0
    @argcheck typeof(parallel) <: Union{Val{:serial}, Val{:threads}}

    if isodd(accuracy_order)
        warn("Central derivative operators support only even orders of accuracy.")
    end

    if derivative_order == 1
        # Exact evaluation of the coefficients, see
        # Beljadid, LeFloch, Mishra, Parés (2017)
        # Schemes with Well-Controlled Dissipation. Hyperbolic Systems in Nonconservative Form.
        # Communications in Computational Physics 21.4, pp. 913-946.
        n = (accuracy_order+1) ÷ 2
        coef = Vector{T}(n)
        for j in 1:n
            tmp = one(Rational{Int128})
            if iseven(j)
                tmp *= -1
            end
            for i in 1:j
                tmp *= (n+1-i) // (n+i)
            end
            coef[j] =  tmp / j
        end
        upper_coef = SVector{n,T}(coef)
        central_coef = zero(T)
        lower_coef = -upper_coef
    elseif derivative_order == 2
        # Exact evaluation of the coefficients, see
        # Beljadid, LeFloch, Mishra, Parés (2017)
        # Schemes with Well-Controlled Dissipation. Hyperbolic Systems in Nonconservative Form.
        # Communications in Computational Physics 21.4, pp. 913-946.
        n = (accuracy_order+1) ÷ 2
        coef = Vector{T}(n)
        for j in 1:n
            tmp = one(Rational{Int128})
            if iseven(j)
                tmp *= -1
            end
            for i in 1:j
                tmp *= (n+1-i) // (n+i)
            end
            coef[j] =  2*tmp / j^2
        end
        upper_coef = SVector{n,T}(coef)
        central_coef = -2*sum(upper_coef)
        lower_coef = upper_coef
    else
        throw(ArgumentError("Derivative order $derivative_order not implemented yet."))
    end
    PeriodicDerivativeCoefficients(lower_coef, central_coef, upper_coef, parallel, derivative_order, accuracy_order)
end


"""
    PeriodicDerivativeOperator{T,StencilWidth,Parallel}

A derivative operator on a uniform periodic grid with grid spacing `Δx`.
"""
struct PeriodicDerivativeOperator{T,LowerOffset,UpperOffset,Parallel} <: AbstractDerivativeOperator{T}
    coefficients::PeriodicDerivativeCoefficients{T,LowerOffset,UpperOffset,Parallel}
    Δx::T
    factor::T

    function PeriodicDerivativeOperator(coefficients::PeriodicDerivativeCoefficients{T,LowerOffset,UpperOffset,Parallel}, Δx::T) where {T,LowerOffset,UpperOffset,Parallel}
        factor = inv(Δx^coefficients.derivative_order)
        new{T,LowerOffset,UpperOffset,Parallel}(coefficients, Δx, factor)
    end
end

derivative_order(D::AbstractDerivativeOperator) = derivative_order(D.coefficients)
accuracy_order(D::AbstractDerivativeOperator) = accuracy_order(D.coefficients)
Base.eltype(D::AbstractDerivativeOperator{T}) where {T} = T


function Base.A_mul_B!(dest, D::PeriodicDerivativeOperator, u)
    mul!(dest, D, u, one(eltype(dest)), zero(eltype(dest)))
end

"""
    mul!(dest::AbstractVector, D::PeriodicDerivativeOperator, u::AbstractVector, α, β)

Compute `α*D*u + β*dest` and store the result in `dest`.
"""
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::PeriodicDerivativeOperator, u::AbstractVector, α, β)
    mul!(dest, D.coefficients, u, α*D.factor, β)
end


"""
    periodic_central_derivative_operator(Δx, derivative_order, accuracy_order, parallel=Val{:serial}())

Create a `PeriodicDerivativeOperator` approximating the `derivative_order`-th
derivative on a uniform grid with spacing `Δx` and up to order of accuracy
`accuracy_order`.
The evaluation of the derivative can be parallised using threads by chosing
`parallel=Val{:threads}())`.
"""
function periodic_central_derivative_operator(Δx, derivative_order, accuracy_order, parallel=Val{:serial}())
    coefficients = periodic_central_derivative_coefficients(derivative_order, accuracy_order, typeof(Δx), parallel)
    PeriodicDerivativeOperator(coefficients, Δx)
end
