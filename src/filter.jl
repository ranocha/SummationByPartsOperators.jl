
"""
    FactorisationWrapper

A small wrapper around a a factorisation `fact`, allowing to represent
multiplication with the inverse of `fact`.
"""
struct FactorisationWrapper{Fact}
    fact::Fact
end

Base.size(fact::FactorisationWrapper) = size(fact.fact)
Base.eltype(fact::FactorisationWrapper) = eltype(fact.fact)
Base.A_mul_B!(dest, fact::FactorisationWrapper, u) = A_ldiv_B!(dest, fact.fact, u)


struct AdaptiveFilter{T<:Real} <: AbstractFilter{T}

end


"""
    ConstantFilter

Represents the action of a modal filter on values in a nodal basis with fixed
strength.
"""
struct ConstantFilter{T<:Real,Nodal2Modal,Modal2Nodal,Tmp,FilterFunction} <: AbstractFilter{T}
    coefficients::Vector{T}
    nodal2modal::Nodal2Modal
    modal2nodal::Modal2Nodal
    tmp::Tmp
    filter::FilterFunction

    function ConstantFilter(coefficients::Vector{T}, nodal2modal::Nodal2Modal, modal2nodal::Modal2Nodal, tmp::Tmp, filter::FilterFunction) where {T<:Real,Nodal2Modal,Modal2Nodal,Tmp,FilterFunction}
        @argcheck size(coefficients) == size(tmp)

        new{T,Nodal2Modal,Modal2Nodal,Tmp,FilterFunction}(coefficients, nodal2modal, modal2nodal, tmp, filter)
    end
end

function (filter::ConstantFilter)(u::AbstractVector)
    @unpack coefficients, nodal2modal, modal2nodal, tmp = filter
    @boundscheck begin
        length(coefficients) == length(tmp)
    end

    A_mul_B!(tmp, nodal2modal, u)
    @inbounds tmp .*= coefficients
    A_mul_B!(u, modal2nodal, tmp)

    nothing
end



doc"
    ExponentialFilter

Represents the exponential filter function `σ(η) = exp(-α*η^p)`.
"
struct ExponentialFilter{T<:Real} <: AbstractFilterFunction
    α::T
    p::Int

    function ExponentialFilter(α::T, p::Int=2) where {T<:Real}
        α < 0 && warn("α should be nonnegative [α = $α].")
        p < 0 && warn("p should be nonnegative [p = $p].")
        new{T}(α, p)
    end
end

function ExponentialFilter(::Type{T}=Float64, p::Int=2) where {T<:Real}
    α = -log(eps(T))
    ExponentialFilter(α, p)
end

function ExponentialFilter(p::Int, ::Type{T}=Float64) where {T<:Real}
    α = -log(eps(T))
    ExponentialFilter(α, p)
end

function set_filter_coefficients!(coefficients::Vector, filter::ExponentialFilter)
    T = eltype(coefficients)
    N = length(coefficients) - 1
    @unpack α, p = filter

    @inbounds for i in Base.OneTo(length(coefficients))
        η = T(i-1) / N
        coefficients[i] = exp(-α*η^p)
    end

    nothing
end
