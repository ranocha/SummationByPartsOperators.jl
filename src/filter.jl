
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
mul!(dest, fact::FactorisationWrapper, u) = ldiv!(dest, fact.fact, u)


#TODO
#struct AdaptiveFilter{T<:Real} <: AbstractFilter{T}
#
#end


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
        @argcheck length(coefficients) == size(tmp,1)

        new{T,Nodal2Modal,Modal2Nodal,Tmp,FilterFunction}(coefficients, nodal2modal, modal2nodal, tmp, filter)
    end
end

(filter::ConstantFilter)(u) = filter(u, filter.tmp)

function (filter::ConstantFilter)(u::AbstractVector, tmp::AbstractVector)
    @unpack coefficients, nodal2modal, modal2nodal = filter
    @boundscheck begin
        length(coefficients) == length(tmp)
    end

    mul!(tmp, nodal2modal, u)
    @inbounds tmp .*= coefficients
    mul!(u, modal2nodal, tmp)

    nothing
end

function (filter::ConstantFilter)(u::AbstractMatrix, tmp::AbstractVector)
    @boundscheck begin
        length(tmp) == size(u,1)
    end

    for j in Base.OneTo(size(u,2))
        v = view(u, :, j)
        @inbounds filter(v, tmp)
    end

    nothing
end

function (filter::ConstantFilter)(u::AbstractMatrix, tmp::AbstractMatrix)
    @unpack coefficients, nodal2modal, modal2nodal = filter
    @boundscheck begin
        size(tmp) == size(u)
    end

    mul!(tmp, nodal2modal, u)
    @inbounds tmp .*= coefficients
    mul!(u, modal2nodal, tmp)

    nothing
end



"""
    ExponentialFilter

Represents the exponential filter function `σ(η) = exp(-α*η^p)`.
"""
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
