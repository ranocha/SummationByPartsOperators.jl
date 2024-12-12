
"""
    LegendreDerivativeOperator{T<:Real}

A derivative operator on a nonperiodic Lobatto-Legendre grid with scalar type
`T` computing the first derivative using a Legendre expansion.
"""
@auto_hash_equals struct LegendreDerivativeOperator{T<:Real} <: AbstractNonperiodicDerivativeOperator{T}
    jac::T
    Δx::T
    grid::Vector{T} # N nodes, including the left and the right boundary
    basis::LobattoLegendre{T}

    function LegendreDerivativeOperator(xmin::T, xmax::T, basis::LobattoLegendre{T}) where {T<:Real}
        grid = map_from_canonical.(basis.nodes, xmin, xmax, basis)
        jac = 2 / (xmax - xmin)
        Δx = inv(jac)

        new{T}(jac, Δx, grid, basis)
    end
end

"""
    LegendreDerivativeOperator(xmin::T, xmax::T, N::Int) where {T<:Real}

Construct the `LegendreDerivativeOperator` on a uniform grid between `xmin` and
`xmax` using `N` nodes and `N-1` Legendre modes.
"""
function LegendreDerivativeOperator(xmin::T, xmax::T, N::Int) where {T<:Real}
    @argcheck N >= 2

    basis = LobattoLegendre(N-1, T)

    LegendreDerivativeOperator(xmin, xmax, basis)
end

"""
    legendre_derivative_operator(xmin::Real, xmax::Real, N::Integer)
    legendre_derivative_operator(; xmin::Real, xmax::Real, N::Integer)

Construct the `LegendreDerivativeOperator` on a uniform grid between `xmin` and
`xmax` using `N` nodes and `N-1` Legendre modes.
"""
function legendre_derivative_operator(xmin::Real, xmax::Real, N::Integer)
    LegendreDerivativeOperator(promote(xmin, xmax)..., N)
end

function legendre_derivative_operator(; xmin::Real, xmax::Real, N::Integer)
    legendre_derivative_operator(xmin, xmax, N)
end

derivative_order(D::LegendreDerivativeOperator) = 1
LinearAlgebra.issymmetric(D::LegendreDerivativeOperator) = false


"""
    LegendreSecondDerivativeOperator{T<:Real}

A derivative operator on a nonperiodic Lobatto-Legendre grid with scalar type
`T` computing the second derivative using a Legendre expansion.
"""
@auto_hash_equals struct LegendreSecondDerivativeOperator{T<:Real} <: AbstractNonperiodicDerivativeOperator{T}
    D2::Matrix{T}
    jac::T
    Δx::T
    grid::Vector{T} # N nodes, including the left and the right boundary
    basis::LobattoLegendre{T}

    function LegendreSecondDerivativeOperator(xmin::T, xmax::T, N::Int) where {T<:Real}
        @argcheck N >= 2

        basis = LobattoLegendre(N-1, T)

        D2 = basis.D^2
        grid = map_from_canonical.(basis.nodes, xmin, xmax, basis)
        jac = 2 / (xmax - xmin)
        Δx = inv(jac)

        new{T}(D2, jac, Δx, grid, basis)
    end
end

"""
    legendre_second_derivative_operator(xmin::Real, xmax::Real, N::Integer)
    legendre_second_derivative_operator(; xmin::Real, xmax::Real, N::Integer)

Construct the `LegendreDerivativeOperator` on a uniform grid between `xmin` and
`xmax` using `N` nodes and `N-1` Legendre modes.
"""
function legendre_second_derivative_operator(xmin::Real, xmax::Real, N::Integer)
  LegendreSecondDerivativeOperator(promote(xmin, xmax)..., N)
end

function legendre_second_derivative_operator(; xmin::Real, xmax::Real, N::Integer)
    legendre_second_derivative_operator(xmin, xmax, N)
end

derivative_order(D::LegendreSecondDerivativeOperator) = 2


integrate(func, u, D::Union{LegendreDerivativeOperator,LegendreSecondDerivativeOperator}) = D.Δx*integrate(func, u, D.basis)
function evaluate_coefficients(u, D::Union{LegendreDerivativeOperator,LegendreSecondDerivativeOperator},
                               npoints=2*size(D,2)+1)
    evaluate_coefficients(u, D.basis, npoints)
end

mass_matrix(D::Union{LegendreDerivativeOperator,LegendreSecondDerivativeOperator}) = Diagonal(D.Δx * D.basis.weights)

Base.eltype(D::Union{LegendreDerivativeOperator{T},LegendreSecondDerivativeOperator{T}}) where {T} = T

function scale_by_mass_matrix!(u::AbstractVector, D::Union{LegendreDerivativeOperator,LegendreSecondDerivativeOperator}, factor=true)
    Base.require_one_based_indexing(u)
    @boundscheck begin
        length(u) == size(D, 2) || throw(DimensionMismatch("sizes of input vector and operator do not match"))
    end
    @unpack Δx, basis = D

    @inbounds @simd for i in eachindex(u, basis.weights)
        u[i] = factor * u[i] * (Δx * basis.weights[i])
    end

    return u
end

function scale_by_inverse_mass_matrix!(u::AbstractVector, D::Union{LegendreDerivativeOperator,LegendreSecondDerivativeOperator}, factor=true)
    Base.require_one_based_indexing(u)
    @boundscheck begin
        length(u) == size(D, 2) || throw(DimensionMismatch("sizes of input vector and operator do not match"))
    end
    @unpack Δx, basis = D

    @inbounds @simd for i in eachindex(u, basis.weights)
        u[i] = factor * u[i] / (Δx * basis.weights[i])
    end

    u
end

function get_weight(D::Union{LegendreDerivativeOperator,LegendreSecondDerivativeOperator}, i::Int)
    @unpack Δx, basis = D
    @unpack weights = basis
    N, _ = size(D)
    @boundscheck begin
        @argcheck 1 <= i <= N
    end
    @inbounds ω = Δx * weights[i]
    ω
end

function Base.show(io::IO, D::LegendreDerivativeOperator)
    if get(io, :compact, false)
        summary(io, D)
    else
        x = grid(D)
        print(io, "First derivative operator {T=", eltype(D), "}")
        print(io, " on ", length(x), " Lobatto Legendre nodes in [", first(x), ", ", last(x), "]")
    end
end

function Base.show(io::IO, D::LegendreSecondDerivativeOperator)
    if get(io, :compact, false)
        summary(io, D)
    else
        x = grid(D)
        print(io, "Second derivative operator {T=", eltype(D), "}")
        print(io, " on ", length(x), " Lobatto Legendre nodes in [", first(x), ", ", last(x), "]")
    end
end


function mul!(dest::AbstractVector, D::LegendreDerivativeOperator, u::AbstractVector, α=true, β=false)
    @unpack jac, basis = D
    N, _ = size(D)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    mul!(dest, basis.D, u, α*jac, β)
end

function mul!(dest::AbstractVector, D::LegendreSecondDerivativeOperator, u::AbstractVector, α=true, β=false)
    @unpack jac, D2 = D
    N, _ = size(D2)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    mul!(dest, D2, u, α*jac^2, β)
end

function derivative_left(D::LegendreSecondDerivativeOperator, u::AbstractVector, ::Val{1})
    @unpack jac, basis = D
    N, _ = size(D)
    @boundscheck begin
        @argcheck N == length(u)
    end

    jac * dot(view(basis.D, 1, :), u)
end

function derivative_right(D::LegendreSecondDerivativeOperator, u::AbstractVector, ::Val{1})
    @unpack jac, basis = D
    N, _ = size(D)
    @boundscheck begin
        @argcheck N == length(u)
    end

    jac * dot(view(basis.D, size(D,1), :), u)
end


function lower_bandwidth(D::Union{LegendreDerivativeOperator,LegendreSecondDerivativeOperator})
    size(D, 1) - 1
end

function upper_bandwidth(D::Union{LegendreDerivativeOperator,LegendreSecondDerivativeOperator})
    size(D, 1) - 1
end

function accuracy_order(D::Union{LegendreDerivativeOperator,LegendreSecondDerivativeOperator})
    size(D, 1) - 1
end


function left_boundary_weight(D::Union{LegendreDerivativeOperator,LegendreSecondDerivativeOperator})
    @inbounds retval = D.Δx * D.basis.weights[1]
    retval
end

function right_boundary_weight(D::Union{LegendreDerivativeOperator,LegendreSecondDerivativeOperator})
    @inbounds retval = D.Δx * D.basis.weights[end]
    retval
end


function ConstantFilter(basis::LobattoLegendre{T}, filter, tmp::Array) where {T}
    Np1 = length(basis.nodes)
    coefficients = Array{T}(undef, Np1)
    set_filter_coefficients!(coefficients, filter)
    modal2nodal = legendre_vandermonde(basis)
    nodal2modal = FactorisationWrapper(factorize(modal2nodal))

    ConstantFilter(coefficients, nodal2modal, modal2nodal, tmp, filter)
end

function ConstantFilter(basis::LobattoLegendre{T}, filter, TmpEltype::DataType=T) where {T}
    Np1 = length(basis.nodes)
    tmp = Array{TmpEltype}(undef, Np1)

    ConstantFilter(basis, filter, tmp)
end

"""
    ConstantFilter(D::LegendreDerivativeOperator, filter, TmpEltype=T)

Create a modal filter with constant parameters adapted to the Legendre
derivative operator `D` with parameters given by the filter function `filter`.
"""
function ConstantFilter(D::LegendreDerivativeOperator{T}, filter, TmpEltype=T) where {T}
    ConstantFilter(D.basis, filter, TmpEltype)
end


# TODO: LegendreSuperSpectralViscosity
