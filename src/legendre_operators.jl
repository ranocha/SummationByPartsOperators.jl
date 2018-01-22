
"""
    LegendreDerivativeOperator{T<:Real, Grid}

A derivative operator on a nonperiodic Lobatto-Legendre grid with scalar type
`T` computing the first derivative using a Legendre expansion.
"""
struct LegendreDerivativeOperator{T<:Real} <: AbstractDerivativeOperator{T}
    jac::T
    Δx::T
    grid::Vector{T} # N nodes, including the left and the right boundary
    basis::LobattoLegendre{T}

    function LegendreDerivativeOperator(xmin::T, xmax::T, basis::LobattoLegendre{T}) where {T<:Real}
        grid = map_from_canonical(basis.nodes, xmin, xmax, basis)
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

function legendre_derivative_operator(xmin::T, xmax::T, N::Int) where {T<:Real}
    LegendreDerivativeOperator(xmin, xmax, N)
end

derivative_order(D::LegendreDerivativeOperator) = 1
Base.issymmetric(D::LegendreDerivativeOperator) = false
integrate(func, u, D::LegendreDerivativeOperator) = D.Δx*integrate(func, u, D.basis)
function evaluate_coefficients(u, D::LegendreDerivativeOperator, 
                               npoints=2*size(D,2)+1)
    evaluate_coefficients(u, D.basis, npoints)
end

function Base.show(io::IO, D::LegendreDerivativeOperator{T}) where {T}
    x = grid(D)
    print(io, "First derivative operator {T=", T, "} \n")
    print(io, "on the Lobatto Legendre nodes in [", first(x), ", ", last(x),
                "] using ", length(x), " nodes. \n")
end


function Base.A_mul_B!(dest::AbstractVector{T}, D::LegendreDerivativeOperator,
                        u::AbstractVector{T}) where {T}
    @unpack jac, basis = D
    N, _ = size(D)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    A_mul_B!(dest, basis.D, u)
    dest .*= jac

    nothing
end


function left_boundary_weight(D::LegendreDerivativeOperator)
    @inbounds retval = D.Δx * D.basis.weights[1]
    retval
end

function right_boundary_weight(D::LegendreDerivativeOperator)
    @inbounds retval = D.Δx * D.basis.weights[end]
    retval
end



# TODO: LegendreSuperSpectralViscosity
