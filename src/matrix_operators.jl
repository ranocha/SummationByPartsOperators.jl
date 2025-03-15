
"""
    MatrixDerivativeOperator{T <: Real}
    MatrixDerivativeOperator(xmin, xmax, nodes, weights, D, accuracy_order, source)

A derivative operator on a nonperiodic grid with scalar type `T` computing a
derivative as matrix vector product. This type is designed to make it easy to
experiment with new operators given in matrix form.

An instance of this type can be constructed by passing the endpoints
`xmin`, `xmax` of the desired grid as well as the `nodes`, `weights`, and the
derivative operator `D::AbstractMatrix` on a reference interval, assuming that the
`nodes` contain the boundary points of the reference interval. `source` is
the source of coefficients and can be `nothing` for experimentation.
"""
@auto_hash_equals struct MatrixDerivativeOperator{T, Dtype <: AbstractMatrix{T}, SourceOfCoefficients} <: AbstractMatrixDerivativeOperator{T}
  grid::Vector{T}
  weights::Vector{T}
  D::Dtype
  accuracy_order::Int
  source::SourceOfCoefficients

  function MatrixDerivativeOperator(xmin::T, xmax::T,
                                    nodes::Vector{T},
                                    weights::Vector{T},
                                    D::AbstractMatrix{T},
                                    accuracy_order::Int,
                                    source::SourceOfCoefficients) where {T <: Real, SourceOfCoefficients}
      # The `nodes`, `weights`, and `D` are given on a reference interval.
      # We need to scale them by the Jacobian to get their values on the
      # given interval.
      jac = (last(nodes) - first(nodes)) / (xmax - xmin)
      grid = (nodes .- first(nodes)) ./ jac .+ xmin
      Δx = inv(jac)

      new{T, typeof(D), SourceOfCoefficients}(grid, Δx * weights, jac * D, accuracy_order, source)
  end
end

derivative_order(D::AbstractMatrixDerivativeOperator) = 1
LinearAlgebra.issymmetric(D::MatrixDerivativeOperator) = false

source_of_coefficients(D::AbstractMatrixDerivativeOperator) = D.source

function integrate(func, u, D::AbstractMatrixDerivativeOperator)
  return integrate(func, u, D.weights)
end

mass_matrix(D::AbstractMatrixDerivativeOperator) = Diagonal(D.weights)

Base.eltype(D::MatrixDerivativeOperator{T}) where {T} = T

function scale_by_mass_matrix!(u::AbstractVector, D::AbstractMatrixDerivativeOperator, factor=true)
    Base.require_one_based_indexing(u)
    @boundscheck begin
        length(u) == size(D, 2) || throw(DimensionMismatch("sizes of input vector and operator do not match"))
    end

    @inbounds @simd for i in eachindex(u, D.weights)
        u[i] = factor * u[i] * D.weights[i]
    end

    return u
end

function scale_by_inverse_mass_matrix!(u::AbstractVector, D::AbstractMatrixDerivativeOperator, factor=true)
    Base.require_one_based_indexing(u)
    @boundscheck begin
        length(u) == size(D, 2) || throw(DimensionMismatch("sizes of input vector and operator do not match"))
    end

    @inbounds @simd for i in eachindex(u, D.weights)
        u[i] = factor * u[i] / D.weights[i]
    end

    u
end

function get_weight(D::AbstractMatrixDerivativeOperator, i::Int)
    @unpack weights = D
    N, _ = size(D)
    @boundscheck begin
        @argcheck 1 <= i <= N
    end
    @inbounds ω = weights[i]
    ω
end

function Base.show(io::IO, D::MatrixDerivativeOperator)
    if get(io, :compact, false)
        summary(io, D)
    else
        x = grid(D)
        print(io, "Matrix-based first-derivative operator {T=", eltype(D), "}")
        print(io, " on ", length(x), " nodes in [", first(x), ", ", last(x), "]")
    end
end


# TODO: Enable different evaluation modes
function mul!(dest::AbstractVector, Dop::MatrixDerivativeOperator, u::AbstractVector, α=true, β=false)
    @unpack D = Dop
    N, _ = size(D)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    mul!(dest, D, u, α, β)
end


function lower_bandwidth(D::AbstractMatrixDerivativeOperator)
    size(D, 1) - 1
end

function upper_bandwidth(D::AbstractMatrixDerivativeOperator)
    size(D, 1) - 1
end

function accuracy_order(D::AbstractMatrixDerivativeOperator)
    D.accuracy_order
end


function left_boundary_weight(D::AbstractMatrixDerivativeOperator)
    @inbounds retval = D.weights[begin]
    retval
end

function right_boundary_weight(D::AbstractMatrixDerivativeOperator)
    @inbounds retval = D.weights[end]
    retval
end
