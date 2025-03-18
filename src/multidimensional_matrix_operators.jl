"""
    MultidimensionalMatrixDerivativeOperator{Dim, T}
    MultidimensionalMatrixDerivativeOperator(nodes::NodesType, boundary_indices::Vector{Int}, normals::Vector{SVector{Dim,T}},
                                             weights::Vector{T}, weights_boundary::Vector{T}, Ds::NTuple{Dim,DType},
                                             accuracy_order::Int, source::SourceOfCoefficients) where {Dim,T<:Real,NodesType,DType<:AbstractMatrix{T}}

Multidimensional operator that represents a first-derivative operator based on matrices.

An instance of this type can be constructed by passing the nodes `nodes` (e.g. a `Vector{SVector}`), a vector of indices `boundary_indices` that
indicates which nodes are on the boundary, the normal vectors `normals` of the boundary nodes, the weights `weights` and `weights_boundary`
of the operator, the derivative matrices `Ds` in each direction, the `accuracy_order` of the operator, and the `source` of coefficients, which
can be `nothing` for experimentation. The lengths of `boundary_indices`, `normals`, and `weights_boundary` should be the same and should coincide
with the number of boundary nodes. Indices in `boundary_indices` can appear multiple times if the boundary nodes are shared between
different boundaries, e.g. corners.

To obtain the derivative operator in a specific direction, use `D[dim]`. The boundary operator in a specific direction can be obtained with
`mass_matrix_boundary(D, dim)` and will be constructed as a mimetic operator based on `weights_boundary`. The mass matrix of the operator
is given by `mass_matrix(D)`.

See also [`TensorProductOperator`](@ref).

References:
- Jason E. Hicken, David C. Del Rey Fernandez, and David W. Zingg (2016)
  Multidimensional Summation-by-Parts Operators: General Theory and Application to Simplex Elements
  SIAM Journal of Scientific Computing 38(4), pp. A1935-A1958, [DOI: 10.1137/15M1038360](https://doi.org/10.1137/15M1038360).
- Jan Glaubitz, Simon-Christian Klein, Jan Nordström, Philipp Öffner (2023)
  Multi-dimensional summation-by-parts operators for general function spaces: Theory and construction
  Journal of Computational Physics 491, 112370, [DOI: 10.1016/j.jcp.2023.112370](https://doi.org/10.1016/j.jcp.2023.112370).
"""
@auto_hash_equals struct MultidimensionalMatrixDerivativeOperator{
    Dim,
    T,
    NodesType,
    DType<:AbstractMatrix{T},
    SourceOfCoefficients,
} <: AbstractMultidimensionalMatrixDerivativeOperator{Dim,T}
    grid::NodesType # length(grid) == N, e.g. Vector{SVector{Dim, T}} or `NodeSet` from KernelInterpolation.jl
    boundary_indices::Vector{Int} # length(boundary_indices) == N_boundary
    normals::Vector{SVector{Dim,T}} # length(normals) == N_boundary
    weights::Vector{T} # length(weights) == N
    weights_boundary::Vector{T} # length(weights_boundary) == N_boundary
    Ds::NTuple{Dim,DType}
    accuracy_order::Int
    source::SourceOfCoefficients

    function MultidimensionalMatrixDerivativeOperator(
        nodes::NodesType,
        boundary_indices::Vector{Int},
        normals::Vector{SVector{Dim,T}},
        weights::Vector{T},
        weights_boundary::Vector{T},
        Ds::NTuple{Dim,DType},
        accuracy_order::Int,
        source::SourceOfCoefficients,
    ) where {Dim,T<:Real,NodesType,DType<:AbstractMatrix{T},SourceOfCoefficients}
        if length(nodes) != length(weights)
            throw(ArgumentError("The number of nodes and weights should be the same"))
        end
        if !(length(boundary_indices) == length(normals) == length(weights_boundary))
            throw(
                ArgumentError(
                    "The number of boundary indices, normals, and boundary weights should be the same",
                ),
            )
        end
        new{Dim,T,NodesType,DType,SourceOfCoefficients}(
            nodes,
            boundary_indices,
            normals,
            weights,
            weights_boundary,
            Ds,
            accuracy_order,
            source,
        )
    end
end

Base.ndims(::AbstractMultidimensionalMatrixDerivativeOperator{Dim}) where {Dim} = Dim
Base.getindex(D::AbstractMultidimensionalMatrixDerivativeOperator, i::Int) = D.Ds[i]
Base.eltype(::AbstractMultidimensionalMatrixDerivativeOperator{Dim,T}) where {Dim,T} = T

restrict_boundary(u, D::AbstractMultidimensionalMatrixDerivativeOperator) =
    u[D.boundary_indices]

"""
    integrate_boundary([func = identity,] u, D::MultidimensionalMatrixDerivativeOperator, dim)

Map the function `func` to the coefficients `u` and integrate along the boundary in direction `dim` with respect to
the surface quadrature rule associated with the [`MultidimensionalMatrixDerivativeOperator`](@ref) `D`.
"""
function integrate_boundary(
    func,
    u,
    D::AbstractMultidimensionalMatrixDerivativeOperator,
    dim,
)
    return integrate(func, restrict_boundary(u, D), weights_boundary_scaled(D, dim))
end

integrate_boundary(u, D::AbstractMultidimensionalMatrixDerivativeOperator, dim) =
    integrate_boundary(identity, u, D, dim)

weights_boundary(D::AbstractMultidimensionalMatrixDerivativeOperator) =
    get_weight_boundary.(Ref(D), 1:length(D.weights_boundary))
weights_boundary_scaled(D::AbstractMultidimensionalMatrixDerivativeOperator, dim::Int) =
    get_weight_boundary_scaled.(Ref(D), Ref(dim), 1:length(D.weights_boundary))

"""
    mass_matrix_boundary(D::AbstractMultidimensionalMatrixDerivativeOperator, dim)

Construct the mass matrix at the boundary of a [`MultidimensionalMatrixDerivativeOperator`](@ref) `D` in direction `dim`.
The boundary mass matrix is constructed to be mimetic, see

- Jan Glaubitz, Simon-Christian Klein, Jan Nordström, Philipp Öffner (2023)
  Multi-dimensional summation-by-parts operators for general function spaces: Theory and construction
  Journal of Computational Physics 491, 112370, [DOI: 10.1016/j.jcp.2023.112370](https://doi.org/10.1016/j.jcp.2023.112370).
"""
function mass_matrix_boundary(D::MultidimensionalMatrixDerivativeOperator, dim::Int)
    if length(unique(D.boundary_indices)) != length(D.boundary_indices)
        throw(ArgumentError("The boundary indices should be unique"))
    end
    b = zeros(eltype(D), length(D.grid))
    b[D.boundary_indices] .= weights_boundary_scaled(D, dim)
    return Diagonal(b)
end

function get_weight_boundary(D::AbstractMultidimensionalMatrixDerivativeOperator, i::Int)
    @unpack weights_boundary = D
    N_boundary = length(weights_boundary)
    @boundscheck begin
        @argcheck 1 <= i <= N_boundary
    end
    @inbounds ω = weights_boundary[i]
    return ω
end

function get_weight_boundary_scaled(
    D::AbstractMultidimensionalMatrixDerivativeOperator,
    dim::Int,
    i::Int,
)
    @unpack normals = D
    ω = get_weight_boundary(D, i)
    return ω * normals[i][dim]
end

function Base.show(io::IO, D::MultidimensionalMatrixDerivativeOperator)
    if get(io, :compact, false)
        summary(io, D)
    else
        x = grid(D)
        print(
            io,
            ndims(D),
            "-dimensional matrix-based first-derivative operator {T=",
            eltype(D),
            "}",
        )
        print(io, " on ", length(x), " nodes")
    end
end
