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
with the number of boundary nodes.

To obtain the derivative operator in a specific direction, use `D[dim]`. The boundary operator in a specific direction can be obtained with
`mass_matrix_boundary(D, dim)` and will be constructed as a mimetic operator based on `weights_boundary`. The mass matrix of the operator
is given by `mass_matrix(D)`.

See also [`tensor_product_operator_2D`](@ref).

References:
- Jason E. Hicken, David C. Del Rey Fernandez, and David W. Zingg (2016)
  Multidimensional Summation-by-Parts Operators: General Theory and Application to Simplex Elements
  SIAM Journal of Scientific Computing 38(4), pp. A1935-A1958, [DOI: 10.1137/15M1038360](https://doi.org/10.1137/15M1038360).
- Jan Glaubitz, Simon-Christian Klein, Jan Nordström, Philipp Öffner (2023)
  Multi-dimensional summation-by-parts operators for general function spaces: Theory and construction
  Journal of Computational Physics 491, 112370, [DOI: 10.1016/j.jcp.2023.112370](https://doi.org/10.1016/j.jcp.2023.112370).
"""
@auto_hash_equals struct MultidimensionalMatrixDerivativeOperator{Dim,T,NodesType,DType<:AbstractMatrix{T},SourceOfCoefficients} <: AbstractMatrixDerivativeOperator{T}
    grid::NodesType # length(grid) == N, e.g. Vector{SVector{Dim, T}} or `NodeSet` from KernelInterpolation.jl
    boundary_indices::Vector{Int} # length(boundary_indices) == N_boundary
    normals::Vector{SVector{Dim,T}} # length(normals) == N_boundary
    weights::Vector{T} # length(weights) == N
    weights_boundary::Vector{T} # length(weights_boundary) == N_boundary
    Ds::NTuple{Dim,DType}
    accuracy_order::Int
    source::SourceOfCoefficients

    function MultidimensionalMatrixDerivativeOperator(nodes::NodesType,
        boundary_indices::Vector{Int},
        normals::Vector{SVector{Dim,T}},
        weights::Vector{T}, weights_boundary::Vector{T},
        Ds::NTuple{Dim,DType}, accuracy_order::Int,
        source::SourceOfCoefficients) where {Dim,T<:Real,NodesType,DType<:AbstractMatrix{T},SourceOfCoefficients}
        if length(nodes) != length(weights)
            throw(ArgumentError("The number of nodes and weights should be the same"))
        end
        if !(length(boundary_indices) == length(normals) == length(weights_boundary))
            throw(ArgumentError("The number of boundary indices, normals, and boundary weights should be the same"))
        end
        new{Dim,T,NodesType,DType,SourceOfCoefficients}(nodes, boundary_indices, normals, weights, weights_boundary, Ds, accuracy_order, source)
    end
end

Base.ndims(::MultidimensionalMatrixDerivativeOperator{Dim}) where {Dim} = Dim
Base.getindex(D::MultidimensionalMatrixDerivativeOperator, i::Int) = D.Ds[i]
Base.eltype(::MultidimensionalMatrixDerivativeOperator{Dim,T}) where {Dim,T} = T

restrict_boundary(u, D::MultidimensionalMatrixDerivativeOperator) = u[D.boundary_indices]

"""
    integrate_boundary([func = identity,] u, D::MultidimensionalMatrixDerivativeOperator, dim)

Map the function `func` to the coefficients `u` and integrate along the boundary in direction `dim` with respect to
the surface quadrature rule associated with the [`MultidimensionalMatrixDerivativeOperator`](@ref) `D`.
"""
function integrate_boundary(func, u, D::MultidimensionalMatrixDerivativeOperator, dim)
    return integrate(func, restrict_boundary(u, D), weights_boundary_scaled(D, dim))
end

integrate_boundary(u, D::MultidimensionalMatrixDerivativeOperator, dim) = integrate_boundary(identity, u, D, dim)

weights_boundary(D::MultidimensionalMatrixDerivativeOperator) = get_weight_boundary.(Ref(D), 1:length(D.weights_boundary))
weights_boundary_scaled(D::MultidimensionalMatrixDerivativeOperator, dim::Int) = get_weight_boundary_scaled.(Ref(D), Ref(dim), 1:length(D.weights_boundary))

"""
    mass_matrix_boundary(D::MultidimensionalMatrixDerivativeOperator, dim)

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

function get_weight_boundary(D::MultidimensionalMatrixDerivativeOperator, i::Int)
    @unpack weights_boundary = D
    N_boundary = length(weights_boundary)
    @boundscheck begin
        @argcheck 1 <= i <= N_boundary
    end
    @inbounds ω = weights_boundary[i]
    return ω
end

function get_weight_boundary_scaled(D::MultidimensionalMatrixDerivativeOperator, dim::Int, i::Int)
    @unpack normals = D
    ω = get_weight_boundary(D, i)
    return ω * normals[i][dim]
end

function Base.show(io::IO, D::MultidimensionalMatrixDerivativeOperator)
    if get(io, :compact, false)
        summary(io, D)
    else
        x = grid(D)
        print(io, ndims(D), "-dimensional matrix-based first-derivative operator {T=", eltype(D), "}")
        print(io, " on ", length(x), " nodes")
    end
end

"""
    SourceOfCoefficientsCombination(sources::SourceOfCoefficients...)

Combine multiple sources of coefficients into a single source.
"""
struct SourceOfCoefficientsCombination{N} <: SourceOfCoefficients
    sources::NTuple{N,SourceOfCoefficients}
end

function SourceOfCoefficientsCombination(sources::SourceOfCoefficients...)
    return SourceOfCoefficientsCombination(sources)
end

function Base.show(io::IO, source::SourceOfCoefficientsCombination)
    println(io, "Combined sources of coefficients with ", length(source.sources), " sources:")
    for s in source.sources
        println(io, "- ", s)
    end
end

"""
    tensor_product_operator_2D(D_x, D_y = D_x)

Create a 2D [`MultidimensionalMatrixDerivativeOperator`](@ref) on a square based on a 1D derivative operators `D_x` and `D_y` using a tensor product structure.
The operator `D_x` is used in the x-direction and `D_y` in the y-direction.

For the construction, see also:
- Sigrun Ortleb (2021)
  Numerical Methods for Fluid Flow:
    High Order SBP Schemes, IMEX Advection-Diffusion Splitting and Positivity Preservation for Production-Destruction-PDEs
  Habilitation thesis, University of Kassel, [DOI: 10.17170/kobra-202301037274](https://doi.org/10.17170/kobra-202301037274),
  Chapter 1.2.4.
- Magnus Svärd, Jan Nordström (2014)
  Review of summation-by-parts schemes for initial-boundary-value problems
  Journal of Computational Physics 268, pp. 17-38, [DOI: 10.1016/j.jcp.2014.02.031](https://doi.org/10.1016/j.jcp.2014.02.031).
"""
function tensor_product_operator_2D(D_x, D_y = D_x)
    T = promote_type(eltype(D_x), eltype(D_y))
    nodes_1D_x = grid(D_x)
    nodes_1D_y = grid(D_y)
    N_x = length(nodes_1D_x)
    N_y = length(nodes_1D_y)
    nodes = SVector.(vec(nodes_1D_x' .* ones(T, N_y)), vec(ones(T, N_x)' .* nodes_1D_y))

    D_1D_x = sparse(D_x)
    M_1D_x = mass_matrix(D_x)
    D_1D_y = sparse(D_y)
    M_1D_y = mass_matrix(D_y)
    M_t = kron(M_1D_x, M_1D_y)
    D_t_x = kron(D_1D_x, Diagonal(I, N_y))
    D_t_y = kron(Diagonal(I, N_x), D_1D_y)

    weights = diag(M_t)
    Ds = (D_t_x, D_t_y)
    # Since the left and right end points are always included, we need to store the corners twice in the boundary indices,
    # once with the weight from the x-direction and once with the weight from the y-direction with corresponding normals.
    N_boundary = 2 * (N_x + N_y)
    boundary_indices = [1:N_y..., # left boundary, N_y
                        1:N_y:((N_x - 1) * N_y + 1)..., # lower boundary, N_x
                        N_y:N_y:(N_x * N_y)..., # upper boundary, N_x
                        ((N_x - 1)* N_y + 1):N_x * N_y...] # right boundary, N_y
    normals = Vector{SVector{2,T}}(undef, N_boundary)
    weights_boundary = Vector{T}(undef, N_boundary)
    for i in eachindex(boundary_indices)
        if i <= N_y # left boundary
            normals[i] = SVector(-1.0, 0.0)
            k = i
            weights_boundary[i] = M_1D_y[k, k]
        elseif i <= N_x + N_y # lower boundary
            normals[i] = SVector(0.0, -1.0)
            k = i - N_y
            weights_boundary[i] = M_1D_x[k, k]
        elseif i <= 2 * N_x + N_y # upper boundary
            normals[i] = SVector(0.0, 1.0)
            k = i - N_x - N_y
            weights_boundary[i] = M_1D_x[k, k]
        else # right boundary
            normals[i] = SVector(1.0, 0.0)
            k = i - 2 * N_x - N_y
            weights_boundary[i] = M_1D_y[k, k]
        end
    end
    acc_order = min(accuracy_order(D_x), accuracy_order(D_y))
    if source_of_coefficients(D_x) == source_of_coefficients(D_y)
        source = source_of_coefficients(D_x)
    else
        source = SourceOfCoefficientsCombination(source_of_coefficients(D_x), source_of_coefficients(D_y))
    end
    return MultidimensionalMatrixDerivativeOperator(nodes, boundary_indices, normals, weights, weights_boundary, Ds,
                                                    acc_order, source)
end