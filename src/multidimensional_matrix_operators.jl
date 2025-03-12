"""
    MultidimensionalMatrixDerivativeOperator{Dim, T}
    MultidimensionalMatrixDerivativeOperator(nodes::NodesType, on_boundary::Vector{Bool}, normals::Vector{SVector{Dim,T}},
                                             weights::Vector{T}, weights_boundary::Vector{T}, Ds::NTuple{Dim,DType},
                                             accuracy_order::Int, source::SourceOfCoefficients) where {Dim,T<:Real,NodesType,DType<:AbstractMatrix{T}}

Multidimensional operator that represents a first-derivative operator based on matrices.

An instance of this type can be constructed by passing the nodes `nodes` (e.g. a `Vector{SVector}`), a vector of booleans `on_boundary` that
indicates whether a node is on the boundary or not, the normal vectors `normals` of the boundary nodes, the weights `weights` and `weights_boundary`
of the operator, the derivative matrices `Ds` in each direction, the `accuracy_order` of the operator, and the `source` of coefficients, which
can be `nothing` for experimentation. The lengths of `normals` and `weights_boundary` should be the same and should coincide with the number of
boundary nodes (i.e. the number if `true`s in `on_boundary`).

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
    on_boundary::Vector{Bool} # length(on_boundary) == N
    normals::Vector{SVector{Dim,T}} # length(normals) == N_boundary < N
    weights::Vector{T} # length(weights) == N
    weights_boundary::Vector{T} # length(weights_boundary) == N_boundary < N
    Ds::NTuple{Dim,DType}
    accuracy_order::Int
    source::SourceOfCoefficients

    function MultidimensionalMatrixDerivativeOperator(nodes::NodesType,
        on_boundary::Vector{Bool},
        normals::Vector{SVector{Dim,T}},
        weights::Vector{T}, weights_boundary::Vector{T},
        Ds::NTuple{Dim,DType}, accuracy_order::Int,
        source::SourceOfCoefficients) where {Dim,T<:Real,NodesType,DType<:AbstractMatrix{T},SourceOfCoefficients}
        new{Dim,T,NodesType,DType,SourceOfCoefficients}(nodes, on_boundary, normals, weights, weights_boundary, Ds, accuracy_order, source)
    end
end

Base.ndims(::MultidimensionalMatrixDerivativeOperator{Dim}) where {Dim} = Dim
Base.getindex(D::MultidimensionalMatrixDerivativeOperator, i::Int) = D.Ds[i]

function integrate_boundary(func, u, D::MultidimensionalMatrixDerivativeOperator, dim)
    return integrate(func, u, weights_boundary_scaled(D, dim))
end

weights_boundary(D::MultidimensionalMatrixDerivativeOperator) = get_weight_boundary.(Ref(D), 1:length(grid(D)))
weights_boundary_scaled(D::MultidimensionalMatrixDerivativeOperator, dim::Int) = get_weight_boundary_scaled.(Ref(D), Ref(dim), 1:length(grid(D)))

"""
    mass_matrix_boundary(D::MultidimensionalMatrixDerivativeOperator, dim)

Construct the mass matrix at the boundary of a `MultidimensionalMatrixDerivativeOperator` `D` in direction `dim`.
The boundary mass matrix is constructed to be mimetix, see

- Jan Glaubitz, Simon-Christian Klein, Jan Nordström, Philipp Öffner (2023)
  Multi-dimensional summation-by-parts operators for general function spaces: Theory and construction
  Journal of Computational Physics 491, 112370, [DOI: 10.1016/j.jcp.2023.112370](https://doi.org/10.1016/j.jcp.2023.112370).
"""
mass_matrix_boundary(D::MultidimensionalMatrixDerivativeOperator, dim::Int) = Diagonal(weights_boundary_scaled(D, dim))

Base.eltype(::MultidimensionalMatrixDerivativeOperator{Dim,T}) where {Dim,T} = T

function get_weight_boundary(D::MultidimensionalMatrixDerivativeOperator, i::Int)
    @unpack weights_boundary, on_boundary = D
    N, _ = size(D)
    @boundscheck begin
        @argcheck 1 <= i <= N
    end
    if !on_boundary[i]
        return zero(eltype(D))
    end
    j = sum(view(on_boundary, 1:i))
    @inbounds ω = weights_boundary[j]
    return ω
end

function get_weight_boundary_scaled(D::MultidimensionalMatrixDerivativeOperator, dim::Int, i::Int)
    @unpack normals, on_boundary = D
    if !on_boundary[i]
        return zero(eltype(D))
    end
    ω = get_weight_boundary(D, i)
    j = sum(view(on_boundary, 1:i))
    ω * normals[j][dim]
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

# TODO: mul! How? Depends on direction

"""
    tensor_product_operator_2D(D)

Create a 2D [`MultidimensionalMatrixDerivativeOperator`](@ref) on a square based on a 1D derivative operator `D` using a tensor product structure.

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
function tensor_product_operator_2D(D)
    nodes_1D = grid(D)
    N = length(nodes_1D)
    nodes = SVector.(vec(nodes_1D' .* ones(N)), vec(ones(N)' .* nodes_1D))
    on_boundary = zeros(Bool, N^2)
    on_boundary[1:N] .= true
    on_boundary[N:N:end-N + 1] .= true
    on_boundary[N + 1:N:end-N + 1] .= true
    on_boundary[end-N+1:end] .= true

    D_1D = Matrix(D)
    M_1D = mass_matrix(D)
    M = kron(M_1D, M_1D)
    D_x = kron(D_1D, Matrix(I, (N, N)))
    D_y = kron(Matrix(I, (N, N)), D_1D)

    weights = diag(M)
    Ds = (sparse(D_x), sparse(D_y))
    normals = Vector{SVector{2,Float64}}(undef, 4 * N - 4)
    # weights_boundary is chosen such that
    # mass_matrix_boundary(D, 1) == Diagonal(kron(B_1D, M_1D)) ( = Q_x + Q_x') and
    # mass_matrix_boundary(D, 2) == Diagonal(kron(M_1D, B_1D)) ( = Q_y + Q_y')
    weights_boundary = Vector{Float64}(undef, 4 * N - 4)
    j = 0
    for i in eachindex(normals)
        if i == 1 # lower left corner
            normals[i] = SVector(-1.0, -1.0)
            weights_boundary[i] = M_1D[1, 1]
        elseif i < N # left boundary
            normals[i] = SVector(-1.0, 0.0)
            weights_boundary[i] = M_1D[i, i]
        elseif i == N # upper left corner
            normals[i] = SVector(-1.0, 1.0)
            weights_boundary[i] = M_1D[N, N]
        elseif i < 3 * N - 3
            if (i - N - 1) % 2 == 0 # lower boundary
                normals[i] = SVector(0.0, -1.0)
                k = i - N + 1 - j
                weights_boundary[i] = M_1D[k, k]
            else # upper boundary
                normals[i] = SVector(0.0, 1.0)
                k = i - N - j
                weights_boundary[i] = M_1D[k, k]
                j += 1
            end
        elseif i == 3 * N - 3 # lower right corner
            normals[i] = SVector(1.0, -1.0)
            weights_boundary[i] = M_1D[1, 1]
        elseif i < 4 * N - 4 # right boundary
            normals[i] = SVector(1.0, 0.0)
            k = i - (3 * N - 4)
            weights_boundary[i] = M_1D[k, k]
        else # i == 4 * N - 4 # upper right corner
            normals[i] = SVector(1.0, 1.0)
            weights_boundary[i] = M_1D[N, N]
        end
    end
    return MultidimensionalMatrixDerivativeOperator(nodes, on_boundary, normals, weights, weights_boundary, Ds,
                                                    accuracy_order(D), source_of_coefficients(D))
end