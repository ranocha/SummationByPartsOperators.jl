# We use an own structure for tensor product operators to be able to specialize some functions.
"""
    TensorProductOperator{Dim,T}

A tensor product operator is a [`MultidimensionalMatrixDerivativeOperator`](@ref) that is the tensor product of one-dimensional derivative operators.

See also [`tensor_product_operator_2D`](@ref).

References:
- Sigrun Ortleb (2021)
  Numerical Methods for Fluid Flow:
    High Order SBP Schemes, IMEX Advection-Diffusion Splitting and Positivity Preservation for Production-Destruction-PDEs
  Habilitation thesis, University of Kassel, [DOI: 10.17170/kobra-202301037274](https://doi.org/10.17170/kobra-202301037274),
  Chapter 1.2.4.
- Magnus Svärd, Jan Nordström (2014)
  Review of summation-by-parts schemes for initial-boundary-value problems
  Journal of Computational Physics 268, pp. 17-38, [DOI: 10.1016/j.jcp.2014.02.031](https://doi.org/10.1016/j.jcp.2014.02.031).
"""
@auto_hash_equals struct TensorProductOperator{Dim,T} <: AbstractMultidimensionalMatrixDerivativeOperator{Dim,T}
    D_multi::MultidimensionalMatrixDerivativeOperator{Dim,T}
    N_x::Int
    N_y::Int
end

# This allows us to treat a `TensorProductOperator` as a `MultidimensionalMatrixDerivativeOperator`.
function Base.getproperty(D::TensorProductOperator, name::Symbol)
    if name in [:D_multi, :N_x, :N_y]
        return getfield(D, name)
    else
        return getproperty(D.D_multi, name)
    end
end

function mass_matrix_boundary(D::TensorProductOperator{2}, dim::Int)
    # boundary_weights contains two weights for each of the four corners. We need to select the correct one based on the dimension.
    boundary_weights = weights_boundary_scaled(D, dim)
    boundary_indices = copy(D.boundary_indices)
    N_x, N_y = D.N_x, D.N_y
    # This exploits how we construct the tensor product operator in `tensor_product_operator_2D`.
    if dim == 1
        corners_x_dir = [N_y + 1, # lower left corner
                         N_x + N_y, # lower right corner
                         N_x + N_y + 1, # upper left corner
                         2 * N_x + N_y] # upper right corner
        deleteat!(boundary_weights, corners_x_dir)
        deleteat!(boundary_indices, corners_x_dir)
    else # dim == 2
        corners_y_dir = [1, # lower left corner
                         N_y,  # upper left corner
                         2 * N_x + N_y + 1, # lower right corner
                         2 * (N_x + N_y)] # upper right corner
        deleteat!(boundary_weights, corners_y_dir)
        deleteat!(boundary_indices, corners_y_dir)
    end
    b = zeros(eltype(D), length(D.grid))
    b[boundary_indices] .= boundary_weights
    return Diagonal(b)
end

function Base.show(io::IO, D::TensorProductOperator)
    if get(io, :compact, false)
        summary(io, D)
    else
        x = grid(D)
        print(io, ndims(D), "-dimensional tensor product operator {T=", eltype(D), "}")
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

Create a 2D [`TensorProductOperator`](@ref) on a square based on a 1D derivative operators `D_x` and `D_y` using a tensor product structure.
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
    D_multi = MultidimensionalMatrixDerivativeOperator(nodes, boundary_indices, normals, weights, weights_boundary, Ds, acc_order, source)
    return TensorProductOperator(D_multi, N_x, N_y)
end
