"""
    SphericalSBPOperators

Spherical SBP operators on a uniform radial grid. `S` and `V` are the scalar
and vector norm matrices, respectively. `G` maps an even scalar field to an
odd radial vector field and `D` maps an odd radial vector field to an even
scalar field. They satisfy
```
S * D + G' * V = B.
```

This is a *pair* of compatible operators, not a subtype of
`AbstractDerivativeOperator`: the scalar and vector fields have distinct norm
matrices. Use [`scalar_mass`](@ref) and [`vector_mass`](@ref) to access them,
and [`gradient_operator`](@ref), [`odd_derivative_operator`](@ref), or
[`divergence_operator`](@ref) to obtain a matrix-like action.
"""
@auto_hash_equals struct SphericalSBPOperators{T <: Real, S <: AbstractMatrix{T},
                                               V <: AbstractMatrix{T},
                                               B <: AbstractMatrix{T},
                                               G <: AbstractMatrix{T},
                                               Godd <: AbstractMatrix{T},
                                               D <: AbstractMatrix{T},
                                               Source <: SourceOfCoefficients,
                                               CartesianSource <: SourceOfCoefficients}
    grid::Vector{T}
    scalar_mass::S
    vector_mass::V
    boundary_matrix::B
    gradient::G
    odd_derivative::Godd
    divergence::D
    accuracy_order::Int
    p::Int
    grid_kind::Symbol
    source::Source
    cartesian_source::CartesianSource
end

grid(operators::SphericalSBPOperators) = operators.grid
accuracy_order(operators::SphericalSBPOperators) = operators.accuracy_order
source_of_coefficients(operators::SphericalSBPOperators) = operators.source
scalar_mass(operators::SphericalSBPOperators) = operators.scalar_mass
vector_mass(operators::SphericalSBPOperators) = operators.vector_mass
boundary_matrix(operators::SphericalSBPOperators) = operators.boundary_matrix

"""
    SphericalOperatorView

Matrix-like view of one action in a [`SphericalSBPOperators`](@ref) pair.
Unlike an ordinary SBP derivative operator, a spherical action maps between
fields with potentially different norms. Use [`scalar_mass`](@ref) and
[`vector_mass`](@ref) on the paired spherical operators instead of
[`mass_matrix`](@ref).
"""
@auto_hash_equals struct SphericalOperatorView{T <: Real, M <: AbstractMatrix{T},
                                               Operators} <:
                         AbstractMatrix{T}
    operators::Operators
    matrix::M
    kind::Symbol

    function SphericalOperatorView(operators::Operators, matrix::M,
                                   kind::Symbol) where {T <: Real, M <: AbstractMatrix{T},
                                                        Operators}
        kind in (:gradient, :odd_derivative, :divergence) ||
            throw(ArgumentError("unknown spherical operator action `$kind`"))
        return new{T, M, Operators}(operators, matrix, kind)
    end
end

Base.eltype(::SphericalOperatorView{T}) where {T} = T
Base.size(operator::SphericalOperatorView) = size(operator.matrix)
function Base.size(operator::SphericalOperatorView, dimension::Int)
    size(operator.matrix, dimension)
end
Base.axes(operator::SphericalOperatorView) = axes(operator.matrix)
function Base.axes(operator::SphericalOperatorView, dimension::Int)
    axes(operator.matrix, dimension)
end
function Base.getindex(operator::SphericalOperatorView, indices...)
    getindex(operator.matrix, indices...)
end
Base.Matrix(operator::SphericalOperatorView) = Matrix(operator.matrix)
SparseArrays.sparse(operator::SphericalOperatorView) = sparse(operator.matrix)
Base.:*(operator::SphericalOperatorView, u::AbstractVector) = operator.matrix * u
Base.:*(operator::SphericalOperatorView, u::AbstractMatrix) = operator.matrix * u

function mul!(dest::AbstractVector, operator::SphericalOperatorView,
              u::AbstractVector, α::Number = true, β::Number = false)
    mul!(dest, operator.matrix, u, α, β)
end

function mul!(dest::AbstractMatrix, operator::SphericalOperatorView,
              u::AbstractMatrix, α::Number = true, β::Number = false)
    mul!(dest, operator.matrix, u, α, β)
end

grid(operator::SphericalOperatorView) = grid(operator.operators)
accuracy_order(operator::SphericalOperatorView) = accuracy_order(operator.operators)
function source_of_coefficients(operator::SphericalOperatorView)
    source_of_coefficients(operator.operators)
end

"""Return the even-scalar to odd-vector spherical gradient action."""
function gradient_operator(operators::SphericalSBPOperators)
    return SphericalOperatorView(operators, operators.gradient, :gradient)
end

"""Return the odd-vector derivative action used for origin regularity."""
function odd_derivative_operator(operators::SphericalSBPOperators)
    return SphericalOperatorView(operators, operators.odd_derivative, :odd_derivative)
end

"""Return the odd-vector to even-scalar spherical divergence action."""
function divergence_operator(operators::SphericalSBPOperators)
    return SphericalOperatorView(operators, operators.divergence, :divergence)
end

function Base.show(io::IO, operators::SphericalSBPOperators)
    if get(io, :compact, false)
        print(io, "SphericalSBPOperators(order=", accuracy_order(operators), ")")
    else
        r = grid(operators)
        print(io, "Spherical SBP operators of order ", accuracy_order(operators),
              " on a ", operators.grid_kind, " grid in [", first(r), ", ", last(r), "]")
    end
end

const _SphericalRational = Rational{BigInt}

"""Solve an exact, possibly underdetermined, linear constraint system."""
function _spherical_solve(A::Matrix{_SphericalRational}, b::Vector{_SphericalRational})
    size(A, 1) == length(b) || throw(DimensionMismatch("constraint dimensions disagree"))
    augmented = hcat(copy(A), copy(b))
    pivots = Int[]
    pivot_row = 1
    for column in axes(A, 2)
        pivot = nothing
        for row in pivot_row:size(A, 1)
            if augmented[row, column] != 0
                pivot = row
                break
            end
        end
        isnothing(pivot) && continue
        augmented[pivot_row, :], augmented[pivot, :] = augmented[pivot, :],
                                                       augmented[pivot_row, :]
        augmented[pivot_row, :] ./= augmented[pivot_row, column]
        for row in axes(augmented, 1)
            row == pivot_row && continue
            augmented[row, :] .-= augmented[row, column] .* augmented[pivot_row, :]
        end
        push!(pivots, column)
        pivot_row += 1
        pivot_row > size(A, 1) && break
    end
    for row in axes(augmented, 1)
        all(iszero, augmented[row, axes(A, 2)]) && !iszero(augmented[row, end]) &&
            throw(ArgumentError("spherical SBP constraints are inconsistent"))
    end
    solution = zeros(_SphericalRational, size(A, 2))
    for (row, column) in enumerate(pivots)
        solution[column] = augmented[row, end]
    end
    return solution
end

"""Create folding maps for an exactly symmetric Cartesian grid."""
function _spherical_folding_maps(x::Vector{T}; origin::Bool) where {T <: Real}
    indices = findall(x -> origin ? x >= zero(T) : x > zero(T), x)
    isempty(indices) && throw(ArgumentError("could not identify the radial half grid"))
    r = x[indices]
    origin && first(r) != zero(T) &&
        throw(ArgumentError("the collocated grid must contain zero"))

    nfull = length(x)
    nhalf = length(r)
    R = sparse(collect(1:nhalf), indices, ones(T, nhalf), nhalf, nfull)
    even_i = Int[]
    even_j = Int[]
    odd_i = Int[]
    odd_j = Int[]
    odd_v = T[]
    for i in eachindex(x)
        j = findfirst(==(abs(x[i])), r)
        isnothing(j) && throw(ArgumentError("the Cartesian grid is not mirror symmetric"))
        push!(even_i, i)
        push!(even_j, j)
        if !iszero(x[i])
            push!(odd_i, i)
            push!(odd_j, j)
            push!(odd_v, x[i] < zero(T) ? -one(T) : one(T))
        end
    end
    Eeven = sparse(even_i, even_j, ones(T, nfull), nfull, nhalf)
    Eodd = sparse(odd_i, odd_j, odd_v, nfull, nhalf)
    return r, R, Eeven, Eodd
end

"""Build and fold the existing Cartesian SBP operator on a canonical grid."""
function _spherical_cartesian_data(source, accuracy_order::Int, N::Int;
                                   staggered::Bool, mode::AbstractExecutionMode)
    radius = staggered ? _SphericalRational(2 * N - 1, 2) : _SphericalRational(N)
    points = staggered ? 2 * N : 2 * N + 1
    cartesian = derivative_operator(source;
                                    derivative_order = 1,
                                    accuracy_order,
                                    xmin = -radius,
                                    xmax = radius,
                                    N = points,
                                    mode)
    x = collect(grid(cartesian))
    Gfull = sparse(Matrix(cartesian))
    Hfull = sparse(mass_matrix(cartesian))
    r, R, Eeven, Eodd = _spherical_folding_maps(x; origin = !staggered)
    G = sparse(R * Gfull * Eeven)
    Godd = sparse(R * Gfull * Eodd)
    Hhalf = sparse((_SphericalRational(1, 2)) * transpose(Eeven) * Hfull * Eeven)
    return (r = r, G = G, Godd = Godd, Hhalf = Hhalf, radius = radius)
end

"""Extract the diagonal of an exactly diagonal sparse matrix."""
function _spherical_diagonal(A::SparseMatrixCSC{T, Ti}) where {T, Ti}
    size(A, 1) == size(A, 2) || throw(DimensionMismatch("mass matrix must be square"))
    diagonal = zeros(T, size(A, 1))
    for (i, j, value) in zip(findnz(A)...)
        i == j || throw(ArgumentError("scalar spherical mass must be diagonal"))
        diagonal[i] = value
    end
    return diagonal
end

"""Recover divergence from the paired SBP identity for a nonsingular scalar norm."""
function _spherical_divergence(S::SparseMatrixCSC{T, Ti}, V::SparseMatrixCSC{T, Ti},
                               G::SparseMatrixCSC{T, Ti}, r::AbstractVector,
                               p::Int) where {T, Ti}
    diagonal = _spherical_diagonal(S)
    B = spzeros(T, length(r), length(r))
    B[end, end] = r[end]^p
    rhs = sparse(B - transpose(G) * V)
    I, J, values = findnz(rhs)
    Dvalues = similar(values)
    for k in eachindex(values)
        iszero(diagonal[I[k]]) && throw(ArgumentError("scalar spherical mass is singular"))
        Dvalues[k] = values[k] / diagonal[I[k]]
    end
    return sparse(I, J, Dvalues, size(rhs)...), B
end

"""Assemble the fourth-order collocated split norms from their exact constraints."""
function _spherical_sbp4(data, p::Int)
    r, G, Hhalf, radius = data.r, data.G, data.Hhalf, data.radius
    n = length(r)
    n > 5 ||
        throw(ArgumentError("fourth-order spherical operators require at least six radial nodes"))
    fixed = _spherical_diagonal(sparse(Hhalf * spdiagm(0 => r .^ p)))
    free_scalar = min(5, n)
    diagonal_slots = collect(2:min(5, n))
    offdiagonal_slots = [(2, 3), (3, 4), (4, 5)]
    offdiagonal_slots = filter(((i, j),) -> j <= n, offdiagonal_slots)
    unknowns = free_scalar + length(diagonal_slots) + length(offdiagonal_slots)
    scalar_index(i) = i
    diagonal_index = Dict(j => free_scalar + k for (k, j) in enumerate(diagonal_slots))
    offdiagonal_index = Dict(pair => free_scalar + length(diagonal_slots) + k
                             for (k, pair) in enumerate(offdiagonal_slots))
    equations = Vector{Vector{_SphericalRational}}()
    targets = _SphericalRational[]
    Gq = sparse(_SphericalRational.(G))
    rq = _SphericalRational.(r)

    for (degree, rows) in ((1, 1:(n - 4)), (3, 1:min(5, n)))
        u = rq .^ degree
        for i in rows
            row = zeros(_SphericalRational, unknowns)
            exact = (p + degree) * rq[i]^(degree - 1)
            constant = i == n ? rq[end]^p * u[end] : zero(_SphericalRational)
            if i <= free_scalar
                row[scalar_index(i)] -= exact
            else
                constant -= fixed[i] * exact
            end
            for j in (free_scalar + 1):n
                constant -= Gq[j, i] * fixed[j] * u[j]
            end
            for j in diagonal_slots
                row[diagonal_index[j]] -= Gq[j, i] * u[j]
            end
            for pair in offdiagonal_slots
                a, b = pair
                row[offdiagonal_index[pair]] -= Gq[a, i] * u[b] + Gq[b, i] * u[a]
            end
            constant -= Gq[1, i] * u[1]
            push!(equations, row)
            push!(targets, -constant)
        end
    end
    quadrature = zeros(_SphericalRational, unknowns)
    quadrature[1:free_scalar] .= one(_SphericalRational)
    push!(equations, quadrature)
    push!(targets, radius^(p + 1) / (p + 1) - sum(fixed[(free_scalar + 1):end]))
    solution = _spherical_solve(reduce(vcat, permutedims.(equations)), targets)

    scalar = _SphericalRational.(fixed)
    for i in 1:free_scalar
        scalar[i] = solution[scalar_index(i)]
    end
    vector = copy(scalar)
    vector[1] = one(_SphericalRational)
    for j in diagonal_slots
        vector[j] = solution[diagonal_index[j]]
    end
    V = spdiagm(0 => vector)
    for pair in offdiagonal_slots
        i, j = pair
        V[i, j] = solution[offdiagonal_index[pair]]
        V[j, i] = V[i, j]
    end
    S = spdiagm(0 => scalar)
    D, B = _spherical_divergence(S, V, Gq, rq, p)
    return (S = S, V = sparse(V), D = D, B = B)
end

"""Assemble the sixth-order collocated split norms from publication constraints."""
function _spherical_sbp6(data, p::Int)
    r, G, radius = data.r, data.G, data.radius
    n = length(r)
    n >= 17 ||
        throw(ArgumentError("sixth-order spherical operators require at least 17 radial nodes"))
    diagonal_slots = collect(2:min(10, n))
    offdiagonal_slots = [(2, 3), (2, 4), (3, 4), (3, 5), (4, 5), (5, 6), (6, 7)]
    scalar_index(i) = i
    diagonal_index = Dict(j => n + k for (k, j) in enumerate(diagonal_slots))
    offdiagonal_index = Dict(pair => n + length(diagonal_slots) + k
                             for (k, pair) in enumerate(offdiagonal_slots))
    unknowns = n + length(diagonal_slots) + length(offdiagonal_slots)
    Gq = sparse(_SphericalRational.(G))
    rq = _SphericalRational.(r)
    equations = Vector{Vector{_SphericalRational}}()
    targets = _SphericalRational[]
    interior = (maximum(diagonal_slots) + 1):n

    for (degree, rows) in ((1, 1:n), (3, 1:(n - 6)), (5, 1:min(3, n)))
        u = rq .^ degree
        for i in rows
            row = zeros(_SphericalRational, unknowns)
            row[scalar_index(i)] -= (p + degree) * rq[i]^(degree - 1)
            for j in interior
                row[scalar_index(j)] -= Gq[j, i] * u[j]
            end
            for j in diagonal_slots
                row[diagonal_index[j]] -= Gq[j, i] * u[j]
            end
            for pair in offdiagonal_slots
                a, b = pair
                row[offdiagonal_index[pair]] -= Gq[a, i] * u[b] + Gq[b, i] * u[a]
            end
            constant = (i == n ? rq[end]^p * u[end] : zero(_SphericalRational)) -
                       Gq[1, i] * u[1]
            push!(equations, row)
            push!(targets, -constant)
        end
    end
    quadrature = zeros(_SphericalRational, unknowns)
    quadrature[1:n] .= one(_SphericalRational)
    push!(equations, quadrature)
    push!(targets, radius^(p + 1) / (p + 1))
    solution = _spherical_solve(reduce(vcat, permutedims.(equations)), targets)

    scalar = solution[1:n]
    vector = copy(scalar)
    vector[1] = one(_SphericalRational)
    for j in diagonal_slots
        vector[j] = solution[diagonal_index[j]]
    end
    V = spdiagm(0 => vector)
    for pair in offdiagonal_slots
        i, j = pair
        V[i, j] = solution[offdiagonal_index[pair]]
        V[j, i] = V[i, j]
    end
    S = spdiagm(0 => scalar)
    D, B = _spherical_divergence(S, sparse(V), Gq, rq, p)
    return (S = S, V = sparse(V), D = D, B = B)
end

"""Scale canonical rational matrices and nodes to the requested physical radius."""
function _spherical_rescale(data, radius::Real, p::Int)
    T = radius isa AbstractFloat ? typeof(radius) : _SphericalRational
    ratio = convert(T, radius) / convert(T, data.radius)
    scale(A, factor) = begin
        I, J, values = findnz(A)
        sparse(I, J, convert.(T, values) .* factor, size(A)...)
    end
    return (r = convert.(T, data.r) .* ratio,
            S = scale(data.S, ratio^(p + 1)),
            V = scale(data.V, ratio^(p + 1)),
            B = scale(data.B, ratio^p),
            G = scale(data.G, inv(ratio)),
            Godd = scale(data.Godd, inv(ratio)),
            D = scale(data.D, inv(ratio)))
end

"""Construct a native staggered spherical SBP pair on the canonical half grid."""
function _spherical_staggered(source, accuracy_order::Int, N::Int, p::Int,
                              mode::AbstractExecutionMode)
    cartesian = _spherical_cartesian_data(source, accuracy_order, N;
                                          staggered = true, mode)
    H = sparse(cartesian.Hhalf * spdiagm(0 => cartesian.r .^ p))
    D, B = _spherical_divergence(H, H, cartesian.G, cartesian.r, p)
    return (r = cartesian.r, S = H, V = H, B = B,
            G = cartesian.G, Godd = cartesian.Godd, D = D,
            radius = cartesian.radius)
end

"""Construct a native collocated spherical SBP pair by exact closure solves."""
function _spherical_collocated(source, accuracy_order::Int, N::Int, p::Int,
                               mode::AbstractExecutionMode)
    cartesian = _spherical_cartesian_data(source, accuracy_order, N;
                                          staggered = false, mode)
    closure = accuracy_order == 4 ? _spherical_sbp4(cartesian, p) :
              _spherical_sbp6(cartesian, p)
    return (r = cartesian.r, S = closure.S, V = closure.V, B = closure.B,
            G = cartesian.G, Godd = cartesian.Godd, D = closure.D,
            radius = cartesian.radius)
end

"""
    spherical_operators(VretinarisSchnetter2026(); accuracy_order, N, R,
                        p=2, grid=:collocated,
                        cartesian_source=MattssonNordström2004(), mode=FastMode())

Construct fourth- or sixth-order spherical SBP operators. `N` is the
number of radial subintervals on a collocated grid and the number of radial
nodes on a staggered grid. The supported radial metric exponents are `p = 1`
and `p = 2`. For `p = 0`, use [`derivative_operator`](@ref) to construct the
usual Cartesian SBP operators instead.
"""
function spherical_operators(::VretinarisSchnetter2026;
                             accuracy_order::Integer,
                             N::Integer,
                             R::Real,
                             p::Integer = 2,
                             grid::Symbol = :collocated,
                             cartesian_source::SourceOfCoefficients = MattssonNordström2004(),
                             mode::AbstractExecutionMode = FastMode())
    @argcheck accuracy_order in (4, 6) ArgumentError("only fourth- and sixth-order spherical operators are implemented")
    @argcheck N > 0 ArgumentError("`N` must be positive")
    @argcheck R > zero(R) ArgumentError("`R` must be positive")
    p == 0 &&
        throw(ArgumentError("`p=0` describes the usual Cartesian operators; use `derivative_operator` instead"))
    @argcheck p in (1, 2) ArgumentError("only `p=1` and `p=2` spherical operators are implemented")
    @argcheck grid in (:collocated, :staggered) ArgumentError("`grid` must be :collocated or :staggered")

    canonical = grid === :collocated ?
                _spherical_collocated(cartesian_source, Int(accuracy_order), Int(N), Int(p),
                                      mode) :
                _spherical_staggered(cartesian_source, Int(accuracy_order), Int(N), Int(p),
                                     mode)
    construction = _spherical_rescale(canonical, R, Int(p))
    return SphericalSBPOperators(construction.r,
                                 construction.S, construction.V, construction.B,
                                 construction.G, construction.Godd, construction.D,
                                 Int(accuracy_order), Int(p), grid,
                                 VretinarisSchnetter2026(), cartesian_source)
end
