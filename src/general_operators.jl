
# general interface functions
"""
    grid(D)

Return the grid associated to a derivative operator `D`.
"""
function grid end

"""
    xmin(D)

Return the left boundary `xmin` of the domain specified when constructing the
derivative operator `D`. Note that this might be different from the leftmost
node of the [`grid`](@ref) of `D` when not all boundary nodes are included,
e.g., for periodic derivative operators.
"""
function xmin end

"""
    xmax(D)

Return the right boundary `xmax` of the domain specified when constructing the
derivative operator `D`. Note that this might be different from the rightmost
node of the [`grid`](@ref) of `D` when not all boundary nodes are included,
e.g., for periodic derivative operators.
"""
function xmax end

"""
    accuracy_order(D)

Return the order of accuracy of a derivative operator `D`. For SBP finite difference
operators, this refers to the interior order of accuracy.
"""
function accuracy_order end

"""
    derivative_order(D)

Return the order of the derivative associated to the derivative operator `D`.
For example, it will return `1` for a first-derivative SBP operator.
"""
function derivative_order end

"""
    source_of_coefficients(D)

Return the source of coefficients of the derivative operator `D`. If you use the
operator `D` for your research, please cite this source in addition to
[`SummationByPartsOperators`](@ref).
"""
source_of_coefficients(D) = SummationByPartsOperators

"""
    left_boundary_weight(D)

Return the left-boundary weight of the (diagonal) mass matrix `M` associated to
the derivative operator `D`.
"""
function left_boundary_weight end

"""
    right_boundary_weight(D)

Return the left-boundary weight of the (diagonal) mass matrix `M` associated to
the derivative operator `D`.
"""
function right_boundary_weight end

function Base.summary(io::IO, D::AbstractDerivativeOperator)
    print(io, nameof(typeof(D)), "(derivative:", derivative_order(D),
          ", accuracy:", accuracy_order(D), ")")
end

"""
    scale_by_mass_matrix!(u, D)

Scale the vector `u` by the mass matrix associated to the
derivative operator `D`.
"""
function scale_by_mass_matrix! end

"""
    scale_by_inverse_mass_matrix!(u, D)

Scale the vector `u` by the inverse of the mass matrix associated to the
derivative operator `D`.
"""
function scale_by_inverse_mass_matrix! end

abstract type AbstractExecutionMode end
"""
    SafeMode()

A safe execution mode relying only on basic functionality of Julia.
"""
struct SafeMode <: AbstractExecutionMode end
"""
    FastMode()

A (probably) faster execution mode that might depend on packages such as
[LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl).
"""
struct FastMode <: AbstractExecutionMode end
"""
    ThreadedMode()

An execution mode using multiple threads and possibly further optimizations, cf.
[`FastMode`](@ref).
"""
struct ThreadedMode <: AbstractExecutionMode end

# TODO: deprecated in v0.5
_parallel_to_mode(::Val{:threads}) = ThreadedMode()
_parallel_to_mode(::Val{:serial}) = FastMode()

function derivative_order(coefficients::AbstractDerivativeCoefficients)
    coefficients.derivative_order
end
accuracy_order(coefficients::AbstractDerivativeCoefficients) = coefficients.accuracy_order
Base.eltype(coefficients::AbstractDerivativeCoefficients{T}) where {T} = T
function LinearAlgebra.issymmetric(coefficients::AbstractDerivativeCoefficients)
    coefficients.symmetric
end

derivative_order(D::AbstractDerivativeOperator) = derivative_order(D.coefficients)
accuracy_order(D::AbstractDerivativeOperator) = accuracy_order(D.coefficients)
Base.eltype(D::AbstractDerivativeOperator{T}) where {T} = T
Base.real(D::AbstractDerivativeOperator) = real(eltype(D))
LinearAlgebra.issymmetric(D::AbstractDerivativeOperator) = issymmetric(D.coefficients)
function Base.size(D::AbstractDerivativeOperator)
    N = length(grid(D))
    (N, N)
end
function Base.size(D::AbstractDerivativeOperator, i::Int)
    if i <= 0
        error("arraysize: dimension out of range")
    elseif i <= 2
        length(grid(D))
    else
        1
    end
end

@inline grid(D::AbstractDerivativeOperator) = D.grid
@inline grid(D::AbstractPeriodicDerivativeOperator) = D.grid_compute

xmin(D::AbstractDerivativeOperator) = first(grid(D))
xmax(D::AbstractDerivativeOperator) = last(grid(D))

Base.@propagate_inbounds function mul!(dest, D::AbstractDerivativeOperator, u)
    mul!(dest, D, u, one(recursive_bottom_eltype(dest)))
end

function Base.:*(D::AbstractDerivativeOperator, u)
    @boundscheck begin
        @argcheck size(D, 1)==size(D, 2)==length(u) DimensionMismatch
    end
    T = typeof(one(eltype(D)) * first(u))
    dest = similar(u, T)
    fill!(dest, zero(eltype(dest)))
    @inbounds mul!(dest, D, u)
    dest
end

function Base.Matrix(D::AbstractDerivativeOperator{T}) where {T}
    v = Array{T}(undef, size(D, 2)...)
    fill!(v, zero(eltype(v)))
    A = Array{T}(undef, size(D)...)
    for i in 1:size(D, 2)
        v[i] = one(T)
        # Using a view here can cause problems with FFT based operators.
        # Since this part is not performance critical, we can also just use copies.
        # mul!(view(A,:,i), D, v)
        A[:, i] .= D * v
        v[i] = zero(T)
    end
    A
end

function SparseArrays.sparse(D::AbstractDerivativeOperator{T}) where {T}
    M, N = size(D)
    rowind = Vector{Int}()
    nzval = Vector{T}()
    colptr = Vector{Int}(undef, N + 1)
    v = fill(zero(T), N)
    dest = Array{T}(undef, M)

    for i in 1:N
        v[i] = one(T)
        mul!(dest, D, v)
        js = findall(!iszero, dest)
        colptr[i] = length(nzval) + 1
        if length(js) > 0
            append!(rowind, js)
            append!(nzval, dest[js])
        end
        v[i] = zero(T)
    end
    colptr[N + 1] = length(nzval) + 1

    return SparseMatrixCSC(M, N, colptr, rowind, nzval)
end

"""
    compute_coefficients(u, D::AbstractDerivativeOperator)

Compute the nodal values of the function `u` at the grid associated to the
derivative operator `D`.
"""
function compute_coefficients(u, D::AbstractDerivativeOperator)
    x = grid(D)
    xmin = first(x)
    xmax = last(x)
    uval = Array{typeof(u((xmin + xmax) / 2))}(undef, size(x)...)
    compute_coefficients!(uval, u, D)
end

"""
    compute_coefficients!(uval, u, D::AbstractDerivativeOperator)

Compute the nodal values of the function `u` at the grid associated to the
derivative operator `D` and stores the result in `uval`.
"""
function compute_coefficients!(uval, u, D::AbstractDerivativeOperator)
    uval .= u.(grid(D))
end

"""
    evaluate_coefficients(u, D::AbstractDerivativeOperator)

Evaluates the nodal coefficients `u` at a grid associated to the derivative
operator `D`.
Returns `xplot, uplot`, where `xplot` contains the nodes and `uplot` the
corresponding values of `u`.
"""
function evaluate_coefficients(u, D::AbstractDerivativeOperator)
    x = grid(D)
    xplot = Array{eltype(x)}(undef, size(x)...)
    uplot = Array{eltype(u)}(undef, size(x)...)

    evaluate_coefficients!(xplot, uplot, u, D)
end

"""
    evaluate_coefficients!(xplot, uplot, u, D::AbstractDerivativeOperator)

Evaluates the nodal coefficients `u` at a grid associated to the derivative
operator `D` and stores the result in `xplot, uplot`.
Returns `xplot, uplot`, where `xplot` contains the nodes and `uplot` the
corresponding values of `u`.
"""
function evaluate_coefficients!(xplot, uplot, u, D::AbstractDerivativeOperator)
    @argcheck size(uplot) == size(xplot)
    @argcheck size(uplot) == size(grid(D))

    xplot .= grid(D)
    uplot .= u

    xplot, uplot
end

"""
    integrate([func = identity,] u, D::AbstractPeriodicDerivativeOperator)

Map the function `func` to the coefficients `u` and integrate with respect to
the quadrature rule associated with the derivative operator `D`.
"""
function integrate(func, u::AbstractVector, D::AbstractPeriodicDerivativeOperator)
    @boundscheck begin
        length(u) == length(grid(D)) ||
            throw(DimensionMismatch("sizes of input vector and operator do not match"))
    end
    @unpack Δx = D

    res = sum(func, u)

    Δx * res
end

"""
    integrate_boundary([func = identity,] u, D::AbstractDerivativeOperator)

Map the function `func` to the coefficients `u` and integrate along the boundary. For classical 1D
operators this is `func(u[end]) - func(u[begin])`. For periodic 1D operators this is zero.
"""
function integrate_boundary(func, u, D::AbstractNonperiodicDerivativeOperator)
    return func(u[end]) - func(u[begin])
end

function integrate_boundary(func, u, D::AbstractPeriodicDerivativeOperator)
    return zero(func(u[begin]))
end

integrate_boundary(u, D) = integrate_boundary(identity, u, D)

"""
    restrict_boundary(u, D::AbstractDerivativeOperator)

Restrict the coefficients `u` to the boundary nodes of the derivative operator `D`.
"""
restrict_boundary(u, D::AbstractNonperiodicDerivativeOperator) = u[[begin, end]]

restrict_boundary(u, D::AbstractPeriodicDerivativeOperator) = eltype(u)[]

"""
    restrict_interior(u, D::AbstractDerivativeOperator)

Restrict the coefficients `u` to the interior nodes of the derivative operator `D`.
"""
restrict_interior(u, D::AbstractNonperiodicDerivativeOperator) = u[2:end-1]

restrict_interior(u, D::AbstractPeriodicDerivativeOperator) = u

"""
    mass_matrix_boundary(D::AbstractDerivativeOperator)

Construct the mass matrix at the boundary of a derivative operator `D`. For classical 1D
non-periodic operators, this is the matrix `Diagonal([-1, 0, ..., 0, 1])`. For periodic 1D
operators this is the zero matrix.
"""
function mass_matrix_boundary(D::AbstractNonperiodicDerivativeOperator)
    T = eltype(D)
    b = zeros(T, length(grid(D)))
    b[begin] = -one(T)
    b[end] = one(T)
    return Diagonal(b)
end

function mass_matrix_boundary(D::AbstractPeriodicDerivativeOperator)
    T = eltype(D)
    return zero(T) * I
end

"""
    LinearlyCombinedDerivativeOperators

Form linear combinations of several derivative operators lazily.
"""
@auto_hash_equals struct LinearlyCombinedDerivativeOperators{T, N,
                                                             Operators <:
                                                             Tuple{Vararg{AbstractDerivativeOperator{T},
                                                                          N}},
                                                             Coefficients <:
                                                             Tuple{Vararg{T, N}}} <:
                         AbstractDerivativeOperator{T}
    operators::Operators
    coefficients::Coefficients

    function LinearlyCombinedDerivativeOperators{T, N, Operators, Coefficients}(operators::Operators,
                                                                                coefficients::Coefficients) where {
                                                                                                                   T,
                                                                                                                   N,
                                                                                                                   Operators <:
                                                                                                                   Tuple{Vararg{AbstractDerivativeOperator{T},
                                                                                                                                N}},
                                                                                                                   Coefficients <:
                                                                                                                   Tuple{Vararg{T,
                                                                                                                                N}}
                                                                                                                   }
        @argcheck all(i -> size(operators[i]) == size(first(operators)),
                      eachindex(operators)) DimensionMismatch
        @argcheck all(i -> grid(operators[i]) ≈ grid(first(operators)),
                      eachindex(operators)) ArgumentError
        new{T, N, Operators, Coefficients}(operators, coefficients)
    end
end

function LinearlyCombinedDerivativeOperators(ops::NTuple{N, AbstractDerivativeOperator{T}},
                                             coefficients::NTuple{N, Number}) where {T, N}
    coefficients = map(c -> convert(T, c), coefficients)
    return LinearlyCombinedDerivativeOperators{T, N, typeof(ops), typeof(coefficients)}(ops,
                                                                                        coefficients)
end

function LinearlyCombinedDerivativeOperators(ops...)
    LinearlyCombinedDerivativeOperators(ops, ntuple(_ -> true, length(ops)))
end

# TODO: deprecated in v0.5.28
Base.@deprecate_binding SumOfDerivativeOperators LinearlyCombinedDerivativeOperators false

Base.size(combi::LinearlyCombinedDerivativeOperators) = size(first(combi.operators))
function Base.size(combi::LinearlyCombinedDerivativeOperators, i::Int)
    size(first(combi.operators), i)
end
function Base.length(::Type{LinearlyCombinedDerivativeOperators{T, N, Operators,
                                                                Coefficients}}) where {T, N,
                                                                                       Operators,
                                                                                       Coefficients
                                                                                       }
    N
end
grid(combi::LinearlyCombinedDerivativeOperators) = grid(first(combi.operators))

function Base.show(io::IO, combi::LinearlyCombinedDerivativeOperators)
    print(io, "Linear combination of ", length(combi.operators), " operators")
    if get(io, :compact, false) == false
        print(io, ":")
        for (D, c) in zip(combi.operators, combi.coefficients)
            print(io, "\n", D)
            print(io, "\nwith coefficient ", c)
        end
    end
end

function Base.:+(D1::AbstractDerivativeOperator, D2::AbstractDerivativeOperator)
    LinearlyCombinedDerivativeOperators(D1, D2)
end

function Base.:+(combi::LinearlyCombinedDerivativeOperators, D::AbstractDerivativeOperator)
    coefficients = (combi.coefficients..., one(eltype(D)))
    LinearlyCombinedDerivativeOperators((combi.operators..., D), coefficients)
end

function Base.:+(D::AbstractDerivativeOperator, combi::LinearlyCombinedDerivativeOperators)
    coefficients = (one(eltype(D)), combi.coefficients...)
    LinearlyCombinedDerivativeOperators((combi.operators..., D), coefficients)
end

function Base.:+(combi1::LinearlyCombinedDerivativeOperators,
                 combi2::LinearlyCombinedDerivativeOperators)
    LinearlyCombinedDerivativeOperators((combi1.operators..., combi2.operators...),
                                        (combi1.coefficients..., combi2.coefficients...))
end

function Base.:-(D1::AbstractDerivativeOperator, D2::AbstractDerivativeOperator)
    T = eltype(D1)
    LinearlyCombinedDerivativeOperators((D1, D2), (one(T), -one(T)))
end

function Base.:-(combi::LinearlyCombinedDerivativeOperators, D::AbstractDerivativeOperator)
    coefficients = (combi.coefficients..., -one(eltype(D)))
    LinearlyCombinedDerivativeOperators((combi.operators..., D), coefficients)
end

function Base.:-(D::AbstractDerivativeOperator, combi::LinearlyCombinedDerivativeOperators)
    coefficients = (one(eltype(D)), map(c -> -c, combi.coefficients)...)
    LinearlyCombinedDerivativeOperators((D, combi.operators...), coefficients)
end

function Base.:-(combi1::LinearlyCombinedDerivativeOperators,
                 combi2::LinearlyCombinedDerivativeOperators)
    coefficients = (combi1.coefficients..., map(c -> -c, combi2.coefficients)...)
    LinearlyCombinedDerivativeOperators((combi1.operators..., combi2.operators...),
                                        coefficients)
end

function Base.:+(D::AbstractDerivativeOperator)
    D
end

function Base.:+(combi::LinearlyCombinedDerivativeOperators)
    combi
end

function Base.:-(D::AbstractDerivativeOperator)
    LinearlyCombinedDerivativeOperators((D,), (-one(eltype(D)),))
end

function Base.:-(combi::LinearlyCombinedDerivativeOperators)
    coefficients = map(c -> -c, combi.coefficients)
    LinearlyCombinedDerivativeOperators(combi.operators, coefficients)
end

function Base.:*(c::Number, D::AbstractDerivativeOperator)
    LinearlyCombinedDerivativeOperators((D,), (c,))
end

function Base.:*(D::AbstractDerivativeOperator, c::Number)
    # TODO: Assume associativity
    c * D
end

function Base.:*(c::Number, combi::LinearlyCombinedDerivativeOperators)
    coefficients = map(ci -> c * ci, combi.coefficients)
    LinearlyCombinedDerivativeOperators(combi.operators, coefficients)
end

function Base.:*(combi::LinearlyCombinedDerivativeOperators, c::Number)
    # TODO: Assume associativity
    c * combi
end

function Base.:\(c::Number, D::AbstractDerivativeOperator)
    if eltype(D) <: AbstractFloat
        factor = inv(c)
    else
        factor = 1 // c
    end
    LinearlyCombinedDerivativeOperators((D,), (factor,))
end

function Base.:/(D::AbstractDerivativeOperator, c::Number)
    # TODO: Assume associativity
    c \ D
end

function Base.:\(c::Number, combi::LinearlyCombinedDerivativeOperators)
    coefficients = map(ci -> c \ ci, combi.coefficients)
    LinearlyCombinedDerivativeOperators(combi.operators, coefficients)
end

function Base.:/(combi::LinearlyCombinedDerivativeOperators, c::Number)
    # TODO: Assume associativity
    c \ combi
end

@unroll function mul!(dest::AbstractVector, combi::LinearlyCombinedDerivativeOperators,
                      u::AbstractVector,
                      α, β)
    @unpack operators, coefficients = combi
    @boundscheck begin
        @argcheck size(first(operators), 2)==length(u) DimensionMismatch
        @argcheck size(first(operators), 1)==length(dest) DimensionMismatch
    end

    @inbounds mul!(dest, operators[1], u, α * coefficients[1], β)
    @unroll for i in 1:length(combi)
        if i != 1
            @inbounds mul!(dest, operators[i], u, α * coefficients[i], one(β))
        end
    end

    nothing
end

@unroll function mul!(dest::AbstractVector, combi::LinearlyCombinedDerivativeOperators,
                      u::AbstractVector,
                      α)
    @unpack operators, coefficients = combi
    @boundscheck begin
        @argcheck size(first(operators), 2)==length(u) DimensionMismatch
        @argcheck size(first(operators), 1)==length(dest) DimensionMismatch
    end

    @inbounds mul!(dest, operators[1], u, α * coefficients[1])
    @unroll for i in 1:length(combi)
        if i != 1
            @inbounds mul!(dest, operators[i], u, α * coefficients[i], one(α))
        end
    end

    nothing
end
