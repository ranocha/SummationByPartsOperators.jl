
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

# LoopVectorization.jl can only vectorize loops over native numbers. For other
# element types, `@turbo` falls back to a scalar loop using `Base.FastMath`
# operations. This fallback can be much slower than plain arithmetic since,
# e.g., `Base.FastMath.mul_fast(α::Float64, u::ForwardDiff.Dual)` promotes `α`
# to a `ForwardDiff.Dual` and evaluates the full product rule.
# Fortunately, vectors of composite element types such as `StaticVector`s or
# `ForwardDiff.Dual`s can be reinterpreted as arrays of native numbers.
# Applying the stencils to such an array allows `@turbo` to vectorize across
# the components, which is usually even faster than a plain scalar loop.
# See https://github.com/ranocha/SummationByPartsOperators.jl/issues/421
#
# All of this is decided by the singleton values below. They must be singletons
# since the kernels are `@generated` functions: their generators run in the
# world age of their own definition and would thus not see methods of
# `reinterpreted_components` added later by package extensions. Dispatching on
# the layout moves the decision to ordinary inference instead.

"""
    SummationByPartsOperators.Components{N, V}()

Values of the element type described by this trait consist of `N` components of
the native number type `V`. See [`reinterpreted_components`](@ref).
"""
struct Components{N, V} end

"""
    SummationByPartsOperators.ArrayLayout{N, Vdest, Vu}()

The kernels access `dest` and `u` as arrays of the native number types `Vdest`
and `Vu` with `N` components per value. See [`reinterpreted_layout`](@ref).
"""
struct ArrayLayout{N, Vdest, Vu} end

"""
    SummationByPartsOperators.reinterpreted_components(::Type{T})

Return `nothing` if values of type `T` must be treated as opaque scalars and
`Components{N, V}()` if `reinterpret(reshape, V, ::AbstractVector{T})` yields an
`N × length` matrix of native numbers `V`.
"""
reinterpreted_components(::Type) = nothing

function reinterpreted_components(::Type{S}) where {N, T, S <: StaticVector{N, T}}
    reinterpreted_components(S, Val{N}(), T)
end

function reinterpreted_components(::Type{C}) where {T, C <: Complex{T}}
    reinterpreted_components(C, Val{2}(), T)
end

# Values of type `T` consist of `N` components of type `V`, which may be
# composite themselves.
@inline function reinterpreted_components(::Type{T}, ::Val{N}, ::Type{V}) where {T, N, V}
    checked_components(T, combine_components(Val{N}(), V, reinterpreted_components(V)))
end

@inline function combine_components(::Val{N}, ::Type{V}, ::Nothing) where {N, V}
    native_components(Val{N}(), V, Val(LoopVectorization.check_type(V)))
end
function combine_components(::Val{N}, ::Type{V}, ::Components{M, W}) where {N, V, M, W}
    Components{N * M, W}()
end
native_components(::Val{N}, ::Type{V}, ::Val{true}) where {N, V} = Components{N, V}()
native_components(::Val{N}, ::Type{V}, ::Val{false}) where {N, V} = nothing

# `reinterpret(reshape, ...)` requires an exact match of the memory layout
checked_components(::Type, ::Nothing) = nothing
@inline function checked_components(::Type{T}, components::Components{N, V}) where {T, N, V}
    matching_components(components, Val(isbitstype(T) && sizeof(T) == N * sizeof(V)))
end
matching_components(components::Components, ::Val{true}) = components
matching_components(::Components, ::Val{false}) = nothing

# Factors multiplying the values must be plain numbers acting on each component
# in the same way - multiplying by, e.g., a `ForwardDiff.Dual` mixes them.
# They must also be native numbers themselves since they end up inside the
# vectorized loops over the components.
is_plain_factor(factor) = false
function is_plain_factor(factor::Real)
    reinterpreted_components(typeof(factor)) === nothing &&
        LoopVectorization.check_type(typeof(factor))
end

all_plain_factors(::Tuple{}) = true
@inline function all_plain_factors(factors::Tuple)
    is_plain_factor(first(factors)) && all_plain_factors(Base.tail(factors))
end

"""
    SummationByPartsOperators.reinterpreted_layout(mode, dest, u, factors...)

Return the [`ArrayLayout`](@ref) the kernels should use to compute
`α * D * u [+ β * dest]` with the given scaling `factors`, or `nothing` if
`dest` and `u` have to be used as they are.
"""
@inline function reinterpreted_layout(mode::AbstractExecutionMode, dest::AbstractVector,
                                      u::AbstractVector, factors...)
    array_layout(reinterpreted_components(eltype(dest)),
                 reinterpreted_components(eltype(u)),
                 Val(all_plain_factors(factors)))
end

# `SafeMode` relies only on basic functionality of Julia. Its plain loops over
# the values are also usually faster than loops over the reinterpreted
# components, which pay off only in combination with `@turbo`.
function reinterpreted_layout(::SafeMode, dest::AbstractVector, u::AbstractVector,
                              factors...)
    nothing
end

array_layout(dest_components, u_components, ::Val{false}) = nothing
array_layout(::Nothing, u_components, ::Val{true}) = nothing
array_layout(dest_components, ::Nothing, ::Val{true}) = nothing
array_layout(::Nothing, ::Nothing, ::Val{true}) = nothing
array_layout(::Components, ::Components, ::Val{true}) = nothing
function array_layout(::Components{N, Vdest}, ::Components{N, Vu},
                      ::Val{true}) where {N, Vdest, Vu}
    ArrayLayout{N, Vdest, Vu}()
end

# A single component per value is best handled as a plain vector of native numbers
function reinterpret_vector_expression(array::Symbol, ::Type{V}) where {V}
    :(reinterpret($V, $array))
end
function reinterpret_matrix_expression(array::Symbol, ::Type{V}) where {V}
    :(reinterpret(reshape, $V, $array))
end

# Stencil coefficients and scaling factors that are not native numbers cannot be
# handled by LoopVectorization.jl. Since the reinterpreted layouts pay off only
# in combination with `@turbo`, the kernels use the values as they are instead.
function layout_for_factors(layout::Type, factor_types)
    all(LoopVectorization.check_type, factor_types) ? layout : Nothing
end

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

"""
    isperiodic(D)

Return `true` if the derivative operator `D` is periodic and `false` otherwise.
"""
isperiodic(D::AbstractDerivativeOperator) = false
isperiodic(D::AbstractPeriodicDerivativeOperator) = true

xmin(D::AbstractDerivativeOperator) = first(grid(D))
xmax(D::AbstractDerivativeOperator) = last(grid(D))

Base.@propagate_inbounds function mul!(dest, D::AbstractDerivativeOperator, u)
    mul!(dest, D, u, one(scaling_eltype(dest)))
end

# The scaling factor used by `mul!(dest, D, u)` should be a plain real number
# instead of, e.g., a `ForwardDiff.Dual` or a `Complex`. Otherwise, multiplying
# by it is unnecessarily expensive and prevents the vectorized kernels from
# being used. A real one is a multiplicative identity just as well. Element
# types we do not know anything about are used as they are since they need not
# support `real`.
@inline function scaling_eltype(dest)
    T = recursive_bottom_eltype(dest)
    native_eltype(T, reinterpreted_components(T))
end
native_eltype(::Type{T}, ::Nothing) where {T} = T
native_eltype(::Type{Complex{T}}, ::Nothing) where {T} = T
native_eltype(::Type{T}, ::Components{N, V}) where {T, N, V} = V

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
function integrate(func::Func, u::AbstractVector,
                   D::AbstractPeriodicDerivativeOperator) where {Func}
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
restrict_interior(u, D::AbstractNonperiodicDerivativeOperator) = u[(begin + 1):(end - 1)]

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
