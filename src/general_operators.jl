
# general interface functions
"""
    grid(D)

Return the grid associated to a derivative operator `D`.
"""
function grid end

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
    print(io, nameof(typeof(D)), "(derivative_order=", derivative_order(D),
              ", accuracy_order=", accuracy_order(D), ")")
end



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


derivative_order(coefficients::AbstractDerivativeCoefficients) = coefficients.derivative_order
accuracy_order(coefficients::AbstractDerivativeCoefficients) = coefficients.accuracy_order
Base.eltype(coefficients::AbstractDerivativeCoefficients{T}) where {T} = T
LinearAlgebra.issymmetric(coefficients::AbstractDerivativeCoefficients) = coefficients.symmetric


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


Base.@propagate_inbounds function mul!(dest, D::AbstractDerivativeOperator, u)
    mul!(dest, D, u, one(recursive_bottom_eltype(dest)))
end

function Base.:*(D::AbstractDerivativeOperator, u)
    @boundscheck begin
        @argcheck size(D,1) == size(D,2) == length(u) DimensionMismatch
    end
    T = typeof(one(eltype(D)) * first(u))
    dest = similar(u, T); fill!(dest, zero(eltype(dest)))
    @inbounds mul!(dest, D, u)
    dest
end


function Base.Matrix(D::AbstractDerivativeOperator{T}) where {T}
    v = Array{T}(undef, size(D, 2)...)
    fill!(v, zero(eltype(v)))
    A = Array{T}(undef, size(D)...)
    for i in 1:size(D,2)
        v[i] = one(T)
        # Using a view here can cause problems with FFT based operators.
        # Since this part is not performance critical, we can also just use copies.
        # mul!(view(A,:,i), D, v)
        A[:,i] .= D * v
        v[i] = zero(T)
    end
    A
end


function SparseArrays.sparse(D::AbstractDerivativeOperator{T}) where {T}
    M, N = size(D)
    rowind = Vector{Int}()
    nzval = Vector{T}()
    colptr = Vector{Int}(undef, N+1)
    v = fill(zero(T), N)
    dest = Array{T}(undef, M)

    for i = 1:N
        v[i] = one(T)
        mul!(dest, D, v)
        js = findall(!iszero, dest)
        colptr[i] = length(nzval)+1
        if length(js) > 0
            append!(rowind, js)
            append!(nzval, dest[js])
        end
        v[i] = zero(T)
    end
    colptr[N+1] = length(nzval)+1

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
    uval = Array{typeof(u((xmin+xmax)/2))}(undef, size(x)...)
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
    integrate(func, u, D::AbstractPeriodicDerivativeOperator)

Map the function `func` to the coefficients `u` and integrate with respect to
the quadrature rule associated with the derivative operator `D`.
"""
function integrate(func, u::AbstractVector, D::AbstractPeriodicDerivativeOperator)
    @boundscheck begin
        length(u) == length(grid(D))
    end
    @unpack Δx = D

    res = sum(func, u)

    Δx * res
end



"""
    SumOfDerivativeOperators

Sum several derivative operators lazily.
"""
struct SumOfDerivativeOperators{T,N,Operators<:Tuple{Vararg{AbstractDerivativeOperator{T},N}}} <: AbstractDerivativeOperator{T}
    operators::Operators

    function SumOfDerivativeOperators(operators::Operators) where {T,N,Operators<:Tuple{Vararg{AbstractDerivativeOperator{T},N}}}
        @argcheck all(i->size(operators[i]) == size(first(operators)), eachindex(operators)) DimensionMismatch
        @argcheck all(i->grid(operators[i]) ≈ grid(first(operators)), eachindex(operators)) ArgumentError
        new{T,N,Operators}(operators)
    end
end

SumOfDerivativeOperators(ops...) = SumOfDerivativeOperators(ops)

Base.size(sum::SumOfDerivativeOperators) = size(first(sum.operators))
Base.size(sum::SumOfDerivativeOperators, i::Int) = size(first(sum.operators), i)
function Base.length(::Type{SumOfDerivativeOperators{T,N,Operators}}) where {T,N,Operators}
    N
end
grid(sum::SumOfDerivativeOperators) = grid(first(sum.operators))

function Base.:+(D1::AbstractDerivativeOperator, D2::AbstractDerivativeOperator)
    SumOfDerivativeOperators(D1, D2)
end

function Base.:+(sum::SumOfDerivativeOperators, D::AbstractDerivativeOperator)
    SumOfDerivativeOperators(sum.operators..., D)
end

function Base.:+(D::AbstractDerivativeOperator, sum::SumOfDerivativeOperators)
    SumOfDerivativeOperators(D, sum.operators...)
end

function Base.:+(sum1::SumOfDerivativeOperators, sum2::SumOfDerivativeOperators)
    SumOfDerivativeOperators(sum1.operators..., sum2.operators...)
end

function Base.show(io::IO, sum::SumOfDerivativeOperators)
    print(io, "Sum of ", length(sum.operators), " operators")
    if get(io, :compact, false) == false
        print(io, ":")
        for D in sum.operators
            print(io, "\n", D)
        end
    end
end

@unroll function mul!(dest::AbstractVector, sum::SumOfDerivativeOperators, u::AbstractVector,
                      α, β)
    @unpack operators = sum
    @boundscheck begin
        @argcheck size(first(operators), 2) == length(u) DimensionMismatch
        @argcheck size(first(operators), 1) == length(dest) DimensionMismatch
    end

    @inbounds mul!(dest, operators[1], u, α, β)
    @unroll for i in 1:length(sum)
        if i != 1
            @inbounds mul!(dest, operators[i], u, α, one(β))
        end
    end

    nothing
end

@unroll function mul!(dest::AbstractVector, sum::SumOfDerivativeOperators, u::AbstractVector,
                      α)
    @unpack operators = sum
    @boundscheck begin
        @argcheck size(first(operators), 2) == length(u) DimensionMismatch
        @argcheck size(first(operators), 1) == length(dest) DimensionMismatch
    end

    @inbounds mul!(dest, operators[1], u, α)
    @unroll for i in 1:length(sum)
        if i != 1
            @inbounds mul!(dest, operators[i], u, α, one(α))
        end
    end

    nothing
end
