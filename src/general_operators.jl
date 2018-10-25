
derivative_order(coefficients::AbstractDerivativeCoefficients) = coefficients.derivative_order
accuracy_order(coefficients::AbstractDerivativeCoefficients) = coefficients.accuracy_order
Base.eltype(coefficients::AbstractDerivativeCoefficients{T}) where {T} = T
LinearAlgebra.issymmetric(coefficients::AbstractDerivativeCoefficients) = coefficients.symmetric


derivative_order(D::AbstractDerivativeOperator) = derivative_order(D.coefficients)
accuracy_order(D::AbstractDerivativeOperator) = accuracy_order(D.coefficients)
Base.eltype(D::AbstractDerivativeOperator{T}) where {T} = T
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
    mul!(dest, D, u, one(eltype(dest)))
end

@noinline function *(D::AbstractDerivativeOperator, u)
    @boundscheck begin
        @argcheck size(D,1) == size(D,2) == length(u) DimensionMismatch
    end
    T = promote_type(eltype(D), eltype(u))
    dest = similar(u, T); fill!(dest, zero(T))
    @inbounds mul!(dest, D, u)
    dest
end


function Base.Matrix(D::AbstractDerivativeOperator{T}) where {T}
    v = Array{T}(undef, size(D, 2)...)
    fill!(v, T(0))
    A = Array{T}(undef, size(D)...)
    for i in 1:size(D,2)
        v[i] = T(1)
        mul!(view(A,:,i), D, v)
        v[i] = T(0)
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
        v[i] = T(1)
        mul!(dest, D, v)
        js = findall(!iszero, dest)
        colptr[i] = length(nzval)+1
        if length(js) > 0
            append!(rowind, js)
            append!(nzval, dest[js])
        end
        v[i] = T(0)
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
    uval = Array{typeof(u((xmin+xmax)/2))}(undef, length(x))
    compute_coefficients!(uval, u, D)
end

"""
    compute_coefficients!(uval::AbstractVector, u, D::AbstractDerivativeOperator)

Compute the nodal values of the function `u` at the grid associated to the
derivative operator `D` and stores the result in `uval`.
"""
function compute_coefficients!(uval::AbstractVector, u, D::AbstractDerivativeOperator)
    uval .= u.(grid(D))
end


"""
    compute_coefficients(u, D::AbstractPeriodicDerivativeOperator)

Compute the nodal values of the function `u` at the grid associated to the
derivative operator `D`.
"""
function compute_coefficients(u, D::AbstractPeriodicDerivativeOperator)
    x = D.grid_compute
    xmin = first(x)
    xmax = last(x)
    uval = Array{typeof(u((xmin+xmax)/2))}(undef, length(x))
    compute_coefficients!(uval, u, D)
end

"""
    compute_coefficients!(uval::AbstractVector, u, D::AbstractPeriodicDerivativeOperator)

Compute the nodal values of the function `u` at the grid associated to the
derivative operator `D` and stores the result in `uval`.
"""
function compute_coefficients!(uval::AbstractVector, u, D::AbstractPeriodicDerivativeOperator)
    uval .= u.(D.grid_compute)
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
    xplot = Array{eltype(x)}(undef, length(x))
    uplot = Array{eltype(u)}(undef, length(x))

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
    @argcheck length(uplot) == length(xplot)
    @argcheck length(uplot) == length(grid(D))

    xplot .= grid(D)
    uplot .= u

    xplot, uplot
end


"""
    evaluate_coefficients(u, D::AbstractPeriodicDerivativeOperator)

Evaluates the nodal coefficients `u` at a grid including both endpoints
associated to the derivative periodic operator `D`.
Returns `xplot, uplot`, where `xplot` contains the equally spaced nodes and
`uplot` the corresponding values of `u`.
"""
function evaluate_coefficients(u, D::AbstractPeriodicDerivativeOperator)
    x = D.grid_evaluate
    xplot = Array{eltype(x)}(undef, length(x))
    uplot = Array{eltype(u)}(undef, length(x))

    evaluate_coefficients!(xplot, uplot, u, D)
end

"""
    evaluate_coefficients!(xplot, uplot, u, D::AbstractPeriodicDerivativeOperator)

Evaluates the nodal coefficients `u` at a grid including both endpoints
associated to the derivative periodic operator `D` and stores the result in
`xplot, uplot`.
Returns `xplot, uplot`, where `xplot` contains the equally spaced nodes and
`uplot` the corresponding values of `u`.
"""
function evaluate_coefficients!(xplot, uplot, u, D::AbstractPeriodicDerivativeOperator)
    @argcheck length(uplot) == length(xplot)
    @argcheck length(uplot) == length(D.grid_evaluate)

    xplot .= D.grid_evaluate
    uplot[1:end-1] = u
    uplot[end] = uplot[1]

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
        @argcheck all(i->grid(operators[i]) == grid(first(operators)), eachindex(operators)) ArgumentError
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
    print(io, "Sum of operators:\n")
    for D in sum.operators
        print(io, D)
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
