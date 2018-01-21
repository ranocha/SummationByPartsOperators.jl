
derivative_order(coefficients::AbstractDerivativeCoefficients) = coefficients.derivative_order
accuracy_order(coefficients::AbstractDerivativeCoefficients) = coefficients.accuracy_order
Base.eltype(coefficients::AbstractDerivativeCoefficients{T}) where {T} = T
Base.issymmetric(coefficients::AbstractDerivativeCoefficients) = coefficients.symmetric


derivative_order(D::AbstractDerivativeOperator) = derivative_order(D.coefficients)
accuracy_order(D::AbstractDerivativeOperator) = accuracy_order(D.coefficients)
Base.eltype(D::AbstractDerivativeOperator{T}) where {T} = T
Base.issymmetric(D::AbstractDerivativeOperator) = issymmetric(D.coefficients)
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


Base.@propagate_inbounds function Base.A_mul_B!(dest, D::AbstractDerivativeOperator, u)
    mul!(dest, D, u, one(eltype(dest)))
end

@noinline function *(D::AbstractDerivativeOperator, u)
    @boundscheck begin
        @argcheck size(D,1) == size(D,2) == length(u) DimensionMismatch
    end
    T = promote_type(eltype(D), eltype(u))
    dest = similar(u, T); fill!(dest, zero(T))
    @inbounds A_mul_B!(dest, D, u)
    dest
end


function Base.full(D::AbstractDerivativeOperator{T}) where {T}
    v = Array{T}(size(D, 2))
    fill!(v, T(0))
    A = Array{T}(size(D)...)
    for i in 1:size(D,2)
        v[i] = T(1)
        A_mul_B!(view(A,:,i), D, v)
        v[i] = T(0)
    end
    A
end


function Base.sparse(D::AbstractDerivativeOperator{T}) where {T}
    M, N = size(D)
    rowind = Vector{Int}()
    nzval = Vector{T}()
    colptr = Vector{Int}(N+1)
    v = Array{T}(N)
    fill!(v, T(0))
    dest = Array{T}(M)

    for i = 1:N
        v[i] = T(1)
        A_mul_B!(dest, D, v)
        js = find(dest)
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
    uval = Array{typeof(u((xmin+xmax)/2))}(length(x))
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
    uval = Array{typeof(u((xmin+xmax)/2))}(length(x))
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
    xplot = Array{eltype(x)}(length(x))
    uplot = Array{eltype(u)}(length(x))

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
    xplot = Array{eltype(x)}(length(x))
    uplot = Array{eltype(u)}(length(x))

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
