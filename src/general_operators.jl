
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


Base.@propagate_inbounds function Base.A_mul_B!(dest, D::AbstractDerivativeOperator, u)
    mul!(dest, D, u, one(eltype(dest)))
end

@noinline function *(D::AbstractDerivativeOperator, u)
    T = promote_type(eltype(D), eltype(u))
    dest = similar(u, T); fill!(dest, zero(T))
    @inbounds A_mul_B!(dest, D, u)
    dest
end


function Base.full(D::AbstractDerivativeOperator{T}) where {T}
    v = zeros(T, size(D, 2))
    A = zeros(T, size(D)...)
    for i in 1:size(D,2)
        v .= zero(T)
        v[i] = one(T)
        A_mul_B!(view(A,:,i), D, v)
    end
    A
end


function Base.sparse(D::AbstractDerivativeOperator{T}) where {T}
    M, N = size(D)
    rowind = Int[]
    nzval = T[]
    colptr = Vector{Int}(N+1)
    u = zeros(T, N)
    dest = zeros(T, N)

    for i = 1:N
        u[i] = one(T)
        A_mul_B!(dest, D, u)
        js = find(dest)
        colptr[i] = length(nzval)+1
        if length(js) > 0
            append!(rowind, js)
            append!(nzval, dest[js])
        end
        u[i] = zero(T)
    end
    colptr[N+1] = length(nzval)+1

    return SparseMatrixCSC(M, N, colptr, rowind, nzval)
end
