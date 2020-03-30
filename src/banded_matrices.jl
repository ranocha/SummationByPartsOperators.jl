import .BandedMatrices: BandedMatrix, isbanded, bandwidth

const BandedDerivativeOperator = Union{DerivativeOperator, DissipationOperator, VarCoefDerivativeOperator, UniformNonperiodicCoupledOperator}

@inline function bandwidth(D::BandedDerivativeOperator, k::Int)
    if k == 1
        lower_bandwidth(D)
     else
        upper_bandwidth(D)
    end
end

isbanded(D::BandedDerivativeOperator) = true

function BandedMatrix(D::BandedDerivativeOperator)
    T = eltype(D)
    l = lower_bandwidth(D)
    u = upper_bandwidth(D)
    n,m = size(D)

    # create uninitialised matrix
    B = BandedMatrix(BandedMatrices.Zeros{T}(n, m), (l, u))

    # B.data[i,:] is the ith band, starting with the upper ones
    # B.data[:,j] is the jth column, entries out of range are undefined
    e = fill(zero(T), length(grid(D)))
    dest = similar(e)
    # left boundary
    for j in 1:max(l,u)
        e[j] = 1
        mul!(dest, D, e)
        e[j] = 0
        B.data[end-j-l+1:end,j] = dest[1:j+l]
    end
    # inner part
    for j in max(l,u)+1:length(e)-max(l,u)
        e[j] = 1
        mul!(dest, D, e)
        e[j] = 0
        B.data[:,j] = dest[j-u:j+l]
    end
    # right boundary
    for j in length(e)-max(l,u)+1:length(e)
        e[j] = 1
        mul!(dest, D, e)
        e[j] = 0
        B.data[1:length(e)+1-j+u,j] = dest[j-u:end]
     end

    B
end
