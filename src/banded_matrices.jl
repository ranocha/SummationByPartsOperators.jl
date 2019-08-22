import .BandedMatrices: BandedMatrix, isbanded, bandwidth

const BandedDerivativeOperator= Union{DerivativeOperator,DissipationOperator,VarCoefDerivativeOperator}

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
    for j in 1:l
        e[j] = 1
        mul!(dest, D, e)
        e[j] = 0
        B.data[end-j-u+1:end,j] = dest[1:j+u]
    end
    # inner part
    for j in l+1:length(e)-u
        e[j] = 1
        mul!(dest, D, e)
        e[j] = 0
        B.data[:,j] = dest[j-l:j+u]
    end
    # right boundary
    for j in length(e)-u+1:length(e)
        e[j] = 1
        mul!(dest, D, e)
        e[j] = 0
        B.data[1:length(e)+1-j+l,j] = dest[j-l:end]
     end

    B
end
