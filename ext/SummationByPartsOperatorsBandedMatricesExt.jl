module SummationByPartsOperatorsBandedMatricesExt

if isdefined(Base, :get_extension)
    import BandedMatrices: BandedMatrices, BandedMatrix, isbanded, bandwidth
else
    import ..BandedMatrices: BandedMatrices, BandedMatrix, isbanded, bandwidth
end

using SummationByPartsOperators:
    DerivativeOperator,
    DissipationOperator,
    VarCoefDerivativeOperator,
    UniformNonperiodicCoupledOperator,
    lower_bandwidth,
    upper_bandwidth,
    grid,
    mul!

const BandedDerivativeOperator = Union{
    DerivativeOperator,
    DissipationOperator,
    VarCoefDerivativeOperator,
    UniformNonperiodicCoupledOperator,
}

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
    n, m = size(D)

    # create uninitialised matrix
    B = BandedMatrix(BandedMatrices.Zeros{T}(n, m), (l, u))

    # B.data[i,:] is the ith band, starting with the upper ones
    # B.data[:,j] is the jth column, entries out of range are undefined
    e = fill(zero(T), length(grid(D)))
    dest = similar(e)
    # left boundary
    for j = 1:max(l, u)
        e[j] = 1
        mul!(dest, D, e)
        e[j] = 0
        B.data[end-j-l+1:end-j-l+min(j + l, length(e)), j] = dest[1:min(j + l, length(e))]
    end
    # inner part
    for j = max(l, u)+1:length(e)-max(l, u)
        e[j] = 1
        mul!(dest, D, e)
        e[j] = 0
        B.data[:, j] = dest[j-u:j+l]
    end
    # right boundary
    for j = length(e)-max(l, u)+1:length(e)
        e[j] = 1
        mul!(dest, D, e)
        e[j] = 0
        # B.data[1:length(e)+1-j+u,j] = dest[j-u:end]
        B.data[1+max(1, j - u)-(j-u):length(e)+1-j+u, j] = dest[max(1, j - u):end]
    end

    B
end

end # module
