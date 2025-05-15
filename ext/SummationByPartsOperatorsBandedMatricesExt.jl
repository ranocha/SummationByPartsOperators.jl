module SummationByPartsOperatorsBandedMatricesExt

if isdefined(Base, :get_extension)
    import BandedMatrices: BandedMatrices, BandedMatrix, isbanded, bandwidth
else
    import ..BandedMatrices: BandedMatrices, BandedMatrix, isbanded, bandwidth
end

using SummationByPartsOperators: DerivativeOperator, DissipationOperator,
                                 VarCoefDerivativeOperator,
                                 UniformNonperiodicCoupledOperator,
                                 lower_bandwidth, upper_bandwidth, grid, mul!

const BandedDerivativeOperator = Union{DerivativeOperator, DissipationOperator,
                                       VarCoefDerivativeOperator,
                                       UniformNonperiodicCoupledOperator}

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
    e = fill(zero(T), m)
    dest = similar(e)
    # left boundary
    for j in 1:max(l, u)
        e[j] = one(T)
        mul!(dest, D, e)
        e[j] = zero(T)
        # Same as
        #   B.data[(end - j - l + 1):(end - j - l + min(j + l, length(e))), j] = dest[1:min(j + l,
        #                                                                                   length(e))]
        # but without allocations
        num_entries = min(j + l, length(e))
        for i in 1:num_entries
            B.data[end - j - l + i, j] = dest[i]
        end
    end
    # inner part
    for j in (max(l, u) + 1):(length(e) - max(l, u))
        e[j] = one(T)
        mul!(dest, D, e)
        e[j] = zero(T)
        # Same as
        #   B.data[:, j] = dest[(j - u):(j + l)]
        # but without allocations
        num_entries = u + l + 1
        for i in 1:num_entries
            B.data[i, j] = dest[j - u - 1 + i]
        end
    end
    # right boundary
    for j in (length(e) - max(l, u) + 1):length(e)
        e[j] = one(T)
        mul!(dest, D, e)
        e[j] = zero(T)
        # B.data[1:length(e)+1-j+u,j] = dest[j-u:end]
        # Same as
        #   B.data[(1 + max(1, j - u) - (j - u)):(length(e) + 1 - j + u), j] = dest[max(1,
        #                                                                               j - u):end]
        # but without allocations
        num_entries = length(e) + 1 - max(1, j - u)
        for i in 1:num_entries
            B.data[max(1, j - u) - (j - u) + i, j] = dest[max(1, j - u) - 1 + i]
        end
    end

    return B
end

end # module
