
"""
    MattssonSvärdNordström2004

Coefficients of the SBP operators given in
  Mattsson, Svärd, Nordström (2004)
  Stable and Accurate Artificial Dissipation.
  Journal of Scientific Computing 21.1, pp. 57-79.
"""
struct MattssonSvärdNordström2004 <: SourceOfCoefficients end

function Base.show(io::IO, ::MattssonSvärdNordström2004)
    print(io,
        "  Mattsson, Svärd, Nordström (2004) \n",
        "  Stable and Accurate Artificial Dissipation. \n",
        "  Journal of Scientific Computing 21.1, pp. 57-79. \n")
end


function dissipation_coefficients(source::MattssonSvärdNordström2004, order::Int, grid, parallel=Val{:serial}())
    T = eltype(grid)
    if order == 2
        left_boundary = (
            # d1
            (
                DerivativeCoefficientRow{T,0,2}(SVector(T(1), T(1) )),
                DerivativeCoefficientRow{T,-1,2}(SVector(T(1), T(-1) )),
            ),
        )
        right_boundary = (
            .- left_boundary[1],
        )
        lower_coef = (
            DerivativeCoefficientRow{T,1,1}(SVector(T(-1), )),
        )
        central_coef = DerivativeCoefficientRow{T,0,2}(SVector(T(1), T(1) ))
        upper_coef = (
            DerivativeCoefficientRow{T,0,1}(SVector(T(-1), )),
        )

        coef = DissipationCoefficients(left_boundary, right_boundary,
                                       lower_coef, central_coef, upper_coef,
                                       parallel, order, 2, source)
        b = ones(grid)
        b[1] = b[end] = T(0)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end

    coef, b
end
