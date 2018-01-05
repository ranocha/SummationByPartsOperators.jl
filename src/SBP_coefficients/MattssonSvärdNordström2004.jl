
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


function dissipation_coefficients(source::MattssonSvärdNordström2004, order::Int, grid, left_weights, right_weights, parallel=Val{:serial}())
    T = eltype(grid)
    if order == 2
        boundary_length = 1
        if length(left_weights) < boundary_length
            left_weights = (left_weights..., ntuple(j->one(T), boundary_length-length(left_weights)))
        end
        if length(right_weights) < boundary_length
            right_weights = (right_weights..., ntuple(j->one(T), boundary_length-length(right_weights)))
        end

        left_boundary = (
            # d1
            left_weights[1] .\ (
                DerivativeCoefficientRow{T,0,2}(SVector(T(1), T(1) )),
                DerivativeCoefficientRow{T,-1,2}(SVector(T(1), T(-1) )),
            ),
        )
        for i in boundary_length+1:length(left_weights)
            di = left_weights[i] .\ (
                ntuple(j->DerivativeCoefficientRow{T,0,1}(SVector(T(0), )), i-2)...,
                DerivativeCoefficientRow{T,1,1}(SVector(T(-1), )),
                DerivativeCoefficientRow{T,0,2}(SVector(T(1), T(1) )),
                DerivativeCoefficientRow{T,0,1}(SVector(T(-1), )),
            )
            left_boundary = (left_boundary..., di)
        end

        right_boundary = (
            # d1
            right_weights[1] .\ (
                DerivativeCoefficientRow{T,0,1}(SVector(T(1), )),
                DerivativeCoefficientRow{T,1,1}(SVector(T(-1), )),
            ),
        )
        for i in boundary_length+1:length(right_weights)
            di = right_weights[i] .\ (
                ntuple(j->DerivativeCoefficientRow{T,0,1}(SVector(T(0), )), i-2)...,
                DerivativeCoefficientRow{T,0,1}(SVector(T(-1), )),
                DerivativeCoefficientRow{T,0,2}(SVector(T(1), T(1) )),
                DerivativeCoefficientRow{T,1,1}(SVector(T(-1), )),
            )
            right_boundary = (right_boundary..., di)
        end

        lower_coef = (
            DerivativeCoefficientRow{T,1,1}(SVector(T(-1), )),
        )
        central_coef = DerivativeCoefficientRow{T,0,2}(SVector(T(1), T(1) ))
        upper_coef = (
            DerivativeCoefficientRow{T,-1,1}(SVector(T(-1), )),
        )

        coef = DissipationCoefficients(left_boundary, right_boundary,
                                       lower_coef, central_coef, upper_coef,
                                       parallel, order, 2, source)
        b = ones(grid)
        b[1] = T(0)
    elseif order == 4
        boundary_length = 3
        if length(left_weights) < boundary_length
            left_weights = (left_weights..., ntuple(j->one(T), boundary_length-length(left_weights)))
        end
        if length(right_weights) < boundary_length
            right_weights = (right_weights..., ntuple(j->one(T), boundary_length-length(right_weights)))
        end

        left_boundary = (
            # d1
            left_weights[1] .\ (
                DerivativeCoefficientRow{T,0,2}(SVector(T(1), T(1) )),
                DerivativeCoefficientRow{T,-1,2}(SVector(T(-2), T(-2) )),
                DerivativeCoefficientRow{T,-2,2}(SVector(T(1), T(1) )),
            ),
            # d2
            left_weights[2] .\ (
                DerivativeCoefficientRow{T,0,2}(SVector(T(-2), T(-2) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(4), T(4), T(1) )),
                DerivativeCoefficientRow{T,-2,3}(SVector(T(-2), T(-2), T(-2) )),
                DerivativeCoefficientRow{T,-1,1}(SVector(T(1), )),
            ),
            # d3
            left_weights[3] .\ (
                DerivativeCoefficientRow{T,0,2}(SVector(T(1), T(1) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(-2), T(-2), T(-2) )),
                DerivativeCoefficientRow{T,-2,4}(SVector(T(1), T(1), T(4), T(1) )),
                DerivativeCoefficientRow{T,-1,2}(SVector(T(-2), T(-2) )),
                DerivativeCoefficientRow{T,-1,1}(SVector(T(1), )),
            ),
        )
        for i in boundary_length+1:length(left_weights)
            di = left_weights[i] .\ (
                ntuple(j->DerivativeCoefficientRow{T,0,1}(SVector(T(0), )), i-boundary_length)...,
                DerivativeCoefficientRow{T,1,1}(SVector(T(1), )),
                DerivativeCoefficientRow{T,0,2}(SVector(T(-2), T(-2) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(1), T(4), T(1) )),
                DerivativeCoefficientRow{T,-1,2}(SVector(T(-2), T(-2) )),
                DerivativeCoefficientRow{T,0,1}(SVector(T(1), )),
            )
            left_boundary = (left_boundary..., di)
        end

        right_boundary = (
            # d1
            right_weights[1] .\ (
                DerivativeCoefficientRow{T,-1,2}(SVector(T(1), T(1) )),
                DerivativeCoefficientRow{T,0,2}(SVector(T(-2), T(-2) )),
                DerivativeCoefficientRow{T,1,2}(SVector(T(1), T(1) )),
            ),
            # d2
            right_weights[2] .\ (
                DerivativeCoefficientRow{T,-1,2}(SVector(T(-2), T(-2) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(1), T(4), T(4) )),
                DerivativeCoefficientRow{T,0,3}(SVector(T(-2), T(-2), T(-2) )),
                DerivativeCoefficientRow{T,-1,1}(SVector(T(1), )),
            ),
            # d3
            right_weights[3] .\ (
                DerivativeCoefficientRow{T,-1,2}(SVector(T(1), T(1) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(-2), T(-2), T(-2) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(1), T(4), T(1), T(1) )),
                DerivativeCoefficientRow{T,0,2}(SVector(T(-2), T(-2) )),
                DerivativeCoefficientRow{T,1,1}(SVector(T(1), )),
            ),
        )
        for i in boundary_length+1:length(right_weights)
            di = right_weights[i] .\ (
                ntuple(j->DerivativeCoefficientRow{T,0,1}(SVector(T(0), )), i-boundary_length)...,
                DerivativeCoefficientRow{T,-1,1}(SVector(T(1), )),
                DerivativeCoefficientRow{T,-1,2}(SVector(T(-2), T(-2) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(1), T(4), T(1) )),
                DerivativeCoefficientRow{T,0,2}(SVector(T(-2), T(-2) )),
                DerivativeCoefficientRow{T,1,1}(SVector(T(1), )),
            )
            right_boundary = (right_boundary..., di)
        end

        lower_coef = (
            DerivativeCoefficientRow{T,0,2}(SVector(T(-2), T(-2) )),
            DerivativeCoefficientRow{T,1,1}(SVector(T(1), )),
        )
        central_coef = DerivativeCoefficientRow{T,-1,3}(SVector(T(1), T(4), T(1) ))
        upper_coef = (
            DerivativeCoefficientRow{T,-1,2}(SVector(T(-2), T(-2) )),
            DerivativeCoefficientRow{T,0,1}(SVector(T(1), )),
        )

        coef = DissipationCoefficients(left_boundary, right_boundary,
                                       lower_coef, central_coef, upper_coef,
                                       parallel, order, 2, source)
        b = ones(grid)
        b[1] = b[end] = T(0)
    elseif order == 6
        boundary_length = 4
        if length(left_weights) < boundary_length
            left_weights = (left_weights..., ntuple(j->one(T), boundary_length-length(left_weights)))
        end
        if length(right_weights) < boundary_length
            right_weights = (right_weights..., ntuple(j->one(T), boundary_length-length(right_weights)))
        end

        left_boundary = (
            # d1
            left_weights[1] .\ (
                DerivativeCoefficientRow{T,0,3}(SVector(T(1), T(1), T(1) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(-3), T(-3), T(-3) )),
                DerivativeCoefficientRow{T,-2,3}(SVector(T(3), T(3), T(3) )),
                DerivativeCoefficientRow{T,-3,3}(SVector(T(-1), T(-1), T(-1) )),
            ),
            # d2
            left_weights[2] .\ (
                DerivativeCoefficientRow{T,0,3}(SVector(T(-3), T(-3), T(-3) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(9), T(9), T(9), T(1) )),
                DerivativeCoefficientRow{T,-2,4}(SVector(T(-9), T(-9), T(-9), T(-3) )),
                DerivativeCoefficientRow{T,-3,4}(SVector(T(3), T(3), T(3), T(3) )),
                DerivativeCoefficientRow{T,-1,1}(SVector(T(-1), )),
            ),
            # d3
            left_weights[3] .\ (
                DerivativeCoefficientRow{T,0,3}(SVector(T(3), T(3), T(3) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(-9), T(-9), T(-9), T(-3) )),
                DerivativeCoefficientRow{T,-2,5}(SVector(T(9), T(9), T(9), T(9), T(1) )),
                DerivativeCoefficientRow{T,-3,5}(SVector(T(-3), T(-3), T(-3), T(-9), T(-3) )),
                DerivativeCoefficientRow{T,-1,2}(SVector(T(3), T(3) )),
                DerivativeCoefficientRow{T,-1,1}(SVector(T(-1), )),
            ),
            # d4
            left_weights[4] .\ (
                DerivativeCoefficientRow{T,0,3}(SVector(T(-1), T(-1), T(-1) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(3), T(3), T(3), T(3) )),
                DerivativeCoefficientRow{T,-2,5}(SVector(T(-3), T(-3), T(-3), T(-9), T(-3) )),
                DerivativeCoefficientRow{T,-3,6}(SVector(T(1), T(1), T(1), T(9), T(9), T(1) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(-3), T(-9), T(-3) )),
                DerivativeCoefficientRow{T,-1,2}(SVector(T(3), T(3) )),
                DerivativeCoefficientRow{T,-1,1}(SVector(T(-1), )),
            )
        )
        for i in boundary_length+1:length(left_weights)
            di = left_weights[i] .\ (
                ntuple(j->DerivativeCoefficientRow{T,0,1}(SVector(T(0), )), i-boundary_length)...,
                DerivativeCoefficientRow{T,2,1}(SVector(T(-1), )),
                DerivativeCoefficientRow{T,1,2}(SVector(T(3), T(3) )),
                DerivativeCoefficientRow{T,0,3}(SVector(T(-3), T(-9), T(-3) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(1), T(9), T(9), T(1) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(-3), T(-9), T(-3) )),
                DerivativeCoefficientRow{T,-1,2}(SVector(T(3), T(3) )),
                DerivativeCoefficientRow{T,-1,1}(SVector(T(-1), )),
            )
            left_boundary = (left_boundary..., di)
        end

        right_boundary = (
            # d1
            right_weights[1] .\ (
                DerivativeCoefficientRow{T,-1,2}(SVector(T(1), T(1) )),
                DerivativeCoefficientRow{T,0,2}(SVector(T(-3), T(-3) )),
                DerivativeCoefficientRow{T,1,2}(SVector(T(3), T(3) )),
                DerivativeCoefficientRow{T,2,2}(SVector(T(-1), T(-1) )),
            ),
            # d2
            right_weights[2] .\ (
                DerivativeCoefficientRow{T,-1,2}(SVector(T(-3), T(-3) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(1), T(9), T(9) )),
                DerivativeCoefficientRow{T,0,3}(SVector(T(-3), T(-9), T(-9) )),
                DerivativeCoefficientRow{T,1,3}(SVector(T(3), T(3), T(3) )),
                DerivativeCoefficientRow{T,1,1}(SVector(T(-1), )),
            ),
            # d3
            right_weights[3] .\ (
                DerivativeCoefficientRow{T,-1,2}(SVector(T(3), T(3) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(-3), T(-9), T(-9) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(1), T(9), T(9), T(9) )),
                DerivativeCoefficientRow{T,0,4}(SVector(T(-3), T(-9), T(-3), T(-3) )),
                DerivativeCoefficientRow{T,1,2}(SVector(T(3), T(3) )),
                DerivativeCoefficientRow{T,2,1}(SVector(T(-1), )),
            ),
            # d4
            right_weights[4] .\ (
                DerivativeCoefficientRow{T,-1,2}(SVector(T(-1), T(-1) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(3), T(3), T(3) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(-3), T(-9), T(-3), T(-3) )),
                DerivativeCoefficientRow{T,-1,5}(SVector(T(1), T(9), T(9), T(1), T(1) )),
                DerivativeCoefficientRow{T,0,3}(SVector(T(-3), T(-9), T(-3) )),
                DerivativeCoefficientRow{T,1,2}(SVector(T(3), T(3) )),
                DerivativeCoefficientRow{T,2,1}(SVector(T(-1), )),
            )
        )
        for i in boundary_length+1:length(right_weights)
            di = right_weights[i] .\ (
                ntuple(j->DerivativeCoefficientRow{T,0,1}(SVector(T(0), )), i-boundary_length)...,
                DerivativeCoefficientRow{T,-1,1}(SVector(T(-1), )),
                DerivativeCoefficientRow{T,-1,2}(SVector(T(3), T(3) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(-3), T(-9), T(-3) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(1), T(9), T(9), T(1) )),
                DerivativeCoefficientRow{T,0,3}(SVector(T(-3), T(-9), T(-3) )),
                DerivativeCoefficientRow{T,1,2}(SVector(T(3), T(3) )),
                DerivativeCoefficientRow{T,2,1}(SVector(T(-1), )),
            )
            right_boundary = (right_boundary..., di)
        end

        lower_coef = (
        DerivativeCoefficientRow{T,0,3}(SVector(T(-3), T(-9), T(-3) )),
        DerivativeCoefficientRow{T,1,2}(SVector(T(3), T(3) )),
        DerivativeCoefficientRow{T,2,1}(SVector(T(-1), )),
        )
        central_coef = DerivativeCoefficientRow{T,-1,4}(SVector(T(1), T(9), T(9), T(1) ))
        upper_coef = (
            DerivativeCoefficientRow{T,-1,3}(SVector(T(-3), T(-9), T(-3) )),
            DerivativeCoefficientRow{T,-1,2}(SVector(T(3), T(3) )),
            DerivativeCoefficientRow{T,-1,1}(SVector(T(-1), )),
        )

        coef = DissipationCoefficients(left_boundary, right_boundary,
                                       lower_coef, central_coef, upper_coef,
                                       parallel, order, 2, source)
        b = ones(grid)
        b[1] = b[2] = b[end] = T(0)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end

    coef, b
end
