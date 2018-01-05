
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
            left_weights = (left_weights..., ntuple(j->one(T), boundary_length-length(left_weights))...)
        end
        if length(right_weights) < boundary_length
            right_weights = (right_weights..., ntuple(j->one(T), boundary_length-length(right_weights))...)
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
            left_weights = (left_weights..., ntuple(j->one(T), boundary_length-length(left_weights))...)
        end
        if length(right_weights) < boundary_length
            right_weights = (right_weights..., ntuple(j->one(T), boundary_length-length(right_weights))...)
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
            left_weights = (left_weights..., ntuple(j->one(T), boundary_length-length(left_weights))...)
        end
        if length(right_weights) < boundary_length
            right_weights = (right_weights..., ntuple(j->one(T), boundary_length-length(right_weights))...)
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
    elseif order == 8
        boundary_length = 5
        if length(left_weights) < boundary_length
            left_weights = (left_weights..., ntuple(j->one(T), boundary_length-length(left_weights))...)
        end
        if length(right_weights) < boundary_length
            right_weights = (right_weights..., ntuple(j->one(T), boundary_length-length(right_weights))...)
        end

        left_boundary = (
            # d1
            left_weights[1] .\ (
                DerivativeCoefficientRow{T,0,3}(SVector(T(1), T(1), T(1) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-2,3}(SVector(T(6), T(6), T(6) )),
                DerivativeCoefficientRow{T,-3,3}(SVector(T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-4,3}(SVector(T(1), T(1), T(1) )),
            ),
            # d2
            left_weights[2] .\ (
                DerivativeCoefficientRow{T,0,3}(SVector(T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(16), T(16), T(16), T(1) )),
                DerivativeCoefficientRow{T,-2,4}(SVector(T(-24), T(-24), T(-24), T(-4) )),
                DerivativeCoefficientRow{T,-3,4}(SVector(T(16), T(16), T(16), T(6) )),
                DerivativeCoefficientRow{T,-4,4}(SVector(T(-4), T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-2,1}(SVector(T(1), )),
            ),
            # d3
            left_weights[3] .\ (
                DerivativeCoefficientRow{T,0,3}(SVector(T(6), T(6), T(6) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(-24), T(-24), T(-24), T(-4) )),
                DerivativeCoefficientRow{T,-2,5}(SVector(T(36), T(36), T(36), T(16), T(1) )),
                DerivativeCoefficientRow{T,-3,5}(SVector(T(-24), T(-24), T(-24), T(-24), T(-4) )),
                DerivativeCoefficientRow{T,-4,5}(SVector(T(6), T(6), T(6), T(16), T(6) )),
                DerivativeCoefficientRow{T,-2,2}(SVector(T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-2,1}(SVector(T(1), )),
            ),
            # d4
            left_weights[4] .\ (
                DerivativeCoefficientRow{T,0,3}(SVector(T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(16), T(16), T(16), T(6) )),
                DerivativeCoefficientRow{T,-2,5}(SVector(T(-24), T(-24), T(-24), T(-24), T(-4) )),
                DerivativeCoefficientRow{T,-3,6}(SVector(T(16), T(16), T(16), T(36), T(16), T(1) )),
                DerivativeCoefficientRow{T,-4,6}(SVector(T(-4), T(-4), T(-4), T(-24), T(-24), T(-4) )),
                DerivativeCoefficientRow{T,-2,3}(SVector(T(6), T(16), T(6) )),
                DerivativeCoefficientRow{T,-2,2}(SVector(T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-2,1}(SVector(T(1), )),
            ),
            # d5
            left_weights[5] .\ (
                DerivativeCoefficientRow{T,0,3}(SVector(T(1), T(1), T(1) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(-4), T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-2,5}(SVector(T(6), T(6), T(6), T(16), T(6) )),
                DerivativeCoefficientRow{T,-3,6}(SVector(T(-4), T(-4), T(-4), T(-24), T(-24), T(-4) )),
                DerivativeCoefficientRow{T,-4,7}(SVector(T(1), T(1), T(1), T(16), T(36), T(16), T(1) )),
                DerivativeCoefficientRow{T,-2,4}(SVector(T(-4), T(-24), T(-24), T(-4) )),
                DerivativeCoefficientRow{T,-2,3}(SVector(T(6), T(16), T(6) )),
                DerivativeCoefficientRow{T,-2,2}(SVector(T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-2,1}(SVector(T(1), )),
            ),
        )
        for i in boundary_length+1:length(left_weights)
            di = left_weights[i] .\ (
                ntuple(j->DerivativeCoefficientRow{T,0,1}(SVector(T(0), )), i-boundary_length)...,
                DerivativeCoefficientRow{T,2,1}(SVector(T(1), )),
                DerivativeCoefficientRow{T,1,2}(SVector(T(-4), T(-4) )),
                DerivativeCoefficientRow{T,0,3}(SVector(T(6), T(16), T(6) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(-4), T(-24), T(-24), T(-4) )),
                DerivativeCoefficientRow{T,-2,5}(SVector(T(1), T(16), T(36), T(16), T(1) )),
                DerivativeCoefficientRow{T,-2,4}(SVector(T(-4), T(-24), T(-24), T(-4) )),
                DerivativeCoefficientRow{T,-2,3}(SVector(T(6), T(16), T(6) )),
                DerivativeCoefficientRow{T,-2,2}(SVector(T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-2,1}(SVector(T(1), )),
            )
            left_boundary = (left_boundary..., di)
        end

        right_boundary = (
            # d1
            right_weights[1] .\ (
                DerivativeCoefficientRow{T,-2,3}(SVector(T(1), T(1), T(1) )),
                DerivativeCoefficientRow{T,-1,3}(SVector(T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,0,3}(SVector(T(6), T(6), T(6) )),
                DerivativeCoefficientRow{T,1,3}(SVector(T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,2,3}(SVector(T(1), T(1), T(1) )),
            ),
            # d2
            right_weights[2] .\ (
                DerivativeCoefficientRow{T,-2,3}(SVector(T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-2,4}(SVector(T(1), T(16), T(16), T(16) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(-4), T(-24), T(-24), T(-24) )),
                DerivativeCoefficientRow{T,0,4}(SVector(T(6), T(16), T(16), T(16) )),
                DerivativeCoefficientRow{T,1,4}(SVector(T(-4), T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,2,1}(SVector(T(1), )),
            ),
            # d3
            right_weights[3] .\ (
                DerivativeCoefficientRow{T,-2,3}(SVector(T(6), T(6), T(6) )),
                DerivativeCoefficientRow{T,-2,4}(SVector(T(-4), T(-24), T(-24), T(-24) )),
                DerivativeCoefficientRow{T,-2,5}(SVector(T(1), T(16), T(36), T(36), T(36) )),
                DerivativeCoefficientRow{T,-1,5}(SVector(T(-4), T(-24), T(-24), T(-24), T(-24) )),
                DerivativeCoefficientRow{T,0,5}(SVector(T(6), T(16), T(6), T(6), T(6) )),
                DerivativeCoefficientRow{T,1,2}(SVector(T(-4), T(-4) )),
                DerivativeCoefficientRow{T,2,1}(SVector(T(1), )),
            ),
            # d4
            right_weights[4] .\ (
                DerivativeCoefficientRow{T,-2,3}(SVector(T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-2,4}(SVector(T(6), T(16), T(16), T(16) )),
                DerivativeCoefficientRow{T,-2,5}(SVector(T(-4), T(-24), T(-24), T(-24), T(-24) )),
                DerivativeCoefficientRow{T,-2,6}(SVector(T(1), T(16), T(36), T(16), T(16), T(16) )),
                DerivativeCoefficientRow{T,-1,6}(SVector(T(-4), T(-24), T(-24), T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,0,3}(SVector(T(6), T(16), T(6) )),
                DerivativeCoefficientRow{T,1,2}(SVector(T(-4), T(-4) )),
                DerivativeCoefficientRow{T,2,1}(SVector(T(1), )),
            ),
            # d5
            right_weights[5] .\ (
                DerivativeCoefficientRow{T,-2,3}(SVector(T(1), T(1), T(1) )),
                DerivativeCoefficientRow{T,-2,4}(SVector(T(-4), T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-2,5}(SVector(T(6), T(16), T(6), T(6), T(6) )),
                DerivativeCoefficientRow{T,-2,6}(SVector(T(-4), T(-24), T(-24), T(-4), T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-2,7}(SVector(T(1), T(16), T(36), T(16), T(1), T(1), T(1) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(-4), T(-24), T(-24), T(-4) )),
                DerivativeCoefficientRow{T,0,3}(SVector(T(6), T(16), T(6) )),
                DerivativeCoefficientRow{T,1,2}(SVector(T(-4), T(-4) )),
                DerivativeCoefficientRow{T,2,1}(SVector(T(1), )),
            ),
        )
        for i in boundary_length+1:length(right_weights)
            di = right_weights[i] .\ (
                ntuple(j->DerivativeCoefficientRow{T,0,1}(SVector(T(0), )), i-boundary_length)...,
                DerivativeCoefficientRow{T,-2,1}(SVector(T(1), )),
                DerivativeCoefficientRow{T,-2,2}(SVector(T(-4), T(-4) )),
                DerivativeCoefficientRow{T,-2,3}(SVector(T(6), T(16), T(6) )),
                DerivativeCoefficientRow{T,-2,4}(SVector(T(-4), T(-24), T(-24), T(-4) )),
                DerivativeCoefficientRow{T,-2,5}(SVector(T(1), T(16), T(36), T(16), T(1) )),
                DerivativeCoefficientRow{T,-1,4}(SVector(T(-4), T(-24), T(-24), T(-4) )),
                DerivativeCoefficientRow{T,0,3}(SVector(T(6), T(16), T(6) )),
                DerivativeCoefficientRow{T,1,2}(SVector(T(-4), T(-4) )),
                DerivativeCoefficientRow{T,2,1}(SVector(T(1), )),
            )
            right_boundary = (right_boundary..., di)
        end

        lower_coef = (
            DerivativeCoefficientRow{T,-1,4}(SVector(T(-4), T(-24), T(-24), T(-4) )),
            DerivativeCoefficientRow{T,0,3}(SVector(T(6), T(16), T(6) )),
            DerivativeCoefficientRow{T,1,2}(SVector(T(-4), T(-4) )),
            DerivativeCoefficientRow{T,2,1}(SVector(T(1), )),
        )
        central_coef = DerivativeCoefficientRow{T,-2,5}(SVector(T(1), T(16), T(36), T(16), T(1) ))
        upper_coef = (
            DerivativeCoefficientRow{T,-2,4}(SVector(T(-4), T(-24), T(-24), T(-4) )),
            DerivativeCoefficientRow{T,-2,3}(SVector(T(6), T(16), T(6) )),
            DerivativeCoefficientRow{T,-2,2}(SVector(T(-4), T(-4) )),
            DerivativeCoefficientRow{T,-2,1}(SVector(T(1), )),
        )

        coef = DissipationCoefficients(left_boundary, right_boundary,
                                       lower_coef, central_coef, upper_coef,
                                       parallel, order, 2, source)
        b = ones(grid)
        b[1] = b[2] = b[end] = b[end-1] = T(0)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end

    coef, b
end


function new_dissipation_coefficients(source::MattssonSvärdNordström2004, order::Int, grid, left_weights, right_weights, parallel=Val{:serial}())
    T = promote_type(eltype(grid), eltype(left_weights), eltype(right_weights))
    if order == 2
        inv_left_weights = one(T) ./ left_weights
        inv_right_weights = one(T) ./ right_weights
        coefficient_cache = MattssonSvärdNordström2004Cache2{T,
            typeof(inv_left_weights),typeof(inv_right_weights)}(inv_left_weights, inv_right_weights)
        b = ones(grid)
        b[1] = T(0)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end


    NewDissipationCoefficients(coefficient_cache, parallel, order, 2, source), b
end



struct MattssonSvärdNordström2004Cache2{T,InvLeftWeights,InvRightWeights} <: AbstractCoefficientCache{T}
    inv_left_weights::InvLeftWeights
    inv_right_weights::InvRightWeights
end

function Base.checkbounds(::Type{Bool}, u::AbstractVector, cache::MattssonSvärdNordström2004Cache2)
    L = length(u)
    (L > 2) && (L > length(cache.inv_left_weights)+length(cache.inv_right_weights))
end

function convolve_boundary_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache2, u::AbstractVector, b::AbstractVector, α)
    #TODO: Modified boundaries...
    #TODO: Not finished...
    @unpack inv_left_weights, inv_right_weights = cache

    @inbounds begin
        dest[1] = α * inv_left_weights[1] * ( (b[1] + b[2]) * u[1] + (b[1] - b[2]) * u[2] )

        dest[end] = α * inv_right_weights[1] * ( b[end] * u[end] - b[end] * u[end-1] )
    end
end

function convolve_interior_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache2, u::AbstractVector, b::AbstractVector, α, parallel)
    for i in (length(cache.inv_left_weights)+1):(length(dest)-length(cache.inv_right_weights))
        convolve_interior_coefficients_loopbody!(dest, i, cache, u, b, α)
    end
end

function convolve_interior_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache2, u::AbstractVector, b::AbstractVector, α, parallel::Val{:threads})
    Threads.@threads for i in (length(cache.inv_left_weights)+1):(length(dest)-length(cache.inv_right_weights))
        convolve_interior_coefficients_loopbody!(dest, i, cache, u, b, α)
    end
end

@inline function convolve_interior_coefficients_loopbody!(dest, i, cache::MattssonSvärdNordström2004Cache2, u, b, α)
    #TODO: Modified boundaries...
    @inbounds begin
        b_i = b[i]
        b_ip1 = b[i+1]

        dest[i] = α * ( -b_i * u[i-1] + (b_i + b_ip1) * u[i] - b_ip1 * u[i+1] )
    end
end
