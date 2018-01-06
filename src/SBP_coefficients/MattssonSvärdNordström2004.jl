
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
        coefficient_cache = MattssonSvärdNordström2004Cache2(inv_left_weights, inv_right_weights)
        b = ones(grid)
        b[1] = T(0)
    elseif order == 4
        inv_left_weights = one(T) ./ left_weights
        inv_right_weights = one(T) ./ right_weights
        coefficient_cache = MattssonSvärdNordström2004Cache4(inv_left_weights, inv_right_weights)
        b = ones(grid)
        b[1] = b[end] = T(0)
    elseif order == 6
        inv_left_weights = one(T) ./ left_weights
        inv_right_weights = one(T) ./ right_weights
        coefficient_cache = MattssonSvärdNordström2004Cache6(inv_left_weights, inv_right_weights)
        b = ones(grid)
        b[1] = b[2] = b[end] = T(0)
    elseif order == 8
        inv_left_weights = one(T) ./ left_weights
        inv_right_weights = one(T) ./ right_weights
        coefficient_cache = MattssonSvärdNordström2004Cache8(inv_left_weights, inv_right_weights)
        b = ones(grid)
        b[1] = b[2] = b[end] = b[end-1] = T(0)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end


    NewDissipationCoefficients(coefficient_cache, parallel, order, 2, source), b
end



struct MattssonSvärdNordström2004Cache2{T,LeftWidth,RightWidth} <: AbstractCoefficientCache{T}
    inv_left_weights::SVector{LeftWidth,T}
    inv_right_weights::SVector{RightWidth,T}

    function MattssonSvärdNordström2004Cache2(inv_left_weights::SVector{LeftWidth,T}, inv_right_weights::SVector{RightWidth,T}) where {T,LeftWidth,RightWidth}
        boundary_length = 2

        if length(inv_left_weights) < boundary_length
            inv_left_weights = SVector(inv_left_weights..., ntuple(j->one(T), boundary_length-length(inv_left_weights))...)
        end
        if length(inv_right_weights) < boundary_length
            inv_right_weights = SVector(inv_right_weights..., ntuple(j->one(T), boundary_length-length(inv_right_weights))...)
        end

        new{T,length(inv_left_weights),length(inv_right_weights)}(inv_left_weights, inv_right_weights)
    end
end

function convolve_boundary_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache2, u::AbstractVector, b::AbstractVector, α)
    @unpack inv_left_weights, inv_right_weights = cache

    @inbounds begin
        dest[1] = α * inv_left_weights[1] * (
                    (b[1] + b[2]) * u[1]
                    + (b[1] - b[2]) * u[2]
                )
        dest[2] = α * inv_left_weights[2] * (
                    (b[1] - b[2]) * u[1]
                    + (b[1] + b[2] + b[3]) * u[2]
                    - b[3] * u[3]
                )
        for i in 3:length(inv_left_weights)
            dest[i] = α * inv_left_weights[i] * (
                        -b[i] * u[i-1]
                        + (b[i] + b[i+1]) * u[i]
                        - b[i+1] * u[i+1]
                    )
        end

        dest[end] = α * inv_right_weights[1] * (
                        b[end] * u[end]
                        - b[end] * u[end-1]
                    )
        dest[end-1] = α * inv_right_weights[2] * (
                        -b[end-1] * u[end-2]
                        + (b[end-1] + b[end]) * u[end-1]
                        - b[end] * u[end]
                    )
        for i in 2:length(inv_right_weights)-1
            dest[end-i] = α * inv_right_weights[i+1] * (
                        -b[end-i] * u[end-(i-1)]
                        + (b[end-i] + b[end-(i+1)]) * u[end-i]
                        - b[end-(i+1)] * u[end-(i+1)]
                    )
        end
    end
end

@inline function convolve_interior_coefficients_loopbody(i, cache::MattssonSvärdNordström2004Cache2, u, b)
    @inbounds begin
        b_i = b[i]
        b_ip1 = b[i+1]

        retval = -b_i * u[i-1] + (b_i + b_ip1) * u[i] - b_ip1 * u[i+1]
    end

    retval
end



struct MattssonSvärdNordström2004Cache4{T,LeftWidth,RightWidth} <: AbstractCoefficientCache{T}
    inv_left_weights::SVector{LeftWidth,T}
    inv_right_weights::SVector{RightWidth,T}

    function MattssonSvärdNordström2004Cache4(inv_left_weights::SVector{LeftWidth,T}, inv_right_weights::SVector{RightWidth,T}) where {T,LeftWidth,RightWidth}
        boundary_length = 3

        if length(inv_left_weights) < boundary_length
            inv_left_weights = SVector(inv_left_weights..., ntuple(j->one(T), boundary_length-length(inv_left_weights))...)
        end
        if length(inv_right_weights) < boundary_length
            inv_right_weights = SVector(inv_right_weights..., ntuple(j->one(T), boundary_length-length(inv_right_weights))...)
        end

        new{T,length(inv_left_weights),length(inv_right_weights)}(inv_left_weights, inv_right_weights)
    end
end

function convolve_boundary_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache4, u::AbstractVector, b::AbstractVector, α)
    @unpack inv_left_weights, inv_right_weights = cache

    @inbounds begin
        dest[1] = α * inv_left_weights[1] * (
                    (b[1] + b[2]) * u[1]
                    - 2*(b[1] + b[2]) * u[2]
                    + (b[1] + b[2]) * u[3]
                )
        dest[2] = α * inv_left_weights[2] * (
                    -2*(b[1] + b[2]) * u[1]
                    + (4*b[1] + 4*b[2] + b[3]) * u[2]
                    -2*(b[1] + b[2] + b[3]) * u[3]
                    + b[3] * u[4]
                )
        dest[3] = α * inv_left_weights[3] * (
                    (b[1] + b[2]) * u[1]
                    -2*(b[1] + b[2] + b[3]) * u[2]
                    + (b[1] + b[2] + 4*b[3] + b[4]) * u[3]
                    -2*(b[3] + b[4]) * u[4]
                    + b[4] * u[5]
                )
        for i in 4:length(inv_left_weights)
            dest[i] = α * inv_left_weights[i] * (
                        b[i-1] * u[i-2]
                        - 2*(b[i-1] + b[i]) * u[i-1]
                        + (b[i-1] + 4*b[i] + b[i+1]) * u[i]
                        - 2*(b[i] + b[i+1]) * u[i+1]
                        + b[i+1] * u[i+2]
                    )
        end

        dest[end] = α * inv_right_weights[1] * (
                        (b[end] + b[end-1]) * u[end]
                        -2*(b[end] + b[end-1]) * u[end-1]
                        + (b[end] + b[end-1]) * u[end-2]
                    )
        dest[end-1] = α * inv_right_weights[2] * (
                        -2*(b[end] + b[end-1]) * u[end]
                        + (b[end-2] + 4*b[end-1] + 4*b[end]) * u[end-1]
                        -2*(b[end-2] + b[end-1] + b[end]) * u[end-2]
                        + b[end-2] * u[end-3]
                    )
        dest[end-2] = α * inv_right_weights[3] * (
                        (b[end-1] + b[end]) * u[end]
                        -2*(b[end-2] + b[end-1] + b[end]) * u[end-1]
                        + (b[end-3] + 4*b[end-2] + b[end-1] + b[end]) * u[end-2]
                        -2*(b[end-3] + b[end-2]) * u[end-3]
                        + b[end-3] * u[end-4]
                    )
        for i in 3:length(inv_right_weights)-1
            dest[end-i] = α * inv_right_weights[i+1] * (
                        b[end-(i-1)] * u[end-(i-2)]
                        - 2*(b[end-(i-1)] + b[i]) * u[end-(i-1)]
                        + (b[end-(i-1)] + 4*b[end-i] + b[end-(i+1)]) * u[end-i]
                        - 2*(b[end-i] + b[end-(i+1)]) * u[end-(i+1)]
                        + b[end-(i+1)] * u[end-(i+2)]
                    )
        end
    end
end

@inline function convolve_interior_coefficients_loopbody(i, cache::MattssonSvärdNordström2004Cache4, u, b)
    @inbounds begin
        b_im1 = b[i-1]
        b_i   = b[i]
        b_ip1 = b[i+1]

        retval = (
                    b_im1 * u[i-2]
                    - 2*(b_im1 + b_i) * u[i-1]
                    + (b_im1 + 4b_i + b_ip1) * u[i]
                    - 2*(b_i + b_ip1) * u[i+1]
                    + b_ip1 * u[i+2]
                )
    end

    retval
end



struct MattssonSvärdNordström2004Cache6{T,LeftWidth,RightWidth} <: AbstractCoefficientCache{T}
    inv_left_weights::SVector{LeftWidth,T}
    inv_right_weights::SVector{RightWidth,T}

    function MattssonSvärdNordström2004Cache6(inv_left_weights::SVector{LeftWidth,T}, inv_right_weights::SVector{RightWidth,T}) where {T,LeftWidth,RightWidth}
        boundary_length = 4

        if length(inv_left_weights) < boundary_length
            inv_left_weights = SVector(inv_left_weights..., ntuple(j->one(T), boundary_length-length(inv_left_weights))...)
        end
        if length(inv_right_weights) < boundary_length
            inv_right_weights = SVector(inv_right_weights..., ntuple(j->one(T), boundary_length-length(inv_right_weights))...)
        end

        new{T,length(inv_left_weights),length(inv_right_weights)}(inv_left_weights, inv_right_weights)
    end
end

function convolve_boundary_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache6, u::AbstractVector, b::AbstractVector, α)
    @unpack inv_left_weights, inv_right_weights = cache

    @inbounds begin
        dest[1] = α * inv_left_weights[1] * (
                    (b[1] + b[2] + b[3]) * u[1]
                    - 3*(b[1] + b[2] + b[3]) * u[2]
                    + 3*(b[1] + b[2] + b[3]) * u[3]
                    - (b[1] + b[2] + b[3]) * u[4]
                )
        dest[2] = α * inv_left_weights[2] * (
                    - 3*(b[1] + b[2] + b[3]) * u[1]
                    + (9*b[1] + 9*b[2] + 9*b[3] + b[4]) * u[2]
                    - 3*(3*b[1] + 3*b[2] + 3*b[3] + b[4]) * u[3]
                    + 3*(b[1] + b[2] + b[3] + b[4]) * u[4]
                    - b[4] * u[5]
                )
        dest[3] = α * inv_left_weights[3] * (
                    3*(b[1] + b[2] + b[3]) * u[1]
                    - 3*(3*b[1] + 3*b[2] + 3*b[3] + b[4]) * u[2]
                    + (9*b[1] + 9*b[2] + 9*b[3] + 9*b[4] + b[5]) * u[3]
                    - 3*(b[1] + b[2] + b[3] + 3*b[4] + b[5]) * u[4]
                    + 3*(b[4] + b[5]) * u[5]
                    - b[5] * u[6]
                )
        dest[4] = α * inv_left_weights[4] * (
                    - (b[1] + b[2] + b[3]) * u[1]
                    + 3*(b[1] + b[2] + b[3] + b[4]) * u[2]
                    - 3*(b[1] + b[2] + b[3] + 3*b[4] + b[5]) * u[3]
                    + (b[1] + b[2] + b[3] + 9*b[4] + 9*b[5] + b[6]) * u[4]
                    - 3*(b[4] + 3*b[5] + b[6]) * u[5]
                    + 3*(b[5] + b[6]) * u[6]
                    - b[6] * u[7]
                )
        for i in 5:length(inv_left_weights)
            dest[i] = α * inv_left_weights[i] * (
                        - b[i-1] * u[i-3]
                        + 3*(b[i-1] + b[i]) * u[i-2]
                        - 3*(b[i-1] + 3*b[i] + b[i+1]) * u[i-1]
                        + (b[i-1] + 9*b[i] + 9*b[i+1] + b[i+2]) * u[i]
                        - 3*(b[i] + 3*b[i+1] + b[i+2]) * u[i+1]
                        + 3*(b[i+1] + b[i+2]) * u[i+2]
                        - b[i+2] * u[i+3]
                    )
        end

        dest[end] = α * inv_right_weights[1] * (
                        (b[end] + b[end-1]) * u[end]
                        - 3*(b[end] + b[end-1]) * u[end-1]
                        + 3*(b[end] + b[end-1]) * u[end-2]
                        - (b[end] + b[end-1]) * u[end-3]
                    )
        dest[end-1] = α * inv_right_weights[2] * (
                        - 3*(b[end] + b[end-1]) * u[end]
                        + (b[end-2] + 9*b[end-1] + 9*b[end]) * u[end-1]
                        - 3*(b[end-2] + 3*b[end-1] + 3*b[end]) * u[end-2]
                        + 3*(b[end-2] + b[end-1] + b[end]) * u[end-3]
                        - b[end-2] * u[end-4]
                    )
        dest[end-2] = α * inv_right_weights[3] * (
                        3*(b[end] + b[end-1]) * u[end]
                        - 3*(b[end-2] + 3*b[end-1] + 3*b[end]) * u[end-1]
                        + (b[end-3] + 9*b[end-2] + 9*b[end-1] + 9*b[end]) * u[end-2]
                        - 3*(b[end-3] + 3*b[end-2] + b[end-1] + b[end]) * u[end-3]
                        + 3*(b[end-3] + b[end-2]) * u[end-4]
                        - b[end-3] * u[end-5]
                    )
        dest[end-3] = α * inv_right_weights[4] * (
                        - (b[end] + b[end-1]) * u[end]
                        + 3*(b[end-2] + b[end-1] + b[end]) * u[end-1]
                        - 3*(b[end-3] + 3*b[end-2] + b[end-1] + b[end]) * u[end-2]
                        + (b[end-4] + 9*b[end-3] + 9*b[end-2] + b[end-1] + b[end]) * u[end-3]
                        - 3*(b[end-4] + 3*b[end-3] + b[end-2]) * u[end-4]
                        + 3*(b[end-4] + b[end-3]) * u[end-5]
                        - b[end-4] * u[end-6]
                    )
        for i in 4:length(inv_right_weights)-1
            dest[end-i] = α * inv_right_weights[i+1] * (
                        - b[end-(i+1)] * u[end-(i+3)]
                        + 3*(b[end-(i+1)] + b[end-i]) * u[end-(i+2)]
                        - 3*(b[end-(i+1)] + 3*b[end-i] + b[end-(i-1)]) * u[end-(i+1)]
                        + (b[end-(i+1)] + 9*b[end-i] + 9*b[end-(i-1)] + b[end-(i-2)]) * u[end-i]
                        - 3*(b[end-i] + 3*b[end-(i-1)] + b[end-(i-2)]) * u[end-(i-1)]
                        + 3*(b[end-(i-1)] + b[end-(i-2)]) * u[end-(i-2)]
                        - b[end-(i-2)] * u[end-(i-3)]
                    )
        end
    end
end

@inline function convolve_interior_coefficients_loopbody(i, cache::MattssonSvärdNordström2004Cache6, u, b)
    @inbounds begin
        b_im1 = b[i-1]
        b_i   = b[i]
        b_ip1 = b[i+1]
        b_ip2 = b[i+2]

        retval = (
                    -b_im1 * u[i-3]
                    + 3*(b_im1 + b_i) * u[i-2]
                    - 3*(b_im1 + 3*b_i + b_ip1) * u[i-1]
                    + (b_im1 + 9*b_i + 9*b_ip1 + b_ip2) * u[i]
                    - 3*(b_i + 3*b_ip1 + b_ip2) * u[i+1]
                    + 3*(b_ip1 + b_ip2) * u[i+2]
                    - b_ip2 * u[i+3]
                )
    end

    retval
end



struct MattssonSvärdNordström2004Cache8{T,LeftWidth,RightWidth} <: AbstractCoefficientCache{T}
    inv_left_weights::SVector{LeftWidth,T}
    inv_right_weights::SVector{RightWidth,T}

    function MattssonSvärdNordström2004Cache8(inv_left_weights::SVector{LeftWidth,T}, inv_right_weights::SVector{RightWidth,T}) where {T,LeftWidth,RightWidth}
        boundary_length = 5

        if length(inv_left_weights) < boundary_length
            inv_left_weights = SVector(inv_left_weights..., ntuple(j->one(T), boundary_length-length(inv_left_weights))...)
        end
        if length(inv_right_weights) < boundary_length
            inv_right_weights = SVector(inv_right_weights..., ntuple(j->one(T), boundary_length-length(inv_right_weights))...)
        end

        new{T,length(inv_left_weights),length(inv_right_weights)}(inv_left_weights, inv_right_weights)
    end
end

function convolve_boundary_coefficients!(dest::AbstractVector, cache::MattssonSvärdNordström2004Cache8, u::AbstractVector, b::AbstractVector, α)
    @unpack inv_left_weights, inv_right_weights = cache

    @inbounds begin
        dest[1] = α * inv_left_weights[1] * (
                    (b[1] + b[2] + b[3]) * u[1]
                    - 4*(b[1] + b[2] + b[3]) * u[2]
                    + 6*(b[1] + b[2] + b[3]) * u[3]
                    - 4*(b[1] + b[2] + b[3]) * u[4]
                    + (b[1] + b[2] + b[3]) * u[5]
                )
        dest[2] = α * inv_left_weights[2] * (
                    - 4*(b[1] + b[2] + b[3]) * u[1]
                    + (16*b[1] + 16*b[2] + 16*b[3] + b[4]) * u[2]
                    - 4*(6*b[1] + 6*b[2] + 6*b[3] + b[4]) * u[3]
                    + (16*b[1] + 16*b[2] + 16*b[3] + 6*b[4]) * u[4]
                    - 4*(b[1] + b[2] + b[3] + b[4]) * u[5]
                    + b[4] * u[6]
                )
        dest[3] = α * inv_left_weights[3] * (
                    6*(b[1] + b[2] + b[3]) * u[1]
                    - 4*(6*b[1] + 6*b[2] + 6*b[3] + b[4]) * u[2]
                    + (36*b[1] + 36*b[2] + 36*b[3] + 16*b[4] + b[5]) * u[3]
                    - 4*(6*b[1] + 6*b[2] + 6*b[3] + 6*b[4] + b[5]) * u[4]
                    + (6*b[1] + 6*b[2] + 6*b[3] + 16*b[4] + 6*b[5]) * u[5]
                    - 4*(b[4] + b[5]) * u[6]
                    + b[5] * u[7]
                )
        dest[4] = α * inv_left_weights[4] * (
                    - 4*(b[1] + b[2] + b[3]) * u[1]
                    + (16*b[1] + 16*b[2] + 16*b[3] + 6*b[4]) * u[2]
                    - 4*(6*b[1] + 6*b[2] + 6*b[3] + 6*b[4] + b[5]) * u[3]
                    + (16*b[1] + 16*b[2] + 16*b[3] + 36*b[4] + 16*b[5] + b[6]) * u[4]
                    - 4*(b[1] + b[2] + b[3] + 6*b[4] + 6*b[5] + b[6]) * u[5]
                    + (6*b[4] + 16*b[5] + 6*b[6]) * u[6]
                    - 4*(b[5] + b[6]) * u[7]
                    + b[6] * u[8]
                )
        dest[5] = α * inv_left_weights[5] * (
                    (b[1] + b[2] + b[3]) * u[1]
                    - 4*(b[1] + b[2] + b[3] + b[4]) * u[2]
                    + (6*b[1] + 6*b[2] + 6*b[3] + 16*b[4] + 6*b[5]) * u[3]
                    - 4*(b[1] + b[2] + b[3] + 6*b[4] + 6*b[5] + b[6]) * u[4]
                    + (b[1] + b[2] + b[3] + 16*b[4] + 36*b[5] + 16*b[6] + b[7]) * u[5]
                    - 4*(b[4] + 6*b[5] + 6*b[6] + b[7]) * u[6]
                    + (6*b[5] + 16*b[6] + 6*b[7]) * u[7]
                    - 4*(b[6] + b[7]) * u[8]
                    + b[7] * u[9]
                )
        for i in 6:length(inv_left_weights)
            dest[i] = α * inv_left_weights[i] * (
                        b[i-2] * u[i-4]
                        - 4*(b[i-2] + b[i-1]) * u[i-3]
                        + (6*b[i-2] + 16*b[i-1] + 6*b[i]) * u[i-2]
                        - 4*(b[i-2] + 6*b[i-1] + 6*b[i] + b[i+1]) * u[i-1]
                        + (b[i-2] + 16*b[i-1] + 36*b[i] + 16*b[i+1] + b[i+2]) * u[i]
                        - 4*(b[i-1] + 6*b[i] + 6*b[i+1] + b[i+2]) * u[i+1]
                        + (6*b[i] + 16*b[i+1] + 6*b[i+2]) * u[i+2]
                        - 4*(b[i+1] + b[i+2]) * u[i+3]
                        + b[i+2] * u[i+4]
                    )
        end

        dest[end] = α * inv_right_weights[1] * (
                        (b[end-2] + b[end-1] + b[end]) * u[end]
                        - 4*(b[end-2] + b[end-1] + b[end]) * u[end-1]
                        + 6*(b[end-2] + b[end-1] + b[end]) * u[end-2]
                        - 4*(b[end-2] + b[end-1] + b[end]) * u[end-3]
                        + (b[end-2] + b[end-1] + b[end]) * u[end-4]
                    )
        dest[end-1] = α * inv_right_weights[2] * (
                        - 4*(b[end-2] + b[end-1] + b[end]) * u[end]
                        + (b[end-3] + 16*b[end-2] + 16*b[end-1] + 16*b[end]) * u[end-1]
                        - 4*(b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-2]
                        + (6*b[end-3] + 16*b[end-2] + 16*b[end-1] + 16*b[end]) * u[end-3]
                        - 4*(b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-4]
                        + b[end-3] * u[end-5]
                    )
        dest[end-2] = α * inv_right_weights[3] * (
                        6*(b[end-2] + b[end-1] + b[end]) * u[end]
                        - 4*(b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-1]
                        + (b[end-4] + 16*b[end-3] + 36*b[end-2] + 36*b[end-1] + 36*b[end]) * u[end-2]
                        - 4*(b[end-4] + 6*b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-3]
                        + (6*b[end-4] + 16*b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-4]
                        - 4*(b[end-4] + b[end-3]) * u[end-5]
                        + b[end-4] * u[end-6]
                    )
        dest[end-3] = α * inv_right_weights[4] * (
                        - 4*(b[end-2] + b[end-1] + b[end]) * u[end]
                        + (6*b[end-3] + 16*b[end-2] + 16*b[end-1] + 16*b[end]) * u[end-1]
                        - 4*(b[end-4] + 6*b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-2]
                        + (b[end-5] + 16*b[end-4] + 36*b[end-3] + 16*b[end-2] + 16*b[end-1] + 16*b[end]) * u[end-3]
                        - 4*(b[end-5] + 6*b[end-4] + 6*b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-4]
                        + (6*b[end-5] + 16*b[end-4] + 6*b[end-3]) * u[end-5]
                        - 4*(b[end-5] + b[end-4]) * u[end-6]
                        + b[end-5] * u[end-7]
                    )
        dest[end-4] = α * inv_right_weights[5] * (
                        (b[end-2] + b[end-1] + b[end]) * u[end]
                        - 4*(b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-1]
                        + (6*b[end-4] + 16*b[end-3] + 6*b[end-2] + 6*b[end-1] + 6*b[end]) * u[end-2]
                        - 4*(b[end-5] + 6*b[end-4] + 6*b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-3]
                        + (b[end-6] + 16*b[end-5] + 36*b[end-4] + 16*b[end-3] + b[end-2] + b[end-1] + b[end]) * u[end-4]
                        - 4*(b[end-6] + 6*b[end-5] + 6*b[end-4] + b[end-3]) * u[end-5]
                        + (6*b[end-6] + 16*b[end-5] + 6*b[end-4]) * u[end-6]
                        - 4*(b[end-6] + b[end-5]) * u[end-7]
                        + b[end-6] * u[end-8]
                    )
        for i in 5:length(inv_right_weights)-1
            dest[end-i] = α * inv_right_weights[i+1] * (
                        b[end-(i+2)] * u[end-(i+4)]
                        - 4*(b[end-(i+2)] + b[end-(i+1)]) * u[end-(i+3)]
                        + (6*b[end-(i+2)] + 16*b[end-(i+1)] + 6*b[end-i]) * u[end-(i+2)]
                        - 4*(b[end-(i+2)] + 6*b[end-(i+1)] + 6*b[end-i] + b[end-(i-1)]) * u[end-(i+1)]
                        + (b[end-(i+2)] + 16*b[end-(i+1)] + 36*b[end-i] + 16*b[end-(i-1)] + b[end-(i-2)]) * u[end-i]
                        - 4*(b[end-(i+1)] + 6*b[end-i] + 6*b[end-(i-1)] + b[end-(i-2)]) * u[end-(i-1)]
                        + (6*b[end-i] + 16*b[end-(i-1)] + 6*b[end-(i-2)]) * u[end-(i-2)]
                        - 4*(b[end-(i-1)] + b[end-(i-2)]) * u[end-(i-3)]
                        + b[end-(i-2)] * u[end-(i-4)]
                    )
        end
    end
end

@inline function convolve_interior_coefficients_loopbody(i, cache::MattssonSvärdNordström2004Cache8, u, b)
    @inbounds begin
        b_im2 = b[i-2]
        b_im1 = b[i-1]
        b_i   = b[i]
        b_ip1 = b[i+1]
        b_ip2 = b[i+2]

        retval = (
                    b_im2 * u[i-4]
                    - 4*(b_im2 + b_im1) * u[i-3]
                    + (6*b_im2 + 16*b_im1 + 6*b_i) * u[i-2]
                    - 4*(b_im2 + 6*b_im1 + 6*b_i + b_ip1) * u[i-1]
                    + (b_im2 + 16*b_im1 + 36*b_i + 16*b_ip1 + b_ip2) * u[i]
                    - 4*(b_im1 + 6*b_i + 6*b_ip1 + b_ip2) * u[i+1]
                    + (6*b_i + 16*b_ip1 + 6*b_ip2) * u[i+2]
                    - 4*(b_ip1 + b_ip2) * u[i+3]
                    + b_ip2 * u[i+4]
                )
    end

    retval
end
