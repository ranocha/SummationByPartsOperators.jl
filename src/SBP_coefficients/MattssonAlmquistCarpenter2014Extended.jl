
"""
    MattssonAlmquistCarpenter2014Extended

Coefficients of the extended SBP operators given in
  Mattsson, Almquist, Carpenter (2014)
  Optimal diagonal-norm SBP operators.
  Journal of Computational Physics 264, pp. 91-111.
"""
struct MattssonAlmquistCarpenter2014Extended <: SourceOfCoefficients end

function Base.show(io::IO, ::MattssonAlmquistCarpenter2014Extended)
    print(io,
        "  Mattsson, Almquist, Carpenter (2014) \n",
        "  Optimal diagonal-norm SBP operators. \n",
        "  Journal of Computational Physics 264, pp. 91-111. \n")
end


function first_derivative_coefficients(source::MattssonAlmquistCarpenter2014Extended, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,3}(SVector(T(-6//5),
                                                    T(7//5),
                                                    T(-1//5) )),
            # d2
            DerivativeCoefficientRow{T,1,3}(SVector(T(-1//2),
                                                    T(0),
                                                    T(1//2) )),
            # d3
            DerivativeCoefficientRow{T,1,4}(SVector(T(1//11),
                                                    T(-7//11),
                                                    T(0),
                                                    T(6//11) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(1//2))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(5//12),
                                T(7//6),
                                T(11//12) )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 1, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end
