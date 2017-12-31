
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
    elseif order == 4
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,6}(SVector(T(-800//511),
                                                    T(1744244//869211),
                                                    T(-312391//1738422),
                                                    T(-100109//289737),
                                                    T(54107//1738422),
                                                    T(6575//124173) )),
            # d2
            DerivativeCoefficientRow{T,1,6}(SVector(T(-1744244//3752595),
                                                    T(0),
                                                    T(238831//750519),
                                                    T(290357//1501038),
                                                    T(-625//35739),
                                                    T(-220357//7505190) )),
            # d3
            DerivativeCoefficientRow{T,1,6}(SVector(T(312391//3511620),
                                                    T(-238831//351162),
                                                    T(0),
                                                    T(17335//25083),
                                                    T(-74827//702324),
                                                    T(3859//585270) )),
            # d4
            DerivativeCoefficientRow{T,1,6}(SVector(T(100109//1109430),
                                                    T(-290357//1331316),
                                                    T(-17335//47547),
                                                    T(0),
                                                    T(342799//665658),
                                                    T(-149959//6656580) )),
            # d5
            DerivativeCoefficientRow{T,1,7}(SVector(T(-54107//5025510),
                                                    T(625//23931),
                                                    T(74827//1005102),
                                                    T(-342799//502551),
                                                    T(0),
                                                    T(1715156//2512755),
                                                    T(-240//2659) )),
            # d6
            DerivativeCoefficientRow{T,1,8}(SVector(T(-6575//392877),
                                                    T(220357//5500278),
                                                    T(-3859//916713),
                                                    T(149959//5500278),
                                                    T(-1715156//2750139),
                                                    T(0),
                                                    T(9600//14551),
                                                    T(-1200//14551) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(2//3), T(-1//12))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(511//1600),
                                T(3971//2880),
                                T(929//1440),
                                T(587//480),
                                T(2659//2880),
                                T(14551//14400) )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 1, order, source)
    elseif order == 6
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,3}(SVector(T(),
                                                    T(),
                                                    T() )),
            # d2
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(3//4), T(-3//20), T(1//60))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(),
                                T(),
                                T() )
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
