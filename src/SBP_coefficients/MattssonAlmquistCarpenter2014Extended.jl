
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
            DerivativeCoefficientRow{T,1,8}(SVector(T(-12700800//7497391),
                                                    T(84512094//37486955),
                                                    T(-12926907//149947820),
                                                    T(-1378838//2044743),
                                                    T(-2874//40093),
                                                    T(2857680//7497391),
                                                    T(-39279943//449843460),
                                                    T(-825552//37486955) )),
            # d2
            DerivativeCoefficientRow{T,1,8}(SVector(T(-84512094//193589725),
                                                    T(0),
                                                    T(56626691//464615340),
                                                    T(97579603//232307670),
                                                    T(3642717//77435890),
                                                    T(-246914711//1161538350),
                                                    T(23380423//464615340),
                                                    T(54432//5531135) )),
            # d3
            DerivativeCoefficientRow{T,1,8}(SVector(T(205189//2106340),
                                                    T(-56626691//79619652),
                                                    T(0),
                                                    T(969659//1153908),
                                                    T(-30244385//79619652),
                                                    T(9737351//44233140),
                                                    T(-508309//6634971),
                                                    T(780937//79619652) )),
            # d4
            DerivativeCoefficientRow{T,1,8}(SVector(T(15167218//136593975),
                                                    T(-97579603//273187950),
                                                    T(-22302157//182125300),
                                                    T(0),
                                                    T(8120941//54637590),
                                                    T(13871884//45531325),
                                                    T(-42123859//546375900),
                                                    T(-2075677//273187950) )),
            # d5
            DerivativeCoefficientRow{T,1,8}(SVector(T(537438//10627085),
                                                    T(-3642717//21254170),
                                                    T(6048877//25505004),
                                                    T(-8120941//12752502),
                                                    T(0),
                                                    T(29205107//63762510),
                                                    T(15263287//127525020),
                                                    T(-243381//4250834) )),
            # d6
            DerivativeCoefficientRow{T,1,9}(SVector(T(-45360//513973),
                                                    T(246914711//971408970),
                                                    T(-9737351//215868660),
                                                    T(-13871884//32380299),
                                                    T(-29205107//194281794),
                                                    T(0),
                                                    T(321139631//647605980),
                                                    T(-49594423//971408970),
                                                    T(6720//513973) )),
            # d7
            DerivativeCoefficientRow{T,1,10}(SVector(T(39279943//1409459100),
                                                    T(-23380423//281891820),
                                                    T(508309//23490985),
                                                    T(42123859//281891820),
                                                    T(-15263287//281891820),
                                                    T(-321139631//469819700),
                                                    T(0),
                                                    T(30841499//40270260),
                                                    T(-108864//671171),
                                                    T(12096//671171) )),
            # d8
            DerivativeCoefficientRow{T,1,11}(SVector(T(825552//128159995),
                                                    T(-381024//25631999),
                                                    T(-780937//307583988),
                                                    T(2075677//153791994),
                                                    T(1216905//51263998),
                                                    T(49594423//768959970),
                                                    T(-215890493//307583988),
                                                    T(0),
                                                    T(19051200//25631999),
                                                    T(-3810240//25631999),
                                                    T(423360//25631999) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(3//4), T(-3//20), T(1//60))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(7497391//25401600),
                                T(1106227//725760),
                                T(105317//403200),
                                T(260179//145152),
                                T(303631//725760),
                                T(513973//403200),
                                T(671171//725760),
                                T(25631999//25401600) )
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
