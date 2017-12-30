
"""
    Mattsson2014

Coefficients of the SBP operators given in
  Mattsson (2014)
  Diagonal-norm summation by parts operators for fiite difference approximations
    of third and fourth derivatives.
  Journal of Computational Physics 274, pp. 432-454.
"""
struct Mattsson2014 <: SourceOfCoefficients end

function Base.show(io::IO, ::Mattsson2014)
    print(io,
        "  Mattsson (2014) \n",
        "  Diagonal-norm summation by parts operators for fiite difference approximations \n",
        "    of third and fourth derivatives. \n",
        "  Journal of Computational Physics 274, pp. 432-454. \n")
end


function first_derivative_coefficients(source::Mattsson2014, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,2}(SVector(T(-1),
                                                    T(1) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(1//2))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector(T(1//2))
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
            DerivativeCoefficientRow{T,1,6}(SVector(T(-50400//35809),
                                                    T(526249//322281),
                                                    T(-75733//1933686),
                                                    T(-50767//322281),
                                                    T(-4417//71618),
                                                    T(31850//966843) )),
            # d2
            DerivativeCoefficientRow{T,1,6}(SVector(T(-526249//1077057),
                                                    T(0),
                                                    T(1421209//3231171),
                                                    T(16657//239346),
                                                    T(-16934//1077057),
                                                    T(-33059//6462342) )),
            # d3
            DerivativeCoefficientRow{T,1,6}(SVector(T(75733//5541372),
                                                    T(-1421209//2770686),
                                                    T(0),
                                                    T(631187//1385343),
                                                    T(400139//5541372),
                                                    T(-8789//307854) )),
            # d4
            DerivativeCoefficientRow{T,1,6}(SVector(T(50767//811962),
                                                    T(-16657//180436),
                                                    T(-631187//1217943),
                                                    T(0),
                                                    T(496403//811962),
                                                    T(-308533//4871772) )),
            # d5
            DerivativeCoefficientRow{T,1,7}(SVector(T(4417//211146),
                                                    T(16934//950157),
                                                    T(-400139//5700942),
                                                    T(-496403//950157),
                                                    T(0),
                                                    T(1805647//2850471),
                                                    T(-2800//35191) )),
            # d6
            DerivativeCoefficientRow{T,1,8}(SVector(T(-31850//2713743),
                                                    T(33059//5427486),
                                                    T(8789//301527),
                                                    T(308533//5427486),
                                                    T(-1805647//2713743),
                                                    T(0),
                                                    T(22400//33503),
                                                    T(-2800//33503) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(2//3), T(-1//12))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(35809//100800),
                                T(13297//11200),
                                T(5701//5600),
                                T(45109//50400),
                                T(35191//33600),
                                T(33503//33600) )
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
            DerivativeCoefficientRow{T,1,8}(SVector(T(-508032//318365),
                                                    T(113221347//55519750),
                                                    T(-3338172//27759875),
                                                    T(-1002751721//2731571700),
                                                    T(-46815789//455261950),
                                                    T(1228638//9105239),
                                                    T(7120579//166559250),
                                                    T(-10874619//350201500) )),
            # d2
            DerivativeCoefficientRow{T,1,8}(SVector(T(-4642075227//10228748530),
                                                    T(0),
                                                    T(69095487//265681780),
                                                    T(32352081//146124979),
                                                    T(28592150//438374937),
                                                    T(-137946309//1461249790),
                                                    T(-720387//44961532),
                                                    T(538846039//30686245590) )),
            # d3
            DerivativeCoefficientRow{T,1,8}(SVector(T(45621684//696580885),
                                                    T(-23031829//36186020),
                                                    T(0),
                                                    T(28368209//39804622),
                                                    T(-9693137//39804622),
                                                    T(3868089//30618940),
                                                    T(-2468403//199023110),
                                                    T(-1686470//139316177) )),
            # d4
            DerivativeCoefficientRow{T,1,8}(SVector(T(1002751721//11591207628),
                                                    T(-32352081//137990567),
                                                    T(-85104627//275981134),
                                                    T(0),
                                                    T(17499453//42458636),
                                                    T(13059537//275981134),
                                                    T(4924918//413971701),
                                                    T(-29088585//1931867938) )),
            # d5
            DerivativeCoefficientRow{T,1,8}(SVector(T(46815789//1188140954),
                                                    T(-28592150//254601633),
                                                    T(29079411//169734422),
                                                    T(-17499453//26112988),
                                                    T(0),
                                                    T(112822635//169734422),
                                                    T(-9119079//84867211),
                                                    T(103152839//7128845724) )),
            # d6
            DerivativeCoefficientRow{T,1,9}(SVector(T(-2047730//52061009),
                                                    T(45982103//371864350),
                                                    T(-3868089//57209900),
                                                    T(-4353179//74372870),
                                                    T(-7521509//14874574),
                                                    T(0),
                                                    T(474569879//743728700),
                                                    T(-138109864//1301525225),
                                                    T(4032//260045) )),
            # d7
            DerivativeCoefficientRow{T,1,10}(SVector(T(-291943739//21305233950),
                                                    T(720387//31216460),
                                                    T(7405209//1014534950),
                                                    T(-4924918//304360485),
                                                    T(9119079//101453495),
                                                    T(-1423709637//2029069900),
                                                    T(0),
                                                    T(5309800707//7101744650),
                                                    T(-108864//709465),
                                                    T(12096//709465) )),
            # d8
            DerivativeCoefficientRow{T,1,11}(SVector(T(10874619//1121684300),
                                                    T(-538846039//21872843850),
                                                    T(1011882//145818959),
                                                    T(5817717//291637918),
                                                    T(-103152839//8749137540),
                                                    T(414329592//3645473975),
                                                    T(-5309800707//7290947950),
                                                    T(0),
                                                    T(762048//1019713),
                                                    T(-762048//5098565),
                                                    T(84672//5098565) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(3//4), T(-3//20), T(1//60))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(318365//1016064),
                                T(145979//103680),
                                T(139177//241920),
                                T(964969//725760),
                                T(593477//725760),
                                T(52009//48384),
                                T(141893//145152),
                                T(1019713//1016064) )
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


function second_derivative_coefficients(source::Mattsson2014, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,3}(SVector(T(1),
                                                    T(-2),
                                                    T(1) )),
        )
        right_boundary = left_boundary
        upper_coef = SVector(T(1))
        central_coef = T(-2)
        lower_coef = upper_coef
        left_weights = SVector(T(1//2))
        right_weights = left_weights
        left_boundary_derivatives = (
            # first derivative
            DerivativeCoefficientRow{T,1,3}(SVector(T(-3//2),
                                                    T(2),
                                                    T(-1//2) )),
        )
        right_boundary_derivatives = (
            -left_boundary_derivatives[1],
        )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 2, order, source)
    elseif order == 4
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,6}(SVector(T(8027765//3867372),
                                                    T(-1690702//322281),
                                                    T(8240267//1933686),
                                                    T(-1030334//966843),
                                                    T(-6817//143236),
                                                    T(21380//966843) )),
            # d2
            DerivativeCoefficientRow{T,1,6}(SVector(T(1030898//1077057),
                                                    T(-23780195//12924684),
                                                    T(2540018//3231171),
                                                    T(26257//239346),
                                                    T(-12268//3231171),
                                                    T(-119459//12924684) )),
            # d3
            DerivativeCoefficientRow{T,1,6}(SVector(T(75467//5541372),
                                                    T(1270009//1385343),
                                                    T(-558115//307854),
                                                    T(1111174//1385343),
                                                    T(551339//5541372),
                                                    T(-8789//461781) )),
            # d4
            DerivativeCoefficientRow{T,1,6}(SVector(T(-61567//1217943),
                                                    T(26257//180436),
                                                    T(1111174//1217943),
                                                    T(-5126635//2435886),
                                                    T(464003//405981),
                                                    T(-222133//4871772) )),
            # d5
            DerivativeCoefficientRow{T,1,7}(SVector(T(-6817//422292),
                                                    T(-12268//2850471),
                                                    T(551339//5700942),
                                                    T(928006//950157),
                                                    T(-25370795//11401884),
                                                    T(3568094//2850471),
                                                    T(-2800//35191) )),
            # d6
            DerivativeCoefficientRow{T,1,8}(SVector(T(21380//2713743),
                                                    T(-119459//10854972),
                                                    T(-17578//904581),
                                                    T(-222133//5427486),
                                                    T(3568094//2713743),
                                                    T(-9063745//3618324),
                                                    T(44800//33503),
                                                    T(-2800//33503) )),
        )
        right_boundary = left_boundary
        upper_coef = SVector(T(4//3), T(-1//12))
        central_coef = T(-5//2)
        lower_coef = upper_coef
        left_weights = SVector( T(35809//100800),
                                T(13297//11200),
                                T(5701//5600),
                                T(45109//50400),
                                T(35191//33600),
                                T(33503//33600) )
        right_weights = left_weights
        left_boundary_derivatives = (
            # first derivative
            DerivativeCoefficientRow{T,1,4}(SVector(T(-11/6),
                                                    T(3),
                                                    T(-3//2),
                                                    T(1//3) )),
        )
        right_boundary_derivatives = (
            -left_boundary_derivatives[1],
        )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 2, order, source)
    elseif order == 6
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d2
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d3
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d4
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d5
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d6
            DerivativeCoefficientRow{T,1,9}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d7
            DerivativeCoefficientRow{T,1,10}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d8
            DerivativeCoefficientRow{T,1,11}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
        )
        right_boundary = left_boundary
        upper_coef = SVector(T(3//2), T(-3//20), T(1//90))
        central_coef = T(-49//18)
        lower_coef = upper_coef
        left_weights = SVector( T(318365//1016064),
                                T(145979//103680),
                                T(139177//241920),
                                T(964969//725760),
                                T(593477//725760),
                                T(52009//48384),
                                T(141893//145152),
                                T(1019713//1016064) )
        right_weights = left_weights
        left_boundary_derivatives = (
            # first derivative
            DerivativeCoefficientRow{T,1,5}(SVector(T(-25/12),
                                                    T(4),
                                                    T(-3),
                                                    T(4//3),
                                                    T(-1//4) )),
        )
        right_boundary_derivatives = (
            -left_boundary_derivatives[1],
        )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 2, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end


function third_derivative_coefficients(source::Mattsson2014, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,4}(SVector(T(1//4),
                                                    T(-5//8),
                                                    T(1//2),
                                                    T(-1//8) )),
            # d2
            DerivativeCoefficientRow{T,1,4}(SVector(T(-11//16),
                                                    T(2),
                                                    T(-31//16),
                                                    T(5//8) )),
            # d3
            DerivativeCoefficientRow{T,1,5}(SVector(T(-1//2),
                                                    T(15//16),
                                                    T(1//8),
                                                    T(-17//16),
                                                    T(1//2) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(-1), T(1//2))
        central_coef = T(0)
        lower_coef = -upper_coef
        left_weights = SVector(T(1//2))
        right_weights = left_weights
        left_boundary_derivatives = (
            # first derivative
            DerivativeCoefficientRow{T,1,3}(SVector(T(-3//2),
                                                    T(2),
                                                    T(-1//2) )),
            # second derivative
            DerivativeCoefficientRow{T,1,3}(SVector(T(1),
                                                    T(-2),
                                                    T(1) )),
        )
        right_boundary_derivatives = (
            -left_boundary_derivatives[1],
            left_boundary_derivatives[2]
        )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 3, order, source)
    elseif order == 4
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,6}(SVector(T(-32200//35809),
                                                    T(188187//71618),
                                                    T(-90183//35809),
                                                    T(27988//35809),
                                                    T(-801//35809),
                                                    T(2205//71618) )),
            # d2
            DerivativeCoefficientRow{T,1,6}(SVector(T(-96329//79782),
                                                    T(50400//13297),
                                                    T(-163046//39891),
                                                    T(63583//39891),
                                                    T(-1337//26594),
                                                    T(-1567//39891) )),
            # d3
            DerivativeCoefficientRow{T,1,6}(SVector(T(-11939//34206),
                                                    T(5923//17103),
                                                    T(6300//5701),
                                                    T(-32543//17103),
                                                    T(29083//34206),
                                                    T(-284//5701) )),
            # d4
            DerivativeCoefficientRow{T,1,7}(SVector(T(5606//45109),
                                                    T(-89949//90218),
                                                    T(72429//45109),
                                                    T(2800//45109),
                                                    T(-77319//45109),
                                                    T(95517//90218),
                                                    T(-6300//45109) )),
            # d5
            DerivativeCoefficientRow{T,1,8}(SVector(T(267//35191),
                                                    T(4011//70382),
                                                    T(-29083//35191),
                                                    T(51546//35191),
                                                    T(0),
                                                    T(-108271//70382),
                                                    T(33600//35191),
                                                    T(-4200//35191) )),
            # d6
            DerivativeCoefficientRow{T,1,9}(SVector(T(-735//67006),
                                                    T(1567//33503),
                                                    T(1704//33503),
                                                    T(-31839//33503),
                                                    T(108271//67006),
                                                    T(0),
                                                    T(-54600//33503),
                                                    T(33600//33503),
                                                    T(-4200//33503) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(-13//8), T(1), T(-1//8))
        central_coef = T(0)
        lower_coef = -upper_coef
        left_weights = SVector( T(35809//100800),
                                T(13297//11200),
                                T(5701//5600),
                                T(45109//50400),
                                T(35191//33600),
                                T(33503//33600) )
        right_weights = left_weights
        left_boundary_derivatives = (
            # first derivative
            DerivativeCoefficientRow{T,1,4}(SVector(T(-11/6),
                                                    T(3),
                                                    T(-3//2),
                                                    T(1//3) )),
            # second derivative
            DerivativeCoefficientRow{T,1,4}(SVector(T(2),
                                                    T(-5),
                                                    T(4),
                                                    T(-1) )),
        )
        right_boundary_derivatives = (
            -left_boundary_derivatives[1],
            left_boundary_derivatives[2]
        )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 3, order, source)
    elseif order == 6
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d2
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d3
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d4
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d5
            DerivativeCoefficientRow{T,1,9}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d6
            DerivativeCoefficientRow{T,1,10}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d7
            DerivativeCoefficientRow{T,1,11}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d8
            DerivativeCoefficientRow{T,1,12}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(-61//30), T(169//120), T(-3//10), T(7//240))
        central_coef = T(0)
        lower_coef = -upper_coef
        left_weights = SVector( T(318365//1016064),
                                T(145979//103680),
                                T(139177//241920),
                                T(964969//725760),
                                T(593477//725760),
                                T(52009//48384),
                                T(141893//145152),
                                T(1019713//1016064) )
        right_weights = left_weights
        left_boundary_derivatives = (
            # first derivative
            DerivativeCoefficientRow{T,1,5}(SVector(T(-25//12),
                                                    T(4),
                                                    T(-3),
                                                    T(4//3),
                                                    T(-1//4) )),
            # second derivative
            DerivativeCoefficientRow{T,1,5}(SVector(T(35//12),
                                                    T(-26//3),
                                                    T(19//2),
                                                    T(-14//3),
                                                    T(11//12) )),
        )
        right_boundary_derivatives = (
            -left_boundary_derivatives[1],
            left_boundary_derivatives[2]
        )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 3, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end


function fourth_derivative_coefficients(source::Mattsson2014, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,4}(SVector(T(8//5),
                                                    T(-24//5),
                                                    T(24//5),
                                                    T(-8//5) )),
            # d2
            DerivativeCoefficientRow{T,1,4}(SVector(T(-2//5),
                                                    T(6//5),
                                                    T(-6//5),
                                                    T(2//5) )),
            # d3
            DerivativeCoefficientRow{T,1,5}(SVector(T(2//5),
                                                    T(-11//5),
                                                    T(21//5),
                                                    T(-17//5),
                                                    T(1) )),
            # d4
            DerivativeCoefficientRow{T,1,6}(SVector(T(1//5),
                                                    T(2//5),
                                                    T(-17//5),
                                                    T(-29//5),
                                                    T(-4),
                                                    T(1) )),
        )
        right_boundary = left_boundary
        upper_coef = SVector(T(-4), T(1))
        central_coef = T(6)
        lower_coef = upper_coef
        left_weights = SVector(T(1//2))
        right_weights = left_weights
        left_boundary_derivatives = (
            # first derivative
            DerivativeCoefficientRow{T,1,3}(SVector(T(-3//2),
                                                    T(2),
                                                    T(-1//2) )),
            # second derivative
            DerivativeCoefficientRow{T,1,3}(SVector(T(1),
                                                    T(-2),
                                                    T(1) )),
            # third derivative
            DerivativeCoefficientRow{T,1,4}(SVector(T(-1),
                                                    T(3),
                                                    T(-3),
                                                    T(1) )),
        )
        right_boundary_derivatives = (
            -left_boundary_derivatives[1],
            left_boundary_derivatives[2],
            -left_boundary_derivatives[3],
        )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 4, order, source)
    elseif order == 4
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,6}(SVector(T(-242219//644562),
                                                    T(881057//644562),
                                                    T(-183673//107427),
                                                    T(220981//322281),
                                                    T(109057//644562),
                                                    T(-29273//214854) )),
            # d2
            DerivativeCoefficientRow{T,1,6}(SVector(T(578657//2154114),
                                                    T(-703457//718038),
                                                    T(1327457//1077057),
                                                    T(-544543//1077057),
                                                    T(-79457//718038),
                                                    T(204257//2154114) )),
            # d3
            DerivativeCoefficientRow{T,1,6}(SVector(T(219527//307854),
                                                    T(-2754943//923562),
                                                    T(2216981//461781),
                                                    T(-559673//153927),
                                                    T(1141057//923562),
                                                    T(-120619//923562) )),
            # d4
            DerivativeCoefficientRow{T,1,7}(SVector(T(69781//811962),
                                                    T(665057//811962),
                                                    T(-584873//135327),
                                                    T(2995381//405981),
                                                    T(-4614143//811962),
                                                    T(172109//90218),
                                                    T(-8400//45109) )),
            # d5
            DerivativeCoefficientRow{T,1,8}(SVector(T(8389//146178),
                                                    T(-79457//633438),
                                                    T(1141057//950157),
                                                    T(-4614143//950157),
                                                    T(557127//70382),
                                                    T(-11293343//1900314),
                                                    T(67200//35191),
                                                    T(-5600//35191) )),
            # d6
            DerivativeCoefficientRow{T,1,9}(SVector(T(-29273//603054),
                                                    T(204257//1809162),
                                                    T(-120619//904581),
                                                    T(172109//100509),
                                                    T(-11293343//1809162),
                                                    T(16787381//1809162),
                                                    T(-218400//33503),
                                                    T(67200//33503),
                                                    T(-5600//33503) )),
        )
        right_boundary = left_boundary
        upper_coef = SVector(T(-13//2), T(2), T(-1//6))
        central_coef = T(28//3)
        lower_coef = upper_coef
        left_weights = SVector( T(35809//100800),
                                T(13297//11200),
                                T(5701//5600),
                                T(45109//50400),
                                T(35191//33600),
                                T(33503//33600) )
        right_weights = left_weights
        left_boundary_derivatives = (
            # first derivative
            DerivativeCoefficientRow{T,1,4}(SVector(T(-11/6),
                                                    T(3),
                                                    T(-3//2),
                                                    T(1//3) )),
            # second derivative
            DerivativeCoefficientRow{T,1,4}(SVector(T(2),
                                                    T(-5),
                                                    T(4),
                                                    T(-1) )),
            # third derivative
            DerivativeCoefficientRow{T,1,4}(SVector(T(-1),
                                                    T(3),
                                                    T(-3),
                                                    T(1) )),
        )
        right_boundary_derivatives = (
            -left_boundary_derivatives[1],
            left_boundary_derivatives[2],
            -left_boundary_derivatives[3],
        )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 4, order, source)
    elseif order == 6
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d2
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d3
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d4
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d5
            DerivativeCoefficientRow{T,1,9}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d6
            DerivativeCoefficientRow{T,1,10}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d7
            DerivativeCoefficientRow{T,1,11}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d8
            DerivativeCoefficientRow{T,1,12}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T(),
                                                    T() )),
        )
        right_boundary = left_boundary
        upper_coef = SVector(T(-122//15), T(169//60), T(-2//5), T(7//240))
        central_coef = T(91//8)
        lower_coef = upper_coef
        left_weights = SVector( T(318365//1016064),
                                T(145979//103680),
                                T(139177//241920),
                                T(964969//725760),
                                T(593477//725760),
                                T(52009//48384),
                                T(141893//145152),
                                T(1019713//1016064) )
        right_weights = left_weights
        left_boundary_derivatives = (
        # first derivative
        DerivativeCoefficientRow{T,1,5}(SVector(T(-25//12),
                                                T(4),
                                                T(-3),
                                                T(4//3),
                                                T(-1//4) )),
        # second derivative
        DerivativeCoefficientRow{T,1,5}(SVector(T(35//12),
                                                T(-26//3),
                                                T(19//2),
                                                T(-14//3),
                                                T(11//12) )),
        # third derivative
        DerivativeCoefficientRow{T,1,5}(SVector(T(-5//2),
                                                T(9),
                                                T(-12),
                                                T(7),
                                                T(-3//2) )),
        )
        right_boundary_derivatives = (
            -left_boundary_derivatives[1],
            left_boundary_derivatives[2],
            -left_boundary_derivatives[3],
        )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 4, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end
