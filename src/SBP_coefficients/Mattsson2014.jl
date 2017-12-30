
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
        "  Diagonal-norm summation by parts operators for fiite difference approximations\n",
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
            DerivativeCoefficientRow{T,1,8}(SVector(T(24055498439//8194715100),
                                                    T(-9925742373//1138154875),
                                                    T(983223468//103468625),
                                                    T(-18528585641//4097357550),
                                                    T(16391451//22207900),
                                                    T(13725804//227630975),
                                                    T(355447739//20486787750),
                                                    T(-24731721//2276309750) )),
            # d2
            DerivativeCoefficientRow{T,1,8}(SVector(T(270821931//300845545),
                                                    T(-42416226217//26302496220),
                                                    T(685962357//1461249790),
                                                    T(2160993//8595587),
                                                    T(59905900//1315124811),
                                                    T(-158509509//2922499580),
                                                    T(-9667431//1461249790),
                                                    T(634102039//92058736770) )),
            # d3
            DerivativeCoefficientRow{T,1,8}(SVector(T(-2479644//63325535),
                                                    T(228654119//199023110),
                                                    T(-16197861//7237204),
                                                    T(24739409//19902311),
                                                    T(-7878737//39804622),
                                                    T(1917829//18093010),
                                                    T(-7508403//398046220),
                                                    T(-41036//10716629) )),
            # d4
            DerivativeCoefficientRow{T,1,8}(SVector(T(-1092927401//17386811442),
                                                    T(36736881//137990567),
                                                    T(74218227//137990567),
                                                    T(-7780367599//4967660412),
                                                    T(210256089//275981134),
                                                    T(1500627//25089194),
                                                    T(1246172//95531931),
                                                    T(-37555785//3863735876) )),
            # d5
            DerivativeCoefficientRow{T,1,8}(SVector(T(-54436269//2376281908),
                                                    T(59905900//763804899),
                                                    T(-23636211//169734422),
                                                    T(210256089//169734422),
                                                    T(-7116321131//3055219596),
                                                    T(8176215//6528247),
                                                    T(-7304679//84867211),
                                                    T(84101639//10693268586) )),
            # d6
            DerivativeCoefficientRow{T,1,9}(SVector(T(4575268//260305045),
                                                    T(-52836503//743728700),
                                                    T(1917829//33805850),
                                                    T(500209//6761170),
                                                    T(545081//572099),
                                                    T(-324760747//148745740),
                                                    T(461969879//371864350),
                                                    T(-11753624//118320475),
                                                    T(2688//260045) )),
            # d7
            DerivativeCoefficientRow{T,1,10}(SVector(T(355447739//63915701850),
                                                    T(-9667431//1014534950),
                                                    T(-22525209//2029069900),
                                                    T(1246172//70237035),
                                                    T(-7304679//101453495),
                                                    T(1385909637//1014534950),
                                                    T(-48284442317//18261629100),
                                                    T(5288632707//3550872325),
                                                    T(-108864//709465),
                                                    T(8064//709465) )),
            # d8
            DerivativeCoefficientRow{T,1,11}(SVector(T(-24731721//7290947950),
                                                    T(634102039//65618531550),
                                                    T(-123108//56084215),
                                                    T(-7511157//583275836),
                                                    T(84101639//13123706310),
                                                    T(-35260872//331406725),
                                                    T(5288632707//3645473975),
                                                    T(-70820489957//26247412620),
                                                    T(1524096//1019713),
                                                    T(-762048//5098565),
                                                    T(56448//5098565) )),
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
            DerivativeCoefficientRow{T,1,8}(SVector(T(-151704//63673),
                                                    T(15270676769//1821047800),
                                                    T(-443349971//41387450),
                                                    T(2063356637//364209560),
                                                    T(-39300617//45526195),
                                                    T(-11473393//364209560),
                                                    T(-38062741//455261950),
                                                    T(40315779//1821047800) )),
            # d2
            DerivativeCoefficientRow{T,1,8}(SVector(T(-13333381409//8182998824),
                                                    T(829440//145979),
                                                    T(-8702160983//1168999832),
                                                    T(1321219979//292249958),
                                                    T(-1463113021//1168999832),
                                                    T(1240729//20874997),
                                                    T(102110955//1168999832),
                                                    T(-50022767//2045749706) )),
            # d3
            DerivativeCoefficientRow{T,1,8}(SVector(T(14062931//75990642),
                                                    T(-1261072297//477655464),
                                                    T(1088640//139177),
                                                    T(-4530616889//477655464),
                                                    T(602572103//119413866),
                                                    T(-116503713//159218488),
                                                    T(-17846623//59706933),
                                                    T(343537955//3343588248) )),
            # d4
            DerivativeCoefficientRow{T,1,8}(SVector(T(661223855//7727471752),
                                                    T(-214194059//275981134),
                                                    T(1209539129//1103924536),
                                                    T(645120//964969),
                                                    T(-2321979501//1103924536),
                                                    T(327603877//275981134),
                                                    T(-175223717//1103924536),
                                                    T(1353613//965933969) )),
            # d5
            DerivativeCoefficientRow{T,1,9}(SVector(T(-91064195//594070477),
                                                    T(632843581//678937688),
                                                    T(-446896583//169734422),
                                                    T(2045223021//678937688),
                                                    T(22680//593477),
                                                    T(-1804641793//678937688),
                                                    T(311038417//169734422),
                                                    T(-1932566239//4752563816),
                                                    T(21168//593477) )),
            # d6
            DerivativeCoefficientRow{T,1,10}(SVector(T(11473393//1249464216),
                                                    T(-8685103//111559305),
                                                    T(116503713//297491480),
                                                    T(-327603877//223118610),
                                                    T(1804641793//892474440),
                                                    T(0),
                                                    T(-1760949511//892474440),
                                                    T(2105883973//1561830270),
                                                    T(-72576//260045),
                                                    T(7056//260045) )),
            # d7
            DerivativeCoefficientRow{T,1,11}(SVector(T(38062741//1420348930),
                                                    T(-20422191//162325592),
                                                    T(17846623//101453495),
                                                    T(175223717//811627960),
                                                    T(-311038417//202906990),
                                                    T(1760949511//811627960),
                                                    T(0),
                                                    T(-1081094773//516490520),
                                                    T(1022112//709465),
                                                    T(-217728//709465),
                                                    T(21168//709465) )),
            # d8
            DerivativeCoefficientRow{T,1,12}(SVector(T(-40315779//5832758360),
                                                    T(50022767//1458189590),
                                                    T(-68707591//1166551672),
                                                    T(-1353613//729094795),
                                                    T(1932566239//5832758360),
                                                    T(-2105883973//1458189590),
                                                    T(1081094773//530250760),
                                                    T(0),
                                                    T(-10329984//5098565),
                                                    T(7154784//5098565),
                                                    T(-1524096//5098565),
                                                    T(148176//5098565) )),
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
            DerivativeCoefficientRow{T,1,8}(SVector(T(37567391168//53948541075),
                                                    T(-95834307667//35965694050),
                                                    T(252350074//65392171),
                                                    T(-58232913019//21579416430),
                                                    T(4040770588//3596569405),
                                                    T(-15248255797//35965694050),
                                                    T(4832196698//53948541075),
                                                    T(134156001//7193138810) )),
            # d2
            DerivativeCoefficientRow{T,1,8}(SVector(T(29125918379//23087746682),
                                                    T(-1255810938848//242421340161),
                                                    T(1289206067431//161614226774),
                                                    T(-431078362378//80807113387),
                                                    T(494586219497//484842680322),
                                                    T(31446420748//80807113387),
                                                    T(-21701585799//161614226774),
                                                    T(334788562//242421340161) )),
            # d3
            DerivativeCoefficientRow{T,1,8}(SVector(T(1308658570//3001630359),
                                                    T(-88210933529//66035867898),
                                                    T(13622370452//11005977983),
                                                    T(-27138341627//66035867898),
                                                    T(23881355534//33017933949),
                                                    T(-26412188989//22011955966),
                                                    T(21399717536//33017933949),
                                                    T(-928716467//9433695414) )),
            # d4
            DerivativeCoefficientRow{T,1,8}(SVector(T(110582060185//457852701306),
                                                    T(-22954806538//76308783551),
                                                    T(-180184675067//152617567102),
                                                    T(678091654628//228926350653),
                                                    T(-378329435643//152617567102),
                                                    T(69519106966//76308783551),
                                                    T(-98928859751//457852701306),
                                                    T(4720003312//76308783551) )),
            # d5
            DerivativeCoefficientRow{T,1,9}(SVector(T(1870177580//46931567683),
                                                    T(-21945155863//281589406098),
                                                    T(45403496174//46931567683),
                                                    T(-384706366203//93863135366),
                                                    T(974238057544//140794703049),
                                                    T(-520477408939//93863135366),
                                                    T(99162460006//46931567683),
                                                    T(-99640101991//281589406098),
                                                    T(21168//593477) )),
            # d6
            DerivativeCoefficientRow{T,1,10}(SVector(T(-15248255797//123384591330),
                                                    T(31446420748//61692295665),
                                                    T(-26412188989//41128197110),
                                                    T(69519106966//61692295665),
                                                    T(-520477408939//123384591330),
                                                    T(155376599432//20564098555),
                                                    T(-772894368601//123384591330),
                                                    T(21159425698//8813185095),
                                                    T(-96768//260045),
                                                    T(7056//260045) )),
            # d7
            DerivativeCoefficientRow{T,1,11}(SVector(T(690313814//24044478315),
                                                    T(-21701585799//112207565470),
                                                    T(21399717536//56103782735),
                                                    T(-98928859751//336622696410),
                                                    T(99162460006//56103782735),
                                                    T(-772894368601//112207565470),
                                                    T(1826861184956//168311348205),
                                                    T(-915425403107//112207565470),
                                                    T(2044224//709465),
                                                    T(-290304//709465),
                                                    T(21168//709465) )),
            # d8
            DerivativeCoefficientRow{T,1,12}(SVector(T(134156001//23039395522),
                                                    T(334788562//172795466415),
                                                    T(-6501015269//115196977610),
                                                    T(4720003312//57598488805),
                                                    T(-99640101991//345590932830),
                                                    T(148115979886//57598488805),
                                                    T(-915425403107//115196977610),
                                                    T(1952118169516//172795466415),
                                                    T(-41319936//5098565),
                                                    T(14309568//5098565),
                                                    T(-2032128//5098565),
                                                    T(148176//5098565) )),
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
