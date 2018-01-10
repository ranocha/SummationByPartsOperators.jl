
"""
MattssonNordström2004

Coefficients of the SBP operators given in
Mattsson, Nordström (2004)
Summation by parts operators for finite difference approximations of second
  derivatives.
Journal of Computational Physics 199, pp.503-540.
"""
struct MattssonNordström2004 <: SourceOfCoefficients end

function Base.show(io::IO, ::MattssonNordström2004)
print(io,
    "  Mattsson, Nordström (2004) \n",
    "  Summation by parts operators for finite difference approximations of second \n",
    "    derivaties. \n",
    "  Journal of Computational Physics 199, pp.503-540. \n")
end


function first_derivative_coefficients(source::MattssonNordström2004, order::Int, T=Float64, parallel=Val{:serial}())
if order == 2
    left_boundary = (
        DerivativeCoefficientRow{T,1,2}(SVector(-one(T), one(T))),
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
        # q1
        DerivativeCoefficientRow{T,1,4}(SVector(T(-24//17),
                                                T(59//34),
                                                T(-4//17),
                                                T(-3//34) )),
        # q2
        DerivativeCoefficientRow{T,1,3}(SVector(T(-1//2),
                                                T(0),
                                                T(1//2))),
        # q3
        DerivativeCoefficientRow{T,1,5}(SVector(T(4//43),
                                                T(-59//86),
                                                T(0),
                                                T(59//86),
                                                T(-4//43))),
        # q4
        DerivativeCoefficientRow{T,1,6}(SVector(T(3//98),
                                                T(0),
                                                T(-59//98),
                                                T(0),
                                                T(32//49),
                                                T(-4//49))),
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(2//3), T(-1//12))
    central_coef = zero(T)
    lower_coef = -upper_coef
    left_weights = SVector( T(17//48),
                            T(59//48),
                            T(43//48),
                            T(49//48) )
    right_weights = left_weights
    left_boundary_derivatives = Tuple{}()
    right_boundary_derivatives = left_boundary_derivatives

    DerivativeCoefficients(left_boundary, right_boundary,
                            left_boundary_derivatives, right_boundary_derivatives,
                            lower_coef, central_coef, upper_coef,
                            left_weights, right_weights, parallel, 1, order, source)
elseif order == 6
    left_boundary = (
        # q1
        DerivativeCoefficientRow{T,1,6}(SVector(T(-21600//13649),
                                                T(104009//54596),
                                                T(30443//81894),
                                                T(-33311//27298),
                                                T(16863//27298),
                                                T(-15025//163788) )),
        # q2
        DerivativeCoefficientRow{T,1,6}(SVector(T(-104009//240260),
                                                T(0),
                                                T(-311//72078),
                                                T(20229//24026),
                                                T(-24337//48052),
                                                T(36661//360390) )),
        # q3
        DerivativeCoefficientRow{T,1,6}(SVector(T(-30443//162660),
                                                T(311//32532),
                                                T(0),
                                                T(-11155//16266),
                                                T(41287//32532),
                                                T(-21999//54220) )),
        # q4 #TODO
        DerivativeCoefficientRow{T,1,7}(SVector(T(33311//107180),
                                                T(-20229//21436),
                                                T(485//1398),
                                                T(0),
                                                T(4147//21436),
                                                T(25427//321540),
                                                T(72//5359) )),
        # q5
        DerivativeCoefficientRow{T,1,8}(SVector(T(-16863//78770),
                                                T(24337//31508),
                                                T(-41287//47262),
                                                T(-4147//15754),
                                                T(0),
                                                T(342523//472620),
                                                T(-1296//7877),
                                                T(144//7877) )),
        # q6
        DerivativeCoefficientRow{T,1,9}(SVector(T(15025//525612),
                                                T(-36661//262806),
                                                T(21999//87602),
                                                T(-25427//262806),
                                                T(-342523//525612),
                                                T(0),
                                                T(32400//43801),
                                                T(-6480//43801),
                                                T(720//43801) ))
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(3//4), T(-3//20), T(1//60))
    central_coef = zero(T)
    lower_coef = -upper_coef
    left_weights = SVector( T(13649//43200),
                            T(12013//8640),
                            T(2711//4320),
                            T(5359//4320),
                            T(7877//8640),
                            T(43801//43200) )
    right_weights = left_weights
    left_boundary_derivatives = Tuple{}()
    right_boundary_derivatives = left_boundary_derivatives

    DerivativeCoefficients(left_boundary, right_boundary,
                            left_boundary_derivatives, right_boundary_derivatives,
                            lower_coef, central_coef, upper_coef,
                            left_weights, right_weights, parallel, 1, order, source)
elseif order == 8
    left_boundary = (
        # q1
        DerivativeCoefficientRow{T,1,8}(SVector(T(-2540160//1498139),
                                                T(5544277//5992556),
                                                T(198794991//29962780),
                                                T(-256916579//17977668),
                                                T(20708767//1498139),
                                                T(-41004357//5992556),
                                                T(27390659//17977668),
                                                T(-2323531//29962780) )),
        # q2
        DerivativeCoefficientRow{T,1,8}(SVector(T(-5544277//31004596),
                                                T(0),
                                                T(-85002381//22146140),
                                                T(49607267//4429228),
                                                T(-165990199//13287684),
                                                T(7655859//1107307),
                                                T(-7568311//4429228),
                                                T(48319961//465068940) )),
        # q3
        DerivativeCoefficientRow{T,1,8}(SVector(T(-66264997//8719620),
                                                T(9444709//415220),
                                                T(0),
                                                T(-20335981//249132),
                                                T(32320879//249132),
                                                T(-35518713//415220),
                                                T(2502774//103805),
                                                T(-3177073//1743924) )),
        # q4
        DerivativeCoefficientRow{T,1,8}(SVector(T(256916579//109619916),
                                                T(-49607267//5219996),
                                                T(61007943//5219996),
                                                T(0),
                                                T(-68748371//5219996),
                                                T(65088123//5219996),
                                                T(-66558305//15659988),
                                                T(3870214//9134993) )),
        # q5
        DerivativeCoefficientRow{T,1,9}(SVector(T(-20708767//2096689),
                                                T(165990199//3594324),
                                                T(-96962637//1198108),
                                                T(68748371//1198108),
                                                T(0),
                                                T(-27294549//1198108),
                                                T(14054993//1198108),
                                                T(-42678199//25160268),
                                                T(-2592//299527) )),
        # q6
        DerivativeCoefficientRow{T,1,10}(SVector(T(13668119//8660148),
                                                 T(-850651//103097),
                                                 T(35518713//2061940),
                                                 T(-21696041//1237164),
                                                 T(9098183//1237164),
                                                 T(0),
                                                 T(-231661//412388),
                                                 T(7120007//43300740),
                                                 T(3072//103097),
                                                 T(-288//103097) )),
         # q7
         DerivativeCoefficientRow{T,1,11}(SVector(T(-27390659//56287644),
                                                  T(7568311//2680364),
                                                  T(-22524966//3350455),
                                                  T(66558305//8041092),
                                                  T(-14054993//2680364),
                                                  T(2084949//2680364),
                                                  T(0),
                                                  T(70710683//93812740),
                                                  T(-145152//670091),
                                                  T(27648//670091),
                                                  T(-2592//670091) )),
          # q8
          DerivativeCoefficientRow{T,1,12}(SVector(T(2323531//102554780),
                                                   T(-48319961//307664340),
                                                   T(9531219//20510956),
                                                   T(-3870214//5127739),
                                                   T(2246221//3238572),
                                                   T(-21360021//102554780),
                                                   T(-70710683//102554780),
                                                   T(0),
                                                   T(4064256//5127739),
                                                   T(-1016064//5127739),
                                                   T(193536//5127739),
                                                   T(-18144//5127739) )),
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(4//5), T(-1//5), T(4//105), T(-1//280))
    central_coef = zero(T)
    lower_coef = -upper_coef
    left_weights = SVector( T(1498139//5080320),
                            T(1107307//725760),
                            T(20761//80640),
                            T(1304999//725760),
                            T(299527//725760),
                            T(103097//80640),
                            T(670091//725760),
                            T(5127739//5080320) )
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


function second_derivative_coefficients(source::MattssonNordström2004, order::Int, T=Float64, parallel=Val{:serial}())
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
        DerivativeCoefficientRow{T,1,3}(SVector(T(-3//2),
                                                T(2),
                                                T(-1//2) )),
    )
    right_boundary_derivatives = (-left_boundary_derivatives[1], )

    DerivativeCoefficients(left_boundary, right_boundary,
                            left_boundary_derivatives, right_boundary_derivatives,
                            lower_coef, central_coef, upper_coef,
                            left_weights, right_weights, parallel, 2, order, source)
elseif order == 4
    left_boundary = (
        #TODO: 
        # d1
        DerivativeCoefficientRow{T,1,4}(SVector(T(2),
                                                T(-5),
                                                T(4),
                                                T(-1) )),
        # d2
        DerivativeCoefficientRow{T,1,3}(SVector(T(1),
                                                T(-2),
                                                T(1) )),
        # d3
        DerivativeCoefficientRow{T,1,5}(SVector(T(-4//43),
                                                T(59//43),
                                                T(-110//43),
                                                T(59//43),
                                                T(-4//43) )),
        # d4
        DerivativeCoefficientRow{T,1,6}(SVector(T(-1//49),
                                                T(0),
                                                T(59//49),
                                                T(-118//49),
                                                T(64//49),
                                                T(-4//49) ))
    )
    right_boundary = left_boundary
    upper_coef = SVector(T(4//3), T(-1//12))
    central_coef = T(-5//2)
    lower_coef = upper_coef
    left_weights = SVector(T(17//48), T(59//48), T(43//48), T(49//48))
    right_weights = left_weights
    left_boundary_derivatives = (
        DerivativeCoefficientRow{T,1,4}(SVector(T(-11/6),
                                                T(3),
                                                T(-3//2),
                                                T(1//3) )),
    )
    right_boundary_derivatives = (-left_boundary_derivatives[1], )

    DerivativeCoefficients(left_boundary, right_boundary,
                            left_boundary_derivatives, right_boundary_derivatives,
                            lower_coef, central_coef, upper_coef,
                            left_weights, right_weights, parallel, 2, order, source)
elseif order == 6
    left_boundary = (
        # d1
        DerivativeCoefficientRow{T,1,6}(SVector(T(114170//40947),
                                                T(-438107//54596),
                                                T(336409//40947),
                                                T(-276997//81894),
                                                T(3747//13649),
                                                T(21035//163788) )),
        # d2
        DerivativeCoefficientRow{T,1,6}(SVector(T(6173//5860),
                                                T(-2066//879),
                                                T(3283//1758),
                                                T(-303//293),
                                                T(2111//3516),
                                                T(-601//4395) )),
        # d3
        DerivativeCoefficientRow{T,1,6}(SVector(T(-52391//81330),
                                                T(134603//32532),
                                                T(-21982//2711),
                                                T(112915//16266),
                                                T(-46969//16266),
                                                T(30409//54220) )),
        # d4
        DerivativeCoefficientRow{T,1,7}(SVector(T(68603//321540),
                                                T(-12423//10718),
                                                T(112915//32154),
                                                T(-75934//16077),
                                                T(53369//21436),
                                                T(-54899//160770),
                                                T(48//5359) )),
        # d5
        DerivativeCoefficientRow{T,1,8}(SVector(T(-7053//39385),
                                                T(86551//94524),
                                                T(-46969//23631),
                                                T(53369//15754),
                                                T(-87904//23631),
                                                T(820271//472620),
                                                T(-1296//7877),
                                                T(96//7877) )),
        # d6
        DerivativeCoefficientRow{T,1,9}(SVector(T(21035//525612),
                                                T(-24641//131403),
                                                T(30409//87602),
                                                T(-54899//131403),
                                                T(820271//525612),
                                                T(-117600//43801),
                                                T(64800//43801),
                                                T(-6480//43801),
                                                T(480//43801) )),
    )
    right_boundary = left_boundary
    upper_coef = SVector(T(3//2), T(-3//20), T(1//90))
    central_coef = T(-49//18)
    lower_coef = upper_coef
    left_weights = SVector( T(13649//43200),
                            T(12013//8640),
                            T(2711//4320),
                            T(5359//4320),
                            T(7877//8640),
                            T(43801//43200) )
    right_weights = left_weights
    left_boundary_derivatives = (
        DerivativeCoefficientRow{T,1,5}(SVector(T(-25//12),
                                                T(4),
                                                T(-3),
                                                T(4//3),
                                                T(-1//4) )),
    )
    right_boundary_derivatives = (-left_boundary_derivatives[1], )

    DerivativeCoefficients(left_boundary, right_boundary,
                            left_boundary_derivatives, right_boundary_derivatives,
                            lower_coef, central_coef, upper_coef,
                            left_weights, right_weights, parallel, 2, order, source)
elseif order == 8
    left_boundary = (
        # d1
        DerivativeCoefficientRow{T,1,7}(SVector(T(4870382994799//1358976868290),
                                                T(-893640087518//75498714905),
                                                T(926594825119//60398971924),
                                                T(-1315109406200//135897686829),
                                                T(39126983272//15099742981),
                                                T(12344491342//75498714905),
                                                T(-451560522577//2717953736580) )),
        # d2
        DerivativeCoefficientRow{T,1,8}(SVector(T(333806012194//390619153855),
                                                T(-154646272029//111605472530),
                                                T(1168338040//33481641759),
                                                T(82699112501//133926567036),
                                                T(-171562838//11160547253),
                                                T(-28244698346//167408208795),
                                                T(11904122576//167408208795),
                                                T(-2598164715//312495323084) )),
        # d3
        DerivativeCoefficientRow{T,1,8}(SVector(T(7838984095//52731029988),
                                                T(1168338040//5649753213),
                                                T(-88747895//144865467),
                                                T(423587231//627750357),
                                                T(-43205598281//22599012852),
                                                T(4876378562//1883251071),
                                                T(-5124426509//3766502142),
                                                T(10496900965//39548272491) )),
        # d4
        DerivativeCoefficientRow{T,1,8}(SVector(T(-94978241528//828644350023),
                                                T(82699112501//157837019052),
                                                T(1270761693//13153084921),
                                                T(-167389605005//118377764289),
                                                T(48242560214//39459254763),
                                                T(-31673996013//52612339684),
                                                T(43556319241//118377764289),
                                                T(-44430275135//552429566682) )),
        # d5
        DerivativeCoefficientRow{T,1,9}(SVector(T(1455067816//21132528431),
                                                T(-171562838//3018932633),
                                                T(-43205598281//36227191596),
                                                T(48242560214//9056797899),
                                                T(-52276055645//6037865266),
                                                T(57521587238//9056797899),
                                                T(-80321706377//36227191596),
                                                T(8078087158//21132528431),
                                                T(-1296//299527) )),
        # d6
        DerivativeCoefficientRow{T,1,10}(SVector(T(10881504334//327321118845),
                                                 T(-28244698346//140280479505),
                                                 T(4876378562//9352031967),
                                                 T(-10557998671//12469375956),
                                                 T(57521587238//28056095901),
                                                 T(-278531401019//93520319670),
                                                 T(73790130002//46760159835),
                                                 T(-137529995233//785570685228),
                                                 T(2048//103097),
                                                 T(-144//103097) )),
         # d7
         DerivativeCoefficientRow{T,1,11}(SVector(T(-135555328849//8509847458140),
                                                  T(11904122576//101307707835),
                                                  T(-5124426509//13507694378),
                                                  T(43556319241//60784624701),
                                                  T(-80321706377//81046166268),
                                                  T(73790130002//33769235945),
                                                  T(-950494905688//303923123505),
                                                  T(239073018673//141830790969),
                                                  T(-145152//670091),
                                                  T(18432//670091),
                                                  T(-1296//670091) )),
          # d8
          DerivativeCoefficientRow{T,1,12}(SVector(T(0),
                                                   T(-2598164715//206729925524),
                                                   T(10496900965//155047444143),
                                                   T(-44430275135//310094888286),
                                                   T(425162482//2720130599),
                                                   T(-137529995233//620189776572),
                                                   T(239073018673//155047444143),
                                                   T(-144648000000//51682481381),
                                                   T(8128512//5127739),
                                                   T(-1016064//5127739),
                                                   T(129024//5127739),
                                                   T(-9072//5127739) )),

    )
    right_boundary = left_boundary
    upper_coef = SVector(T(8//5), T(-1//5), T(8//315), T(-1//560))
    central_coef = T(-205//72)
    lower_coef = upper_coef
    left_weights = SVector( T(1498139//5080320),
                            T(1107307//725760),
                            T(20761//80640),
                            T(1304999//725760),
                            T(299527//725760),
                            T(103097//80640),
                            T(670091//725760),
                            T(5127739//5080320) )
    right_weights = left_weights
    left_boundary_derivatives = (
        DerivativeCoefficientRow{T,1,7}(SVector(T(-4723//2100),
                                                T(839//175),
                                                T(-157//35),
                                                T(278//105),
                                                T(-103//140),
                                                T(-1//175),
                                                T(6//175) )),
    )
    right_boundary_derivatives = (-left_boundary_derivatives[1], )

    DerivativeCoefficients(left_boundary, right_boundary,
                            left_boundary_derivatives, right_boundary_derivatives,
                            lower_coef, central_coef, upper_coef,
                            left_weights, right_weights, parallel, 2, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end
