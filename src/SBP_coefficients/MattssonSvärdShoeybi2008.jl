
"""
    MattssonSvärdShoeybi2008

Coefficients of the SBP operators given in
  Mattsson, Svärd, Shoeybi (2008)
  Stable and accurate schemes for the compressible Navier-Stokes equations.
  Journal of Computational Physics 227, pp. 2293-2316.
"""
struct MattssonSvärdShoeybi2008 end

function Base.show(io::IO, ::MattssonSvärdShoeybi2008)
    print(io,
        "  Mattsson, Svärd, Shoeybi (2008) \n",
        "  Stable and accurate schemes for the compressible Navier-Stokes equations. \n",
        "  Journal of Computational Physics 227, pp. 2293-2316. \n")
end


function first_derivative_coefficients(source::MattssonSvärdShoeybi2008, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            DerivativeCoefficientRow{T,1,2}(SVector(-one(T), one(T))),
        )
        right_boundary = (
            -left_boundary[1],
        )
        upper_coef = SVector(T(1//2))
        central_coef = zero(T)
        lower_coef = -upper_coef

        DerivativeCoefficients(left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel, 1, order, source)
    elseif order == 4
        left_boundary = (
            DerivativeCoefficientRow{T,1,4}(SVector(T(-24//17), T(59//34), T(-4//17), T(-3//34))),
            DerivativeCoefficientRow{T,1,3}(SVector(T(-1//2), T(0), T(1//2))),
            DerivativeCoefficientRow{T,1,5}(SVector(T(4//43), T(-59//86), T(0), T(59//86), T(-4//43))),
            DerivativeCoefficientRow{T,1,6}(SVector(T(3//98), T(0), T(-59//98), T(0), T(32//49), T(-4//49))),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(2//3), T(-1//12))
        central_coef = zero(T)
        lower_coef = -upper_coef

        DerivativeCoefficients(left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel, 1, order, source)
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
            # q4
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

        DerivativeCoefficients(left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel, 1, order, source)
    elseif order == 8
        x1 =
        left_boundary = (
            # q1
            DerivativeCoefficientRow{T,1,8}(SVector(T(-2540160//1498139),
                                                    T(515174618457408334671//81327545493318772580),
                                                    T(-319653445960068136431//16265509098663754516),
                                                    T(1756838164798071830975//48796527295991263548),
                                                    T(-132855587210457907869//4066377274665938629),
                                                    T(205031990442925032261//16265509098663754516),
                                                    T(-147365687//1707878460),
                                                    T(-13621478277275954493//16265509098663754516) )),
            # q2
            DerivativeCoefficientRow{T,1,8}(SVector(T(-515174618457408334671//420776658856749814780),
                                                    T(0),
                                                    T(335158710375082986831//30055475632624986770),
                                                    T(-81614753349680015895//3005547563262498677),
                                                    T(996306452293020120805//36066570759149984124),
                                                    T(-338011221795904141851//30055475632624986770),
                                                    T(1109310//21038833),
                                                    T(503120655918727839389//631164988285124722170) )),
            # q3
            DerivativeCoefficientRow{T,1,8}(SVector(T(106551148653356045477//4733507987139058764),
                                                    T(-111719570125027662277//1690538566835378130),
                                                    T(0),
                                                    T(66611314247587471205//338107713367075626),
                                                    T(-177186717090009550283//676215426734151252),
                                                    T(13307564278128858831//112702571122358542),
                                                    T(-354462//1972295),
                                                    T(-7342453837193297117//788917997856509794) )),
            # q4
            DerivativeCoefficientRow{T,1,8}(SVector(T(-1756838164798071830975//297539771191584439476),
                                                    T(81614753349680015895//3542140133233148089),
                                                    T(-199833942742762413615//7084280266466296178),
                                                    T(0),
                                                    T(334165099234925485449//14168560532932592356),
                                                    T(-97726720448760690009//7084280266466296178),
                                                    T(-2934266//74384943),
                                                    T(65502786250599341781//49589961865264073246) )),
            # q5
            DerivativeCoefficientRow{T,1,9}(SVector(T(132855587210457907869//5691012984537517679),
                                                    T(-996306452293020120805//9756022259207173164),
                                                    T(531560151270028650849//3252007419735724388),
                                                    T(-334165099234925485449//3252007419735724388),
                                                    T(0),
                                                    T(68552420813742056079//3252007419735724388),
                                                    T(-2343561//22764052),
                                                    T(-197397151320811747355//68292155814450212148),
                                                    T(-2592//299527) )),
            # q6
            DerivativeCoefficientRow{T,1,10}(SVector(T(-9763428116329763441//3358016562304840404),
                                                     T(112670407265301380617//8395041405762101010),
                                                     T(-13307564278128858831//559669427050806734),
                                                     T(32575573482920230003//1679008281152420202),
                                                     T(-22850806937914018693//3358016562304840404),
                                                     T(0),
                                                     T(5346432//9794215),
                                                     T(1388463154122240//14728142817126493),
                                                     T(3072//103097),
                                                     T(-288//103097) )),
             # q7
             DerivativeCoefficientRow{T,1,11}(SVector(T(21052241//763903740),
                                                      T(-1109310//12731729),
                                                      T(3190158//63658645),
                                                      T(2934266//38195187),
                                                      T(2343561//50926916),
                                                      T(-48117888//63658645),
                                                      T(0),
                                                      T(52309152//63658645),
                                                      T(-145152//670091),
                                                      T(27648//670091),
                                                      T(-2592//670091) )),
              # q8
              DerivativeCoefficientRow{T,1,12}(SVector(T(13621478277275954493//55672594705880416916),
                                                       T(-503120655918727839389//417544460294103126870),
                                                       T(66082084534739674053//27836297352940208458),
                                                       T(-65502786250599341781//27836297352940208458),
                                                       T(197397151320811747355//167017784117641250748),
                                                       T(-87473178709701120//732534140866847591),
                                                       T(-366164064//487135205),
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

        DerivativeCoefficients(left_boundary, right_boundary, lower_coef, central_coef, upper_coef, parallel, 1, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end
