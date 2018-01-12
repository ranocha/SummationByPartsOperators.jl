
"""
    Mattsson2012

Coefficients of the SBP operators given in
  Mattsson (2012)
  Summation by Parts Operators for Finite Difference Approximations of
    Second-Derivatives with Variable Coefficients.
  Journal of Scientific Computing 51, pp. 650-682.
"""
struct Mattsson2012 <: SourceOfCoefficients end

function Base.show(io::IO, ::Mattsson2012)
    print(io,
        "  Mattsson (2012) \n",
        "  Summation by Parts Operators for Finite Difference Approximations of\n",
        "    Second-Derivatives with Variable Coefficients. \n",
        "  Journal of Scientific Computing 51, pp. 650-682. \n",
        "See also (first derivatives) \n",
        "  Mattsson, Nordström (2004) \n",
        "  Summation by parts operators for finite difference approximations of second \n",
        "    derivaties. \n",
        "  Journal of Computational Physics 199, pp.503-540. \n")
end


@inline function first_derivative_coefficients(source::Mattsson2012, order::Int, T=Float64, parallel=Val{:serial}())
    first_derivative_coefficients(MattssonNordström2004(), order, T, parallel)
end

@inline function second_derivative_coefficients(source::Mattsson2012, order::Int, T=Float64, parallel=Val{:serial}())
    second_derivative_coefficients(MattssonNordström2004(), order, T, parallel)
end


function var_coef_derivative_coefficients(source::Mattsson2012, derivative_order::Int, accuracy_order::Int, grid, parallel=Val{:serial}())
    @argcheck derivative_order == 2
    T = eltype(grid)
    if accuracy_order == 2
        coefficient_cache = Mattsson2012Cache2(T)
        left_weights = SVector(T(1//2))
        right_weights = left_weights
    elseif accuracy_order == 4
        coefficient_cache = Mattsson2012Cache4(T)
        left_weights = SVector( T(17//48),
                                T(59//48),
                                T(43//48),
                                T(49//48) )
        right_weights = left_weights
    elseif accuracy_order == 6
        coefficient_cache = Mattsson2012Cache6(T)
        left_weights = SVector( T(13649//43200),
                                T(12013//8640),
                                T(2711//4320),
                                T(5359//4320),
                                T(7877//8640),
                                T(43801//43200) )
        right_weights = left_weights
    elseif accuracy_order == 8
        coefficient_cache = Mattsson2012Cache8(T)
        left_weights = SVector( T(1498139//5080320),
                                T(1107307//725760),
                                T(20761//80640),
                                T(1304999//725760),
                                T(299527//725760),
                                T(103097//80640),
                                T(670091//725760),
                                T(5127739//5080320) )
        right_weights = left_weights
    else
        throw(ArgumentError("Order of accuracy $accuracy_order not implemented/derived."))
    end

    VarCoefDerivativeCoefficients(coefficient_cache, left_weights, right_weights,
                                  parallel, derivative_order, accuracy_order, source)
end



struct Mattsson2012Cache2{T} <: AbstractCoefficientCache{T}
    half::T

    function Mattsson2012Cache2(::Type{T}) where {T}
        half = one(T) / 2

        new{T}(half)
    end
end

lower_bandwidth(cache::Mattsson2012Cache2) = 2
upper_bandwidth(cache::Mattsson2012Cache2) = 2
Base.checkbounds(::Type{Bool}, u::AbstractVector, ::Mattsson2012Cache2) = length(u) > 2
left_length(::Mattsson2012Cache2) = 1
right_length(::Mattsson2012Cache2) = 1

function convolve_boundary_coefficients!(dest::AbstractVector, cache::Mattsson2012Cache2,
                                         u::AbstractVector, b::AbstractVector, α)
    @inbounds begin
        dest[  1] = α * (
                        (2*b[1] - b[2]) * u[1]
                        + (-3*b[1] + b[2]) * u[2]
                        + b[1] * u[3]
                    )

        dest[end] = α * (
                        (2*b[end] - b[end-1]) * u[end]
                        + (-3*b[end] + b[end-1]) * u[end-1]
                        + b[end] * u[end-2]
                    )
    end
end

function convolve_boundary_coefficients!(dest::AbstractVector, cache::Mattsson2012Cache2,
                                         u::AbstractVector, b::AbstractVector, α, β)
    @inbounds begin
        dest[  1] = α * (
                        (2*b[1] - b[2]) * u[1]
                        + (-3*b[1] + b[2]) * u[2]
                        + b[1] * u[3]
                    ) + β*dest[1]

        dest[end] = α * (
                        (2*b[end] - b[end-1]) * u[end]
                        + (-3*b[end] + b[end-1]) * u[end-1]
                        + b[end] * u[end-2]
                    ) + β*dest[end]
    end
end

@inline function convolve_interior_coefficients_loopbody(i, cache::Mattsson2012Cache2, u, b)
    @unpack half = cache
    @inbounds begin
        b_im1 = b[i-1]
        b_i   = b[i]
        b_ip1 = b[i+1]

        retval = half * (
                    (b_im1 + b_i) * u[i-1]
                    - (b_im1 + 2*b_i + b_ip1) * u[i]
                    + (b_i + b_ip1) * u[i+1]
                )
    end

    retval
end



struct Mattsson2012Cache4{T} <: AbstractCoefficientCache{T}
    d111::T
    d112::T
    d113::T
    d114::T
    d121::T
    d123::T
    d124::T
    d131::T
    d132::T
    d133::T
    d134::T
    d141::T
    d143::T
    d144::T
    d153::T
    d154::T
    d163::T
    d164::T
    d211::T
    d213::T
    d214::T
    d221::T
    d223::T
    d224::T
    d231::T
    d233::T
    d234::T
    d241::T
    d243::T
    d244::T
    d253::T
    d254::T
    d263::T
    d264::T
    d311::T
    d312::T
    d313::T
    d314::T
    d321::T
    d323::T
    d324::T
    d331::T
    d332::T
    d333::T
    d334::T
    d335::T
    d341::T
    d343::T
    d344::T
    d345::T
    d353::T
    d354::T
    d355::T
    d363::T
    d364::T
    d365::T
    d411::T
    d413::T
    d414::T
    d421::T
    d423::T
    d424::T
    d431::T
    d433::T
    d434::T
    d435::T
    d441::T
    d443::T
    d444::T
    d445::T
    d446::T
    d453::T
    d454::T
    d455::T
    d456::T
    d463::T
    d464::T
    d465::T
    d466::T
    d513::T
    d514::T
    d523::T
    d524::T
    d533::T
    d534::T
    d535::T
    d543::T
    d544::T
    d545::T
    d546::T
    d553::T
    d554::T
    d555::T
    d556::T
    d557::T
    d563::T
    d564::T
    d565::T
    d566::T
    d567::T
    d575::T
    d576::T
    d577::T
    d613::T
    d614::T
    d623::T
    d624::T
    d633::T
    d634::T
    d635::T
    d643::T
    d644::T
    d645::T
    d646::T
    d653::T
    d654::T
    d655::T
    d656::T
    d657::T
    d663::T
    d664::T
    d665::T
    d666::T
    d667::T
    d668::T
    d675::T
    d676::T
    d677::T
    d678::T
    d686::T
    d687::T
    d688::T
    dim2m2::T
    dim2m1::T
    dim2_0::T
    dim1m2::T
    dim1m1::T
    dim1_0::T
    dim1p1::T
    di_0m2::T
    di_0m1::T
    di_0_0::T
    di_0p1::T
    di_0p2::T


    function Mattsson2012Cache4(::Type{T}) where {T}
        d111 = T(920//289)
        d112 = T(-59//68)
        d113 = T(-81031200387//366633756146)
        d114 = T(-69462376031//733267512292)
        d121 = T(-1740//289)
        d123 = T(6025413881//7482321554)
        d124 = T(1612249989//7482321554)
        d131 = T(1128//289)
        d132 = T(59//68)
        d133 = T(-6251815797//8526366422)
        d134 = T(-639954015//17052732844)
        d141 = T(-308//289)
        d143 = T(1244724001//7482321554)
        d144 = T(-752806667//7482321554)
        d153 = T(-148737261//10783345769)
        d154 = T(148737261//10783345769)
        d163 = T(-3//833)
        d164 = T(3//833)

        d211 = T(12/17)
        d213 = T(102125659//440136562)
        d214 = T(27326271//440136562)
        d221 = T(-59//68)
        d223 = T(-156920047993625//159775733917868)
        d224 = T(-12001237118451//79887866958934)
        d231 = T(2//17)
        d233 = T(1489556735319//1857857371138)
        d234 = T(149729180391//1857857371138)
        d241 = T(3//68)
        d243 = T(-13235456910147//159775733917868)
        d244 = T(3093263736297//79887866958934)
        d253 = T(67535018271//2349643145851)
        d254 = T(-67535018271//2349643145851)
        d263 = T(441//181507)
        d264 = T(-441//181507)

        d311 = T(-96//731)
        d312 = T(59//172)
        d313 = T(-6251815797//21566691538)
        d314 = T(-639954015//43133383076)
        d321 = T(118//731)
        d323 = T(87883847383821//79887866958934)
        d324 = T(8834021643069//79887866958934)
        d331 = T(-16//731)
        d332 = T(-59//172)
        d333 = T(-1134866646907639536627//727679167377258785038)
        d334 = T(-13777050223300597//23487032885926596)
        d335 = T(-26254//557679)
        d341 = T(-6//731)
        d343 = T(14509020271326561681//14850595252597118062)
        d344 = T(17220493277981//79887866958934)
        d345 = T(1500708//7993399)
        d353 = T(-4841930283098652915//21402328452272317207)
        d354 = T(31597236232005//115132514146699)
        d355 = T(-26254//185893)
        d363 = T(-2318724711//1653303156799)
        d364 = T(960119//1147305747)
        d365 = T(13564//23980197)

        d411 = T(-36//833)
        d413 = T(1244724001//21566691538)
        d414 = T(-752806667//21566691538)
        d421 = T(177//3332)
        d423 = T(-780891957698673//7829010961975532)
        d424 = T(3724542049827//79887866958934)
        d431 = T(-6//833)
        d433 = T(14509020271326561681//16922771334354855466)
        d434 = T(2460070468283//13005001597966)
        d435 = T(1500708//9108757)
        d441 = T(-9//3332)
        d443 = T(-217407431400324796377//207908333536359652868)
        d444 = T(-1950062198436997//3914505480987766)
        d445 = T(-7476412//9108757)
        d446 = T(-2//49)
        d453 = T(4959271814984644613//21402328452272317207)
        d454 = T(47996144728947//115132514146699)
        d455 = T(4502124//9108757)
        d456 = T(8//49)
        d463 = T(-2258420001//1653303156799)
        d464 = T(-1063649//8893843)
        d465 = T(1473580//9108757)
        d466 = T(-6//49)

        d513 = T(-49579087//10149031312)
        d514 = T(49579087//10149031312)
        d523 = T(1328188692663//37594290333616)
        d524 = T(-1328188692663//37594290333616)
        d533 = T(-1613976761032884305//7963657098519931984)
        d534 = T(10532412077335//42840005263888)
        d535 = T(-564461//4461432)
        d543 = T(4959271814984644613//20965546238960637264)
        d544 = T(15998714909649//37594290333616)
        d545 = T(375177//743572)
        d546 = T(1//6)
        d553 = T(-8386761355510099813//128413970713633903242)
        d554 = T(-2224717261773437//2763180339520776)
        d555 = T(-280535//371786)
        d556 = T(-5//6)
        d557 = T(-1//24)
        d563 = T(13091810925//13226425254392)
        d564 = T(35039615//213452232)
        d565 = T(1118749//2230716)
        d566 = T(1//2)
        d567 = T(1//6)
        d575 = T(-1//8)
        d576 = T(1//6)
        d577 = T(-1//8)

        d613 = T(-1//784)
        d614 = T(1//784)
        d623 = T(8673//2904112)
        d624 = T(-8673//2904112)
        d633 = T(-33235054191//26452850508784)
        d634 = T(960119//1280713392)
        d635 = T(3391//6692148)
        d643 = T(-752806667//539854092016)
        d644 = T(-1063649//8712336)
        d645 = T(368395//2230716)
        d646 = T(-1//8)
        d653 = T(13091810925//13226425254392)
        d654 = T(35039615//213452232)
        d655 = T(1118749//2230716)
        d656 = T(1//2)
        d657 = T(1//6)
        d663 = T(-660204843//13226425254392)
        d664 = T(-3290636//80044587)
        d665 = T(-5580181//6692148)
        d666 = T(-3//4)
        d667 = T(-5//6)
        d668 = T(-1//24)
        d675 = T(1//6)
        d676 = T(1//2)
        d677 = T(1//2)
        d678 = T(1//6)
        d686 = T(-1//8)
        d687 = T(1//6)
        d688 = T(-1//8)

        dim2m2 = T(-1//8)
        dim2m1 = T(1//6)
        dim2_0 = T(-1//8)
        dim1m2 = T(1//6)
        dim1m1 = T(1//2)
        dim1_0 = T(1//2)
        dim1p1 = T(1//6)
        di_0m2 = T(-1//24)
        di_0m1 = T(-5//6)
        di_0_0 = T(-3//4)
        di_0p1 = T(-5//6)
        di_0p2 = T(-1//24)

        new{T}( d111, d112, d113, d114, d121, d123, d124, d131, d132, d133, d134,
                d141, d143, d144, d153, d154, d163, d164,
                d211, d213, d214, d221, d223, d224, d231, d233, d234, d241, d243, d244,
                d253, d254, d263, d264,
                d311, d312, d313, d314, d321, d323, d324, d331, d332, d333, d334, d335,
                d341, d343, d344, d345, d353, d354, d355, d363, d364, d365,
                d411, d413, d414, d421, d423, d424, d431, d433, d434, d435,
                d441, d443, d444, d445, d446, d453, d454, d455, d456,
                d463, d464, d465, d466,
                d513, d514, d523, d524, d533, d534, d535, d543, d544, d545, d546,
                d553, d554, d555, d556, d557, d563, d564, d565, d566, d567,
                d575, d576, d577,
                d613, d614, d623, d624, d633, d634, d635, d643, d644, d645, d646,
                d653, d654, d655, d656, d657, d663, d664, d665, d666, d667, d668,
                d675, d676, d677, d678, d686, d687, d688,
                dim2m2, dim2m1, dim2_0, dim1m2, dim1m1, dim1_0, dim1p1,
                di_0m2, di_0m1, di_0_0, di_0p1, di_0p2)
    end
end

lower_bandwidth(cache::Mattsson2012Cache4) = 3
upper_bandwidth(cache::Mattsson2012Cache4) = 3
Base.checkbounds(::Type{Bool}, u::AbstractVector, ::Mattsson2012Cache4) = length(u) > 8
left_length(::Mattsson2012Cache4) = 6
right_length(::Mattsson2012Cache4) = 6

function convolve_boundary_coefficients!(dest::AbstractVector, cache::Mattsson2012Cache4,
                                         u::AbstractVector, b::AbstractVector, α)
    @unpack d111, d112, d113, d114, d121, d123, d124, d131, d132, d133, d134,
            d141, d143, d144, d153, d154, d163, d164,
            d211, d213, d214, d221, d223, d224, d231, d233, d234, d241, d243, d244,
            d253, d254, d263, d264,
            d311, d312, d313, d314, d321, d323, d324, d331, d332, d333, d334, d335,
            d341, d343, d344, d345, d353, d354, d355, d363, d364, d365,
            d411, d413, d414, d421, d423, d424, d431, d433, d434, d435,
            d441, d443, d444, d445, d446, d453, d454, d455, d456,
            d463, d464, d465, d466,
            d513, d514, d523, d524, d533, d534, d535, d543, d544, d545, d546,
            d553, d554, d555, d556, d557, d563, d564, d565, d566, d567,
            d575, d576, d577,
            d613, d614, d623, d624, d633, d634, d635, d643, d644, d645, d646,
            d653, d654, d655, d656, d657, d663, d664, d665, d666, d667, d668,
            d675, d676, d677, d678, d686, d687, d688, = cache
    @inbounds begin
        #b1 = b[1]
        #b2 = b[2]
        #b3 = b[3]
        #b4 = b[4]
        #b5 = b[5]
        #b6 = b[6]
        #b7 = b[7]
        #b8 = b[8]

        dest[  1] = α * (
                          (d111*b[1] + d112*b[2] + d113*b[3] + d114*b[4]) * u[1]
                        + (d121*b[1] + d123*b[3] + d124*b[4]) * u[2]
                        + (d131*b[1] + d132*b[2] + d133*b[3] + d134*b[4]) * u[3]
                        + (d141*b[1] + d143*b[3] + d144*b[4]) * u[4]
                        + (d153*b[3] + d154*b[4]) * u[5]
                        + (d163*b[3] + d164*b[4]) * u[6]
                    )
        dest[  2] = α * (
                          (d211*b[1] + d213*b[3] + d214*b[4]) * u[1]
                        + (d221*b[1] + d223*b[3] + d224*b[4]) * u[2]
                        + (d231*b[1] + d233*b[3] + d234*b[4]) * u[3]
                        + (d241*b[1] + d243*b[3] + d244*b[4]) * u[4]
                        + (d253*b[3] + d254*b[4]) * u[5]
                        + (d263*b[3] + d264*b[4]) * u[6]
                    )
        dest[  3] = α * (
                          (d311*b[1] + d312*b[2] + d313*b[3] + d314*b[4]) * u[1]
                        + (d321*b[1] + d323*b[3] + d324*b[4]) * u[2]
                        + (d331*b[1] + d332*b[2] + d333*b[3] + d334*b[4] + d335*b[5]) * u[3]
                        + (d341*b[1] + d343*b[3] + d344*b[4] + d345*b[5]) * u[4]
                        + (d353*b[3] + d354*b[4] + d355*b[5]) * u[5]
                        + (d363*b[3] + d364*b[4] + d365*b[5]) * u[6]
                    )
        dest[  4] = α * (
                          (d411*b[1] + d413*b[3] + d414) * u[1]
                        + (d421*b[1] + d423*b[3] + d424) * u[2]
                        + (d431*b[1] + d433*b[3] + d434*b[4] + d435*b[5]) * u[3]
                        + (d441*b[1] + d443*b[3] + d444*b[4] + d445*b[5] + d446*b[6]) * u[4]
                        + (d453*b[3] + d454*b[4] + d455*b[5] + d456*b[6]) * u[5]
                        + (d463*b[3] + d464*b[4] + d465*b[5] + d466*b[6]) * u[6]
                    )
        dest[  5] = α * (
                          (d513*b[3] + d514) * u[1]
                        + (d523*b[3] + d524) * u[2]
                        + (d533*b[3] + d534*b[4] + d535*b[5]) * u[3]
                        + (d543*b[3] + d544*b[4] + d545*b[5] + d546*b[6]) * u[4]
                        + (d553*b[3] + d554*b[4] + d555*b[5] + d556*b[6] + d557*b[7]) * u[5]
                        + (d563*b[3] + d564*b[4] + d565*b[5] + d566*b[6] + d567*b[7]) * u[6]
                        + (d575*b[5] + d576*b[6] + d577*b[7]) * u[7]
                    )
        dest[  6] = α * (
                          (d613*b[3] + d614) * u[1]
                        + (d623*b[3] + d624) * u[2]
                        + (d633*b[3] + d634*b[4] + d635*b[5]) * u[3]
                        + (d643*b[3] + d644*b[4] + d645*b[5] + d646*b[6]) * u[4]
                        + (d653*b[3] + d654*b[4] + d655*b[5] + d656*b[6] + d657*b[7]) * u[5]
                        + (d663*b[3] + d664*b[4] + d665*b[5] + d666*b[6] + d667*b[7] + d668*b[8]) * u[6]
                        + (d675*b[5] + d676*b[6] + d677*b[7] + d678) * u[7]
                        + (d686*b[6] + d687*b[7] + d688*b[8]) * u[8]
                    )


        dest[end] = α * (
                          (d111*b[end] + d112*b[end-1] + d113*b[end-2] + d114*b[end-3]) * u[end]
                        + (d121*b[end] + d123*b[end-2] + d124*b[end-3]) * u[end-1]
                        + (d131*b[end] + d132*b[end-1] + d133*b[end-2] + d134*b[end-3]) * u[end-2]
                        + (d141*b[end] + d143*b[end-2] + d144*b[end-3]) * u[end-3]
                        + (d153*b[end-2] + d154*b[end-3]) * u[end-4]
                        + (d163*b[end-2] + d164*b[end-3]) * u[end-5]
                    )
        dest[end-1] = α * (
                          (d211*b[end] + d213*b[end-2] + d214*b[end-3]) * u[end]
                        + (d221*b[end] + d223*b[end-2] + d224*b[end-3]) * u[end-1]
                        + (d231*b[end] + d233*b[end-2] + d234*b[end-3]) * u[end-2]
                        + (d241*b[end] + d243*b[end-2] + d244*b[end-3]) * u[end-3]
                        + (d253*b[end-2] + d254*b[end-3]) * u[end-4]
                        + (d263*b[end-2] + d264*b[end-3]) * u[end-5]
                    )
        dest[end-2] = α * (
                          (d311*b[end] + d312*b[end-1] + d313*b[end-2] + d314*b[end-3]) * u[end]
                        + (d321*b[end] + d323*b[end-2] + d324*b[end-3]) * u[end-1]
                        + (d331*b[end] + d332*b[end-1] + d333*b[end-2] + d334*b[end-3] + d335*b[end-4]) * u[end-2]
                        + (d341*b[end] + d343*b[end-2] + d344*b[end-3] + d345*b[end-4]) * u[end-3]
                        + (d353*b[end-2] + d354*b[end-3] + d355*b[end-4]) * u[end-4]
                        + (d363*b[end-2] + d364*b[end-3] + d365*b[end-4]) * u[end-5]
                    )
        dest[end-3] = α * (
                          (d411*b[end] + d413*b[end-2] + d414) * u[end]
                        + (d421*b[end] + d423*b[end-2] + d424) * u[end-1]
                        + (d431*b[end] + d433*b[end-2] + d434*b[end-3] + d435*b[end-4]) * u[end-2]
                        + (d441*b[end] + d443*b[end-2] + d444*b[end-3] + d445*b[end-4] + d446*b[end-5]) * u[end-3]
                        + (d453*b[end-2] + d454*b[end-3] + d455*b[end-4] + d456*b[end-5]) * u[end-4]
                        + (d463*b[end-2] + d464*b[end-3] + d465*b[end-4] + d466*b[end-5]) * u[end-5]
                    )
        dest[end-4] = α * (
                          (d513*b[end-2] + d514) * u[end]
                        + (d523*b[end-2] + d524) * u[end-1]
                        + (d533*b[end-2] + d534*b[end-3] + d535*b[end-4]) * u[end-2]
                        + (d543*b[end-2] + d544*b[end-3] + d545*b[end-4] + d546*b[end-5]) * u[end-3]
                        + (d553*b[end-2] + d554*b[end-3] + d555*b[end-4] + d556*b[end-5] + d557*b[end-6]) * u[end-4]
                        + (d563*b[end-2] + d564*b[end-3] + d565*b[end-4] + d566*b[end-5] + d567*b[end-6]) * u[end-5]
                        + (d575*b[end-4] + d576*b[end-5] + d577*b[end-6]) * u[end-6]
                    )
        dest[end-5] = α * (
                          (d613*b[end-2] + d614) * u[end]
                        + (d623*b[end-2] + d624) * u[end-1]
                        + (d633*b[end-2] + d634*b[end-3] + d635*b[end-4]) * u[end-2]
                        + (d643*b[end-2] + d644*b[end-3] + d645*b[end-4] + d646*b[end-5]) * u[end-3]
                        + (d653*b[end-2] + d654*b[end-3] + d655*b[end-4] + d656*b[end-5] + d657*b[end-6]) * u[end-4]
                        + (d663*b[end-2] + d664*b[end-3] + d665*b[end-4] + d666*b[end-5] + d667*b[end-6] + d668*b[end-7]) * u[end-5]
                        + (d675*b[end-4] + d676*b[end-5] + d677*b[end-6] + d678) * u[end-6]
                        + (d686*b[end-5] + d687*b[end-6] + d688*b[end-7]) * u[end-7]
                    )
    end
end

function convolve_boundary_coefficients!(dest::AbstractVector, cache::Mattsson2012Cache4,
                                         u::AbstractVector, b::AbstractVector, α, β)
    @unpack d111, d112, d113, d114, d121, d123, d124, d131, d132, d133, d134,
            d141, d143, d144, d153, d154, d163, d164,
            d211, d213, d214, d221, d223, d224, d231, d233, d234, d241, d243, d244,
            d253, d254, d263, d264,
            d311, d312, d313, d314, d321, d323, d324, d331, d332, d333, d334, d335,
            d341, d343, d344, d345, d353, d354, d355, d363, d364, d365,
            d411, d413, d414, d421, d423, d424, d431, d433, d434, d435,
            d441, d443, d444, d445, d446, d453, d454, d455, d456,
            d463, d464, d465, d466,
            d513, d514, d523, d524, d533, d534, d535, d543, d544, d545, d546,
            d553, d554, d555, d556, d557, d563, d564, d565, d566, d567,
            d575, d576, d577,
            d613, d614, d623, d624, d633, d634, d635, d643, d644, d645, d646,
            d653, d654, d655, d656, d657, d663, d664, d665, d666, d667, d668,
            d675, d676, d677, d678, d686, d687, d688, = cache
    @inbounds begin
        #b1 = b[1]
        #b2 = b[2]
        #b3 = b[3]
        #b4 = b[4]
        #b5 = b[5]
        #b6 = b[6]
        #b7 = b[7]
        #b8 = b[8]

        dest[  1] = α * (
                          (d111*b[1] + d112*b[2] + d113*b[3] + d114*b[4]) * u[1]
                        + (d121*b[1] + d123*b[3] + d124*b[4]) * u[2]
                        + (d131*b[1] + d132*b[2] + d133*b[3] + d134*b[4]) * u[3]
                        + (d141*b[1] + d143*b[3] + d144*b[4]) * u[4]
                        + (d153*b[3] + d154*b[4]) * u[5]
                        + (d163*b[3] + d164*b[4]) * u[6]
                    ) + β*dest[1]
        dest[  2] = α * (
                          (d211*b[1] + d213*b[3] + d214*b[4]) * u[1]
                        + (d221*b[1] + d223*b[3] + d224*b[4]) * u[2]
                        + (d231*b[1] + d233*b[3] + d234*b[4]) * u[3]
                        + (d241*b[1] + d243*b[3] + d244*b[4]) * u[4]
                        + (d253*b[3] + d254*b[4]) * u[5]
                        + (d263*b[3] + d264*b[4]) * u[6]
                    ) + β*dest[2]
        dest[  3] = α * (
                          (d311*b[1] + d312*b[2] + d313*b[3] + d314*b[4]) * u[1]
                        + (d321*b[1] + d323*b[3] + d324*b[4]) * u[2]
                        + (d331*b[1] + d332*b[2] + d333*b[3] + d334*b[4] + d335*b[5]) * u[3]
                        + (d341*b[1] + d343*b[3] + d344*b[4] + d345*b[5]) * u[4]
                        + (d353*b[3] + d354*b[4] + d355*b[5]) * u[5]
                        + (d363*b[3] + d364*b[4] + d365*b[5]) * u[6]
                    ) + β*dest[3]
        dest[  4] = α * (
                          (d411*b[1] + d413*b[3] + d414) * u[1]
                        + (d421*b[1] + d423*b[3] + d424) * u[2]
                        + (d431*b[1] + d433*b[3] + d434*b[4] + d435*b[5]) * u[3]
                        + (d441*b[1] + d443*b[3] + d444*b[4] + d445*b[5] + d446*b[6]) * u[4]
                        + (d453*b[3] + d454*b[4] + d455*b[5] + d456*b[6]) * u[5]
                        + (d463*b[3] + d464*b[4] + d465*b[5] + d466*b[6]) * u[6]
                    ) + β*dest[4]
        dest[  5] = α * (
                          (d513*b[3] + d514) * u[1]
                        + (d523*b[3] + d524) * u[2]
                        + (d533*b[3] + d534*b[4] + d535*b[5]) * u[3]
                        + (d543*b[3] + d544*b[4] + d545*b[5] + d546*b[6]) * u[4]
                        + (d553*b[3] + d554*b[4] + d555*b[5] + d556*b[6] + d557*b[7]) * u[5]
                        + (d563*b[3] + d564*b[4] + d565*b[5] + d566*b[6] + d567*b[7]) * u[6]
                        + (d575*b[5] + d576*b[6] + d577*b[7]) * u[7]
                    ) + β*dest[5]
        dest[  6] = α * (
                          (d613*b[3] + d614) * u[1]
                        + (d623*b[3] + d624) * u[2]
                        + (d633*b[3] + d634*b[4] + d635*b[5]) * u[3]
                        + (d643*b[3] + d644*b[4] + d645*b[5] + d646*b[6]) * u[4]
                        + (d653*b[3] + d654*b[4] + d655*b[5] + d656*b[6] + d657*b[7]) * u[5]
                        + (d663*b[3] + d664*b[4] + d665*b[5] + d666*b[6] + d667*b[7] + d668*b[8]) * u[6]
                        + (d675*b[5] + d676*b[6] + d677*b[7] + d678) * u[7]
                        + (d686*b[6] + d687*b[7] + d688*b[8]) * u[8]
                    ) + β*dest[6]


        dest[end] = α * (
                          (d111*b[end] + d112*b[end-1] + d113*b[end-2] + d114*b[end-3]) * u[end]
                        + (d121*b[end] + d123*b[end-2] + d124*b[end-3]) * u[end-1]
                        + (d131*b[end] + d132*b[end-1] + d133*b[end-2] + d134*b[end-3]) * u[end-2]
                        + (d141*b[end] + d143*b[end-2] + d144*b[end-3]) * u[end-3]
                        + (d153*b[end-2] + d154*b[end-3]) * u[end-4]
                        + (d163*b[end-2] + d164*b[end-3]) * u[end-5]
                    ) + β*dest[end]
        dest[end-1] = α * (
                          (d211*b[end] + d213*b[end-2] + d214*b[end-3]) * u[end]
                        + (d221*b[end] + d223*b[end-2] + d224*b[end-3]) * u[end-1]
                        + (d231*b[end] + d233*b[end-2] + d234*b[end-3]) * u[end-2]
                        + (d241*b[end] + d243*b[end-2] + d244*b[end-3]) * u[end-3]
                        + (d253*b[end-2] + d254*b[end-3]) * u[end-4]
                        + (d263*b[end-2] + d264*b[end-3]) * u[end-5]
                    ) + β*dest[end-1]
        dest[end-2] = α * (
                          (d311*b[end] + d312*b[end-1] + d313*b[end-2] + d314*b[end-3]) * u[end]
                        + (d321*b[end] + d323*b[end-2] + d324*b[end-3]) * u[end-1]
                        + (d331*b[end] + d332*b[end-1] + d333*b[end-2] + d334*b[end-3] + d335*b[end-4]) * u[end-2]
                        + (d341*b[end] + d343*b[end-2] + d344*b[end-3] + d345*b[end-4]) * u[end-3]
                        + (d353*b[end-2] + d354*b[end-3] + d355*b[end-4]) * u[end-4]
                        + (d363*b[end-2] + d364*b[end-3] + d365*b[end-4]) * u[end-5]
                    ) + β*dest[end-2]
        dest[end-3] = α * (
                          (d411*b[end] + d413*b[end-2] + d414) * u[end]
                        + (d421*b[end] + d423*b[end-2] + d424) * u[end-1]
                        + (d431*b[end] + d433*b[end-2] + d434*b[end-3] + d435*b[end-4]) * u[end-2]
                        + (d441*b[end] + d443*b[end-2] + d444*b[end-3] + d445*b[end-4] + d446*b[end-5]) * u[end-3]
                        + (d453*b[end-2] + d454*b[end-3] + d455*b[end-4] + d456*b[end-5]) * u[end-4]
                        + (d463*b[end-2] + d464*b[end-3] + d465*b[end-4] + d466*b[end-5]) * u[end-5]
                    ) + β*dest[end-3]
        dest[end-4] = α * (
                          (d513*b[end-2] + d514) * u[end]
                        + (d523*b[end-2] + d524) * u[end-1]
                        + (d533*b[end-2] + d534*b[end-3] + d535*b[end-4]) * u[end-2]
                        + (d543*b[end-2] + d544*b[end-3] + d545*b[end-4] + d546*b[end-5]) * u[end-3]
                        + (d553*b[end-2] + d554*b[end-3] + d555*b[end-4] + d556*b[end-5] + d557*b[end-6]) * u[end-4]
                        + (d563*b[end-2] + d564*b[end-3] + d565*b[end-4] + d566*b[end-5] + d567*b[end-6]) * u[end-5]
                        + (d575*b[end-4] + d576*b[end-5] + d577*b[end-6]) * u[end-6]
                    ) + β*dest[end-4]
        dest[end-5] = α * (
                          (d613*b[end-2] + d614) * u[end]
                        + (d623*b[end-2] + d624) * u[end-1]
                        + (d633*b[end-2] + d634*b[end-3] + d635*b[end-4]) * u[end-2]
                        + (d643*b[end-2] + d644*b[end-3] + d645*b[end-4] + d646*b[end-5]) * u[end-3]
                        + (d653*b[end-2] + d654*b[end-3] + d655*b[end-4] + d656*b[end-5] + d657*b[end-6]) * u[end-4]
                        + (d663*b[end-2] + d664*b[end-3] + d665*b[end-4] + d666*b[end-5] + d667*b[end-6] + d668*b[end-7]) * u[end-5]
                        + (d675*b[end-4] + d676*b[end-5] + d677*b[end-6] + d678) * u[end-6]
                        + (d686*b[end-5] + d687*b[end-6] + d688*b[end-7]) * u[end-7]
                    ) + β*dest[end-5]
    end
end

@inline function convolve_interior_coefficients_loopbody(i, cache::Mattsson2012Cache4, u, b)
    @unpack dim2m2, dim2m1, dim2_0, dim1m2, dim1m1, dim1_0, dim1p1,
            di_0m2, di_0m1, di_0_0, di_0p1, di_0p2 = cache

    @inbounds begin
        b_im2 = b[i-2]
        b_im1 = b[i-1]
        b_i   = b[i]
        b_ip1 = b[i+1]
        b_ip2 = b[i+2]

        retval = (
                      (dim2m2*b_im2 + dim2m1*b_im1 + dim2_0*b_i) * u[i-2]
                    + (dim1m2*b_im2 + dim1m1*b_im1 + dim1_0*b_i + dim1p1*b_ip1) * u[i-1]
                    + (di_0m2*b_im2 + di_0m1*b_im1 + di_0_0*b_i + di_0p1*b_ip1 + di_0p2*b_ip2) * u[i]
                    + (dim1m2*b_ip2 + dim1m1*b_ip1 + dim1_0*b_i + dim1p1*b_im1) * u[i+1]
                    + (dim2m2*b_ip2 + dim2m1*b_ip1 + dim2_0*b_i) * u[i+2]
                )
    end

    retval
end
