
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



struct Mattsson2012Cache6{T} <: AbstractCoefficientCache{T}
    # boundary coefficients
    d010101::T
    d010102::T
    d010103::T
    d010104::T
    d010105::T
    d010106::T
    d010107::T
    d010201::T
    d010203::T
    d010204::T
    d010205::T
    d010206::T
    d010207::T
    d010301::T
    d010302::T
    d010304::T
    d010305::T
    d010306::T
    d010307::T
    d010401::T
    d010402::T
    d010403::T
    d010405::T
    d010406::T
    d010407::T
    d010501::T
    d010502::T
    d010503::T
    d010504::T
    d010505::T
    d010506::T
    d010507::T
    d010601::T
    d010602::T
    d010603::T
    d010604::T
    d010605::T
    d010606::T
    d010607::T
    d010704::T
    d010705::T
    d010706::T
    d010707::T
    d010805::T
    d010806::T
    d010807::T
    d010905::T
    d010906::T
    d010907::T
    d020101::T
    d020103::T
    d020104::T
    d020105::T
    d020106::T
    d020107::T
    d020201::T
    d020203::T
    d020204::T
    d020205::T
    d020206::T
    d020207::T
    d020301::T
    d020304::T
    d020305::T
    d020306::T
    d020307::T
    d020401::T
    d020403::T
    d020405::T
    d020406::T
    d020407::T
    d020501::T
    d020503::T
    d020504::T
    d020505::T
    d020506::T
    d020507::T
    d020601::T
    d020603::T
    d020604::T
    d020605::T
    d020606::T
    d020607::T
    d020704::T
    d020705::T
    d020706::T
    d020707::T
    d020805::T
    d020806::T
    d020807::T
    d020905::T
    d020906::T
    d020907::T
    d030101::T
    d030102::T
    d030104::T
    d030105::T
    d030106::T
    d030107::T
    d030201::T
    d030204::T
    d030205::T
    d030206::T
    d030207::T
    d030301::T
    d030302::T
    d030304::T
    d030305::T
    d030306::T
    d030307::T
    d030401::T
    d030402::T
    d030405::T
    d030406::T
    d030407::T
    d030501::T
    d030502::T
    d030504::T
    d030505::T
    d030506::T
    d030507::T
    d030601::T
    d030602::T
    d030604::T
    d030605::T
    d030606::T
    d030607::T
    d030704::T
    d030705::T
    d030706::T
    d030707::T
    d030805::T
    d030806::T
    d030807::T
    d030905::T
    d030906::T
    d030907::T
    d040101::T
    d040102::T
    d040103::T
    d040105::T
    d040106::T
    d040107::T
    d040201::T
    d040203::T
    d040205::T
    d040206::T
    d040207::T
    d040301::T
    d040302::T
    d040305::T
    d040306::T
    d040307::T
    d040401::T
    d040402::T
    d040403::T
    d040405::T
    d040406::T
    d040407::T
    d040501::T
    d040502::T
    d040503::T
    d040505::T
    d040506::T
    d040507::T
    d040601::T
    d040602::T
    d040603::T
    d040605::T
    d040606::T
    d040607::T
    d040705::T
    d040706::T
    d040707::T
    d040805::T
    d040806::T
    d040807::T
    d040905::T
    d040906::T
    d040907::T
    d050101::T
    d050102::T
    d050103::T
    d050104::T
    d050105::T
    d050106::T
    d050107::T
    d050201::T
    d050203::T
    d050204::T
    d050205::T
    d050206::T
    d050207::T
    d050301::T
    d050302::T
    d050304::T
    d050305::T
    d050306::T
    d050307::T
    d050401::T
    d050402::T
    d050403::T
    d050405::T
    d050406::T
    d050407::T
    d050501::T
    d050502::T
    d050503::T
    d050504::T
    d050505::T
    d050506::T
    d050507::T
    d050508::T
    d050601::T
    d050602::T
    d050603::T
    d050604::T
    d050605::T
    d050606::T
    d050607::T
    d050608::T
    d050704::T
    d050705::T
    d050706::T
    d050707::T
    d050708::T
    d050805::T
    d050806::T
    d050807::T
    d050808::T
    d050905::T
    d050906::T
    d050907::T
    d050908::T
    d060101::T
    d060102::T
    d060103::T
    d060104::T
    d060105::T
    d060106::T
    d060107::T
    d060201::T
    d060203::T
    d060204::T
    d060205::T
    d060206::T
    d060207::T
    d060301::T
    d060302::T
    d060304::T
    d060305::T
    d060306::T
    d060307::T
    d060401::T
    d060402::T
    d060403::T
    d060405::T
    d060406::T
    d060407::T
    d060501::T
    d060502::T
    d060503::T
    d060504::T
    d060505::T
    d060506::T
    d060507::T
    d060508::T
    d060601::T
    d060602::T
    d060603::T
    d060604::T
    d060605::T
    d060606::T
    d060607::T
    d060608::T
    d060609::T
    d060704::T
    d060705::T
    d060706::T
    d060707::T
    d060708::T
    d060709::T
    d060805::T
    d060806::T
    d060807::T
    d060808::T
    d060809::T
    d060905::T
    d060906::T
    d060907::T
    d060908::T
    d060909::T
    d070104::T
    d070105::T
    d070106::T
    d070107::T
    d070204::T
    d070205::T
    d070206::T
    d070207::T
    d070304::T
    d070305::T
    d070306::T
    d070307::T
    d070405::T
    d070406::T
    d070407::T
    d070504::T
    d070505::T
    d070506::T
    d070507::T
    d070508::T
    d070604::T
    d070605::T
    d070606::T
    d070607::T
    d070608::T
    d070609::T
    d070704::T
    d070705::T
    d070706::T
    d070707::T
    d070708::T
    d070709::T
    d070710::T
    d070805::T
    d070806::T
    d070807::T
    d070808::T
    d070809::T
    d070810::T
    d070905::T
    d070906::T
    d070907::T
    d070908::T
    d070909::T
    d070910::T
    d071010::T
    d071007::T
    d071008::T
    d071009::T
    d080105::T
    d080106::T
    d080107::T
    d080205::T
    d080206::T
    d080207::T
    d080305::T
    d080306::T
    d080307::T
    d080405::T
    d080406::T
    d080407::T
    d080505::T
    d080506::T
    d080507::T
    d080508::T
    d080605::T
    d080606::T
    d080607::T
    d080608::T
    d080609::T
    d080705::T
    d080706::T
    d080707::T
    d080708::T
    d080709::T
    d080710::T
    d080805::T
    d080806::T
    d080807::T
    d080808::T
    d080809::T
    d080810::T
    d080811::T
    d080910::T
    d080911::T
    d080905::T
    d080906::T
    d080907::T
    d080908::T
    d080909::T
    d081010::T
    d081011::T
    d081007::T
    d081008::T
    d081009::T
    d081110::T
    d081111::T
    d081108::T
    d081109::T
    d090105::T
    d090106::T
    d090107::T
    d090205::T
    d090206::T
    d090207::T
    d090305::T
    d090306::T
    d090307::T
    d090405::T
    d090406::T
    d090407::T
    d090505::T
    d090506::T
    d090507::T
    d090508::T
    d090605::T
    d090606::T
    d090607::T
    d090608::T
    d090609::T
    d090705::T
    d090706::T
    d090707::T
    d090708::T
    d090709::T
    d090710::T
    d090810::T
    d090811::T
    d090805::T
    d090806::T
    d090807::T
    d090808::T
    d090809::T
    d090910::T
    d090911::T
    d090912::T
    d090905::T
    d090906::T
    d090907::T
    d090908::T
    d090909::T
    d091010::T
    d091011::T
    d091012::T
    d091007::T
    d091008::T
    d091009::T
    d091110::T
    d091111::T
    d091112::T
    d091108::T
    d091109::T
    d091210::T
    d091211::T
    d091212::T
    d091209::T
    # interior coefficients
    di_im3_im3::T
    di_im3_im2::T
    di_im3_im1::T
    di_im3_i_0::T
    di_im2_im3::T
    di_im2_im2::T
    di_im2_im1::T
    di_im2_i_0::T
    di_im2_ip1::T
    di_im1_im3::T
    di_im1_im2::T
    di_im1_im1::T
    di_im1_i_0::T
    di_im1_ip1::T
    di_im1_ip2::T
    di_i_0_im3::T
    di_i_0_im2::T
    di_i_0_im1::T
    di_i_0_i_0::T
    di_i_0_ip1::T
    di_i_0_ip2::T
    di_i_0_ip3::T


    function Mattsson2012Cache6(::Type{T}) where {T}
        d010101 = T(3488579749432192834647//853062500000000000000)
        d010102 = T(-801487464472260200583//853062500000000000000)
        d010103 = T(-107511269246879479263//10663281250000000000000)
        d010104 = T(-220374545749778513709//4265312500000000000000)
        d010105 = T(-853281605905492073679//8530625000000000000000)
        d010106 = T(-213837620491087107723//2132656250000000000000)
        d010107 = T(-170023197633151711707//1706125000000000000000)
        d010201 = T(-201373469593852175253//21326562500000000000)
        d010203 = T(76831939003273506963//853062500000000000000)
        d010204 = T(557284028177155340673//4265312500000000000000)
        d010205 = T(15036639677497517481//34122500000000000000)
        d010206 = T(161429938957512238491//426531250000000000000)
        d010207 = T(322452294262920020943//853062500000000000000)
        d010301 = T(79089576406944196937499//8530625000000000000000)
        d010302 = T(498278548626520783581//853062500000000000000)
        d010304 = T(117832404972021093369//853062500000000000000)
        d010305 = T(-656739334946087764443//853062500000000000000)
        d010306 = T(-10679357746376204247//21326562500000000000)
        d010307 = T(-108158529039623075763//213265625000000000000)
        d010401 = T(-4208037713880665367429//853062500000000000000)
        d010402 = T(448942151996204821731//853062500000000000000)
        d010403 = T(-731843359645141250379//8530625000000000000000)
        d010405 = T(517494470211478775457//853062500000000000000)
        d010406 = T(259338954156716742561//1066328125000000000000)
        d010407 = T(2219288476454558345493//8530625000000000000000)
        d010501 = T(8160588935334555205509//8530625000000000000000)
        d010502 = T(-119881726609766657967//853062500000000000000)
        d010503 = T(55294733554706945223//17061250000000000000000)
        d010504 = T(-399788170847073400887//1706125000000000000000)
        d010505 = T(-1940936268931002346329//8530625000000000000000)
        d010506 = T(-167619451814797080879//4265312500000000000000)
        d010507 = T(-186214004609139292959//4265312500000000000000)
        d010601 = T(493802139746843285229//8530625000000000000000)
        d010602 = T(-2585150954069874676389//85306250000000000000000)
        d010603 = T(2188561823255629155033//853062500000000000000000)
        d010604 = T(248023538587423789821//10663281250000000000000)
        d010605 = T(-28717252642545488919//853062500000000000000)
        d010606 = T(22206943772706685203//533164062500000000000)
        d010607 = T(285921424643231725821//4265312500000000000000)
        d010704 = T(-258104956047683075199//42653125000000000000000)
        d010705 = T(1098395195992870482747//8530625000000000000000)
        d010706 = T(-35631224488991227509//853062500000000000000)
        d010707 = T(-690461959893421592601//8530625000000000000000)
        d010805 = T(-10527656880173303907//213265625000000000000)
        d010806 = T(350253944495903046573//17061250000000000000000)
        d010807 = T(2459793029589806329827//85306250000000000000000)
        d010905 = T(373902084477215020521//85306250000000000000000)
        d010906 = T(-1510375598090063026011//853062500000000000000000)
        d010907 = T(-2228645246682087179091//853062500000000000000000)
        d020101 = T(68626530406147824747//93851562500000000000)
        d020103 = T(76831939003273506963//3754062500000000000000)
        d020104 = T(557284028177155340673//18770312500000000000000)
        d020105 = T(15036639677497517481//150162500000000000000)
        d020106 = T(161429938957512238491//1877031250000000000000)
        d020107 = T(322452294262920020943//3754062500000000000000)
        d020201 = T(-176354841210075129861//187703125000000000000)
        d020203 = T(-137268095064696670149//750812500000000000000)
        d020204 = T(-56370482716998380913//750812500000000000000)
        d020205 = T(-1801528565678670496743//3754062500000000000000)
        d020206 = T(-1264091227125142349097//3754062500000000000000)
        d020207 = T(-1262632160752875968427//3754062500000000000000)
        d020301 = T(2454680772897965443323//37540625000000000000000)
        d020304 = T(-297975054556298634831//3754062500000000000000)
        d020305 = T(34840733714933009949//37540625000000000000)
        d020306 = T(1792693417158462092517//3754062500000000000000)
        d020307 = T(893156552640702834867//1877031250000000000000)
        d020401 = T(390630248413213274991//1877031250000000000000)
        d020403 = T(40859775068590002177//234628906250000000000)
        d020405 = T(-320670968857568377479//375406250000000000000)
        d020406 = T(-8345949222454037577//29328613281250000000)
        d020407 = T(-1042213208874066000189//3754062500000000000000)
        d020501 = T(-181244910846701627703//3754062500000000000000)
        d020503 = T(-12348710120802100743//1877031250000000000000)
        d020504 = T(505492958611028165859//3754062500000000000000)
        d020505 = T(16040554255356295239//45781250000000000000)
        d020506 = T(356226810890357123127//3754062500000000000000)
        d020507 = T(282829404314370997287//3754062500000000000000)
        d020601 = T(-317240276569659297369//18770312500000000000000)
        d020603 = T(-488761113390649738857//93851562500000000000000)
        d020604 = T(-501762329025753451149//37540625000000000000000)
        d020605 = T(2508196356163442567061//37540625000000000000000)
        d020606 = T(-33034858301307807009//375406250000000000000)
        d020607 = T(-300650486517799659153//3754062500000000000000)
        d020704 = T(652696839870332524731//187703125000000000000000)
        d020705 = T(-623884634825711508837//3754062500000000000000)
        d020706 = T(291809018392922157759//3754062500000000000000)
        d020707 = T(159510839817691350279//1877031250000000000000)
        d020805 = T(2259398187288932945499//37540625000000000000000)
        d020806 = T(-138155607003388094829//4692578125000000000000)
        d020807 = T(-1154153331261828186867//37540625000000000000000)
        d020905 = T(-1995137763512691244167//375406250000000000000000)
        d020906 = T(482833592914092724359//187703125000000000000000)
        d020907 = T(1029470577684505795449//375406250000000000000000)
        d030101 = T(-1910423593055803062501//16943750000000000000000)
        d030102 = T(498278548626520783581//1694375000000000000000)
        d030104 = T(117832404972021093369//1694375000000000000000)
        d030105 = T(-656739334946087764443//1694375000000000000000)
        d030106 = T(-10679357746376204247//42359375000000000000)
        d030107 = T(-108158529039623075763//423593750000000000000)
        d030201 = T(2454680772897965443323//16943750000000000000000)
        d030204 = T(-297975054556298634831//1694375000000000000000)
        d030205 = T(34840733714933009949//16943750000000000000)
        d030206 = T(1792693417158462092517//1694375000000000000000)
        d030207 = T(893156552640702834867//847187500000000000000)
        d030301 = T(-854166754863477964053//84718750000000000000000)
        d030302 = T(-12391036567731380097//67775000000000000000)
        d030304 = T(-12600797986220544081//67775000000000000000)
        d030305 = T(-186746229559267472511//42359375000000000000)
        d030306 = T(-144574293144710309067//84718750000000000000)
        d030307 = T(-68388393894725681361//42359375000000000000)
        d030401 = T(-271858870874447628291//8471875000000000000000)
        d030402 = T(-279103858550339328459//1694375000000000000000)
        d030405 = T(786569837710970535561//169437500000000000000)
        d030406 = T(1182033991370105405073//847187500000000000000)
        d030407 = T(373137687908395871031//338875000000000000000)
        d030501 = T(1261372794725374974441//169437500000000000000000)
        d030502 = T(372647704361091231753//8471875000000000000000)
        d030504 = T(106881676125526978623//338875000000000000000)
        d030505 = T(-334577950118615991033//169437500000000000000)
        d030506 = T(-1458806041360826924121//1694375000000000000000)
        d030507 = T(-709570273566786721371//1694375000000000000000)
        d030601 = T(88313266813781942253//33887500000000000000000)
        d030602 = T(1607168324488480097319//169437500000000000000000)
        d030604 = T(-265232174995477154619//8471875000000000000000)
        d030605 = T(-436946401706271291081//847187500000000000000)
        d030606 = T(1258069507171795162953//1694375000000000000000)
        d030607 = T(299331817327134476019//847187500000000000000)
        d030704 = T(172508170140646016049//21179687500000000000000)
        d030705 = T(72304355210959183671//84718750000000000000)
        d030706 = T(-112331309776528689009//211796875000000000000)
        d030707 = T(-561237279618205842633//1694375000000000000000)
        d030805 = T(-246284303508129090897//847187500000000000000)
        d030806 = T(286150328153140980393//1694375000000000000000)
        d030807 = T(206418278863117201401//1694375000000000000000)
        d030905 = T(429705761095573155117//16943750000000000000000)
        d030906 = T(-310806786266228232639//21179687500000000000000)
        d030907 = T(-1810603320825905690031//169437500000000000000000)
        d040101 = T(-608037713880665367429//3349375000000000000000)
        d040102 = T(448942151996204821731//3349375000000000000000)
        d040103 = T(-731843359645141250379//33493750000000000000000)
        d040105 = T(517494470211478775457//3349375000000000000000)
        d040106 = T(259338954156716742561//4186718750000000000000)
        d040107 = T(2219288476454558345493//33493750000000000000000)
        d040201 = T(390630248413213274991//1674687500000000000000)
        d040203 = T(40859775068590002177//209335937500000000000)
        d040205 = T(-320670968857568377479//334937500000000000000)
        d040206 = T(-8345949222454037577//26166992187500000000)
        d040207 = T(-1042213208874066000189//3349375000000000000000)
        d040301 = T(-271858870874447628291//16746875000000000000000)
        d040302 = T(-279103858550339328459//3349375000000000000000)
        d040305 = T(786569837710970535561//334937500000000000000)
        d040306 = T(1182033991370105405073//1674687500000000000000)
        d040307 = T(373137687908395871031//669875000000000000000)
        d040401 = T(-1730510939516540486409//33493750000000000000000)
        d040402 = T(-2514687562477268986203//33493750000000000000000)
        d040403 = T(-622719258651382291857//3349375000000000000000)
        d040405 = T(-996148883236603487481//334937500000000000000)
        d040406 = T(-160724295722827197963//167468750000000000000)
        d040407 = T(-798184738829995630119//1674687500000000000000)
        d040501 = T(200730935600856557361//16746875000000000000000)
        d040502 = T(335750480917671703167//16746875000000000000000)
        d040503 = T(2352490805082173630517//334937500000000000000000)
        d040505 = T(99271705008091699449//83734375000000000000)
        d040506 = T(341147116090392535539//334937500000000000000)
        d040507 = T(825745697417879278263//3349375000000000000000)
        d040601 = T(43918431439409625861//10466796875000000000000)
        d040602 = T(144803666183270647203//33493750000000000000000)
        d040603 = T(1862228546763464577501//334937500000000000000000)
        d040605 = T(92240505254357022933//145625000000000000000)
        d040606 = T(-673451844049760660721//837343750000000000000)
        d040607 = T(-618609961464408342027//3349375000000000000000)
        d040705 = T(-452067580444567202901//837343750000000000000)
        d040706 = T(135357678613947454353//334937500000000000000)
        d040707 = T(484693535638794268071//3349375000000000000000)
        d040805 = T(136546361557543230651//837343750000000000000)
        d040806 = T(-196211950838593897311//1674687500000000000000)
        d040807 = T(-48050482672807852491//1046679687500000000000)
        d040905 = T(-453092720052490978833//33493750000000000000000)
        d040906 = T(324115996979071035099//33493750000000000000000)
        d040907 = T(257953446146839887441//66987500000000000000000)
        d050101 = T(1410588935334555205509//24615625000000000000000)
        d050102 = T(-119881726609766657967//2461562500000000000000)
        d050103 = T(55294733554706945223//49231250000000000000000)
        d050104 = T(-399788170847073400887//4923125000000000000000)
        d050105 = T(-1940936268931002346329//24615625000000000000000)
        d050106 = T(-167619451814797080879//12307812500000000000000)
        d050107 = T(-186214004609139292959//12307812500000000000000)
        d050201 = T(-181244910846701627703//2461562500000000000000)
        d050203 = T(-12348710120802100743//1230781250000000000000)
        d050204 = T(505492958611028165859//2461562500000000000000)
        d050205 = T(657662724469608104799//1230781250000000000000)
        d050206 = T(356226810890357123127//2461562500000000000000)
        d050207 = T(282829404314370997287//2461562500000000000000)
        d050301 = T(1261372794725374974441//246156250000000000000000)
        d050302 = T(372647704361091231753//12307812500000000000000)
        d050304 = T(106881676125526978623//492312500000000000000)
        d050305 = T(-334577950118615991033//246156250000000000000)
        d050306 = T(-1458806041360826924121//2461562500000000000000)
        d050307 = T(-709570273566786721371//2461562500000000000000)
        d050401 = T(200730935600856557361//12307812500000000000000)
        d050402 = T(335750480917671703167//12307812500000000000000)
        d050403 = T(2352490805082173630517//246156250000000000000000)
        d050405 = T(99271705008091699449//61539062500000000000)
        d050406 = T(341147116090392535539//246156250000000000000)
        d050407 = T(825745697417879278263//2461562500000000000000)
        d050501 = T(-232838218973940769701//61539062500000000000000)
        d050502 = T(-224139943105689130929//30769531250000000000000)
        d050503 = T(-177743434496680594359//492312500000000000000000)
        d050504 = T(-226646215261412211//615390625000000000)
        d050505 = T(-479865994921127926953//615390625000000000000)
        d050506 = T(-56603159898633713853//24615625000000000000)
        d050507 = T(-12513545302588885017//49231250000000000000)
        d050508 = T(-329822591275886287131//49231250000000000000000)
        d050601 = T(-326036899522132327647//246156250000000000000000)
        d050602 = T(-15466859990123194161//9846250000000000000000)
        d050603 = T(-703507314469479108471//2461562500000000000000000)
        d050604 = T(224973525144485388879//6153906250000000000000)
        d050605 = T(-141191833817194468299//1230781250000000000000)
        d050606 = T(45509748300015149583//30769531250000000000)
        d050607 = T(189556314573523456143//615390625000000000000)
        d050608 = T(183661295637943143567//6153906250000000000000)
        d050704 = T(-2341180282727167008777//246156250000000000000000)
        d050705 = T(327366260424611050743//1230781250000000000000)
        d050706 = T(-468739821065965712199//2461562500000000000000)
        d050707 = T(-211817059786609416531//1230781250000000000000)
        d050708 = T(-719733886913829430701//12307812500000000000000)
        d050805 = T(-2186512419883225793031//24615625000000000000000)
        d050806 = T(58263866866465575363//615390625000000000000)
        d050807 = T(-145737487465433959149//4923125000000000000000)
        d050808 = T(884645182551772574253//24615625000000000000000)
        d050905 = T(1944584727437707973811//246156250000000000000000)
        d050906 = T(-2014756060823330869791//246156250000000000000000)
        d050907 = T(274105362206317914549//307695312500000000000000)
        d050908 = T(-1491129563794314356511//2461562500000000000000000)
        d060101 = T(493802139746843285229//27375625000000000000000)
        d060102 = T(-2585150954069874676389//273756250000000000000000)
        d060103 = T(2188561823255629155033//2737562500000000000000000)
        d060104 = T(248023538587423789821//34219531250000000000000)
        d060105 = T(-28717252642545488919//2737562500000000000000)
        d060106 = T(22206943772706685203//1710976562500000000000)
        d060107 = T(285921424643231725821//13687812500000000000000)
        d060201 = T(-317240276569659297369//13687812500000000000000)
        d060203 = T(-488761113390649738857//68439062500000000000000)
        d060204 = T(-501762329025753451149//27375625000000000000000)
        d060205 = T(2508196356163442567061//27375625000000000000000)
        d060206 = T(-33034858301307807009//273756250000000000000)
        d060207 = T(-300650486517799659153//2737562500000000000000)
        d060301 = T(88313266813781942253//54751250000000000000000)
        d060302 = T(1607168324488480097319//273756250000000000000000)
        d060304 = T(-265232174995477154619//13687812500000000000000)
        d060305 = T(-436946401706271291081//1368781250000000000000)
        d060306 = T(1258069507171795162953//2737562500000000000000)
        d060307 = T(299331817327134476019//1368781250000000000000)
        d060401 = T(43918431439409625861//8554882812500000000000)
        d060402 = T(144803666183270647203//27375625000000000000000)
        d060403 = T(1862228546763464577501//273756250000000000000000)
        d060405 = T(2121531620850211527459//2737562500000000000000)
        d060406 = T(-673451844049760660721//684390625000000000000)
        d060407 = T(-618609961464408342027//2737562500000000000000)
        d060501 = T(-326036899522132327647//273756250000000000000000)
        d060502 = T(-15466859990123194161//10950250000000000000000)
        d060503 = T(-703507314469479108471//2737562500000000000000000)
        d060504 = T(224973525144485388879//6843906250000000000000)
        d060505 = T(-141191833817194468299//1368781250000000000000)
        d060506 = T(45509748300015149583//34219531250000000000)
        d060507 = T(189556314573523456143//684390625000000000000)
        d060508 = T(183661295637943143567//6843906250000000000000)
        d060601 = T(-45654042673252926459//109502500000000000000000)
        d060602 = T(-52114082811395024331//171097656250000000000000)
        d060603 = T(-278447720397403133727//1368781250000000000000000)
        d060604 = T(-893252719134261429309//273756250000000000000000)
        d060605 = T(-786786347214675245133//1368781250000000000000)
        d060606 = T(-321397759402711123797//273756250000000000000)
        d060607 = T(-135914189968679886117//171097656250000000000)
        d060608 = T(-361358073020709029691//2737562500000000000000)
        d060609 = T(-240//43801)
        d060704 = T(1161950974112665795851//1368781250000000000000000)
        d060705 = T(471980171341954347111//2737562500000000000000)
        d060706 = T(1806733914744940873497//2737562500000000000000)
        d060707 = T(855674902433815903251//2737562500000000000000)
        d060708 = T(16915742190621270891//54751250000000000000)
        d060709 = T(1080//43801)
        d060805 = T(-905758683757205069349//27375625000000000000000)
        d060806 = T(-146460866162013185679//684390625000000000000)
        d060807 = T(18055548120889645587//54751250000000000000)
        d060808 = T(-496358073020709029691//2737562500000000000000)
        d060809 = T(-2160//43801)
        d060905 = T(240542934301893680487//136878125000000000000000)
        d060906 = T(165725243866011988647//5475125000000000000000)
        d060907 = T(-78637998874221125361//2737562500000000000000)
        d060908 = T(-153838704362056856433//6843906250000000000000)
        d060909 = T(206250000000000000003//6843906250000000000000)
        d070104 = T(-9559442816580854637//5000000000000000000000)
        d070105 = T(40681303555291499361//1000000000000000000000)
        d070106 = T(-1319674981073749167//100000000000000000000)
        d070107 = T(-25572665181237836763//1000000000000000000000)
        d070204 = T(24173957032234537953//5000000000000000000000)
        d070205 = T(-23106838326878204031//100000000000000000000)
        d070206 = T(10807741421960079917//100000000000000000000)
        d070207 = T(5907808882136716677//50000000000000000000)
        d070304 = T(6389191486690593187//1250000000000000000000)
        d070305 = T(2677939081887377173//5000000000000000000)
        d070306 = T(-4160418880612173667//12500000000000000000)
        d070307 = T(-20786565911785401579//100000000000000000000)
        d070405 = T(-16743243720169155663//25000000000000000000)
        d070406 = T(5013247356072127939//10000000000000000000)
        d070407 = T(17951612431066454373//100000000000000000000)
        d070504 = T(-86710380841746926251//10000000000000000000000)
        d070505 = T(12124676312022631509//50000000000000000000)
        d070506 = T(-17360734113554285637//100000000000000000000)
        d070507 = T(-7845076288392941353//50000000000000000000)
        d070508 = T(-26656810626438127063//500000000000000000000)
        d070604 = T(43035221263432066513//50000000000000000000000)
        d070605 = T(17480747086739049893//100000000000000000000)
        d070606 = T(66916070916479291611//100000000000000000000)
        d070607 = T(31691663053104292713//100000000000000000000)
        d070608 = T(626508970023010033//2000000000000000000)
        d070609 = T(1//40)
        d070704 = T(-2239223735771599179//10000000000000000000000)
        d070705 = T(-6377188927154783369//50000000000000000000)
        d070706 = T(-5058497419648040823//5000000000000000000)
        d070707 = T(-48231775430312815001//100000000000000000000)
        d070708 = T(-3879526910069030099//4000000000000000000)
        d070709 = T(-1//8)
        d070710 = T(-1//180)
        d070805 = T(37841139730330129499//1000000000000000000000)
        d070806 = T(1873473053209267101//6250000000000000000)
        d070807 = T(4989358584308526473//12500000000000000000)
        d070808 = T(876508970023010033//2000000000000000000)
        d070809 = T(3//10)
        d070810 = T(1//40)
        d070905 = T(-3836682690370322953//1250000000000000000000)
        d070906 = T(-46981462180226839339//1000000000000000000000)
        d070907 = T(-17163557041460064817//100000000000000000000)
        d070908 = T(29668637874712374587//100000000000000000000)
        d070909 = T(-7//40)
        d070910 = T(-1//20)
        d071010 = T(11//360)
        d071007 = T(11//360)
        d071008 = T(-1//40)
        d071009 = T(-1//40)
        d080105 = T(-389913217784196441//25000000000000000000)
        d080106 = T(12972368314663075799//2000000000000000000000)
        d080107 = T(91103445540363197401//10000000000000000000000)
        d080205 = T(83681414344034553537//1000000000000000000000)
        d080206 = T(-5116874333458818327//125000000000000000000)
        d080207 = T(-42746419676364006921//1000000000000000000000)
        d080305 = T(-9121640870671447811//50000000000000000000)
        d080306 = T(10598160301968184459//100000000000000000000)
        d080307 = T(7645121439374711163//100000000000000000000)
        d080405 = T(5057272650279378913//25000000000000000000)
        d080406 = T(-7267109290318292493//50000000000000000000)
        d080407 = T(-1779647506400290833//31250000000000000000)
        d080505 = T(-80981941477156510853//1000000000000000000000)
        d080506 = T(2157920995054280569//25000000000000000000)
        d080507 = T(-5397684720941998487//200000000000000000000)
        d080508 = T(32764636390806391639//1000000000000000000000)
        d080605 = T(-33546617916933521087//1000000000000000000000)
        d080606 = T(-5424476524519006877//25000000000000000000)
        d080607 = T(668724004477394281//2000000000000000000)
        d080608 = T(-18383632334100334433//100000000000000000000)
        d080609 = T(-1//20)
        d080705 = T(37841139730330129499//1000000000000000000000)
        d080706 = T(1873473053209267101//6250000000000000000)
        d080707 = T(4989358584308526473//12500000000000000000)
        d080708 = T(876508970023010033//2000000000000000000)
        d080709 = T(3//10)
        d080710 = T(1//40)
        d080805 = T(-6151644713584022277//500000000000000000000)
        d080806 = T(-473459011858359333//4000000000000000000)
        d080807 = T(-47052559491139716671//50000000000000000000)
        d080808 = T(-11398948689042289109//20000000000000000000)
        d080809 = T(-19//20)
        d080810 = T(-1//8)
        d080811 = T(-1//180)
        d080910 = T(3//10)
        d080911 = T(1//40)
        d080905 = T(5238674302575254013//5000000000000000000000)
        d080906 = T(5770169731679790849//250000000000000000000)
        d080907 = T(7466562634437873743//25000000000000000000)
        d080908 = T(5340113510440635451//12500000000000000000)
        d080909 = T(17//40)
        d081010 = T(-7//40)
        d081011 = T(-1//20)
        d081007 = T(-1//20)
        d081008 = T(-7//40)
        d081009 = T(3//10)
        d081110 = T(-1//40)
        d081111 = T(11//360)
        d081108 = T(11//360)
        d081109 = T(-1//40)
        d090105 = T(13848225351007963723//10000000000000000000000)
        d090106 = T(-55939836966298630593//100000000000000000000000)
        d090107 = T(-82542416543781006633//100000000000000000000000)
        d090205 = T(-73893991241210786821//10000000000000000000000)
        d090206 = T(17882725663484915717//5000000000000000000000)
        d090207 = T(38128539914240955387//10000000000000000000000)
        d090305 = T(15915028188724931671//1000000000000000000000)
        d090306 = T(-11511362454304749357//1250000000000000000000)
        d090307 = T(-67059382252811321853//10000000000000000000000)
        d090405 = T(-16781211853795962179//1000000000000000000000)
        d090406 = T(12004296184410038337//1000000000000000000000)
        d090407 = T(9553831338771847683//2000000000000000000000)
        d090505 = T(72021656571766961993//10000000000000000000000)
        d090506 = T(-74620594845308550733//10000000000000000000000)
        d090507 = T(10152050452085848687//12500000000000000000000)
        d090508 = T(-55227020881270902093//100000000000000000000000)
        d090605 = T(8908997566736802981//5000000000000000000000)
        d090606 = T(6137971995037481061//200000000000000000000)
        d090607 = T(-2912518476823004643//100000000000000000000)
        d090608 = T(-5697729791187290979//250000000000000000000)
        d090609 = T(7638888888888888889//250000000000000000000)
        d090705 = T(-3836682690370322953//1250000000000000000000)
        d090706 = T(-46981462180226839339//1000000000000000000000)
        d090707 = T(-17163557041460064817//100000000000000000000)
        d090708 = T(29668637874712374587//100000000000000000000)
        d090709 = T(-7//40)
        d090710 = T(-1//20)
        d090810 = T(3//10)
        d090811 = T(1//40)
        d090805 = T(5238674302575254013//5000000000000000000000)
        d090806 = T(5770169731679790849//250000000000000000000)
        d090807 = T(7466562634437873743//25000000000000000000)
        d090808 = T(5340113510440635451//12500000000000000000)
        d090809 = T(17//40)
        d090910 = T(-19//20)
        d090911 = T(-1//8)
        d090912 = T(-1//180)
        d090905 = T(-91593624651536418269//1000000000000000000000000)
        d090906 = T(-5139370221149109977//1000000000000000000000)
        d090907 = T(-6238616075047110007//50000000000000000000)
        d090908 = T(-47527613510440635451//50000000000000000000)
        d090909 = T(-56111111111111111111//100000000000000000000)
        d091010 = T(17//40)
        d091011 = T(3//10)
        d091012 = T(1//40)
        d091007 = T(1//40)
        d091008 = T(3//10)
        d091009 = T(17//40)
        d091110 = T(3//10)
        d091111 = T(-7//40)
        d091112 = T(-1//20)
        d091108 = T(-1//20)
        d091109 = T(-7//40)
        d091210 = T(-1//40)
        d091211 = T(-1//40)
        d091212 = T(11//360)
        d091209 = T(11//360)

        di_im3_im3 = T(11//360)
        di_im3_im2 = T(-1//40)
        di_im3_im1 = T(-1//40)
        di_im3_i_0 = T(11//360)
        di_im2_im3 = T(-1//20)
        di_im2_im2 = T(-7//40)
        di_im2_im1 = T(3//10)
        di_im2_i_0 = T(-7//40)
        di_im2_ip1 = T(-1//20)
        di_im1_im3 = T(1//40)
        di_im1_im2 = T(3//10)
        di_im1_im1 = T(17//40)
        di_im1_i_0 = T(17//40)
        di_im1_ip1 = T(3//10)
        di_im1_ip2 = T(1//40)
        di_i_0_im3 = T(-1//180)
        di_i_0_im2 = T(-1//8)
        di_i_0_im1 = T(-19//20)
        di_i_0_i_0 = T(-101//180)
        di_i_0_ip1 = T(-19//20)
        di_i_0_ip2 = T(-1//8)
        di_i_0_ip3 = T(-1//180)

        new{T}( d010101, d010102, d010103, d010104, d010105, d010106, d010107, 
        d010201, d010203, d010204, d010205, d010206, d010207, 
        d010301, d010302, d010304, d010305, d010306, d010307, 
        d010401, d010402, d010403, d010405, d010406, d010407, 
        d010501, d010502, d010503, d010504, d010505, d010506, d010507, 
        d010601, d010602, d010603, d010604, d010605, d010606, d010607, 
        d010704, d010705, d010706, d010707, 
        d010805, d010806, d010807, 
        d010905, d010906, d010907, 
        d020101, d020103, d020104, d020105, d020106, d020107, 
        d020201, d020203, d020204, d020205, d020206, d020207, 
        d020301, d020304, d020305, d020306, d020307, 
        d020401, d020403, d020405, d020406, d020407, 
        d020501, d020503, d020504, d020505, d020506, d020507, 
        d020601, d020603, d020604, d020605, d020606, d020607, 
        d020704, d020705, d020706, d020707, 
        d020805, d020806, d020807, 
        d020905, d020906, d020907, 
        d030101, d030102, d030104, d030105, d030106, d030107, 
        d030201, d030204, d030205, d030206, d030207, 
        d030301, d030302, d030304, d030305, d030306, d030307, 
        d030401, d030402, d030405, d030406, d030407, 
        d030501, d030502, d030504, d030505, d030506, d030507, 
        d030601, d030602, d030604, d030605, d030606, d030607, 
        d030704, d030705, d030706, d030707, 
        d030805, d030806, d030807, 
        d030905, d030906, d030907, 
        d040101, d040102, d040103, d040105, d040106, d040107, 
        d040201, d040203, d040205, d040206, d040207, 
        d040301, d040302, d040305, d040306, d040307, 
        d040401, d040402, d040403, d040405, d040406, d040407, 
        d040501, d040502, d040503, d040505, d040506, d040507, 
        d040601, d040602, d040603, d040605, d040606, d040607, 
        d040705, d040706, d040707, 
        d040805, d040806, d040807, 
        d040905, d040906, d040907, 
        d050101, d050102, d050103, d050104, d050105, d050106, d050107, 
        d050201, d050203, d050204, d050205, d050206, d050207, 
        d050301, d050302, d050304, d050305, d050306, d050307, 
        d050401, d050402, d050403, d050405, d050406, d050407, 
        d050501, d050502, d050503, d050504, d050505, d050506, d050507, d050508, 
        d050601, d050602, d050603, d050604, d050605, d050606, d050607, d050608, 
        d050704, d050705, d050706, d050707, d050708, 
        d050805, d050806, d050807, d050808, 
        d050905, d050906, d050907, d050908, 
        d060101, d060102, d060103, d060104, d060105, d060106, d060107, 
        d060201, d060203, d060204, d060205, d060206, d060207, 
        d060301, d060302, d060304, d060305, d060306, d060307, 
        d060401, d060402, d060403, d060405, d060406, d060407, 
        d060501, d060502, d060503, d060504, d060505, d060506, d060507, d060508, 
        d060601, d060602, d060603, d060604, d060605, d060606, d060607, d060608, d060609, 
        d060704, d060705, d060706, d060707, d060708, d060709, 
        d060805, d060806, d060807, d060808, d060809, 
        d060905, d060906, d060907, d060908, d060909, 
        d070104, d070105, d070106, d070107, 
        d070204, d070205, d070206, d070207, 
        d070304, d070305, d070306, d070307, 
        d070405, d070406, d070407, 
        d070504, d070505, d070506, d070507, d070508, 
        d070604, d070605, d070606, d070607, d070608, d070609,
        d070704, d070705, d070706, d070707, d070708, d070709, d070710,
        d070805, d070806, d070807, d070808, d070809, d070810,
        d070905, d070906, d070907, d070908, d070909, d070910,
        d071010, d071007, d071008, d071009,
        d080105, d080106, d080107, 
        d080205, d080206, d080207, 
        d080305, d080306, d080307, 
        d080405, d080406, d080407, 
        d080505, d080506, d080507, d080508, 
        d080605, d080606, d080607, d080608, d080609,
        d080705, d080706, d080707, d080708, d080709, d080710,
        d080805, d080806, d080807, d080808, d080809, d080810, d080811,
        d080910, d080911, d080905, d080906, d080907, d080908, d080909, 
        d081010, d081011, d081007, d081008, d081009, 
        d081110, d081111, d081108, d081109,
        d090105, d090106, d090107, 
        d090205, d090206, d090207, 
        d090305, d090306, d090307, 
        d090405, d090406, d090407, 
        d090505, d090506, d090507, d090508, 
        d090605, d090606, d090607, d090608, d090609, 
        d090705, d090706, d090707, d090708, d090709, d090710,
        d090810, d090805, d090806, d090807, d090808, d090809, d090811,
        d090910, d090911, d090912, d090905, d090906, d090907, d090908, d090909, 
        d091010, d091011, d091012, d091007, d091008, d091009, 
        d091110, d091111, d091112, d091108, d091109, 
        d091210, d091211, d091212, d091209,
        di_im3_im3, di_im3_im2, di_im3_im1, di_im3_i_0, 
        di_im2_im3, di_im2_im2, di_im2_im1, di_im2_i_0, di_im2_ip1,
        di_im1_im3, di_im1_im2, di_im1_im1, di_im1_i_0, di_im1_ip1, di_im1_ip2,
        di_i_0_im3, di_i_0_im2, di_i_0_im1, di_i_0_i_0, di_i_0_ip1, di_i_0_ip2, di_i_0_ip3
        )
    end
end

lower_bandwidth(cache::Mattsson2012Cache6) = 8
upper_bandwidth(cache::Mattsson2012Cache6) = 8
Base.checkbounds(::Type{Bool}, u::AbstractVector, ::Mattsson2012Cache6) = length(u) > 12
left_length(::Mattsson2012Cache6) = 9
right_length(::Mattsson2012Cache6) = 9

function convolve_boundary_coefficients!(dest::AbstractVector, cache::Mattsson2012Cache6,
                                         u::AbstractVector, b::AbstractVector, α)
    @unpack d010101, d010102, d010103, d010104, d010105, d010106, d010107, 
    d010201, d010203, d010204, d010205, d010206, d010207, 
    d010301, d010302, d010304, d010305, d010306, d010307, 
    d010401, d010402, d010403, d010405, d010406, d010407, 
    d010501, d010502, d010503, d010504, d010505, d010506, d010507, 
    d010601, d010602, d010603, d010604, d010605, d010606, d010607, 
    d010704, d010705, d010706, d010707, 
    d010805, d010806, d010807, 
    d010905, d010906, d010907, 
    d020101, d020103, d020104, d020105, d020106, d020107, 
    d020201, d020203, d020204, d020205, d020206, d020207, 
    d020301, d020304, d020305, d020306, d020307, 
    d020401, d020403, d020405, d020406, d020407, 
    d020501, d020503, d020504, d020505, d020506, d020507, 
    d020601, d020603, d020604, d020605, d020606, d020607, 
    d020704, d020705, d020706, d020707, 
    d020805, d020806, d020807, 
    d020905, d020906, d020907, 
    d030101, d030102, d030104, d030105, d030106, d030107, 
    d030201, d030204, d030205, d030206, d030207, 
    d030301, d030302, d030304, d030305, d030306, d030307, 
    d030401, d030402, d030405, d030406, d030407, 
    d030501, d030502, d030504, d030505, d030506, d030507, 
    d030601, d030602, d030604, d030605, d030606, d030607, 
    d030704, d030705, d030706, d030707, 
    d030805, d030806, d030807, 
    d030905, d030906, d030907, 
    d040101, d040102, d040103, d040105, d040106, d040107, 
    d040201, d040203, d040205, d040206, d040207, 
    d040301, d040302, d040305, d040306, d040307, 
    d040401, d040402, d040403, d040405, d040406, d040407, 
    d040501, d040502, d040503, d040505, d040506, d040507, 
    d040601, d040602, d040603, d040605, d040606, d040607, 
    d040705, d040706, d040707, 
    d040805, d040806, d040807, 
    d040905, d040906, d040907, 
    d050101, d050102, d050103, d050104, d050105, d050106, d050107, 
    d050201, d050203, d050204, d050205, d050206, d050207, 
    d050301, d050302, d050304, d050305, d050306, d050307, 
    d050401, d050402, d050403, d050405, d050406, d050407, 
    d050501, d050502, d050503, d050504, d050505, d050506, d050507, d050508, 
    d050601, d050602, d050603, d050604, d050605, d050606, d050607, d050608, 
    d050704, d050705, d050706, d050707, d050708, 
    d050805, d050806, d050807, d050808, 
    d050905, d050906, d050907, d050908, 
    d060101, d060102, d060103, d060104, d060105, d060106, d060107, 
    d060201, d060203, d060204, d060205, d060206, d060207, 
    d060301, d060302, d060304, d060305, d060306, d060307, 
    d060401, d060402, d060403, d060405, d060406, d060407, 
    d060501, d060502, d060503, d060504, d060505, d060506, d060507, d060508, 
    d060601, d060602, d060603, d060604, d060605, d060606, d060607, d060608, d060609, 
    d060704, d060705, d060706, d060707, d060708, d060709, 
    d060805, d060806, d060807, d060808, d060809, 
    d060905, d060906, d060907, d060908, d060909, 
    d070104, d070105, d070106, d070107, 
    d070204, d070205, d070206, d070207, 
    d070304, d070305, d070306, d070307, 
    d070405, d070406, d070407, 
    d070504, d070505, d070506, d070507, d070508, 
    d070604, d070605, d070606, d070607, d070608, d070609,
    d070704, d070705, d070706, d070707, d070708, d070709, d070710,
    d070805, d070806, d070807, d070808, d070809, d070810,
    d070905, d070906, d070907, d070908, d070909, d070910,
    d071010, d071007, d071008, d071009,
    d080105, d080106, d080107, 
    d080205, d080206, d080207, 
    d080305, d080306, d080307, 
    d080405, d080406, d080407, 
    d080505, d080506, d080507, d080508, 
    d080605, d080606, d080607, d080608, d080609,
    d080705, d080706, d080707, d080708, d080709, d080710,
    d080805, d080806, d080807, d080808, d080809, d080810, d080811,
    d080910, d080911, d080905, d080906, d080907, d080908, d080909, 
    d081010, d081011, d081007, d081008, d081009, 
    d081110, d081111, d081108, d081109,
    d090105, d090106, d090107, 
    d090205, d090206, d090207, 
    d090305, d090306, d090307, 
    d090405, d090406, d090407, 
    d090505, d090506, d090507, d090508, 
    d090605, d090606, d090607, d090608, d090609, 
    d090705, d090706, d090707, d090708, d090709, d090710,
    d090810, d090805, d090806, d090807, d090808, d090809, d090811,
    d090910, d090911, d090912, d090905, d090906, d090907, d090908, d090909, 
    d091010, d091011, d091012, d091007, d091008, d091009, 
    d091110, d091111, d091112, d091108, d091109, 
    d091210, d091211, d091212, d091209 = cache

    @inbounds begin
        #b1 = b[1]
        #b2 = b[2]
        #b3 = b[3]
        #b4 = b[4]
        #b5 = b[5]
        #b6 = b[6]
        #b7 = b[7]
        #b8 = b[8]
        #b9 = b[9]

        dest[  1] = α * (
                          (d010101*b[ 1] + d010102*b[ 2] + d010103*b[ 3] + d010104*b[ 4] + d010105*b[ 5] + d010106*b[ 6] + d010107*b[ 7]) * u[ 1]
                        + (d010201*b[ 1] + d010203*b[ 3] + d010204*b[ 4] + d010205*b[ 5] + d010206*b[ 6] + d010207*b[ 7]) * u[ 2]
                        + (d010301*b[ 1] + d010302*b[ 2] + d010304*b[ 4] + d010305*b[ 5] + d010306*b[ 6] + d010307*b[ 7]) * u[ 3]
                        + (d010401*b[ 1] + d010402*b[ 2] + d010403*b[ 3] + d010405*b[ 5] + d010406*b[ 6] + d010407*b[ 7]) * u[ 4]
                        + (d010501*b[ 1] + d010502*b[ 2] + d010503*b[ 3] + d010504*b[ 4] + d010505*b[ 5] + d010506*b[ 6] + d010507*b[ 7]) * u[ 5]
                        + (d010601*b[ 1] + d010602*b[ 2] + d010603*b[ 3] + d010604*b[ 4] + d010605*b[ 5] + d010606*b[ 6] + d010607*b[ 7]) * u[ 6]
                        + (d010704*b[ 4] + d010705*b[ 5] + d010706*b[ 6] + d010707*b[ 7]) * u[ 7]
                        + (d010805*b[ 5] + d010806*b[ 6] + d010807*b[ 7]) * u[ 8]
                        + (d010905*b[ 5] + d010906*b[ 6] + d010907*b[ 7]) * u[ 9]
        )
        dest[  2] = α * (
                          (d020101*b[ 1] + d020103*b[ 3] + d020104*b[ 4] + d020105*b[ 5] + d020106*b[ 6] + d020107*b[ 7]) * u[ 1]
                        + (d020201*b[ 1] + d020203*b[ 3] + d020204*b[ 4] + d020205*b[ 5] + d020206*b[ 6] + d020207*b[ 7]) * u[ 2]
                        + (d020301*b[ 1] + d020304*b[ 4] + d020305*b[ 5] + d020306*b[ 6] + d020307*b[ 7]) * u[ 3]
                        + (d020401*b[ 1] + d020403*b[ 3] + d020405*b[ 5] + d020406*b[ 6] + d020407*b[ 7]) * u[ 4]
                        + (d020501*b[ 1] + d020503*b[ 3] + d020504*b[ 4] + d020505*b[ 5] + d020506*b[ 6] + d020507*b[ 7]) * u[ 5]
                        + (d020601*b[ 1] + d020603*b[ 3] + d020604*b[ 4] + d020605*b[ 5] + d020606*b[ 6] + d020607*b[ 7]) * u[ 6]
                        + (d020704*b[ 4] + d020705*b[ 5] + d020706*b[ 6] + d020707*b[ 7]) * u[ 7]
                        + (d020805*b[ 5] + d020806*b[ 6] + d020807*b[ 7]) * u[ 8]
                        + (d020905*b[ 5] + d020906*b[ 6] + d020907*b[ 7]) * u[ 9]
        )
        dest[  3] = α * (
                          (d030101*b[ 1] + d030102*b[ 2] + d030104*b[ 4] + d030105*b[ 5] + d030106*b[ 6] + d030107*b[ 7]) * u[ 1]
                        + (d030201*b[ 1] + d030204*b[ 4] + d030205*b[ 5] + d030206*b[ 6] + d030207*b[ 7]) * u[ 2]
                        + (d030301*b[ 1] + d030302*b[ 2] + d030304*b[ 4] + d030305*b[ 5] + d030306*b[ 6] + d030307*b[ 7]) * u[ 3]
                        + (d030401*b[ 1] + d030402*b[ 2] + d030405*b[ 5] + d030406*b[ 6] + d030407*b[ 7]) * u[ 4]
                        + (d030501*b[ 1] + d030502*b[ 2] + d030504*b[ 4] + d030505*b[ 5] + d030506*b[ 6] + d030507*b[ 7]) * u[ 5]
                        + (d030601*b[ 1] + d030602*b[ 2] + d030604*b[ 4] + d030605*b[ 5] + d030606*b[ 6] + d030607*b[ 7]) * u[ 6]
                        + (d030704*b[ 4] + d030705*b[ 5] + d030706*b[ 6] + d030707*b[ 7]) * u[ 7]
                        + (d030805*b[ 5] + d030806*b[ 6] + d030807*b[ 7]) * u[ 8]
                        + (d030905*b[ 5] + d030906*b[ 6] + d030907*b[ 7]) * u[ 9]
        )
        dest[  4] = α * (
                          (d040101*b[ 1] + d040102*b[ 2] + d040103*b[ 3] + d040105*b[ 5] + d040106*b[ 6] + d040107*b[ 7]) * u[ 1]
                        + (d040201*b[ 1] + d040203*b[ 3] + d040205*b[ 5] + d040206*b[ 6] + d040207*b[ 7]) * u[ 2]
                        + (d040301*b[ 1] + d040302*b[ 2] + d040305*b[ 5] + d040306*b[ 6] + d040307*b[ 7]) * u[ 3]
                        + (d040401*b[ 1] + d040402*b[ 2] + d040403*b[ 3] + d040405*b[ 5] + d040406*b[ 6] + d040407*b[ 7]) * u[ 4]
                        + (d040501*b[ 1] + d040502*b[ 2] + d040503*b[ 3] + d040505*b[ 5] + d040506*b[ 6] + d040507*b[ 7]) * u[ 5]
                        + (d040601*b[ 1] + d040602*b[ 2] + d040603*b[ 3] + d040605*b[ 5] + d040606*b[ 6] + d040607*b[ 7]) * u[ 6]
                        + (d040705*b[ 5] + d040706*b[ 6] + d040707*b[ 7]) * u[ 7]
                        + (d040805*b[ 5] + d040806*b[ 6] + d040807*b[ 7]) * u[ 8]
                        + (d040905*b[ 5] + d040906*b[ 6] + d040907*b[ 7]) * u[ 9]
        )
        dest[  5] = α * (
                          (d050101*b[ 1] + d050102*b[ 2] + d050103*b[ 3] + d050104*b[ 4] + d050105*b[ 5] + d050106*b[ 6] + d050107*b[ 7]) * u[ 1]
                        + (d050201*b[ 1] + d050203*b[ 3] + d050204*b[ 4] + d050205*b[ 5] + d050206*b[ 6] + d050207*b[ 7]) * u[ 2]
                        + (d050301*b[ 1] + d050302*b[ 2] + d050304*b[ 4] + d050305*b[ 5] + d050306*b[ 6] + d050307*b[ 7]) * u[ 3]
                        + (d050401*b[ 1] + d050402*b[ 2] + d050403*b[ 3] + d050405*b[ 5] + d050406*b[ 6] + d050407*b[ 7]) * u[ 4]
                        + (d050501*b[ 1] + d050502*b[ 2] + d050503*b[ 3] + d050504*b[ 4] + d050505*b[ 5] + d050506*b[ 6] + d050507*b[ 7] + d050508*b[ 8]) * u[ 5]
                        + (d050601*b[ 1] + d050602*b[ 2] + d050603*b[ 3] + d050604*b[ 4] + d050605*b[ 5] + d050606*b[ 6] + d050607*b[ 7] + d050608*b[ 8]) * u[ 6]
                        + (d050704*b[ 4] + d050705*b[ 5] + d050706*b[ 6] + d050707*b[ 7] + d050708*b[ 8]) * u[ 7]
                        + (d050805*b[ 5] + d050806*b[ 6] + d050807*b[ 7] + d050808*b[ 8]) * u[ 8]
                        + (d050905*b[ 5] + d050906*b[ 6] + d050907*b[ 7] + d050908*b[ 8]) * u[ 9]
        )
        dest[  6] = α * (
                          (d060101*b[ 1] + d060102*b[ 2] + d060103*b[ 3] + d060104*b[ 4] + d060105*b[ 5] + d060106*b[ 6] + d060107*b[ 7]) * u[ 1]
                        + (d060201*b[ 1] + d060203*b[ 3] + d060204*b[ 4] + d060205*b[ 5] + d060206*b[ 6] + d060207*b[ 7]) * u[ 2]
                        + (d060301*b[ 1] + d060302*b[ 2] + d060304*b[ 4] + d060305*b[ 5] + d060306*b[ 6] + d060307*b[ 7]) * u[ 3]
                        + (d060401*b[ 1] + d060402*b[ 2] + d060403*b[ 3] + d060405*b[ 5] + d060406*b[ 6] + d060407*b[ 7]) * u[ 4]
                        + (d060501*b[ 1] + d060502*b[ 2] + d060503*b[ 3] + d060504*b[ 4] + d060505*b[ 5] + d060506*b[ 6] + d060507*b[ 7] + d060508*b[ 8]) * u[ 5]
                        + (d060601*b[ 1] + d060602*b[ 2] + d060603*b[ 3] + d060604*b[ 4] + d060605*b[ 5] + d060606*b[ 6] + d060607*b[ 7] + d060608*b[ 8] + d060609*b[ 9]) * u[ 6]
                        + (d060704*b[ 4] + d060705*b[ 5] + d060706*b[ 6] + d060707*b[ 7] + d060708*b[ 8] + d060709*b[ 9]) * u[ 7]
                        + (d060805*b[ 5] + d060806*b[ 6] + d060807*b[ 7] + d060808*b[ 8] + d060809*b[ 9]) * u[ 8]
                        + (d060905*b[ 5] + d060906*b[ 6] + d060907*b[ 7] + d060908*b[ 8] + d060909*b[ 9]) * u[ 9]
        )
        dest[  7] = α * (
                          (d070104*b[ 4] + d070105*b[ 5] + d070106*b[ 6] + d070107*b[ 7]) * u[ 1]
                        + (d070204*b[ 4] + d070205*b[ 5] + d070206*b[ 6] + d070207*b[ 7]) * u[ 2]
                        + (d070304*b[ 4] + d070305*b[ 5] + d070306*b[ 6] + d070307*b[ 7]) * u[ 3]
                        + (d070405*b[ 5] + d070406*b[ 6] + d070407*b[ 7]) * u[ 4]
                        + (d070504*b[ 4] + d070505*b[ 5] + d070506*b[ 6] + d070507*b[ 7] + d070508*b[ 8]) * u[ 5]
                        + (d070604*b[ 4] + d070605*b[ 5] + d070606*b[ 6] + d070607*b[ 7] + d070608*b[ 8] + d070609*b[ 9]) * u[ 6]
                        + (d070704*b[ 4] + d070705*b[ 5] + d070706*b[ 6] + d070707*b[ 7] + d070708*b[ 8] + d070709*b[ 9] + d070710*b[10]) * u[ 7]
                        + (d070805*b[ 5] + d070806*b[ 6] + d070807*b[ 7] + d070808*b[ 8] + d070809*b[ 9] + d070810*b[10]) * u[ 8]
                        + (d070905*b[ 5] + d070906*b[ 6] + d070907*b[ 7] + d070908*b[ 8] + d070909*b[ 9] + d070910*b[10]) * u[ 9]
                        + (d071010*b[10] + d071007*b[ 7] + d071008*b[ 8] + d071009*b[ 9]) * u[10]
        )
        dest[  8] = α * (
                          (d080105*b[ 5] + d080106*b[ 6] + d080107*b[ 7]) * u[ 1]
                        + (d080205*b[ 5] + d080206*b[ 6] + d080207*b[ 7]) * u[ 2]
                        + (d080305*b[ 5] + d080306*b[ 6] + d080307*b[ 7]) * u[ 3]
                        + (d080405*b[ 5] + d080406*b[ 6] + d080407*b[ 7]) * u[ 4]
                        + (d080505*b[ 5] + d080506*b[ 6] + d080507*b[ 7] + d080508*b[ 8]) * u[ 5]
                        + (d080605*b[ 5] + d080606*b[ 6] + d080607*b[ 7] + d080608*b[ 8] + d080609*b[ 9]) * u[ 6]
                        + (d080705*b[ 5] + d080706*b[ 6] + d080707*b[ 7] + d080708*b[ 8] + d080709*b[ 9] + d080710*b[10]) * u[ 7]
                        + (d080805*b[ 5] + d080806*b[ 6] + d080807*b[ 7] + d080808*b[ 8] + d080809*b[ 9] + d080810*b[10] + d080811*b[11]) * u[ 8]
                        + (d080910*b[10] + d080911*b[11] + d080905*b[ 5] + d080906*b[ 6] + d080907*b[ 7] + d080908*b[ 8] + d080909*b[ 9]) * u[ 9]
                        + (d081010*b[10] + d081011*b[11] + d081007*b[ 7] + d081008*b[ 8] + d081009*b[ 9]) * u[10]
                        + (d081110*b[10] + d081111*b[11] + d081108*b[ 8] + d081109*b[ 9]) * u[11]
        )
        dest[  9] = α * (
                          (d090105*b[ 5] + d090106*b[ 6] + d090107*b[ 7]) * u[ 1]
                        + (d090205*b[ 5] + d090206*b[ 6] + d090207*b[ 7]) * u[ 2]
                        + (d090305*b[ 5] + d090306*b[ 6] + d090307*b[ 7]) * u[ 3]
                        + (d090405*b[ 5] + d090406*b[ 6] + d090407*b[ 7]) * u[ 4]
                        + (d090505*b[ 5] + d090506*b[ 6] + d090507*b[ 7] + d090508*b[ 8]) * u[ 5]
                        + (d090605*b[ 5] + d090606*b[ 6] + d090607*b[ 7] + d090608*b[ 8] + d090609*b[ 9]) * u[ 6]
                        + (d090705*b[ 5] + d090706*b[ 6] + d090707*b[ 7] + d090708*b[ 8] + d090709*b[ 9] + d090710*b[10]) * u[ 7]
                        + (d090810*b[10] + d090811*b[11] + d090805*b[ 5] + d090806*b[ 6] + d090807*b[ 7] + d090808*b[ 8] + d090809*b[ 9]) * u[ 8]
                        + (d090910*b[10] + d090911*b[11] + d090912*b[12] + d090905*b[ 5] + d090906*b[ 6] + d090907*b[ 7] + d090908*b[ 8] + d090909*b[ 9]) * u[ 9]
                        + (d091010*b[10] + d091011*b[11] + d091012*b[12] + d091007*b[ 7] + d091008*b[ 8] + d091009*b[ 9]) * u[10]
                        + (d091110*b[10] + d091111*b[11] + d091112*b[12] + d091108*b[ 8] + d091109*b[ 9]) * u[11]
                        + (d091210*b[10] + d091211*b[11] + d091212*b[12] + d091209*b[ 9]) * u[12]
        )

        
        dest[end] = α * (
                          (d010101*b[end] + d010102*b[end-1] + d010103*b[end-2] + d010104*b[end-3] + d010105*b[end-4] + d010106*b[end-5] + d010107*b[end-6]) * u[end]
                        + (d010201*b[end] + d010203*b[end-2] + d010204*b[end-3] + d010205*b[end-4] + d010206*b[end-5] + d010207*b[end-6]) * u[end-1]
                        + (d010301*b[end] + d010302*b[end-1] + d010304*b[end-3] + d010305*b[end-4] + d010306*b[end-5] + d010307*b[end-6]) * u[end-2]
                        + (d010401*b[end] + d010402*b[end-1] + d010403*b[end-2] + d010405*b[end-4] + d010406*b[end-5] + d010407*b[end-6]) * u[end-3]
                        + (d010501*b[end] + d010502*b[end-1] + d010503*b[end-2] + d010504*b[end-3] + d010505*b[end-4] + d010506*b[end-5] + d010507*b[end-6]) * u[end-4]
                        + (d010601*b[end] + d010602*b[end-1] + d010603*b[end-2] + d010604*b[end-3] + d010605*b[end-4] + d010606*b[end-5] + d010607*b[end-6]) * u[end-5]
                        + (d010704*b[end-3] + d010705*b[end-4] + d010706*b[end-5] + d010707*b[end-6]) * u[end-6]
                        + (d010805*b[end-4] + d010806*b[end-5] + d010807*b[end-6]) * u[end-7]
                        + (d010905*b[end-4] + d010906*b[end-5] + d010907*b[end-6]) * u[end-8]
        )
        dest[end-1] = α * (
                          (d020101*b[end] + d020103*b[end-2] + d020104*b[end-3] + d020105*b[end-4] + d020106*b[end-5] + d020107*b[end-6]) * u[end]
                        + (d020201*b[end] + d020203*b[end-2] + d020204*b[end-3] + d020205*b[end-4] + d020206*b[end-5] + d020207*b[end-6]) * u[end-1]
                        + (d020301*b[end] + d020304*b[end-3] + d020305*b[end-4] + d020306*b[end-5] + d020307*b[end-6]) * u[end-2]
                        + (d020401*b[end] + d020403*b[end-2] + d020405*b[end-4] + d020406*b[end-5] + d020407*b[end-6]) * u[end-3]
                        + (d020501*b[end] + d020503*b[end-2] + d020504*b[end-3] + d020505*b[end-4] + d020506*b[end-5] + d020507*b[end-6]) * u[end-4]
                        + (d020601*b[end] + d020603*b[end-2] + d020604*b[end-3] + d020605*b[end-4] + d020606*b[end-5] + d020607*b[end-6]) * u[end-5]
                        + (d020704*b[end-3] + d020705*b[end-4] + d020706*b[end-5] + d020707*b[end-6]) * u[end-6]
                        + (d020805*b[end-4] + d020806*b[end-5] + d020807*b[end-6]) * u[end-7]
                        + (d020905*b[end-4] + d020906*b[end-5] + d020907*b[end-6]) * u[end-8]
        )
        dest[end-2] = α * (
                          (d030101*b[end] + d030102*b[end-1] + d030104*b[end-3] + d030105*b[end-4] + d030106*b[end-5] + d030107*b[end-6]) * u[end]
                        + (d030201*b[end] + d030204*b[end-3] + d030205*b[end-4] + d030206*b[end-5] + d030207*b[end-6]) * u[end-1]
                        + (d030301*b[end] + d030302*b[end-1] + d030304*b[end-3] + d030305*b[end-4] + d030306*b[end-5] + d030307*b[end-6]) * u[end-2]
                        + (d030401*b[end] + d030402*b[end-1] + d030405*b[end-4] + d030406*b[end-5] + d030407*b[end-6]) * u[end-3]
                        + (d030501*b[end] + d030502*b[end-1] + d030504*b[end-3] + d030505*b[end-4] + d030506*b[end-5] + d030507*b[end-6]) * u[end-4]
                        + (d030601*b[end] + d030602*b[end-1] + d030604*b[end-3] + d030605*b[end-4] + d030606*b[end-5] + d030607*b[end-6]) * u[end-5]
                        + (d030704*b[end-3] + d030705*b[end-4] + d030706*b[end-5] + d030707*b[end-6]) * u[end-6]
                        + (d030805*b[end-4] + d030806*b[end-5] + d030807*b[end-6]) * u[end-7]
                        + (d030905*b[end-4] + d030906*b[end-5] + d030907*b[end-6]) * u[end-8]
        )
        dest[end-3] = α * (
                          (d040101*b[end] + d040102*b[end-1] + d040103*b[end-2] + d040105*b[end-4] + d040106*b[end-5] + d040107*b[end-6]) * u[end]
                        + (d040201*b[end] + d040203*b[end-2] + d040205*b[end-4] + d040206*b[end-5] + d040207*b[end-6]) * u[end-1]
                        + (d040301*b[end] + d040302*b[end-1] + d040305*b[end-4] + d040306*b[end-5] + d040307*b[end-6]) * u[end-2]
                        + (d040401*b[end] + d040402*b[end-1] + d040403*b[end-2] + d040405*b[end-4] + d040406*b[end-5] + d040407*b[end-6]) * u[end-3]
                        + (d040501*b[end] + d040502*b[end-1] + d040503*b[end-2] + d040505*b[end-4] + d040506*b[end-5] + d040507*b[end-6]) * u[end-4]
                        + (d040601*b[end] + d040602*b[end-1] + d040603*b[end-2] + d040605*b[end-4] + d040606*b[end-5] + d040607*b[end-6]) * u[end-5]
                        + (d040705*b[end-4] + d040706*b[end-5] + d040707*b[end-6]) * u[end-6]
                        + (d040805*b[end-4] + d040806*b[end-5] + d040807*b[end-6]) * u[end-7]
                        + (d040905*b[end-4] + d040906*b[end-5] + d040907*b[end-6]) * u[end-8]
        )
        dest[end-4] = α * (
                          (d050101*b[end] + d050102*b[end-1] + d050103*b[end-2] + d050104*b[end-3] + d050105*b[end-4] + d050106*b[end-5] + d050107*b[end-6]) * u[end]
                        + (d050201*b[end] + d050203*b[end-2] + d050204*b[end-3] + d050205*b[end-4] + d050206*b[end-5] + d050207*b[end-6]) * u[end-1]
                        + (d050301*b[end] + d050302*b[end-1] + d050304*b[end-3] + d050305*b[end-4] + d050306*b[end-5] + d050307*b[end-6]) * u[end-2]
                        + (d050401*b[end] + d050402*b[end-1] + d050403*b[end-2] + d050405*b[end-4] + d050406*b[end-5] + d050407*b[end-6]) * u[end-3]
                        + (d050501*b[end] + d050502*b[end-1] + d050503*b[end-2] + d050504*b[end-3] + d050505*b[end-4] + d050506*b[end-5] + d050507*b[end-6] + d050508*b[end-7]) * u[end-4]
                        + (d050601*b[end] + d050602*b[end-1] + d050603*b[end-2] + d050604*b[end-3] + d050605*b[end-4] + d050606*b[end-5] + d050607*b[end-6] + d050608*b[end-7]) * u[end-5]
                        + (d050704*b[end-3] + d050705*b[end-4] + d050706*b[end-5] + d050707*b[end-6] + d050708*b[end-7]) * u[end-6]
                        + (d050805*b[end-4] + d050806*b[end-5] + d050807*b[end-6] + d050808*b[end-7]) * u[end-7]
                        + (d050905*b[end-4] + d050906*b[end-5] + d050907*b[end-6] + d050908*b[end-7]) * u[end-8]
        )
        dest[end-5] = α * (
                          (d060101*b[end] + d060102*b[end-1] + d060103*b[end-2] + d060104*b[end-3] + d060105*b[end-4] + d060106*b[end-5] + d060107*b[end-6]) * u[end]
                        + (d060201*b[end] + d060203*b[end-2] + d060204*b[end-3] + d060205*b[end-4] + d060206*b[end-5] + d060207*b[end-6]) * u[end-1]
                        + (d060301*b[end] + d060302*b[end-1] + d060304*b[end-3] + d060305*b[end-4] + d060306*b[end-5] + d060307*b[end-6]) * u[end-2]
                        + (d060401*b[end] + d060402*b[end-1] + d060403*b[end-2] + d060405*b[end-4] + d060406*b[end-5] + d060407*b[end-6]) * u[end-3]
                        + (d060501*b[end] + d060502*b[end-1] + d060503*b[end-2] + d060504*b[end-3] + d060505*b[end-4] + d060506*b[end-5] + d060507*b[end-6] + d060508*b[end-7]) * u[end-4]
                        + (d060601*b[end] + d060602*b[end-1] + d060603*b[end-2] + d060604*b[end-3] + d060605*b[end-4] + d060606*b[end-5] + d060607*b[end-6] + d060608*b[end-7] + d060609*b[end-8]) * u[end-5]
                        + (d060704*b[end-3] + d060705*b[end-4] + d060706*b[end-5] + d060707*b[end-6] + d060708*b[end-7] + d060709*b[end-8]) * u[end-6]
                        + (d060805*b[end-4] + d060806*b[end-5] + d060807*b[end-6] + d060808*b[end-7] + d060809*b[end-8]) * u[end-7]
                        + (d060905*b[end-4] + d060906*b[end-5] + d060907*b[end-6] + d060908*b[end-7] + d060909*b[end-8]) * u[end-8]
        )
        dest[end-6] = α * (
                          (d070104*b[end-3] + d070105*b[end-4] + d070106*b[end-5] + d070107*b[end-6]) * u[end]
                        + (d070204*b[end-3] + d070205*b[end-4] + d070206*b[end-5] + d070207*b[end-6]) * u[end-1]
                        + (d070304*b[end-3] + d070305*b[end-4] + d070306*b[end-5] + d070307*b[end-6]) * u[end-2]
                        + (d070405*b[end-4] + d070406*b[end-5] + d070407*b[end-6]) * u[end-3]
                        + (d070504*b[end-3] + d070505*b[end-4] + d070506*b[end-5] + d070507*b[end-6] + d070508*b[end-7]) * u[end-4]
                        + (d070604*b[end-3] + d070605*b[end-4] + d070606*b[end-5] + d070607*b[end-6] + d070608*b[end-7] + d070609*b[end-8]) * u[end-5]
                        + (d070704*b[end-3] + d070705*b[end-4] + d070706*b[end-5] + d070707*b[end-6] + d070708*b[end-7] + d070709*b[end-8] + d070710*b[end-9]) * u[end-6]
                        + (d070805*b[end-4] + d070806*b[end-5] + d070807*b[end-6] + d070808*b[end-7] + d070809*b[end-8] + d070810*b[end-9]) * u[end-7]
                        + (d070905*b[end-4] + d070906*b[end-5] + d070907*b[end-6] + d070908*b[end-7] + d070909*b[end-8] + d070910*b[end-9]) * u[end-8]
                        + (d071010*b[end-9] + d071007*b[end-6] + d071008*b[end-7] + d071009*b[end-8]) * u[end-9]
        )
        dest[end-7] = α * (
                          (d080105*b[end-4] + d080106*b[end-5] + d080107*b[end-6]) * u[end]
                        + (d080205*b[end-4] + d080206*b[end-5] + d080207*b[end-6]) * u[end-1]
                        + (d080305*b[end-4] + d080306*b[end-5] + d080307*b[end-6]) * u[end-2]
                        + (d080405*b[end-4] + d080406*b[end-5] + d080407*b[end-6]) * u[end-3]
                        + (d080505*b[end-4] + d080506*b[end-5] + d080507*b[end-6] + d080508*b[end-7]) * u[end-4]
                        + (d080605*b[end-4] + d080606*b[end-5] + d080607*b[end-6] + d080608*b[end-7] + d080609*b[end-8]) * u[end-5]
                        + (d080705*b[end-4] + d080706*b[end-5] + d080707*b[end-6] + d080708*b[end-7] + d080709*b[end-8] + d080710*b[end-9]) * u[end-6]
                        + (d080805*b[end-4] + d080806*b[end-5] + d080807*b[end-6] + d080808*b[end-7] + d080809*b[end-8] + d080810*b[end-9] + d080811*b[end-10]) * u[end-7]
                        + (d080910*b[end-9] + d080911*b[end-10] + d080905*b[end-4] + d080906*b[end-5] + d080907*b[end-6] + d080908*b[end-7] + d080909*b[end-8]) * u[end-8]
                        + (d081010*b[end-9] + d081011*b[end-10] + d081007*b[end-6] + d081008*b[end-7] + d081009*b[end-8]) * u[end-9]
                        + (d081110*b[end-9] + d081111*b[end-10] + d081108*b[end-7] + d081109*b[end-8]) * u[end-10]
        )
        dest[end-8] = α * (
                          (d090105*b[end-4] + d090106*b[end-5] + d090107*b[end-6]) * u[end]
                        + (d090205*b[end-4] + d090206*b[end-5] + d090207*b[end-6]) * u[end-1]
                        + (d090305*b[end-4] + d090306*b[end-5] + d090307*b[end-6]) * u[end-2]
                        + (d090405*b[end-4] + d090406*b[end-5] + d090407*b[end-6]) * u[end-3]
                        + (d090505*b[end-4] + d090506*b[end-5] + d090507*b[end-6] + d090508*b[end-7]) * u[end-4]
                        + (d090605*b[end-4] + d090606*b[end-5] + d090607*b[end-6] + d090608*b[end-7] + d090609*b[end-8]) * u[end-5]
                        + (d090705*b[end-4] + d090706*b[end-5] + d090707*b[end-6] + d090708*b[end-7] + d090709*b[end-8] + d090710*b[end-9]) * u[end-6]
                        + (d090810*b[end-9] + d090811*b[end-10] + d090805*b[end-4] + d090806*b[end-5] + d090807*b[end-6] + d090808*b[end-7] + d090809*b[end-8]) * u[end-7]
                        + (d090910*b[end-9] + d090911*b[end-10] + d090912*b[end-11] + d090905*b[end-4] + d090906*b[end-5] + d090907*b[end-6] + d090908*b[end-7] + d090909*b[end-8]) * u[end-8]
                        + (d091010*b[end-9] + d091011*b[end-10] + d091012*b[end-11] + d091007*b[end-6] + d091008*b[end-7] + d091009*b[end-8]) * u[end-9]
                        + (d091110*b[end-9] + d091111*b[end-10] + d091112*b[end-11] + d091108*b[end-7] + d091109*b[end-8]) * u[end-10]
                        + (d091210*b[end-9] + d091211*b[end-10] + d091212*b[end-11] + d091209*b[end-8]) * u[end-11]
        )
    end
end

function convolve_boundary_coefficients!(dest::AbstractVector, cache::Mattsson2012Cache6,
                                         u::AbstractVector, b::AbstractVector, α, β)
    @unpack d010101, d010102, d010103, d010104, d010105, d010106, d010107, 
    d010201, d010203, d010204, d010205, d010206, d010207, 
    d010301, d010302, d010304, d010305, d010306, d010307, 
    d010401, d010402, d010403, d010405, d010406, d010407, 
    d010501, d010502, d010503, d010504, d010505, d010506, d010507, 
    d010601, d010602, d010603, d010604, d010605, d010606, d010607, 
    d010704, d010705, d010706, d010707, 
    d010805, d010806, d010807, 
    d010905, d010906, d010907, 
    d020101, d020103, d020104, d020105, d020106, d020107, 
    d020201, d020203, d020204, d020205, d020206, d020207, 
    d020301, d020304, d020305, d020306, d020307, 
    d020401, d020403, d020405, d020406, d020407, 
    d020501, d020503, d020504, d020505, d020506, d020507, 
    d020601, d020603, d020604, d020605, d020606, d020607, 
    d020704, d020705, d020706, d020707, 
    d020805, d020806, d020807, 
    d020905, d020906, d020907, 
    d030101, d030102, d030104, d030105, d030106, d030107, 
    d030201, d030204, d030205, d030206, d030207, 
    d030301, d030302, d030304, d030305, d030306, d030307, 
    d030401, d030402, d030405, d030406, d030407, 
    d030501, d030502, d030504, d030505, d030506, d030507, 
    d030601, d030602, d030604, d030605, d030606, d030607, 
    d030704, d030705, d030706, d030707, 
    d030805, d030806, d030807, 
    d030905, d030906, d030907, 
    d040101, d040102, d040103, d040105, d040106, d040107, 
    d040201, d040203, d040205, d040206, d040207, 
    d040301, d040302, d040305, d040306, d040307, 
    d040401, d040402, d040403, d040405, d040406, d040407, 
    d040501, d040502, d040503, d040505, d040506, d040507, 
    d040601, d040602, d040603, d040605, d040606, d040607, 
    d040705, d040706, d040707, 
    d040805, d040806, d040807, 
    d040905, d040906, d040907, 
    d050101, d050102, d050103, d050104, d050105, d050106, d050107, 
    d050201, d050203, d050204, d050205, d050206, d050207, 
    d050301, d050302, d050304, d050305, d050306, d050307, 
    d050401, d050402, d050403, d050405, d050406, d050407, 
    d050501, d050502, d050503, d050504, d050505, d050506, d050507, d050508, 
    d050601, d050602, d050603, d050604, d050605, d050606, d050607, d050608, 
    d050704, d050705, d050706, d050707, d050708, 
    d050805, d050806, d050807, d050808, 
    d050905, d050906, d050907, d050908, 
    d060101, d060102, d060103, d060104, d060105, d060106, d060107, 
    d060201, d060203, d060204, d060205, d060206, d060207, 
    d060301, d060302, d060304, d060305, d060306, d060307, 
    d060401, d060402, d060403, d060405, d060406, d060407, 
    d060501, d060502, d060503, d060504, d060505, d060506, d060507, d060508, 
    d060601, d060602, d060603, d060604, d060605, d060606, d060607, d060608, d060609, 
    d060704, d060705, d060706, d060707, d060708, d060709, 
    d060805, d060806, d060807, d060808, d060809, 
    d060905, d060906, d060907, d060908, d060909, 
    d070104, d070105, d070106, d070107, 
    d070204, d070205, d070206, d070207, 
    d070304, d070305, d070306, d070307, 
    d070405, d070406, d070407, 
    d070504, d070505, d070506, d070507, d070508, 
    d070604, d070605, d070606, d070607, d070608, d070609,
    d070704, d070705, d070706, d070707, d070708, d070709, d070710,
    d070805, d070806, d070807, d070808, d070809, d070810,
    d070905, d070906, d070907, d070908, d070909, d070910,
    d071010, d071007, d071008, d071009,
    d080105, d080106, d080107, 
    d080205, d080206, d080207, 
    d080305, d080306, d080307, 
    d080405, d080406, d080407, 
    d080505, d080506, d080507, d080508, 
    d080605, d080606, d080607, d080608, d080609,
    d080705, d080706, d080707, d080708, d080709, d080710,
    d080805, d080806, d080807, d080808, d080809, d080810, d080811,
    d080910, d080911, d080905, d080906, d080907, d080908, d080909, 
    d081010, d081011, d081007, d081008, d081009, 
    d081110, d081111, d081108, d081109,
    d090105, d090106, d090107, 
    d090205, d090206, d090207, 
    d090305, d090306, d090307, 
    d090405, d090406, d090407, 
    d090505, d090506, d090507, d090508, 
    d090605, d090606, d090607, d090608, d090609, 
    d090705, d090706, d090707, d090708, d090709, d090710,
    d090810, d090805, d090806, d090807, d090808, d090809, d090811,
    d090910, d090911, d090912, d090905, d090906, d090907, d090908, d090909, 
    d091010, d091011, d091012, d091007, d091008, d091009, 
    d091110, d091111, d091112, d091108, d091109, 
    d091210, d091211, d091212, d091209 = cache

    @inbounds begin
        #b1 = b[1]
        #b2 = b[2]
        #b3 = b[3]
        #b4 = b[4]
        #b5 = b[5]
        #b6 = b[6]
        #b7 = b[7]
        #b8 = b[8]
        #b9 = b[9]

        dest[  1] = α * (
                          (d010101*b[ 1] + d010102*b[ 2] + d010103*b[ 3] + d010104*b[ 4] + d010105*b[ 5] + d010106*b[ 6] + d010107*b[ 7]) * u[ 1]
                        + (d010201*b[ 1] + d010203*b[ 3] + d010204*b[ 4] + d010205*b[ 5] + d010206*b[ 6] + d010207*b[ 7]) * u[ 2]
                        + (d010301*b[ 1] + d010302*b[ 2] + d010304*b[ 4] + d010305*b[ 5] + d010306*b[ 6] + d010307*b[ 7]) * u[ 3]
                        + (d010401*b[ 1] + d010402*b[ 2] + d010403*b[ 3] + d010405*b[ 5] + d010406*b[ 6] + d010407*b[ 7]) * u[ 4]
                        + (d010501*b[ 1] + d010502*b[ 2] + d010503*b[ 3] + d010504*b[ 4] + d010505*b[ 5] + d010506*b[ 6] + d010507*b[ 7]) * u[ 5]
                        + (d010601*b[ 1] + d010602*b[ 2] + d010603*b[ 3] + d010604*b[ 4] + d010605*b[ 5] + d010606*b[ 6] + d010607*b[ 7]) * u[ 6]
                        + (d010704*b[ 4] + d010705*b[ 5] + d010706*b[ 6] + d010707*b[ 7]) * u[ 7]
                        + (d010805*b[ 5] + d010806*b[ 6] + d010807*b[ 7]) * u[ 8]
                        + (d010905*b[ 5] + d010906*b[ 6] + d010907*b[ 7]) * u[ 9]
        ) + β*dest[1]
        dest[  2] = α * (
                          (d020101*b[ 1] + d020103*b[ 3] + d020104*b[ 4] + d020105*b[ 5] + d020106*b[ 6] + d020107*b[ 7]) * u[ 1]
                        + (d020201*b[ 1] + d020203*b[ 3] + d020204*b[ 4] + d020205*b[ 5] + d020206*b[ 6] + d020207*b[ 7]) * u[ 2]
                        + (d020301*b[ 1] + d020304*b[ 4] + d020305*b[ 5] + d020306*b[ 6] + d020307*b[ 7]) * u[ 3]
                        + (d020401*b[ 1] + d020403*b[ 3] + d020405*b[ 5] + d020406*b[ 6] + d020407*b[ 7]) * u[ 4]
                        + (d020501*b[ 1] + d020503*b[ 3] + d020504*b[ 4] + d020505*b[ 5] + d020506*b[ 6] + d020507*b[ 7]) * u[ 5]
                        + (d020601*b[ 1] + d020603*b[ 3] + d020604*b[ 4] + d020605*b[ 5] + d020606*b[ 6] + d020607*b[ 7]) * u[ 6]
                        + (d020704*b[ 4] + d020705*b[ 5] + d020706*b[ 6] + d020707*b[ 7]) * u[ 7]
                        + (d020805*b[ 5] + d020806*b[ 6] + d020807*b[ 7]) * u[ 8]
                        + (d020905*b[ 5] + d020906*b[ 6] + d020907*b[ 7]) * u[ 9]
        ) + β*dest[2]
        dest[  3] = α * (
                          (d030101*b[ 1] + d030102*b[ 2] + d030104*b[ 4] + d030105*b[ 5] + d030106*b[ 6] + d030107*b[ 7]) * u[ 1]
                        + (d030201*b[ 1] + d030204*b[ 4] + d030205*b[ 5] + d030206*b[ 6] + d030207*b[ 7]) * u[ 2]
                        + (d030301*b[ 1] + d030302*b[ 2] + d030304*b[ 4] + d030305*b[ 5] + d030306*b[ 6] + d030307*b[ 7]) * u[ 3]
                        + (d030401*b[ 1] + d030402*b[ 2] + d030405*b[ 5] + d030406*b[ 6] + d030407*b[ 7]) * u[ 4]
                        + (d030501*b[ 1] + d030502*b[ 2] + d030504*b[ 4] + d030505*b[ 5] + d030506*b[ 6] + d030507*b[ 7]) * u[ 5]
                        + (d030601*b[ 1] + d030602*b[ 2] + d030604*b[ 4] + d030605*b[ 5] + d030606*b[ 6] + d030607*b[ 7]) * u[ 6]
                        + (d030704*b[ 4] + d030705*b[ 5] + d030706*b[ 6] + d030707*b[ 7]) * u[ 7]
                        + (d030805*b[ 5] + d030806*b[ 6] + d030807*b[ 7]) * u[ 8]
                        + (d030905*b[ 5] + d030906*b[ 6] + d030907*b[ 7]) * u[ 9]
        ) + β*dest[3]
        dest[  4] = α * (
                          (d040101*b[ 1] + d040102*b[ 2] + d040103*b[ 3] + d040105*b[ 5] + d040106*b[ 6] + d040107*b[ 7]) * u[ 1]
                        + (d040201*b[ 1] + d040203*b[ 3] + d040205*b[ 5] + d040206*b[ 6] + d040207*b[ 7]) * u[ 2]
                        + (d040301*b[ 1] + d040302*b[ 2] + d040305*b[ 5] + d040306*b[ 6] + d040307*b[ 7]) * u[ 3]
                        + (d040401*b[ 1] + d040402*b[ 2] + d040403*b[ 3] + d040405*b[ 5] + d040406*b[ 6] + d040407*b[ 7]) * u[ 4]
                        + (d040501*b[ 1] + d040502*b[ 2] + d040503*b[ 3] + d040505*b[ 5] + d040506*b[ 6] + d040507*b[ 7]) * u[ 5]
                        + (d040601*b[ 1] + d040602*b[ 2] + d040603*b[ 3] + d040605*b[ 5] + d040606*b[ 6] + d040607*b[ 7]) * u[ 6]
                        + (d040705*b[ 5] + d040706*b[ 6] + d040707*b[ 7]) * u[ 7]
                        + (d040805*b[ 5] + d040806*b[ 6] + d040807*b[ 7]) * u[ 8]
                        + (d040905*b[ 5] + d040906*b[ 6] + d040907*b[ 7]) * u[ 9]
        ) + β*dest[4]
        dest[  5] = α * (
                          (d050101*b[ 1] + d050102*b[ 2] + d050103*b[ 3] + d050104*b[ 4] + d050105*b[ 5] + d050106*b[ 6] + d050107*b[ 7]) * u[ 1]
                        + (d050201*b[ 1] + d050203*b[ 3] + d050204*b[ 4] + d050205*b[ 5] + d050206*b[ 6] + d050207*b[ 7]) * u[ 2]
                        + (d050301*b[ 1] + d050302*b[ 2] + d050304*b[ 4] + d050305*b[ 5] + d050306*b[ 6] + d050307*b[ 7]) * u[ 3]
                        + (d050401*b[ 1] + d050402*b[ 2] + d050403*b[ 3] + d050405*b[ 5] + d050406*b[ 6] + d050407*b[ 7]) * u[ 4]
                        + (d050501*b[ 1] + d050502*b[ 2] + d050503*b[ 3] + d050504*b[ 4] + d050505*b[ 5] + d050506*b[ 6] + d050507*b[ 7] + d050508*b[ 8]) * u[ 5]
                        + (d050601*b[ 1] + d050602*b[ 2] + d050603*b[ 3] + d050604*b[ 4] + d050605*b[ 5] + d050606*b[ 6] + d050607*b[ 7] + d050608*b[ 8]) * u[ 6]
                        + (d050704*b[ 4] + d050705*b[ 5] + d050706*b[ 6] + d050707*b[ 7] + d050708*b[ 8]) * u[ 7]
                        + (d050805*b[ 5] + d050806*b[ 6] + d050807*b[ 7] + d050808*b[ 8]) * u[ 8]
                        + (d050905*b[ 5] + d050906*b[ 6] + d050907*b[ 7] + d050908*b[ 8]) * u[ 9]
        ) + β*dest[5]
        dest[  6] = α * (
                          (d060101*b[ 1] + d060102*b[ 2] + d060103*b[ 3] + d060104*b[ 4] + d060105*b[ 5] + d060106*b[ 6] + d060107*b[ 7]) * u[ 1]
                        + (d060201*b[ 1] + d060203*b[ 3] + d060204*b[ 4] + d060205*b[ 5] + d060206*b[ 6] + d060207*b[ 7]) * u[ 2]
                        + (d060301*b[ 1] + d060302*b[ 2] + d060304*b[ 4] + d060305*b[ 5] + d060306*b[ 6] + d060307*b[ 7]) * u[ 3]
                        + (d060401*b[ 1] + d060402*b[ 2] + d060403*b[ 3] + d060405*b[ 5] + d060406*b[ 6] + d060407*b[ 7]) * u[ 4]
                        + (d060501*b[ 1] + d060502*b[ 2] + d060503*b[ 3] + d060504*b[ 4] + d060505*b[ 5] + d060506*b[ 6] + d060507*b[ 7] + d060508*b[ 8]) * u[ 5]
                        + (d060601*b[ 1] + d060602*b[ 2] + d060603*b[ 3] + d060604*b[ 4] + d060605*b[ 5] + d060606*b[ 6] + d060607*b[ 7] + d060608*b[ 8] + d060609*b[ 9]) * u[ 6]
                        + (d060704*b[ 4] + d060705*b[ 5] + d060706*b[ 6] + d060707*b[ 7] + d060708*b[ 8] + d060709*b[ 9]) * u[ 7]
                        + (d060805*b[ 5] + d060806*b[ 6] + d060807*b[ 7] + d060808*b[ 8] + d060809*b[ 9]) * u[ 8]
                        + (d060905*b[ 5] + d060906*b[ 6] + d060907*b[ 7] + d060908*b[ 8] + d060909*b[ 9]) * u[ 9]
        ) + β*dest[6]
        dest[  7] = α * (
                          (d070104*b[ 4] + d070105*b[ 5] + d070106*b[ 6] + d070107*b[ 7]) * u[ 1]
                        + (d070204*b[ 4] + d070205*b[ 5] + d070206*b[ 6] + d070207*b[ 7]) * u[ 2]
                        + (d070304*b[ 4] + d070305*b[ 5] + d070306*b[ 6] + d070307*b[ 7]) * u[ 3]
                        + (d070405*b[ 5] + d070406*b[ 6] + d070407*b[ 7]) * u[ 4]
                        + (d070504*b[ 4] + d070505*b[ 5] + d070506*b[ 6] + d070507*b[ 7] + d070508*b[ 8]) * u[ 5]
                        + (d070604*b[ 4] + d070605*b[ 5] + d070606*b[ 6] + d070607*b[ 7] + d070608*b[ 8] + d070609*b[ 9]) * u[ 6]
                        + (d070704*b[ 4] + d070705*b[ 5] + d070706*b[ 6] + d070707*b[ 7] + d070708*b[ 8] + d070709*b[ 9] + d070710*b[10]) * u[ 7]
                        + (d070805*b[ 5] + d070806*b[ 6] + d070807*b[ 7] + d070808*b[ 8] + d070809*b[ 9] + d070810*b[10]) * u[ 8]
                        + (d070905*b[ 5] + d070906*b[ 6] + d070907*b[ 7] + d070908*b[ 8] + d070909*b[ 9] + d070910*b[10]) * u[ 9]
                        + (d071010*b[10] + d071007*b[ 7] + d071008*b[ 8] + d071009*b[ 9]) * u[10]
        ) + β*dest[7]
        dest[  8] = α * (
                          (d080105*b[ 5] + d080106*b[ 6] + d080107*b[ 7]) * u[ 1]
                        + (d080205*b[ 5] + d080206*b[ 6] + d080207*b[ 7]) * u[ 2]
                        + (d080305*b[ 5] + d080306*b[ 6] + d080307*b[ 7]) * u[ 3]
                        + (d080405*b[ 5] + d080406*b[ 6] + d080407*b[ 7]) * u[ 4]
                        + (d080505*b[ 5] + d080506*b[ 6] + d080507*b[ 7] + d080508*b[ 8]) * u[ 5]
                        + (d080605*b[ 5] + d080606*b[ 6] + d080607*b[ 7] + d080608*b[ 8] + d080609*b[ 9]) * u[ 6]
                        + (d080705*b[ 5] + d080706*b[ 6] + d080707*b[ 7] + d080708*b[ 8] + d080709*b[ 9] + d080710*b[10]) * u[ 7]
                        + (d080805*b[ 5] + d080806*b[ 6] + d080807*b[ 7] + d080808*b[ 8] + d080809*b[ 9] + d080810*b[10] + d080811*b[11]) * u[ 8]
                        + (d080910*b[10] + d080911*b[11] + d080905*b[ 5] + d080906*b[ 6] + d080907*b[ 7] + d080908*b[ 8] + d080909*b[ 9]) * u[ 9]
                        + (d081010*b[10] + d081011*b[11] + d081007*b[ 7] + d081008*b[ 8] + d081009*b[ 9]) * u[10]
                        + (d081110*b[10] + d081111*b[11] + d081108*b[ 8] + d081109*b[ 9]) * u[11]
        ) + β*dest[8]
        dest[  9] = α * (
                          (d090105*b[ 5] + d090106*b[ 6] + d090107*b[ 7]) * u[ 1]
                        + (d090205*b[ 5] + d090206*b[ 6] + d090207*b[ 7]) * u[ 2]
                        + (d090305*b[ 5] + d090306*b[ 6] + d090307*b[ 7]) * u[ 3]
                        + (d090405*b[ 5] + d090406*b[ 6] + d090407*b[ 7]) * u[ 4]
                        + (d090505*b[ 5] + d090506*b[ 6] + d090507*b[ 7] + d090508*b[ 8]) * u[ 5]
                        + (d090605*b[ 5] + d090606*b[ 6] + d090607*b[ 7] + d090608*b[ 8] + d090609*b[ 9]) * u[ 6]
                        + (d090705*b[ 5] + d090706*b[ 6] + d090707*b[ 7] + d090708*b[ 8] + d090709*b[ 9] + d090710*b[10]) * u[ 7]
                        + (d090810*b[10] + d090811*b[11] + d090805*b[ 5] + d090806*b[ 6] + d090807*b[ 7] + d090808*b[ 8] + d090809*b[ 9]) * u[ 8]
                        + (d090910*b[10] + d090911*b[11] + d090912*b[12] + d090905*b[ 5] + d090906*b[ 6] + d090907*b[ 7] + d090908*b[ 8] + d090909*b[ 9]) * u[ 9]
                        + (d091010*b[10] + d091011*b[11] + d091012*b[12] + d091007*b[ 7] + d091008*b[ 8] + d091009*b[ 9]) * u[10]
                        + (d091110*b[10] + d091111*b[11] + d091112*b[12] + d091108*b[ 8] + d091109*b[ 9]) * u[11]
                        + (d091210*b[10] + d091211*b[11] + d091212*b[12] + d091209*b[ 9]) * u[12]
        ) + β*dest[9]

        
        dest[end] = α * (
                          (d010101*b[end] + d010102*b[end-1] + d010103*b[end-2] + d010104*b[end-3] + d010105*b[end-4] + d010106*b[end-5] + d010107*b[end-6]) * u[end]
                        + (d010201*b[end] + d010203*b[end-2] + d010204*b[end-3] + d010205*b[end-4] + d010206*b[end-5] + d010207*b[end-6]) * u[end-1]
                        + (d010301*b[end] + d010302*b[end-1] + d010304*b[end-3] + d010305*b[end-4] + d010306*b[end-5] + d010307*b[end-6]) * u[end-2]
                        + (d010401*b[end] + d010402*b[end-1] + d010403*b[end-2] + d010405*b[end-4] + d010406*b[end-5] + d010407*b[end-6]) * u[end-3]
                        + (d010501*b[end] + d010502*b[end-1] + d010503*b[end-2] + d010504*b[end-3] + d010505*b[end-4] + d010506*b[end-5] + d010507*b[end-6]) * u[end-4]
                        + (d010601*b[end] + d010602*b[end-1] + d010603*b[end-2] + d010604*b[end-3] + d010605*b[end-4] + d010606*b[end-5] + d010607*b[end-6]) * u[end-5]
                        + (d010704*b[end-3] + d010705*b[end-4] + d010706*b[end-5] + d010707*b[end-6]) * u[end-6]
                        + (d010805*b[end-4] + d010806*b[end-5] + d010807*b[end-6]) * u[end-7]
                        + (d010905*b[end-4] + d010906*b[end-5] + d010907*b[end-6]) * u[end-8]
        ) + β*dest[end]
        dest[end-1] = α * (
                          (d020101*b[end] + d020103*b[end-2] + d020104*b[end-3] + d020105*b[end-4] + d020106*b[end-5] + d020107*b[end-6]) * u[end]
                        + (d020201*b[end] + d020203*b[end-2] + d020204*b[end-3] + d020205*b[end-4] + d020206*b[end-5] + d020207*b[end-6]) * u[end-1]
                        + (d020301*b[end] + d020304*b[end-3] + d020305*b[end-4] + d020306*b[end-5] + d020307*b[end-6]) * u[end-2]
                        + (d020401*b[end] + d020403*b[end-2] + d020405*b[end-4] + d020406*b[end-5] + d020407*b[end-6]) * u[end-3]
                        + (d020501*b[end] + d020503*b[end-2] + d020504*b[end-3] + d020505*b[end-4] + d020506*b[end-5] + d020507*b[end-6]) * u[end-4]
                        + (d020601*b[end] + d020603*b[end-2] + d020604*b[end-3] + d020605*b[end-4] + d020606*b[end-5] + d020607*b[end-6]) * u[end-5]
                        + (d020704*b[end-3] + d020705*b[end-4] + d020706*b[end-5] + d020707*b[end-6]) * u[end-6]
                        + (d020805*b[end-4] + d020806*b[end-5] + d020807*b[end-6]) * u[end-7]
                        + (d020905*b[end-4] + d020906*b[end-5] + d020907*b[end-6]) * u[end-8]
        ) + β*dest[end-1]
        dest[end-2] = α * (
                          (d030101*b[end] + d030102*b[end-1] + d030104*b[end-3] + d030105*b[end-4] + d030106*b[end-5] + d030107*b[end-6]) * u[end]
                        + (d030201*b[end] + d030204*b[end-3] + d030205*b[end-4] + d030206*b[end-5] + d030207*b[end-6]) * u[end-1]
                        + (d030301*b[end] + d030302*b[end-1] + d030304*b[end-3] + d030305*b[end-4] + d030306*b[end-5] + d030307*b[end-6]) * u[end-2]
                        + (d030401*b[end] + d030402*b[end-1] + d030405*b[end-4] + d030406*b[end-5] + d030407*b[end-6]) * u[end-3]
                        + (d030501*b[end] + d030502*b[end-1] + d030504*b[end-3] + d030505*b[end-4] + d030506*b[end-5] + d030507*b[end-6]) * u[end-4]
                        + (d030601*b[end] + d030602*b[end-1] + d030604*b[end-3] + d030605*b[end-4] + d030606*b[end-5] + d030607*b[end-6]) * u[end-5]
                        + (d030704*b[end-3] + d030705*b[end-4] + d030706*b[end-5] + d030707*b[end-6]) * u[end-6]
                        + (d030805*b[end-4] + d030806*b[end-5] + d030807*b[end-6]) * u[end-7]
                        + (d030905*b[end-4] + d030906*b[end-5] + d030907*b[end-6]) * u[end-8]
        ) + β*dest[end-2]
        dest[end-3] = α * (
                          (d040101*b[end] + d040102*b[end-1] + d040103*b[end-2] + d040105*b[end-4] + d040106*b[end-5] + d040107*b[end-6]) * u[end]
                        + (d040201*b[end] + d040203*b[end-2] + d040205*b[end-4] + d040206*b[end-5] + d040207*b[end-6]) * u[end-1]
                        + (d040301*b[end] + d040302*b[end-1] + d040305*b[end-4] + d040306*b[end-5] + d040307*b[end-6]) * u[end-2]
                        + (d040401*b[end] + d040402*b[end-1] + d040403*b[end-2] + d040405*b[end-4] + d040406*b[end-5] + d040407*b[end-6]) * u[end-3]
                        + (d040501*b[end] + d040502*b[end-1] + d040503*b[end-2] + d040505*b[end-4] + d040506*b[end-5] + d040507*b[end-6]) * u[end-4]
                        + (d040601*b[end] + d040602*b[end-1] + d040603*b[end-2] + d040605*b[end-4] + d040606*b[end-5] + d040607*b[end-6]) * u[end-5]
                        + (d040705*b[end-4] + d040706*b[end-5] + d040707*b[end-6]) * u[end-6]
                        + (d040805*b[end-4] + d040806*b[end-5] + d040807*b[end-6]) * u[end-7]
                        + (d040905*b[end-4] + d040906*b[end-5] + d040907*b[end-6]) * u[end-8]
        ) + β*dest[end-3]
        dest[end-4] = α * (
                          (d050101*b[end] + d050102*b[end-1] + d050103*b[end-2] + d050104*b[end-3] + d050105*b[end-4] + d050106*b[end-5] + d050107*b[end-6]) * u[end]
                        + (d050201*b[end] + d050203*b[end-2] + d050204*b[end-3] + d050205*b[end-4] + d050206*b[end-5] + d050207*b[end-6]) * u[end-1]
                        + (d050301*b[end] + d050302*b[end-1] + d050304*b[end-3] + d050305*b[end-4] + d050306*b[end-5] + d050307*b[end-6]) * u[end-2]
                        + (d050401*b[end] + d050402*b[end-1] + d050403*b[end-2] + d050405*b[end-4] + d050406*b[end-5] + d050407*b[end-6]) * u[end-3]
                        + (d050501*b[end] + d050502*b[end-1] + d050503*b[end-2] + d050504*b[end-3] + d050505*b[end-4] + d050506*b[end-5] + d050507*b[end-6] + d050508*b[end-7]) * u[end-4]
                        + (d050601*b[end] + d050602*b[end-1] + d050603*b[end-2] + d050604*b[end-3] + d050605*b[end-4] + d050606*b[end-5] + d050607*b[end-6] + d050608*b[end-7]) * u[end-5]
                        + (d050704*b[end-3] + d050705*b[end-4] + d050706*b[end-5] + d050707*b[end-6] + d050708*b[end-7]) * u[end-6]
                        + (d050805*b[end-4] + d050806*b[end-5] + d050807*b[end-6] + d050808*b[end-7]) * u[end-7]
                        + (d050905*b[end-4] + d050906*b[end-5] + d050907*b[end-6] + d050908*b[end-7]) * u[end-8]
        ) + β*dest[end-4]
        dest[end-5] = α * (
                          (d060101*b[end] + d060102*b[end-1] + d060103*b[end-2] + d060104*b[end-3] + d060105*b[end-4] + d060106*b[end-5] + d060107*b[end-6]) * u[end]
                        + (d060201*b[end] + d060203*b[end-2] + d060204*b[end-3] + d060205*b[end-4] + d060206*b[end-5] + d060207*b[end-6]) * u[end-1]
                        + (d060301*b[end] + d060302*b[end-1] + d060304*b[end-3] + d060305*b[end-4] + d060306*b[end-5] + d060307*b[end-6]) * u[end-2]
                        + (d060401*b[end] + d060402*b[end-1] + d060403*b[end-2] + d060405*b[end-4] + d060406*b[end-5] + d060407*b[end-6]) * u[end-3]
                        + (d060501*b[end] + d060502*b[end-1] + d060503*b[end-2] + d060504*b[end-3] + d060505*b[end-4] + d060506*b[end-5] + d060507*b[end-6] + d060508*b[end-7]) * u[end-4]
                        + (d060601*b[end] + d060602*b[end-1] + d060603*b[end-2] + d060604*b[end-3] + d060605*b[end-4] + d060606*b[end-5] + d060607*b[end-6] + d060608*b[end-7] + d060609*b[end-8]) * u[end-5]
                        + (d060704*b[end-3] + d060705*b[end-4] + d060706*b[end-5] + d060707*b[end-6] + d060708*b[end-7] + d060709*b[end-8]) * u[end-6]
                        + (d060805*b[end-4] + d060806*b[end-5] + d060807*b[end-6] + d060808*b[end-7] + d060809*b[end-8]) * u[end-7]
                        + (d060905*b[end-4] + d060906*b[end-5] + d060907*b[end-6] + d060908*b[end-7] + d060909*b[end-8]) * u[end-8]
        ) + β*dest[end-5]
        dest[end-6] = α * (
                          (d070104*b[end-3] + d070105*b[end-4] + d070106*b[end-5] + d070107*b[end-6]) * u[end]
                        + (d070204*b[end-3] + d070205*b[end-4] + d070206*b[end-5] + d070207*b[end-6]) * u[end-1]
                        + (d070304*b[end-3] + d070305*b[end-4] + d070306*b[end-5] + d070307*b[end-6]) * u[end-2]
                        + (d070405*b[end-4] + d070406*b[end-5] + d070407*b[end-6]) * u[end-3]
                        + (d070504*b[end-3] + d070505*b[end-4] + d070506*b[end-5] + d070507*b[end-6] + d070508*b[end-7]) * u[end-4]
                        + (d070604*b[end-3] + d070605*b[end-4] + d070606*b[end-5] + d070607*b[end-6] + d070608*b[end-7] + d070609*b[end-8]) * u[end-5]
                        + (d070704*b[end-3] + d070705*b[end-4] + d070706*b[end-5] + d070707*b[end-6] + d070708*b[end-7] + d070709*b[end-8] + d070710*b[end-9]) * u[end-6]
                        + (d070805*b[end-4] + d070806*b[end-5] + d070807*b[end-6] + d070808*b[end-7] + d070809*b[end-8] + d070810*b[end-9]) * u[end-7]
                        + (d070905*b[end-4] + d070906*b[end-5] + d070907*b[end-6] + d070908*b[end-7] + d070909*b[end-8] + d070910*b[end-9]) * u[end-8]
                        + (d071010*b[end-9] + d071007*b[end-6] + d071008*b[end-7] + d071009*b[end-8]) * u[end-9]
        ) + β*dest[end-6]
        dest[end-7] = α * (
                          (d080105*b[end-4] + d080106*b[end-5] + d080107*b[end-6]) * u[end]
                        + (d080205*b[end-4] + d080206*b[end-5] + d080207*b[end-6]) * u[end-1]
                        + (d080305*b[end-4] + d080306*b[end-5] + d080307*b[end-6]) * u[end-2]
                        + (d080405*b[end-4] + d080406*b[end-5] + d080407*b[end-6]) * u[end-3]
                        + (d080505*b[end-4] + d080506*b[end-5] + d080507*b[end-6] + d080508*b[end-7]) * u[end-4]
                        + (d080605*b[end-4] + d080606*b[end-5] + d080607*b[end-6] + d080608*b[end-7] + d080609*b[end-8]) * u[end-5]
                        + (d080705*b[end-4] + d080706*b[end-5] + d080707*b[end-6] + d080708*b[end-7] + d080709*b[end-8] + d080710*b[end-9]) * u[end-6]
                        + (d080805*b[end-4] + d080806*b[end-5] + d080807*b[end-6] + d080808*b[end-7] + d080809*b[end-8] + d080810*b[end-9] + d080811*b[end-10]) * u[end-7]
                        + (d080910*b[end-9] + d080911*b[end-10] + d080905*b[end-4] + d080906*b[end-5] + d080907*b[end-6] + d080908*b[end-7] + d080909*b[end-8]) * u[end-8]
                        + (d081010*b[end-9] + d081011*b[end-10] + d081007*b[end-6] + d081008*b[end-7] + d081009*b[end-8]) * u[end-9]
                        + (d081110*b[end-9] + d081111*b[end-10] + d081108*b[end-7] + d081109*b[end-8]) * u[end-10]
        ) + β*dest[end-7]
        dest[end-8] = α * (
                          (d090105*b[end-4] + d090106*b[end-5] + d090107*b[end-6]) * u[end]
                        + (d090205*b[end-4] + d090206*b[end-5] + d090207*b[end-6]) * u[end-1]
                        + (d090305*b[end-4] + d090306*b[end-5] + d090307*b[end-6]) * u[end-2]
                        + (d090405*b[end-4] + d090406*b[end-5] + d090407*b[end-6]) * u[end-3]
                        + (d090505*b[end-4] + d090506*b[end-5] + d090507*b[end-6] + d090508*b[end-7]) * u[end-4]
                        + (d090605*b[end-4] + d090606*b[end-5] + d090607*b[end-6] + d090608*b[end-7] + d090609*b[end-8]) * u[end-5]
                        + (d090705*b[end-4] + d090706*b[end-5] + d090707*b[end-6] + d090708*b[end-7] + d090709*b[end-8] + d090710*b[end-9]) * u[end-6]
                        + (d090810*b[end-9] + d090811*b[end-10] + d090805*b[end-4] + d090806*b[end-5] + d090807*b[end-6] + d090808*b[end-7] + d090809*b[end-8]) * u[end-7]
                        + (d090910*b[end-9] + d090911*b[end-10] + d090912*b[end-11] + d090905*b[end-4] + d090906*b[end-5] + d090907*b[end-6] + d090908*b[end-7] + d090909*b[end-8]) * u[end-8]
                        + (d091010*b[end-9] + d091011*b[end-10] + d091012*b[end-11] + d091007*b[end-6] + d091008*b[end-7] + d091009*b[end-8]) * u[end-9]
                        + (d091110*b[end-9] + d091111*b[end-10] + d091112*b[end-11] + d091108*b[end-7] + d091109*b[end-8]) * u[end-10]
                        + (d091210*b[end-9] + d091211*b[end-10] + d091212*b[end-11] + d091209*b[end-8]) * u[end-11]
        ) + β*dest[end-8]
    end
end

@inline function convolve_interior_coefficients_loopbody(i, cache::Mattsson2012Cache6, u, b)
    @unpack di_im3_im3, di_im3_im2, di_im3_im1, di_im3_i_0, 
    di_im2_im3, di_im2_im2, di_im2_im1, di_im2_i_0, di_im2_ip1,
    di_im1_im3, di_im1_im2, di_im1_im1, di_im1_i_0, di_im1_ip1, di_im1_ip2,
    di_i_0_im3, di_i_0_im2, di_i_0_im1, di_i_0_i_0, di_i_0_ip1, di_i_0_ip2, di_i_0_ip3 = cache

    @inbounds begin
        b_im3 = b[i-3]
        b_im2 = b[i-2]
        b_im1 = b[i-1]
        b_i_0 = b[i]
        b_ip1 = b[i+1]
        b_ip2 = b[i+2]
        b_ip3 = b[i+3]

        retval = (
                      (di_im3_im3*b_im3 + di_im3_im2*b_im2 + di_im3_im1*b_im1 + di_im3_i_0*b_i_0) * u[i-3]
                    + (di_im2_im3*b_im3 + di_im2_im2*b_im2 + di_im2_im1*b_im1 + di_im2_i_0*b_i_0 + di_im2_ip1*b_ip1) * u[i-2]
                    + (di_im1_im3*b_im3 + di_im1_im2*b_im2 + di_im1_im1*b_im1 + di_im1_i_0*b_i_0 + di_im1_ip1*b_ip1 + di_im1_ip2*b_ip2) * u[i-1]
                    + (di_i_0_im3*b_im3 + di_i_0_im2*b_im2 + di_i_0_im1*b_im1 + di_i_0_i_0*b_i_0 + di_i_0_ip1*b_ip1 + di_i_0_ip2*b_ip2 + di_i_0_ip3*b_ip3) * u[i]
                    + (di_im1_im3*b_ip3 + di_im1_im2*b_ip2 + di_im1_im1*b_ip1 + di_im1_i_0*b_i_0 + di_im1_ip1*b_im1 + di_im1_ip2*b_im2) * u[i+1]
                    + (di_im2_im3*b_ip3 + di_im2_im2*b_ip2 + di_im2_im1*b_ip1 + di_im2_i_0*b_i_0 + di_im2_ip1*b_im1) * u[i+2]
                    + (di_im3_im3*b_ip3 + di_im3_im2*b_ip2 + di_im3_im1*b_ip1 + di_im3_i_0*b_i_0) * u[i+3]
                )
    end

    retval
end
