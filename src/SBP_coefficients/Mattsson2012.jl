
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
    elseif accuracy_order == 6
        coefficient_cache = Mattsson2012Cache6(T)
    elseif accuracy_order == 8
        coefficient_cache = Mattsson2012Cache8(T)
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

lower_bandwidth(cache::Mattsson2012Cache2) = 1
upper_bandwidth(cache::Mattsson2012Cache2) = 1
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


    function Mattsson2012Cache4(::Type{T}) where {T}
        d111 = T(920//289)
        d112 = T(-59//68)
        d113 = T(-81031200387//366633756146)
        d114 = T(-69462376031//733267512292)
        d121 = T(-1740//289)
        d123 = T(6025413881//7482321554)
        d124 = T(1612249989//7482321554)
        d131 = T(1128//289)
        d132 = T(59*b_2/68)
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
        d231 = T(1//17)
        d233 = T(1489556735319//1857857371138)
        d234 = T(149729180391//1857857371138)
        d241 = T(3//68)
        d243 = T(-13235456910147//159775733917868)
        d244 = T(3093263736297//79887866958934)
        d253 = T(67535018271//2349643145851)
        d263 = T(-67535018271//2349643145851)
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
        d451 = T(-9//3332)
        d453 = T(-217407431400324796377//207908333536359652868)
        d454 = T(-1950062198436997//3914505480987766)
        d455 = T(-7476412//9108757)
        d456 = T(-2//49)
        d463 = T(4959271814984644613//21402328452272317207)
        d464 = T(47996144728947//115132514146699)
        d465 = T(4502124//9108757)
        d466 = T(8//49)
        d473 = T(-2258420001//1653303156799)
        d474 = T(-1063649//8893843)
        d475 = T(1473580//9108757)
        d476 = T(-6//49)

        d513 = T(-49579087//10149031312)
        d514 = T(49579087//10149031312)
        d523 = T(1328188692663//37594290333616)
        d534 = T(-1328188692663//37594290333616)
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

        new{T}(d111, d112, d113, d114, d121, d123, d123, d131, d132, d133, d134,
                d141, d143, d144, d153, d154, d163, d164,
                d211, d213, d214, d221, d223, d224, d231, d233, d234, d241, d243, d244,
                d253, d254, d263, d264,
                d311, d312 d313, d314, d321, d323, d323, d331, d332, d33, d334, d335,
                d341, d343, d344, d345, d353, d354, d355, d363, d364, d365,
                d411, d413, d414, d421, d423, d424, d431, d433, d434, d435,
                d451, d453, d44, d455, d456, d463, d464, d465, d466,
                d473, d474, d475, d476,
                d513, d514, d523, d524, d533, d534, d535, d543, d544, d545, d546,
                d553, d554, d555, d556, d557, d563, d564, d565, d566, d567,
                d575, d576, d577)
    end
end

lower_bandwidth(cache::Mattsson2012Cache4) = 1
upper_bandwidth(cache::Mattsson2012Cache4) = 1
Base.checkbounds(::Type{Bool}, u::AbstractVector, ::Mattsson2012Cache4) = length(u) > 2
left_length(::Mattsson2012Cache4) = 1
right_length(::Mattsson2012Cache4) = 1

function convolve_boundary_coefficients!(dest::AbstractVector, cache::Mattsson2012Cache4,
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

function convolve_boundary_coefficients!(dest::AbstractVector, cache::Mattsson2012Cache4,
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

@inline function convolve_interior_coefficients_loopbody(i, cache::Mattsson2012Cache4, u, b)
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
