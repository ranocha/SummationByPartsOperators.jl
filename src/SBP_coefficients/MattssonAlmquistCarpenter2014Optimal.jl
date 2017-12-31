
"""
    MattssonAlmquistCarpenter2014Optimal

Coefficients of the optimal SBP operators with nonuniform grid given in
  Mattsson, Almquist, Carpenter (2014)
  Optimal diagonal-norm SBP operators.
  Journal of Computational Physics 264, pp. 91-111.
"""
struct MattssonAlmquistCarpenter2014Optimal <: SourceOfCoefficients end

function Base.show(io::IO, ::MattssonAlmquistCarpenter2014Optimal)
    print(io,
        "  Mattsson, Almquist, Carpenter (2014) \n",
        "  Optimal diagonal-norm SBP operators. \n",
        "  Journal of Computational Physics 264, pp. 91-111. \n")
end



struct MattssonAlmquistCarpenter2014OptimalGrid{T,Grid} <: AbstractArray{T,1}
    xmin::T
    xmax::T
    d1::T
    d2::T
    d3::T
    uniform_grid::Grid
end

function MattssonAlmquistCarpenter2014OptimalGrid{T}(xmin::T, xmax::T, d1::T, d2::T, d3::T, N::Int)
    @argcheck xmin < xmax
    @argcheck N > 6
    @argcheck d1 > 0
    @argcheck d2 > 0
    @argcheck d3 > 0

    uniform_grid = linspace(xmin+d1+d2+d3, xmax-d1-d2-d3, N-6)
    MattssonAlmquistCarpenter2014OptimalGrid{T,typeof(uniform_grid)}(xmin, xmax, d1, d2, d3, uniform_grid)
end

function MattssonAlmquistCarpenter2014OptimalGrid(xmin, xmax, d1, d2, d3, N)
    MattssonAlmquistCarpenter2014OptimalGrid(promote(xmin, xmax, d1, d2, d3)..., N)
end

function Base.length(grid::MattssonAlmquistCarpenter2014OptimalGrid)
    length(grid.uniform_grid) + 6
end

function Base.getindex(grid::MattssonAlmquistCarpenter2014OptimalGrid, i::Int)
    N = length(grid)
    @boundscheck begin
        @argcheck i > 0
        @argcheck i <= N
    end

    if i == 1
        grid.xmin
    elseif i == 2
        grid.xmin + grid.d1
    elseif i == 3
        grid.xmin + grid.d1 + grid.d2
    elseif i == N
        grid.xmax
    elseif i == N-1
        grid.xmax - grid.d1
    elseif i == N-2
        grid.xmax - grid.d1 - grid.d2
    else
        @inbounds x = getindex(grid.uniform_grid, i-3)
        x
    end
end

Base.size(grid::MattssonAlmquistCarpenter2014OptimalGrid) = (length(grid),)
Base.step(grid::MattssonAlmquistCarpenter2014OptimalGrid) = step(grid.uniform_grid)


function construct_grid(::MattssonAlmquistCarpenter2014Optimal, accuracy_order, xmin, xmax, N)
    @argcheck N > 6
    T = promote_type(typeof(xmin), typeof(xmax))

    if accuracy_order == 2
        d1 = T(78866488858096586513//Int128(10)^20)
        d2 = T(95915098594220826013//Int128(10)^20)
        d3 = T(1)
    elseif accuracy_order == 4
        d1 = T(72181367003646814327//Int128(10)^20)
        d2 = T(13409118421582217252//Int128(10)^19)
        d3 = T(12898797485951900258//Int128(10)^19)
    elseif accuracy_order == 6
        d1 = T(51670081689316731234//Int128(10)^20)
        d2 = T(98190527037374634269//Int128(10)^20)
        d3 = T(10868393364992957832//Int128(10)^19)
    elseif accuracy_order == 8
        d1 = T(41669687672575697416//Int128(10)^20)
        d2 = T(78703773886730090312//Int128(10)^20)
        d3 = T(92685925671601406028//Int128(10)^20)
    else
        throw(ArgumentError("Accuracy order $accuracy_order not implemented/derived."))
    end

    h = (xmax-xmin) / (2*(d1+d2+d3) + N - 7)

    MattssonAlmquistCarpenter2014OptimalGrid(xmin, xmax, d1*h, d2*h, d3*h, N)
end


function first_derivative_coefficients(source::MattssonAlmquistCarpenter2014Optimal, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,3}(SVector(T(-50000000000000000000//33743097329453577701),
                                                    T(55932483188770411252//33743097329453577701),
                                                    T(-5932483188770411252//33743097329453577701) )),
            # d2
            DerivativeCoefficientRow{T,1,3}(SVector(T(-13983120797192602813//24439920504708372824),
                                                    T(0),
                                                    T(13983120797192602813//24439920504708372824) )),
            # d3
            DerivativeCoefficientRow{T,1,4}(SVector(T(2966241594385205626//46639404052015171765),
                                                    T(-27966241594385205626//46639404052015171765),
                                                    T(0),
                                                    T(5000000000000000000//9327880810403034353) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(1//2))
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
    elseif order == 4
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,5}(SVector(T(-50000000000000000000//21427296612044126417),
                                                    T(66884898686930380508//21427296612044126417),
                                                    T(-25171531878753856238//21427296612044126417),
                                                    T(10997619816825822803//21427296612044126417),
                                                    T(-13554933125011735296//107136483060220632085) )),
            # d2
            DerivativeCoefficientRow{T,1,5}(SVector(T(-16721224671732595127//28093989712218483700),
                                                    T(0),
                                                    T(92214436948640491071//112375958848873934800),
                                                    T(-16206184326771260201//56187979424436967400),
                                                    T(17707075979581024571//280939897122184837000) )),
            # d3
            DerivativeCoefficientRow{T,1,5}(SVector(T(12585765939376928119//71722939624706300000),
                                                    T(-92214436948640491071//143445879249412600000),
                                                    T(0),
                                                    T(1636075617843355867//2868917584988252000),
                                                    T(-14760875822281158529//143445879249412600000) )),
            # d4
            DerivativeCoefficientRow{T,1,6}(SVector(T(-10997619816825822803//109173230217361308360),
                                                    T(16206184326771260201//54586615108680654180),
                                                    T(-8180378089216779335//10917323021736130836),
                                                    T(0),
                                                    T(17180591347196107273//27293307554340327090),
                                                    T(-625000000000000000//8187992266302098127) )),
            # d5
            DerivativeCoefficientRow{T,1,7}(SVector(T(2259155520835289216//82365134292746683125),
                                                    T(-17707075979581024571//247095402878240049375),
                                                    T(14760875822281158529//98838161151296019750),
                                                    T(-34361182694392214546//49419080575648009875),
                                                    T(0),
                                                    T(800000000000000000//1186057933815552237),
                                                    T(-100000000000000000//1186057933815552237) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(2//3), T(-1//12))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(21427296612044126417//100000000000000000000),
                                T(280939897122184837//250000000000000000),
                                T(717229396247063//500000000000000),
                                T(2729330755434032709//2500000000000000000),
                                T(395352644605184079//400000000000000000) )
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
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
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
    elseif order == 8
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,8}(SVector(T(),
                                                    T(),
                                                    T(),
                                                    T() )),
            # d2
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(4//5), T(-1//5), T(4//105), T(-1//280))
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
