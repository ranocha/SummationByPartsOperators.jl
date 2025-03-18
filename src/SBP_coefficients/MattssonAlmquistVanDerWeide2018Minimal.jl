
"""
    MattssonAlmquistVanDerWeide2018Minimal()

Coefficients of the optimized SBP operators with nonuniform grid given in
- Mattsson, Almquist, van der Weide (2018)
  Boundary optimized diagonal-norm SBP operators.
  Journal of Computational Physics 374, pp. 1261-1266.
"""
struct MattssonAlmquistVanDerWeide2018Minimal <: SourceOfCoefficients end

function Base.show(io::IO, source::MattssonAlmquistVanDerWeide2018Minimal)
    if get(io, :compact, false)
        summary(io, source)
    else
        print(
            io,
            "Mattsson, Almquist, van der Weide (2018) \n",
            "  Boundary optimized diagonal-norm SBP operators ('Minimal'). \n",
            "  Journal of Computational Physics 374, pp. 1261-1266.",
        )
    end
end



struct BoundaryAdaptedGrid{T,M,Grid} <: AbstractArray{T,1}
    xmin::T
    xmax::T
    xstart::SVector{M,T}
    uniform_grid::Grid
    Δx::T
end

function BoundaryAdaptedGrid(xmin::T, xmax::T, _xstart::SVector{M,T}, N::Int) where {M,T}
    @argcheck xmin < xmax
    @argcheck N > 2M

    Δx = (xmax - xmin) / (2 * _xstart[end] + N + 1 - 2M)
    xstart = Δx .* _xstart
    if N - 2M == 1 # This is to avoid an error if starting and end points are not the same due to rounding errors
        uniform_grid = range(xmax - xstart[end] - Δx, xmax - xstart[end] - Δx, length = 1)
    else
        uniform_grid =
            range(xmin + xstart[end] + Δx, xmax - xstart[end] - Δx, length = N - 2M)
    end
    BoundaryAdaptedGrid{T,M,typeof(uniform_grid)}(xmin, xmax, xstart, uniform_grid, Δx)
end

function BoundaryAdaptedGrid(xmin, xmax, xstart, N)
    T = promote_type(typeof(xmin), typeof(xmax), eltype(xstart))
    BoundaryAdaptedGrid(T(xmin), T(xmax), T.(xstart), N)
end

function Base.length(grid::BoundaryAdaptedGrid)
    length(grid.uniform_grid) + 2 * length(grid.xstart)
end

function Base.getindex(grid::BoundaryAdaptedGrid, i::Int)
    N = length(grid)
    M = length(grid.xstart)
    @boundscheck begin
        @argcheck i > 0
        @argcheck i <= N
    end

    if i <= M
        grid.xmin + grid.xstart[i]
    elseif i >= N - M + 1
        grid.xmax - grid.xstart[1+N-i]
    else
        @inbounds x = getindex(grid.uniform_grid, i - M)
        x
    end
end

Base.size(grid::BoundaryAdaptedGrid) = (length(grid),)
Base.step(grid::BoundaryAdaptedGrid) = grid.Δx


function construct_grid(
    ::MattssonAlmquistVanDerWeide2018Minimal,
    accuracy_order,
    xmin,
    xmax,
    N,
)
    @argcheck N > 6
    T = promote_type(typeof(xmin), typeof(xmax))

    if accuracy_order == 4
        xstart = SVector(
            T(0.0000000000000e+00),
            T(7.7122987842562e-01),
            T(1.7712298784256e+00),
            T(2.7712298784256e+00),
        )
    elseif accuracy_order == 6
        xstart = SVector(
            T(0.0000000000000e+00),
            T(4.0842950991998e-01),
            T(1.1968523189207e+00),
            T(2.1968523189207e+00),
            T(3.1968523189207e+00),
            T(4.1968523189207e+00),
        )
    elseif accuracy_order == 8
        xstart = SVector(
            T(0.0000000000000e+00),
            T(4.9439570885261e-01),
            T(1.4051531374839e+00),
            T(2.4051531374839e+00),
            T(3.4051531374839e+00),
            T(4.4051531374839e+00),
            T(5.4051531374839e+00),
        )
    else
        throw(ArgumentError("Accuracy order $accuracy_order not implemented/derived."))
    end

    BoundaryAdaptedGrid(xmin, xmax, xstart, N)
end


function first_derivative_coefficients(
    source::MattssonAlmquistVanDerWeide2018Minimal,
    order::Int,
    T = Float64,
    mode = FastMode(),
)
    if order == 4
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,3}(
                SVector(
                    T(-1.8612097181862932),
                    T(2.2966302628676907),
                    T(-0.4354205446813973),
                ),
            ),
            # d2
            DerivativeCoefficientRow{T,1,4}(
                SVector(
                    T(-0.611186527373936),
                    T(0),
                    T(0.6937383659498486),
                    T(-0.08255183857591565),
                ),
            ),
            # d3
            DerivativeCoefficientRow{T,1,5}(
                SVector(
                    T(0.11778272127063323),
                    T(-0.7051567887047331),
                    T(0),
                    T(0.6712846484961179),
                    T(-0.08391058106201474),
                ),
            ),
        )
        right_boundary = .-left_boundary
        upper_coef = SVector(T(2 // 3), T(-1 // 12))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights =
            SVector(T(2.6864248295847e-01), T(1.0094667153500e+00), T(9.9312068011715e-01))
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(
            left_boundary,
            right_boundary,
            left_boundary_derivatives,
            right_boundary_derivatives,
            lower_coef,
            central_coef,
            upper_coef,
            left_weights,
            right_weights,
            mode,
            1,
            order,
            source,
        )
    elseif order == 6
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,5}(
                SVector(
                    T(-3.924566448353279),
                    T(4.962014957077398),
                    T(-1.2881968205658438),
                    T(0.28645730739095715),
                    T(-0.03570899554923223),
                ),
            ),
            # d2
            DerivativeCoefficientRow{T,1,5}(
                SVector(
                    T(-1.022587534557978),
                    T(0),
                    T(1.3023996740158794),
                    T(-0.3357380168806994),
                    T(0.05592587742279807),
                ),
            ),
            # d3
            DerivativeCoefficientRow{T,1,6}(
                SVector(
                    T(0.17402334848464412),
                    T(-0.8537429915394141),
                    T(0),
                    T(0.8419418935035105),
                    T(-0.179894658468737),
                    T(0.017672408019993497),
                ),
            ),
            # d4
            DerivativeCoefficientRow{T,1,7}(
                SVector(
                    T(-0.036159059808532096),
                    T(0.20564365765539105),
                    T(-0.7867088686169533),
                    T(0),
                    T(0.749328780370319),
                    T(-0.1486175733002554),
                    T(0.01651306370002838),
                ),
            ),
            # d5
            DerivativeCoefficientRow{T,1,8}(
                SVector(
                    T(0.004554664962414826),
                    T(-0.03461379228073425),
                    T(0.16985243300141972),
                    T(-0.7571710333787599),
                    T(0),
                    T(0.7508648039541874),
                    T(-0.15017296079083747),
                    T(0.016685884532315277),
                ),
            ),
        )
        right_boundary = .-left_boundary
        upper_coef = SVector(T(3 // 4), T(-3 // 20), T(1 // 60))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector(
            T(0.12740260779883),
            T(0.61820981002054),
            T(0.94308973897679),
            T(1.0093019060199),
            T(0.99884825610465),
        )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(
            left_boundary,
            right_boundary,
            left_boundary_derivatives,
            right_boundary_derivatives,
            lower_coef,
            central_coef,
            upper_coef,
            left_weights,
            right_weights,
            mode,
            1,
            order,
            source,
        )
    elseif order == 8
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,6}(
                SVector(
                    T(-3.4425783018277825),
                    T(4.592216499078388),
                    T(-1.5780326084197085),
                    T(0.5114508898159292),
                    T(-0.08276413469602352),
                    T(-0.00029234395080269453),
                ),
            ),
            # d2
            DerivativeCoefficientRow{T,1,6}(
                SVector(
                    T(-0.8677229176949711),
                    T(0),
                    T(1.1480053882687404),
                    T(-0.34675581554193646),
                    T(0.07194779962833377),
                    T(-0.005474454660166564),
                ),
            ),
            # d3
            DerivativeCoefficientRow{T,1,7}(
                SVector(
                    T(0.2312364265180344),
                    T(-0.8902776878396409),
                    T(0),
                    T(0.8364384829720679),
                    T(-0.21345201833077546),
                    T(0.039658060554022004),
                    T(-0.0036032638737076134),
                ),
            ),
            # d4
            DerivativeCoefficientRow{T,1,8}(
                SVector(
                    T(-0.07428870791937443),
                    T(0.26655320361241747),
                    T(-0.8291108452061701),
                    T(0),
                    T(0.8031154600461715),
                    T(-0.200795518743896),
                    T(0.038098105611977996),
                    T(-0.0035716974011229742),
                ),
            ),
            # d5
            DerivativeCoefficientRow{T,1,9}(
                SVector(
                    T(0.012018140909062053),
                    T(-0.0552909326936426),
                    T(0.21152178487239182),
                    T(-0.8028866425212018),
                    T(0),
                    T(0.8000791459238572),
                    T(-0.19995806770071226),
                    T(0.03808725099061186),
                    T(-0.0035706797803698614),
                ),
            ),
            # d6
            DerivativeCoefficientRow{T,1,10}(
                SVector(
                    T(4.246147652870495e-5),
                    T(0.004208071691358718),
                    T(-0.03930901628651875),
                    T(0.2007872492592992),
                    T(-0.800274203791225),
                    T(0),
                    T(0.800027268482503),
                    T(-0.20000681712062576),
                    T(0.03809653659440491),
                    T(-0.00357155030572546),
                ),
            ),
        )
        right_boundary = .-left_boundary
        upper_coef = SVector(T(4 // 5), T(-1 // 5), T(4 // 105), T(-1 // 280))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector(
            T(0.14523997892351),
            T(0.76864793350174),
            T(0.99116487068535),
            T(0.99992473335107),
            T(1.0002097054636),
            T(0.99996591555866),
        )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(
            left_boundary,
            right_boundary,
            left_boundary_derivatives,
            right_boundary_derivatives,
            lower_coef,
            central_coef,
            upper_coef,
            left_weights,
            right_weights,
            mode,
            1,
            order,
            source,
        )
        # elseif order == 8
        #   left_boundary = (
        #       # d1
        #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T() )),
        #       # d2
        #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T() )),
        #       # d3
        #       DerivativeCoefficientRow{T,1,7}(SVector(T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T() )),
        #       # d4
        #       DerivativeCoefficientRow{T,1,8}(SVector(T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T() )),
        #       # d5
        #       DerivativeCoefficientRow{T,1,9}(SVector(T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T() )),
        #       # d6
        #       DerivativeCoefficientRow{T,1,10}(SVector(T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T(),
        #                                               T() )),
        #   )
        #   right_boundary = .- left_boundary
        #   upper_coef = SVector(T(4//5), T(-1//5), T(4//105), T(-1//280))
        #   central_coef = zero(T)
        #   lower_coef = -upper_coef
        #   left_weights = SVector( T(),
        #                           T(),
        #                           T(),
        #                           T(),
        #                           T(),
        #                           T() )
        #   right_weights = left_weights
        #   left_boundary_derivatives = Tuple{}()
        #   right_boundary_derivatives = left_boundary_derivatives

        #   DerivativeCoefficients(left_boundary, right_boundary,
        #                           left_boundary_derivatives, right_boundary_derivatives,
        #                           lower_coef, central_coef, upper_coef,
        #                           left_weights, right_weights, mode, 1, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end
