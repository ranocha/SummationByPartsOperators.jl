
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
    print(io,
        "Mattsson, Almquist, van der Weide (2018) \n",
        "  Boundary optimized diagonal-norm SBP operators ('Minimal'). \n",
        "  Journal of Computational Physics 374, pp. 1261-1266.")
  end
end



struct BoundaryAdaptedGrid{T,M,Grid} <: AbstractArray{T,1}
  xmin::T
  xmax::T
  xstart::SVector{M,T}
  uniform_grid::Grid
end

function BoundaryAdaptedGrid(xmin::T, xmax::T, _xstart::SVector{M,T}, N::Int) where {M,T}
  @argcheck xmin < xmax
  @argcheck N > 2M

  Δx = (xmax - xmin) / (2*_xstart[end] + N + 1 - 2M)
  xstart = Δx .* _xstart
  uniform_grid = range(xmin + xstart[end] + Δx, xmax - xstart[end] - Δx, length=N-2M)
  BoundaryAdaptedGrid{T,M,typeof(uniform_grid)}(xmin, xmax, xstart, uniform_grid)
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
  elseif i >= N-M+1
    grid.xmax - grid.xstart[1+N-i]
  else
    @inbounds x = getindex(grid.uniform_grid, i-M)
    x
  end
end

Base.size(grid::BoundaryAdaptedGrid) = (length(grid),)
Base.step(grid::BoundaryAdaptedGrid) = step(grid.uniform_grid)


function construct_grid(::MattssonAlmquistVanDerWeide2018Minimal, accuracy_order, xmin, xmax, N)
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
  elseif accuracy_order == 10
    xstart = SVector(
      T(0.0000000000000e+00),
      T(5.8556160757529e-01),
      T(1.7473267488572e+00),
      T(3.0000000000000e+00),
      T(4.0000000000000e+00),
      T(5.0000000000000e+00),
      T(6.0000000000000e+00),
      T(7.0000000000000e+00),
      T(8.0000000000000e+00),
    )
  else
    throw(ArgumentError("Accuracy order $accuracy_order not implemented/derived."))
  end

  BoundaryAdaptedGrid(xmin, xmax, xstart, N)
end


function first_derivative_coefficients(source::MattssonAlmquistVanDerWeide2018Minimal,
                                       order::Int, T=Float64, mode=FastMode())
  if order == 4
    left_boundary = (
        # d1
        DerivativeCoefficientRow{T,1,3}(SVector(T(-1.8612097181862932),
                                                T(2.2966302628676907),
                                                T(-0.4354205446813973) )),
        # d2
        DerivativeCoefficientRow{T,1,4}(SVector(T(-0.611186527373936),
                                                T(0),
                                                T(0.6937383659498486),
                                                T(-0.08255183857591565) )),
        # d3
        DerivativeCoefficientRow{T,1,5}(SVector(T(0.11778272127063323),
                                                T(-0.7051567887047331),
                                                T(0),
                                                T(0.6712846484961179),
                                                T(-0.08391058106201474) )),
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(2//3), T(-1//12))
    central_coef = zero(T)
    lower_coef = -upper_coef
    left_weights = SVector( T(2.6864248295847e-01),
                            T(1.0094667153500e+00),
                            T(9.9312068011715e-01) )
    right_weights = left_weights
    left_boundary_derivatives = Tuple{}()
    right_boundary_derivatives = left_boundary_derivatives

    DerivativeCoefficients(left_boundary, right_boundary,
                            left_boundary_derivatives, right_boundary_derivatives,
                            lower_coef, central_coef, upper_coef,
                            left_weights, right_weights, mode, 1, order, source)
  elseif order == 6
    left_boundary = (
        # d1
        DerivativeCoefficientRow{T,1,5}(SVector(T(-3.924566448353279),
                                                T(4.962014957077398),
                                                T(-1.2881968205658438),
                                                T(0.28645730739095715),
                                                T(-0.03570899554923223) )),
        # d2
        DerivativeCoefficientRow{T,1,5}(SVector(T(-1.022587534557978),
                                                T(0),
                                                T(1.3023996740158794),
                                                T(-0.3357380168806994),
                                                T(0.05592587742279807) )),
        # d3
        DerivativeCoefficientRow{T,1,6}(SVector(T(0.17402334848464412),
                                                T(-0.8537429915394141),
                                                T(0),
                                                T(0.8419418935035105),
                                                T(-0.179894658468737),
                                                T(0.017672408019993497) )),
        # d4
        DerivativeCoefficientRow{T,1,7}(SVector(T(-0.036159059808532096),
                                                T(0.20564365765539105),
                                                T(-0.7867088686169533),
                                                T(0),
                                                T(0.749328780370319),
                                                T(-0.1486175733002554),
                                                T(0.01651306370002838) )),
        # d5
        DerivativeCoefficientRow{T,1,8}(SVector(T(0.004554664962414826),
                                                T(-0.03461379228073425),
                                                T(0.16985243300141972),
                                                T(-0.7571710333787599),
                                                T(0),
                                                T(0.7508648039541874),
                                                T(-0.15017296079083747),
                                                T(0.016685884532315277) )),
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(3//4), T(-3//20), T(1//60))
    central_coef = zero(T)
    lower_coef = -upper_coef
    left_weights = SVector( T(0.12740260779883),
                            T(0.61820981002054),
                            T(0.94308973897679),
                            T(1.0093019060199),
                            T(0.99884825610465) )
    right_weights = left_weights
    left_boundary_derivatives = Tuple{}()
    right_boundary_derivatives = left_boundary_derivatives

    DerivativeCoefficients(left_boundary, right_boundary,
                            left_boundary_derivatives, right_boundary_derivatives,
                            lower_coef, central_coef, upper_coef,
                            left_weights, right_weights, mode, 1, order, source)
  elseif order == 8
    left_boundary = (
        # d1
        DerivativeCoefficientRow{T,1,6}(SVector(T(-3.4425783018277825),
                                                T(4.592216499078388),
                                                T(-1.5780326084197085),
                                                T(0.5114508898159292),
                                                T(-0.08276413469602352),
                                                T(-0.00029234395080269453) )),
        # d2
        DerivativeCoefficientRow{T,1,6}(SVector(T(-0.8677229176949711),
                                                T(0),
                                                T(1.1480053882687404),
                                                T(-0.34675581554193646),
                                                T(0.07194779962833377),
                                                T(-0.005474454660166564) )),
        # d3
        DerivativeCoefficientRow{T,1,7}(SVector(T(0.2312364265180344),
                                                T(-0.8902776878396409),
                                                T(0),
                                                T(0.8364384829720679),
                                                T(-0.21345201833077546),
                                                T(0.039658060554022004),
                                                T(-0.0036032638737076134) )),
        # d4
        DerivativeCoefficientRow{T,1,8}(SVector(T(-0.07428870791937443),
                                                T(0.26655320361241747),
                                                T(-0.8291108452061701),
                                                T(0),
                                                T(0.8031154600461715),
                                                T(-0.200795518743896),
                                                T(0.038098105611977996),
                                                T(-0.0035716974011229742) )),
        # d5
        DerivativeCoefficientRow{T,1,9}(SVector(T(0.012018140909062053),
                                                T(-0.0552909326936426),
                                                T(0.21152178487239182),
                                                T(-0.8028866425212018),
                                                T(0),
                                                T(0.8000791459238572),
                                                T(-0.19995806770071226),
                                                T(0.03808725099061186),
                                                T(-0.0035706797803698614)  )),
        # d6
        DerivativeCoefficientRow{T,1,10}(SVector(T(4.246147652870495e-5),
                                                T(0.004208071691358718),
                                                T(-0.03930901628651875),
                                                T(0.2007872492592992),
                                                T(-0.800274203791225),
                                                T(0),
                                                T(0.800027268482503),
                                                T(-0.20000681712062576),
                                                T(0.03809653659440491),
                                                T(-0.00357155030572546) )),
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(4//5), T(-1//5), T(4//105), T(-1//280))
    central_coef = zero(T)
    lower_coef = -upper_coef
    left_weights = SVector( T(0.14523997892351),
                            T(0.76864793350174),
                            T(0.99116487068535),
                            T(0.99992473335107),
                            T(1.0002097054636),
                            T(0.99996591555866) )
    right_weights = left_weights
    left_boundary_derivatives = Tuple{}()
    right_boundary_derivatives = left_boundary_derivatives

    DerivativeCoefficients(left_boundary, right_boundary,
                            left_boundary_derivatives, right_boundary_derivatives,
                            lower_coef, central_coef, upper_coef,
                            left_weights, right_weights, mode, 1, order, source)
  elseif order == 10
    left_boundary = (
        # d1
        DerivativeCoefficientRow{T,1,5}(SVector(T(0.4990746995535701),
                                                T(0.15241969704789585),
                                                T(-0.2751201290253085),
                                                T(0.10484051310768731),
                                                T(-0.012410293169942119) )),
        # d2
        DerivativeCoefficientRow{T,1,8}(SVector(T(-0.7189620019231621),
                                                T(0),
                                                T(0.981928653008664),
                                                T(-0.29102809324190426),
                                                T(-0.054434706146276914),
                                                T(0.1261228081640961),
                                                T(-0.0498454816488874),
                                                T(0.006218821787470101) )),
        # d3
        DerivativeCoefficientRow{T,1,8}(SVector(T(0.19321345045852464),
                                                T(-0.7056320348443407),
                                                T(0),
                                                T(0.60271222136483),
                                                T(-0.018486246612649363),
                                                T(-0.11746485231744051),
                                                T(0.05327862103411152),
                                                T(-0.0076211590830434816) )),
        # d4
        DerivativeCoefficientRow{T,1,8}(SVector(T(-0.07456912923348757),
                                                T(0.24366416251747966),
                                                T(-0.7022122267291899),
                                                T(0),
                                                T(0.5545701660854921),
                                                T(0.013207155223809193),
                                                T(-0.04190893546685442),
                                                T(0.006539460055673846),
                                                T(0.0007093475470824396) )),
        # d5
        DerivativeCoefficientRow{T,1,9}(SVector(T(-0.026359591690994846),
                                                T(0.05275172853187895),
                                                T(0.024929363975844176),
                                                T(-0.6418899400631001),
                                                T(0),
                                                T(0.7182069637035886),
                                                T(-0.1622764380264411),
                                                T(0.04407984796588429),
                                                T(-0.010262972170289429),
                                                T(0.0008210377736231522)  )),
        # d6
        DerivativeCoefficientRow{T,1,10}(SVector(T(0.04561276027376562),
                                                T(-0.11717119178343229),
                                                T(0.1518576925123369),
                                                T(-0.014654795161874228),
                                                T(-0.6885191359475415),
                                                T(0),
                                                T(0.8005320825781802),
                                                T(-0.22763821979760565),
                                                T(0.059032449597837176),
                                                T(-0.009838741599639595),
                                                T(0.0007870993279711656) )),
        # d7
        DerivativeCoefficientRow{T,1,10}(SVector(T(-0.01755120093038144),
                                                T(0.04675916256969975),
                                                T(-0.06954974589403921),
                                                T(0.04695596952022851),
                                                T(0.15708531311233767),
                                                T(-0.8083369920825813),
                                                T(0),
                                                T(0.8326013727190911),
                                                T(-0.2384319818363953),
                                                T(0.059607995459097826),
                                                T(-0.009934665909849703),
                                                T(0.0007947732727879743) )),
        # d8
        DerivativeCoefficientRow{T,1,10}(SVector(T(0.0020744138839504703),
                                                T(-0.005824849493654949),
                                                T(0.009933430738978299),
                                                T(-0.0073157986339953684),
                                                T(-0.042604539226373306),
                                                T(0.22950627385958502),
                                                T(-0.8313287339797591),
                                                T(0),
                                                T(0.8332363773741089),
                                                T(-0.23806753639260542),
                                                T(0.059516884098150355),
                                                T(-0.009919480683025124),
                                                T(0.0007935584546420079) )),
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(5//6), T(-5//21), T(5//84), T(-5//504), T(1//1260))
    central_coef = zero(T)
    lower_coef = -upper_coef
  else
    throw(ArgumentError("Order $order not implemented/derived."))
  end
end
