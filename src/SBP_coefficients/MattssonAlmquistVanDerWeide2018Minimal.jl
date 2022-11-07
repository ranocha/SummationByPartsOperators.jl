
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
  elseif accuracy_order == 12
    xstart = SVector(
      T(0.0000000000000e+00),
      T(4.6552112904489e-01),
      T(1.4647984306493e+00),
      T(2.7620429464763e+00),
      T(4.0000000000000e+00),
      T(5.0000000000000e+00),
      T(6.0000000000000e+00),
      T(7.0000000000000e+00),
      T(8.0000000000000e+00),
      T(9.0000000000000e+00),
      T(1.0000000000000e+01)
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
    left_weights = SVector( T(0.16717213975289),
                            T(0.93675739171278),
                            T(1.3035532379753),
                            T(1.1188461804303),
                            T(0.9666434592266),
                            T(1.0083235564392),
                            T(0.99858767377362),
                            T(1.0001163606893) )
    right_weights = left_weights
    left_boundary_derivatives = Tuple{}()
    right_boundary_derivatives = left_boundary_derivatives

    DerivativeCoefficients(left_boundary, right_boundary,
                            left_boundary_derivatives, right_boundary_derivatives,
                            lower_coef, central_coef, upper_coef,
                            left_weights, right_weights, mode, 1, order, source)
  elseif order == 12
    left_boundary = (
        # d1
        DerivativeCoefficientRow{T,1,10}(SVector(T(-7.684270470438172),
                                                T(5.194807555458745),
                                                T(-2.0579295430730933),
                                                T(1.0796638584520148),
                                                T(-0.41550889074234293),
                                                T(-0.09126595764586296),
                                                T(0.202101651484316),
                                                T(-0.07536518127096767),
                                                T(0.003292629927279952),
                                                T(0.0023391126289983594) )),
        # d2
        DerivativeCoefficientRow{T,1,8}(SVector(T(-0.8878088485050056),
                                                T(0),
                                                T(1.250283162562752),
                                                T(-0.5424654450963853),
                                                T(0.2028021004571883),
                                                T(0.034847787734123656),
                                                T(-0.0891304560108277),
                                                T(0.03471204698149003),
                                                T(-0.0024142417515724767),
                                                T(-0.0008261063717634195) )),
        # d3
        DerivativeCoefficientRow{T,1,8}(SVector(T(0.22346937002606296),
                                                T(-0.7944121537113519),
                                                T(0),
                                                T(0.7879098659795528),
                                                T(0.2534158843840873),
                                                T(-0.014215698959584674),
                                                T(0.08139548515731004),
                                                T(-0.035492265967731106),
                                                T(0.004461798639218733),
                                                T(0.00029948322062261124) )),
        # d4
        DerivativeCoefficientRow{T,1,8}(SVector(T(-0.10532369420351033),
                                                T(0.3096418068834778),
                                                T(-0.7078260723043035),
                                                T(0),
                                                T(0.6099618580432923),
                                                T(-0.06298823739309845),
                                                T(-0.08036793017152877),
                                                T(0.04580031261656808),
                                                T(-0.00945925508773524),
                                                T(0.0005612116168350629) )),
        # d5
        DerivativeCoefficientRow{T,1,9}(SVector(T(0.049373250385308945),
                                                T(-0.14100478191293775),
                                                T(0.2773050181486685),
                                                T(-0.7429790311918253),
                                                T(0),
                                                T(0.6221830719327455),
                                                T(-0.045737501160932524),
                                                T(-0.029361570393079246),
                                                T(0.011203816364579948),
                                                T(-0.0008175732267816134),
                                                T(-0.00016469894575332322)  )),
        # d6
        DerivativeCoefficientRow{T,1,10}(SVector(T(0.012172895383088677),
                                                T(-0.027196328696641285),
                                                T(0.017460866406897767),
                                                T(0.08612059273501356),
                                                T(-0.6983801181665034),
                                                T(0),
                                                T(0.7536732748892866),
                                                T(-0.18003027136972943),
                                                T(0.046996413958795336),
                                                T(-0.013294572041918046),
                                                T(0.002662116073478288),
                                                T(-0.00018486917176933042) )),
        # d7
        DerivativeCoefficientRow{T,1,10}(SVector(T(-0.02613877695266672),
                                                T(0.06745148621069272),
                                                T(-0.09694562474753163),
                                                T(0.10655176894259812),
                                                T(0.04978246944553149),
                                                T(-0.7308250779753221),
                                                T(0),
                                                T(0.8167927121559334),
                                                T(-0.2468944631824365),
                                                T(0.07557056593974017),
                                                T(-0.017747207053940086),
                                                T(0.0025814119351186032),
                                                T(-0.00017926471771657437) )),
        # d8
        DerivativeCoefficientRow{T,1,10}(SVector(T(0.009820060996139302),
                                                T(-0.02646510714247867),
                                                T(0.04258823695499128),
                                                T(-0.061175053786101385),
                                                T(0.03219669689820139),
                                                T(0.1758749143714943),
                                                T(-0.8228864117869611),
                                                T(0),
                                                T(0.85314316830222),
                                                T(-0.2671018982641268),
                                                T(0.07946493544275793),
                                                T(-0.01787961047462061),
                                                T(0.002600670614490316),
                                                T(-0.00018060212600627668) )),
        # d9
        DerivativeCoefficientRow{T,1,10}(SVector(T(-0.0004284166502310784),
                                                T(0.0018380366654838602),
                                                T(-0.005346208397613871),
                                                T(0.012616614614814827),
                                                T(-0.01226811967154658),
                                                T(-0.0458461705888513),
                                                T(0.24838157124468524),
                                                T(-0.8519260536187364),
                                                T(0),
                                                T(0.8568762094346083),
                                                T(-0.2678115448464915),
                                                T(0.07935156884340752),
                                                T(-0.017854102989766765),
                                                T(0.0025969604348752115),
                                                T(-0.00018034447464411663) )),
        # d10
        DerivativeCoefficientRow{T,1,10}(SVector(T(-0.0003044061228920383),
                                                T(0.0006290544170151834),
                                                T(-0.0003589113909262337),
                                                T(-0.0007486716402445799),
                                                T(0.0008954008719030808),
                                                T(0.012971539942668298),
                                                T(-0.07603954755787204),
                                                T(0.2667692618970865),
                                                T(-0.8570317576203546),
                                                T(0),
                                                T(0.8571525136112116),
                                                T(-0.2678601605035305),
                                                T(0.07936597348253019),
                                                T(-0.017857344033569367),
                                                T(0.002597431859428317),
                                                T(-0.00018037721246030451) )),
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(6//7), T(-15//56), T(5//63), T(-1//56), T(1//385), T(-1//5544))
    central_coef = zero(T)
    lower_coef = -upper_coef
    left_weights = SVector( T(0.1301359711175),
                            T(0.7614604507902),
                            T(1.1984222247012),
                            T(1.3340123109301),
                            T(1.0951811473364),
                            T(0.9756909637713),
                            T(1.0061945410831),
                            T(0.99874339446564),
                            T(1.0001702615573),
                            T(0.99998873424721) )
    right_weights = left_weights
    left_boundary_derivatives = Tuple{}()
    right_boundary_derivatives = left_boundary_derivatives

    DerivativeCoefficients(left_boundary, right_boundary,
                            left_boundary_derivatives, right_boundary_derivatives,
                            lower_coef, central_coef, upper_coef,
                            left_weights, right_weights, mode, 1, order, source)
  else
    throw(ArgumentError("Order $order not implemented/derived."))
  end
end
