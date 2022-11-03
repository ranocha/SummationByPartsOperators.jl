
"""
    MattssonAlmquistVanDerWeide2018Accurate()

Coefficients of the optimized SBP operators with nonuniform grid given in
- Mattsson, Almquist, van der Weide (2018)
  Boundary optimized diagonal-norm SBP operators.
  Journal of Computational Physics 374, pp. 1261-1266.
"""
struct MattssonAlmquistVanDerWeide2018Accurate <: SourceOfCoefficients end

function Base.show(io::IO, source::MattssonAlmquistVanDerWeide2018Accurate)
  if get(io, :compact, false)
    summary(io, source)
  else
    print(io,
        "Mattsson, Almquist, van der Weide (2018) \n",
        "  Boundary optimized diagonal-norm SBP operators ('Accurate'). \n",
        "  Journal of Computational Physics 374, pp. 1261-1266.")
  end
end


function construct_grid(::MattssonAlmquistVanDerWeide2018Accurate, accuracy_order, xmin, xmax, N)
  @argcheck N > 6
  T = promote_type(typeof(xmin), typeof(xmax))

  if accuracy_order == 4
    xstart = SVector(
      T(0.0000000000000e+00),
      T(6.8764546205559e-01),
      T(1.8022115125776e+00),
      T(2.8022115125776e+00),
      T(3.8022115125776e+00),
    )
  elseif accuracy_order == 6
    xstart = SVector(
      T(0.0000000000000e+00),
      T(4.4090263368623e-01),
      T(1.2855984345073e+00),
      T(2.2638953951239e+00),
      T(3.2638953951239e+00),
      T(4.2638953951239e+00),
      T(5.2638953951239e+00),
    )
  elseif accuracy_order == 8
    xstart = SVector(
      T(0.0000000000000e+00),
      T(3.8118550247622e-01),
      T(1.1899550868338e+00),
      T(2.2476300175641e+00),
      T(3.3192851303204e+00),
      T(4.3192851303204e+00),
      T(5.3192851303204e+00),
      T(6.3192851303204e+00),
      T(7.3192851303204e+00),
    )
  elseif accuracy_order == 10
    xstart = SVector(
      T(0.0000000000000e+00),
      T(3.5902433622052e-01),
      T(1.1436659188355e+00),
      T(2.2144895894456e+00),
      T(3.3682742337736e+00),
      T(4.4309689056870e+00),
      T(5.4309689056870e+00),
      T(6.4309689056870e+00),
      T(7.4309689056870e+00),
      T(8.4309689056870e+00),
      T(9.4309689056870e+00),
    )
  else
    throw(ArgumentError("Accuracy order $accuracy_order not implemented/derived."))
  end

  BoundaryAdaptedGrid(xmin, xmax, xstart, N)
end


function first_derivative_coefficients(source::MattssonAlmquistVanDerWeide2018Accurate,
                                       order::Int, T=Float64, mode=FastMode())
  if order == 4
    left_boundary = (
        # d1
        DerivativeCoefficientRow{T,1,4}(SVector(T(-2.3518634632279443),
                                                T(3.0858932129093573),
                                                T(-0.9349061508864488),
                                                T(0.20087640120503583) )),
        # d2
        DerivativeCoefficientRow{T,1,4}(SVector(T(-0.6394095807755721),
                                                T(0),
                                                T(0.7917608885016885),
                                                T(-0.15235130772611635) )),
        # d3
        DerivativeCoefficientRow{T,1,5}(SVector(T(0.18446061657605967),
                                                T(-0.7539307161467195),
                                                T(0),
                                                T(0.6468087329936288),
                                                T(-0.07733863342297213) )),
        # d4
        DerivativeCoefficientRow{T,1,6}(SVector(T(-0.043308971097943676),
                                                T(0.1585246807786971),
                                                T(-0.7067880256335334),
                                                T(0),
                                                T(0.676082646803181),
                                                T(-0.08451033085039762) )),
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(2//3), T(-1//12))
    central_coef = zero(T)
    lower_coef = -upper_coef
    left_weights = SVector( T(0.21259737557798),
                            T(1.0260290400758),
                            T(1.0775123588954),
                            T(0.98607273802835) )
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
        DerivativeCoefficientRow{T,1,6}(SVector(T(-3.837232862086773),
                                                T(5.068376175016296),
                                                T(-1.6963756420892036),
                                                T(0.5851295073874213),
                                                T(-0.13275449453716695),
                                                T(0.012857316309426029) )),
        # d2
        DerivativeCoefficientRow{T,1,6}(SVector(T(-0.9591958116075092),
                                                T(0),
                                                T(1.2687130518352143),
                                                T(-0.3860731957925361),
                                                T(0.08345276954506707),
                                                T(-0.006896813980236088) )),
        # d3
        DerivativeCoefficientRow{T,1,6}(SVector(T(0.23226893958753303),
                                                T(-0.9178972823223088),
                                                T(0),
                                                T(0.8585728949642006),
                                                T(-0.19840158580452047),
                                                T(0.0254570335750956) )),
        # d4
        DerivativeCoefficientRow{T,1,7}(SVector(T(-0.07693308450770664),
                                                T(0.2682207441654631),
                                                T(-0.8244592785032373),
                                                T(0),
                                                T(0.7749305933808132),
                                                T(-0.15857634324473752),
                                                T(0.016817368709408833) )),
        # d5
        DerivativeCoefficientRow{T,1,8}(SVector(T(0.017248604775379643),
                                                T(-0.057293724956320936),
                                                T(0.1882698460378496),
                                                T(-0.7657841952075574),
                                                T(0),
                                                T(0.7505104744873311),
                                                T(-0.14956988077877112),
                                                T(0.01661887564208568) )),
        # d6
        DerivativeCoefficientRow{T,1,9}(SVector(T(-0.0016761725437960326),
                                                T(0.004750928272020083),
                                                T(-0.02423857479897864),
                                                T(0.15723370478311835),
                                                T(-0.753044107168037),
                                                T(0),
                                                T(0.7503740531217644),
                                                T(-0.15007481062435288),
                                                T(0.016674978958261432) )),
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(3//4), T(-3//20), T(1//60))
    central_coef = zero(T)
    lower_coef = -upper_coef
    left_weights = SVector( T(0.13030223027124),
                            T(0.68851501587715),
                            T(0.95166202564389),
                            T(0.99103890475697),
                            T(1.0028757074552),
                            T(0.99950151111941) )
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
        DerivativeCoefficientRow{T,1,8}(SVector(T(-4.647545021331372),
                                                T(6.254178662563343),
                                                T(-2.4139100510438665),
                                                T(1.256639509571857),
                                                T(-0.647667696659897),
                                                T(0.24570663393330786),
                                                T(-0.05204535767697368),
                                                T(0.004643320643692943) )),
        # d2
        DerivativeCoefficientRow{T,1,8}(SVector(T(-1.0868211647678514),
                                                T(0),
                                                T(1.5195364183857292),
                                                T(-0.6543667982790691),
                                                T(0.3128620695707949),
                                                T(-0.11086808101811908),
                                                T(0.021234895991755504),
                                                T(-0.0015773398832359302) )),
        # d3
        DerivativeCoefficientRow{T,1,8}(SVector(T(0.26780878369504263),
                                                T(-0.9701235413358851),
                                                T(0),
                                                T(0.9726229663976735),
                                                T(-0.3684397832422662),
                                                T(0.1161859354443521),
                                                T(-0.018907619882971544),
                                                T(0.0008532589240538837) )),
        # d4
        DerivativeCoefficientRow{T,1,8}(SVector(T(-0.12264219453174283),
                                                T(0.3675044946198741),
                                                T(-0.8555984502334516),
                                                T(0),
                                                T(0.7955264157053945),
                                                T(-0.22405034810099697),
                                                T(0.04290097854980221),
                                                T(-0.0036408960088885555) )),
        # d5
        DerivativeCoefficientRow{T,1,9}(SVector(T(0.0680142412945559),
                                                T(-0.18906569319897437),
                                                T(0.34874693977149884),
                                                T(-0.855998538940948),
                                                T(0),
                                                T(0.7918634438230676),
                                                T(-0.19783080393205543),
                                                T(0.03775653807587011),
                                                T(-0.0034861268930073425)  )),
        # d6
        DerivativeCoefficientRow{T,1,10}(SVector(T(-0.026557903380662235),
                                                T(0.06895974248171861),
                                                T(-0.11319486586973267),
                                                T(0.24813802673084387),
                                                T(-0.8150412273910146),
                                                T(0),
                                                T(0.8048396208993126),
                                                T(-0.20182899393713088),
                                                T(0.038273766032179475),
                                                T(-0.0035881655655168633) )),
        # d7
        DerivativeCoefficientRow{T,1,11}(SVector(T(0.005594622888681003),
                                                T(-0.013135637494343391),
                                                T(0.0183198514986364),
                                                T(-0.047252720307139756),
                                                T(0.20250470348441063),
                                                T(-0.8004261411357496),
                                                T(0),
                                                T(0.7997353219580023),
                                                T(-0.19983539676215634),
                                                T(0.03806388509755359),
                                                T(-0.0035684892278956487) )),
        # d8
        DerivativeCoefficientRow{T,1,12}(SVector(T(-0.000499585189812272),
                                                T(0.0009766036909553027),
                                                T(-0.000827480914165871),
                                                T(0.0040138390353474905),
                                                T(-0.038683470124513085),
                                                T(0.20090351330033354),
                                                T(-0.8004576105033867),
                                                T(0),
                                                T(0.8000635199885898),
                                                T(-0.20001587999714746),
                                                T(0.038098262856599514),
                                                T(-0.0035717121428062043) )),
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(4//5), T(-1//5), T(4//105), T(-1//280))
    central_coef = zero(T)
    lower_coef = -upper_coef
    left_weights = SVector( T(0.1075836807831),
                            T(0.61909685107891),
                            T(0.96971176519117),
                            T(1.1023441350947),
                            T(1.0244688965833),
                            T(0.99533550116831),
                            T(1.0008236941028),
                            T(0.99992060631812) )
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
        DerivativeCoefficientRow{T,1,10}(SVector(T(-55),
                                                T(6.7548747038001995),
                                                T(-2.6691978151545994),
                                                T(1.4438714982129999),
                                                T(-0.7727367375076),
                                                T(0.25570078343005),
                                                T(0.042808774693299),
                                                T(-0.082902108933389),
                                                T(0.032031176427907995),
                                                T(-0.0044502749689555995) )),
        # d2
        DerivativeCoefficientRow{T,1,10}(SVector(T(-1.145265719198745),
                                                T(0),
                                                T(1.6131685230292827),
                                                T(-0.7195954106367713),
                                                T(0.3651840332042441),
                                                T(-0.12197141381091767),
                                                T(-0.01399423475053903),
                                                T(0.0337917497680833),
                                                T(-0.013128541778312972),
                                                T(0.0018110141736784769) )),
        # d3
        DerivativeCoefficientRow{T,1,10}(SVector(T(0.27900977459917853),
                                                T(-0.9945564383180177),
                                                T(0),
                                                T(1.004253824704819),
                                                T(-0.4116222831605484),
                                                T(0.13904608960593332),
                                                T(8.487839079429473e-5),
                                                T(-0.024930033516808243),
                                                T(0.010097590982069497),
                                                T(-0.0013834032874125409) )),
        # d4
        DerivativeCoefficientRow{T,1,10}(SVector(T(-0.12555079634350264),
                                                T(0.36905436758383703),
                                                T(-0.8354024892044419),
                                                T(0),
                                                T(0.7960759131488463),
                                                T(-0.2481823909255492),
                                                T(0.035998033476551373),
                                                T(0.01500488078728063),
                                                T(-0.008204065636465682),
                                                T(0.0012065471134400296) )),
        # d5
        DerivativeCoefficientRow{T,1,10}(SVector(T(0.06879174149957458),
                                                T(-0.1917465511011161),
                                                T(0.3505621014987283),
                                                T(-0.8150200626888124),
                                                T(0),
                                                T(0.7435190828161148),
                                                T(-0.18326819871996888),
                                                T(0.027802276682660717),
                                                T(-0.00018667744812003666),
                                                T(-0.0004537125390638899) )),
        # d6
        DerivativeCoefficientRow{T,1,11}(SVector(T(-0.02525933759067232),
                                                T(0.07106552895046016),
                                                T(-0.1314044342558119),
                                                T(0.28194859492566227),
                                                T(-0.8250443110814595),
                                                T(0),
                                                T(0.8006147214131956),
                                                T(-0.2161330773414056),
                                                T(0.052333431144975136),
                                                T(-0.008905122105668357),
                                                T(0.0007840059407332415) )),
        # d7
        DerivativeCoefficientRow{T,1,12}(SVector(T(-0.004286144167448658),
                                                T(0.008264073453404727),
                                                T(-8.130035783403134e-5),
                                                T(-0.04144974144338818),
                                                T(0.2061182017256204),
                                                T(-0.8114609971447189),
                                                T(0),
                                                T(0.8288973712486651),
                                                T(-0.23611473523103743),
                                                T(0.05925148483279516),
                                                T(-0.009932840126144218),
                                                T(0.0007946272100915355) )),
        # d8
        DerivativeCoefficientRow{T,1,13}(SVector(T(0.008289486953575544),
                                                T(-0.01992892123103868),
                                                T(0.023847687855928824),
                                                T(-0.01725455228860255),
                                                T(-0.031227534064274177),
                                                T(0.2187728437910129),
                                                T(-0.8278065503298393),
                                                T(0),
                                                T(0.8329375462609859),
                                                T(-0.23802243145944782),
                                                T(0.05951861162798522),
                                                T(-0.009919768604664269),
                                                T(0.0007935814883731396) )),
        # d9
        DerivativeCoefficientRow{T,1,14}(SVector(T(-0.003203103055575049),
                                                T(0.007743290435329053),
                                                T(-0.009660000293183757),
                                                T(0.00943487651309554),
                                                T(0.00020969357970332722),
                                                T(-0.05297699654295238),
                                                T(0.23582352986415386),
                                                T(-0.8330064950072007),
                                                T(0),
                                                T(0.8333327624136899),
                                                T(-0.23809415379332086),
                                                T(0.059523538448329215),
                                                T(-0.009920589741388267),
                                                T(0.0007936471793110595) )),
        # d10
        DerivativeCoefficientRow{T,1,15}(SVector(T(0.0004450281134593904),
                                                T(-0.0010681530559586649),
                                                T(0.0013234597796542534),
                                                T(-0.0013875669869698303),
                                                T(0.0005096554110994863),
                                                T(0.0090146855416246),
                                                T(-0.0591787601986842),
                                                T(0.23804354829738644),
                                                T(-0.8333377120321315),
                                                T(0),
                                                T(0.8333344878759047),
                                                T(-0.23809556796454703),
                                                T(0.05952389199113576),
                                                T(-0.009920648665189359),
                                                T(0.0007936518932151467) )),
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(5//6), T(-5//21), T(5//84), T(-5//504), T(1//1260))
    central_coef = zero(T)
    lower_coef = -upper_coef
  else
    throw(ArgumentError("Order $order not implemented/derived."))
  end
end
