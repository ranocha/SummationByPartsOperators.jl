
"""
    MattssonAlmquistVanDerWeide2018Accurate

Coefficients of the optimized SBP operators with nonuniform grid given in
  Mattsson, Almquist, van der Weide (2018)
  Boundary optimized diagonal-norm SBP operators.
  Journal of Computational Physics 374, pp. 1261-1266.
"""
struct MattssonAlmquistVanDerWeide2018Accurate <: SourceOfCoefficients end

function Base.show(io::IO, ::MattssonAlmquistVanDerWeide2018Accurate)
  print(io,
      "  Mattsson, Almquist, van der Weide (2018) \n",
      "  Boundary optimized diagonal-norm SBP operators. \n",
      "  Journal of Computational Physics 374, pp. 1261-1266. \n")
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
  else
    throw(ArgumentError("Accuracy order $accuracy_order not implemented/derived."))
  end

  BoundaryAdaptedGrid(xmin, xmax, xstart, N)
end


function first_derivative_coefficients(source::MattssonAlmquistVanDerWeide2018Accurate, order::Int, T=Float64, parallel=Val{:serial}())
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
                            left_weights, right_weights, parallel, 1, order, source)
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
                                                T(-0.19840158580452047),
                                                T(0.0254570335750956) )),
        # d4
        DerivativeCoefficientRow{T,1,7}(SVector(T(-0.07693308450770664),
                                                T(0.2682207441654631),
                                                T(0.1905185095473035),
                                                T(0),
                                                T(0.7749305933808132),
                                                T(-0.15857634324473752),
                                                T(0.016817368709408833) )),
        # d5
        DerivativeCoefficientRow{T,1,8}(SVector(T(0.017248604775379643),
                                                T(-0.057293724956320936),
                                                T(-0.024157023606080553),
                                                T(-0.7657841952075574),
                                                T(0),
                                                T(0.7505104744873311),
                                                T(-0.14956988077877112),
                                                T(0.01661887564208568) )),
        # d6
        DerivativeCoefficientRow{T,1,9}(SVector(T(0),
                                                T(0),
                                                T(-0.016674978958261432),
                                                T(0.15007481062435288),
                                                T(-0.7503740531217644),
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
                            left_weights, right_weights, parallel, 1, order, source)
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
                            left_weights, right_weights, parallel, 1, order, source)
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
  #                           left_weights, right_weights, parallel, 1, order, source)
  else
    throw(ArgumentError("Order $order not implemented/derived."))
  end
end
