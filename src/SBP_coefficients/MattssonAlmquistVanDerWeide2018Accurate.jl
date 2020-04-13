
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
  # elseif order == 6
  #   left_boundary = (
  #       # d1
  #       DerivativeCoefficientRow{T,1,5}(SVector(T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T() )),
  #       # d2
  #       DerivativeCoefficientRow{T,1,5}(SVector(T(),
  #                                               T(0),
  #                                               T(),
  #                                               T(),
  #                                               T() )),
  #       # d3
  #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
  #                                               T(),
  #                                               T(0),
  #                                               T(),
  #                                               T(),
  #                                               T() )),
  #       # d4
  #       DerivativeCoefficientRow{T,1,7}(SVector(T(),
  #                                               T(),
  #                                               T(),
  #                                               T(0),
  #                                               T(),
  #                                               T(),
  #                                               T() )),
  #       # d5
  #       DerivativeCoefficientRow{T,1,8}(SVector(T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T(0),
  #                                               T(),
  #                                               T(),
  #                                               T() )),
  #   )
  #   right_boundary = .- left_boundary
  #   upper_coef = SVector(T(3//4), T(-3//20), T(1//60))
  #   central_coef = zero(T)
  #   lower_coef = -upper_coef
  #   left_weights = SVector( T(),
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
  #                                               T(0),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T() )),
  #       # d3
  #       DerivativeCoefficientRow{T,1,7}(SVector(T(),
  #                                               T(),
  #                                               T(0),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T() )),
  #       # d4
  #       DerivativeCoefficientRow{T,1,8}(SVector(T(),
  #                                               T(),
  #                                               T(),
  #                                               T(0),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T() )),
  #       # d5
  #       DerivativeCoefficientRow{T,1,9}(SVector(T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T(0),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T()  )),
  #       # d6
  #       DerivativeCoefficientRow{T,1,10}(SVector(T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T(0),
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
