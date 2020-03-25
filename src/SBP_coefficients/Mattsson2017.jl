
"""
    Mattsson2017

Coefficients of the upwind SBP operators given in
  Mattsson (2017)
  Diagonal-norm upwind SBP operators.
  Journal of Computational Physics 335, pp.283-310.
"""
struct Mattsson2017 <: SourceOfCoefficients
  kind::Symbol

  function Mattsson2017(kind::Symbol)
    if (kind !== :plus) && (kind !== :minus) && (kind !== :central)
      throw(ArgumentError("The only choices are :plus, :minus, and :central, not :$kind."))
    end
    new(kind)
  end
end

function Base.show(io::IO, source::Mattsson2017)
print(io,
    "  Upwind coefficients (", source.kind, ") of \n",
    "  Mattsson (2017) \n",
    "  Diagonal-norm upwind SBP operators. \n",
    "  Journal of Computational Physics 335, pp.283-310. \n")
end


function first_derivative_coefficients(source::Mattsson2017, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary_plus = (
            DerivativeCoefficientRow{T,1,3}(SVector(T(-3),
                                                    T(5),
                                                    T(-2), )),
            DerivativeCoefficientRow{T,1,4}(SVector(T(-1//5),
                                                    T(-1),
                                                    T(8//5),
                                                    T(-2//5), )),
        )
        right_boundary_plus = (
            DerivativeCoefficientRow{T,1,2}(SVector(T(1),
                                                    T(-1), )),
            DerivativeCoefficientRow{T,1,2}(SVector(T(1),
                                                    T(-1), )),
        )
        upper_coef_plus = SVector( T(2),
                                   T(-1//2), )
        central_coef_plus = T(-3//2)
        lower_coef_plus = SVector{0,T}()
        left_weights = SVector( T(1//4),
                                T(5//4), )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        left_boundary_minus = (
            DerivativeCoefficientRow{T,1,2}(SVector(T(-1),
                                                    T(1), )),
            DerivativeCoefficientRow{T,1,2}(SVector(T(-1),
                                                    T(1), )),
        )
        right_boundary_minus = (
            DerivativeCoefficientRow{T,1,3}(SVector(T(3),
                                                    T(-5),
                                                    T(2), )),
            DerivativeCoefficientRow{T,1,4}(SVector(T(1//5),
                                                    T(1),
                                                    T(-8//5),
                                                    T(2//5), )),
        )
        upper_coef_minus     = .- lower_coef_plus
        central_coef_minus   = .- central_coef_plus
        lower_coef_minus     = .- upper_coef_plus

        left_boundary_central  = (left_boundary_plus  .+ left_boundary_minus)  ./ 2
        right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
        upper_coef_central     = widening_plus(upper_coef_plus, upper_coef_minus) / 2
        central_coef_central   = (central_coef_plus   + central_coef_minus) / 2
        lower_coef_central     = widening_plus(lower_coef_plus, lower_coef_minus) / 2

        if source.kind === :plus
          left_boundary  = left_boundary_plus
          right_boundary = right_boundary_plus
          upper_coef     = upper_coef_plus
          central_coef   = central_coef_plus
          lower_coef     = lower_coef_plus
        elseif source.kind === :minus
          left_boundary  = left_boundary_minus
          right_boundary = right_boundary_minus
          upper_coef     = upper_coef_minus
          central_coef   = central_coef_minus
          lower_coef     = lower_coef_minus
        elseif source.kind === :central
          left_boundary  = left_boundary_central
          right_boundary = right_boundary_central
          upper_coef     = upper_coef_central
          central_coef   = central_coef_central
          lower_coef     = lower_coef_central
        end
        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 1, order, source)
    elseif order == 3
      left_boundary_plus = (
          DerivativeCoefficientRow{T,1,3}(SVector(T(-7//5),
                                                  T(9//5),
                                                  T(-2//5), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(-5//13),
                                                  T(-5//13),
                                                  T(12//13),
                                                  T(-2//13), )),
      )
      right_boundary_plus = (
          DerivativeCoefficientRow{T,1,2}(SVector(T(1),
                                                  T(-1), )),
          DerivativeCoefficientRow{T,1,3}(SVector(T(9//13),
                                                  T(-5//13),
                                                  T(-4//13), )),
      )
      upper_coef_plus = SVector( T(1),
                                 T(-1//6), )
      central_coef_plus = T(-1//2)
      lower_coef_plus = SVector( T(-1//3), )
      left_weights = SVector( T(5//12),
                              T(13//12), )
      right_weights = left_weights
      left_boundary_derivatives = Tuple{}()
      right_boundary_derivatives = left_boundary_derivatives

      left_boundary_minus = (
          DerivativeCoefficientRow{T,1,2}(SVector(T(-1),
                                                  T(1), )),
          DerivativeCoefficientRow{T,1,3}(SVector(T(-9//13),
                                                  T(5//13),
                                                  T(4//13), )),
      )
      right_boundary_minus = (
          DerivativeCoefficientRow{T,1,3}(SVector(T(7//5),
                                                  T(-9//5),
                                                  T(2//5), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(5//13),
                                                  T(5//13),
                                                  T(-12//13),
                                                  T(2//13), )),
      )
      upper_coef_minus     = .- lower_coef_plus
      central_coef_minus   = .- central_coef_plus
      lower_coef_minus     = .- upper_coef_plus

      left_boundary_central  = (left_boundary_plus  .+ left_boundary_minus)  ./ 2
      right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
      upper_coef_central     = widening_plus(upper_coef_plus, upper_coef_minus) / 2
      central_coef_central   = (central_coef_plus   + central_coef_minus) / 2
      lower_coef_central     = widening_plus(lower_coef_plus, lower_coef_minus) / 2

      if source.kind === :plus
        left_boundary  = left_boundary_plus
        right_boundary = right_boundary_plus
        upper_coef     = upper_coef_plus
        central_coef   = central_coef_plus
        lower_coef     = lower_coef_plus
      elseif source.kind === :minus
        left_boundary  = left_boundary_minus
        right_boundary = right_boundary_minus
        upper_coef     = upper_coef_minus
        central_coef   = central_coef_minus
        lower_coef     = lower_coef_minus
      elseif source.kind === :central
        left_boundary  = left_boundary_central
        right_boundary = right_boundary_central
        upper_coef     = upper_coef_central
        central_coef   = central_coef_central
        lower_coef     = lower_coef_central
      end
      DerivativeCoefficients(left_boundary, right_boundary,
                              left_boundary_derivatives, right_boundary_derivatives,
                              lower_coef, central_coef, upper_coef,
                              left_weights, right_weights, parallel, 1, order, source)
    elseif order == 4
      left_boundary_plus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(-75//49),
                                                  T(205//98),
                                                  T(-29//49),
                                                  T(3//98), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(-169//366),
                                                  T(-11//61),
                                                  T(99//122),
                                                  T(-43//183),
                                                  T(4//61), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(11//123),
                                                  T(-39//82),
                                                  T(-29//41),
                                                  T(389//246),
                                                  T(-24//41),
                                                  T(4//41), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(9//298),
                                                  T(-11//149),
                                                  T(-65//298),
                                                  T(-117//149),
                                                  T(216//149),
                                                  T(-72//149),
                                                  T(12//149), )),
      )
      right_boundary_plus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(69//49),
                                                  T(-169//98),
                                                  T(11//49),
                                                  T(9//98), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(205//366),
                                                  T(-11//61),
                                                  T(-39//122),
                                                  T(-11//183), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(-29//123),
                                                  T(99//82),
                                                  T(-29//41),
                                                  T(-65//246), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(3//298),
                                                  T(-43//149),
                                                  T(389//298),
                                                  T(-117//149),
                                                  T(-36//149), )),
      )
      upper_coef_plus = SVector( T(3//2),
                                 T(-1//2),
                                 T(1//12), )
      central_coef_plus = T(-5//6)
      lower_coef_plus = SVector( T(-1//4), )
      left_weights = SVector( T(49//144),
                              T(61//48),
                              T(41//48),
                              T(149//144), )
      right_weights = left_weights
      left_boundary_derivatives = Tuple{}()
      right_boundary_derivatives = left_boundary_derivatives

      left_boundary_minus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(-69//49),
                                                  T(169//98),
                                                  T(-11//49),
                                                  T(-9//98), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(-205//366),
                                                  T(11//61),
                                                  T(39//122),
                                                  T(11//183), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(29//123),
                                                  T(-99//82),
                                                  T(29//41),
                                                  T(65//246), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(-3//298),
                                                  T(43//149),
                                                  T(-389//298),
                                                  T(117//149),
                                                  T(36//149), )),
      )
      right_boundary_minus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(75//49),
                                                  T(-205//98),
                                                  T(29//49),
                                                  T(-3//98), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(169//366),
                                                  T(11//61),
                                                  T(-99//122),
                                                  T(43//183),
                                                  T(-4//61), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-11//123),
                                                  T(39//82),
                                                  T(29//41),
                                                  T(-389//246),
                                                  T(24//41),
                                                  T(-4//41), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(-9//298),
                                                  T(11//149),
                                                  T(65//298),
                                                  T(117//149),
                                                  T(-216//149),
                                                  T(72//149),
                                                  T(-12//149), )),
      )
      upper_coef_minus     = .- lower_coef_plus
      central_coef_minus   = .- central_coef_plus
      lower_coef_minus     = .- upper_coef_plus

      left_boundary_central  = (left_boundary_plus  .+ left_boundary_minus)  ./ 2
      right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
      upper_coef_central     = widening_plus(upper_coef_plus, upper_coef_minus) / 2
      central_coef_central   = (central_coef_plus   + central_coef_minus) / 2
      lower_coef_central     = widening_plus(lower_coef_plus, lower_coef_minus) / 2

      if source.kind === :plus
        left_boundary  = left_boundary_plus
        right_boundary = right_boundary_plus
        upper_coef     = upper_coef_plus
        central_coef   = central_coef_plus
        lower_coef     = lower_coef_plus
      elseif source.kind === :minus
        left_boundary  = left_boundary_minus
        right_boundary = right_boundary_minus
        upper_coef     = upper_coef_minus
        central_coef   = central_coef_minus
        lower_coef     = lower_coef_minus
      elseif source.kind === :central
        left_boundary  = left_boundary_central
        right_boundary = right_boundary_central
        upper_coef     = upper_coef_central
        central_coef   = central_coef_central
        lower_coef     = lower_coef_central
      end
      DerivativeCoefficients(left_boundary, right_boundary,
                              left_boundary_derivatives, right_boundary_derivatives,
                              lower_coef, central_coef, upper_coef,
                              left_weights, right_weights, parallel, 1, order, source)
    elseif order == 5
      left_boundary_plus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(-366//251),
                                                  T(941//502),
                                                  T(-94//251),
                                                  T(-21//502), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(-869//1794),
                                                  T(-22//299),
                                                  T(375//598),
                                                  T(-86//897),
                                                  T(8//299), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(58//633),
                                                  T(-255//422),
                                                  T(-58//211),
                                                  T(1309//1266),
                                                  T(-60//211),
                                                  T(8//211), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(45//1478),
                                                  T(-22//739),
                                                  T(-661//1478),
                                                  T(-234//739),
                                                  T(720//739),
                                                  T(-180//739),
                                                  T(24//739), )),
      )
      right_boundary_plus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(354//251),
                                                  T(-869//502),
                                                  T(58//251),
                                                  T(45//502), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(941//1794),
                                                  T(-22//299),
                                                  T(-255//598),
                                                  T(-22//897), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(-94//633),
                                                  T(375//422),
                                                  T(-58//211),
                                                  T(-661//1266),
                                                  T(12//211), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-21//1478),
                                                  T(-86//739),
                                                  T(1309//1478),
                                                  T(-234//739),
                                                  T(-360//739),
                                                  T(36//739), )),
      )
      upper_coef_plus = SVector( T(1),
                                 T(-1//4),
                                 T(1//30), )
      central_coef_plus = T(-1//3)
      lower_coef_plus = SVector( T(-1//2),
                                 T(1//20), )
      left_weights = SVector( T(251//720),
                              T(299//240),
                              T(211//240),
                              T(739//720), )
      right_weights = left_weights
      left_boundary_derivatives = Tuple{}()
      right_boundary_derivatives = left_boundary_derivatives

      left_boundary_minus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(-354//251),
                                                  T(869//502),
                                                  T(-58//251),
                                                  T(-45//502), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(-941//1794),
                                                  T(22//299),
                                                  T(255//598),
                                                  T(22//897), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(94//633),
                                                  T(-375//422),
                                                  T(58//211),
                                                  T(661//1266),
                                                  T(-12//211), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(21//1478),
                                                  T(86//739),
                                                  T(-1309//1478),
                                                  T(234//739),
                                                  T(360//739),
                                                  T(-36//739), )),
      )
      right_boundary_minus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(366//251),
                                                  T(-941//502),
                                                  T(94//251),
                                                  T(21//502), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(869//1794),
                                                  T(22//299),
                                                  T(-375//598),
                                                  T(86//897),
                                                  T(-8//299), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-58//633),
                                                  T(255//422),
                                                  T(58//211),
                                                  T(-1309//1266),
                                                  T(60//211),
                                                  T(-8//211), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(-45//1478),
                                                  T(22//739),
                                                  T(661//1478),
                                                  T(234//739),
                                                  T(-720//739),
                                                  T(180//739),
                                                  T(-24//739), )),
      )
      upper_coef_minus     = .- lower_coef_plus
      central_coef_minus   = .- central_coef_plus
      lower_coef_minus     = .- upper_coef_plus

      left_boundary_central  = (left_boundary_plus  .+ left_boundary_minus)  ./ 2
      right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
      upper_coef_central     = widening_plus(upper_coef_plus, upper_coef_minus) / 2
      central_coef_central   = (central_coef_plus   + central_coef_minus) / 2
      lower_coef_central     = widening_plus(lower_coef_plus, lower_coef_minus) / 2

      if source.kind === :plus
        left_boundary  = left_boundary_plus
        right_boundary = right_boundary_plus
        upper_coef     = upper_coef_plus
        central_coef   = central_coef_plus
        lower_coef     = lower_coef_plus
      elseif source.kind === :minus
        left_boundary  = left_boundary_minus
        right_boundary = right_boundary_minus
        upper_coef     = upper_coef_minus
        central_coef   = central_coef_minus
        lower_coef     = lower_coef_minus
      elseif source.kind === :central
        left_boundary  = left_boundary_central
        right_boundary = right_boundary_central
        upper_coef     = upper_coef_central
        central_coef   = central_coef_central
        lower_coef     = lower_coef_central
      end
      DerivativeCoefficients(left_boundary, right_boundary,
                              left_boundary_derivatives, right_boundary_derivatives,
                              lower_coef, central_coef, upper_coef,
                              left_weights, right_weights, parallel, 1, order, source)
    # elseif order == 4
    #   left_boundary_plus = (
    #       DerivativeCoefficientRow{T,1,4}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,5}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,7}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #   )
    #   right_boundary_plus = (
    #       DerivativeCoefficientRow{T,1,4}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,4}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,4}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,5}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #   )
    #   upper_coef_plus = SVector( T(),
    #                              T(),
    #                              T(), )
    #   central_coef_plus = T()
    #   lower_coef_plus = SVector( T(), )
    #   left_weights = SVector( T(),
    #                           T(),
    #                           T(),
    #                           T(), )
    #   right_weights = left_weights
    #   left_boundary_derivatives = Tuple{}()
    #   right_boundary_derivatives = left_boundary_derivatives

    #   left_boundary_minus = (
    #       DerivativeCoefficientRow{T,1,4}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,4}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,4}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,5}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #   )
    #   right_boundary_minus = (
    #       DerivativeCoefficientRow{T,1,4}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,5}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,7}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #   )
    #   upper_coef_minus     = .- lower_coef_plus
    #   central_coef_minus   = .- central_coef_plus
    #   lower_coef_minus     = .- upper_coef_plus

    #   left_boundary_central  = (left_boundary_plus  .+ left_boundary_minus)  ./ 2
    #   right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
    #   upper_coef_central     = widening_plus(upper_coef_plus, upper_coef_minus) / 2
    #   central_coef_central   = (central_coef_plus   + central_coef_minus) / 2
    #   lower_coef_central     = widening_plus(lower_coef_plus, lower_coef_minus) / 2

    #   if source.kind === :plus
    #     left_boundary  = left_boundary_plus
    #     right_boundary = right_boundary_plus
    #     upper_coef     = upper_coef_plus
    #     central_coef   = central_coef_plus
    #     lower_coef     = lower_coef_plus
    #   elseif source.kind === :minus
    #     left_boundary  = left_boundary_minus
    #     right_boundary = right_boundary_minus
    #     upper_coef     = upper_coef_minus
    #     central_coef   = central_coef_minus
    #     lower_coef     = lower_coef_minus
    #   elseif source.kind === :central
    #     left_boundary  = left_boundary_central
    #     right_boundary = right_boundary_central
    #     upper_coef     = upper_coef_central
    #     central_coef   = central_coef_central
    #     lower_coef     = lower_coef_central
    #   end
    #   DerivativeCoefficients(left_boundary, right_boundary,
    #                           left_boundary_derivatives, right_boundary_derivatives,
    #                           lower_coef, central_coef, upper_coef,
    #                           left_weights, right_weights, parallel, 1, order, source)
  else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end


# function second_derivative_coefficients(source::Mattsson2017, order::Int, T=Float64, parallel=Val{:serial}())
#     if order == 2
#         left_boundary = (
#             # d1
#             DerivativeCoefficientRow{T,1,3}(SVector(T(1),
#                                                     T(-2),
#                                                     T(1) )),
#         )
#         right_boundary = left_boundary
#         upper_coef = SVector(T(1))
#         central_coef = T(-2)
#         lower_coef = upper_coef
#         left_weights = SVector(T(1//2))
#         right_weights = left_weights
#         left_boundary_derivatives = (
#             DerivativeCoefficientRow{T,1,3}(SVector(T(-3//2),
#                                                     T(2),
#                                                     T(-1//2) )),
#         )
#         right_boundary_derivatives = (-left_boundary_derivatives[1], )

#         DerivativeCoefficients(left_boundary, right_boundary,
#                                 left_boundary_derivatives, right_boundary_derivatives,
#                                 lower_coef, central_coef, upper_coef,
#                                 left_weights, right_weights, parallel, 2, order, source)
#     else
#         throw(ArgumentError("Order $order not implemented/derived."))
#     end
# end
