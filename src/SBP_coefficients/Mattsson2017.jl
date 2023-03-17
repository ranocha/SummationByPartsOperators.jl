
"""
    Mattsson2017(version::Symbol)

Coefficients of the upwind SBP operators given in
- Mattsson (2017)
  Diagonal-norm upwind SBP operators.
  Journal of Computational Physics 335, pp. 283-310.

You can choose between the different versions `:central`, `:plus`, and `:minus`.
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
  if get(io, :compact, false)
    summary(io, source)
  else
    print(io,
        "Mattsson (2017) \n",
        "  Diagonal-norm upwind SBP operators. \n",
        "  Journal of Computational Physics 335, pp. 283-310. \n",
        "  (upwind coefficients ", source.kind, ")")
  end
end


function first_derivative_coefficients(source::Mattsson2017, order::Int,
                                       T=Float64, mode=FastMode())
    if order == 1
        left_boundary_plus = (
            DerivativeCoefficientRow{T,1,2}(SVector(T(-1),    # -2 for continuously coupled = same
                                                    T(1), )), #  2 for continuously coupled = same
        )
        right_boundary_plus = (
            DerivativeCoefficientRow{T,1,2}(SVector(T(0),
                                                    T(0), )),
        )
        upper_coef_plus = SVector( T(1), )
        central_coef_plus = T(-1)
        lower_coef_plus = SVector{0,T}()
        left_weights = SVector( T(1), ) # 1//2 for continuously coupled = same
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        left_boundary_minus = (
            DerivativeCoefficientRow{T,1,2}(SVector(T(0),
                                                    T(0), )),
        )
        right_boundary_minus = (
            DerivativeCoefficientRow{T,1,2}(SVector(T(1),      #  2 for continuously coupled = same
                                                    T(-1), )), # -2 for continuously coupled = same
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
                               left_weights, right_weights, mode, 1, order, source)
    elseif order == 2
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
                                left_weights, right_weights, mode, 1, order, source)
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
                              left_weights, right_weights, mode, 1, order, source)
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
                              left_weights, right_weights, mode, 1, order, source)
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
                              left_weights, right_weights, mode, 1, order, source)
    elseif order == 6
      left_boundary_plus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(-58148100//36496453),
                                                  T(1146190567//547446795),
                                                  T(-14369571//52137790),
                                                  T(-55265831//182482265),
                                                  T(26269819//1094893590),
                                                  T(9858004//182482265), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-1116490567//2422752675),
                                                  T(-954612//32303369),
                                                  T(190538869//484550535),
                                                  T(102705469//969101070),
                                                  T(4964892//161516845),
                                                  T(-191689861//4845505350), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(9869571//102452500),
                                                  T(-135385429//215150250),
                                                  T(-2198412//7171675),
                                                  T(45137333//35858375),
                                                  T(-253641811//430300500),
                                                  T(70665929//358583750),
                                                  T(-72//2675), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(66965831//723199750),
                                                  T(-208765789//867839700),
                                                  T(-17623253//72319975),
                                                  T(-657684//2066285),
                                                  T(410905829//433919850),
                                                  T(-477953317//1446399500),
                                                  T(576//5395),
                                                  T(-72//5395), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(-49219819//3153258150),
                                                  T(3519588//105108605),
                                                  T(26422771//630651630),
                                                  T(-141938309//315325815),
                                                  T(-12476988//21021721),
                                                  T(2217185207//1576629075),
                                                  T(-4320//7841),
                                                  T(1152//7841),
                                                  T(-144//7841), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(-9498004//587634985),
                                                  T(142906261//3525809910),
                                                  T(-3137129//587634985),
                                                  T(-29884283//1175269970),
                                                  T(-630168407//1762904955),
                                                  T(-9609300//16789571),
                                                  T(57600//43837),
                                                  T(-21600//43837),
                                                  T(5760//43837),
                                                  T(-720//43837), )),
      )
      right_boundary_plus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(57671100//36496453),
                                                  T(-1116490567//547446795),
                                                  T(9869571//52137790),
                                                  T(66965831//182482265),
                                                  T(-49219819//1094893590),
                                                  T(-9498004//182482265), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(1146190567//2422752675),
                                                  T(-954612//32303369),
                                                  T(-135385429//484550535),
                                                  T(-208765789//969101070),
                                                  T(3519588//161516845),
                                                  T(142906261//4845505350), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-14369571//102452500),
                                                  T(190538869//215150250),
                                                  T(-2198412//7171675),
                                                  T(-17623253//35858375),
                                                  T(26422771//430300500),
                                                  T(-3137129//358583750), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-55265831//723199750),
                                                  T(102705469//867839700),
                                                  T(45137333//72319975),
                                                  T(-657684//2066285),
                                                  T(-141938309//433919850),
                                                  T(-2298791//111261500), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(26269819//3153258150),
                                                  T(4964892//105108605),
                                                  T(-253641811//630651630),
                                                  T(410905829//315325815),
                                                  T(-12476988//21021721),
                                                  T(-630168407//1576629075),
                                                  T(288//7841), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(9858004//587634985),
                                                  T(-191689861//3525809910),
                                                  T(70665929//587634985),
                                                  T(-477953317//1175269970),
                                                  T(2217185207//1762904955),
                                                  T(-9609300//16789571),
                                                  T(-17280//43837),
                                                  T(1440//43837), )),
      )
      upper_coef_plus = SVector( T(4//3),
                                 T(-1//2),
                                 T(2//15),
                                 T(-1//60), )
      central_coef_plus = T(-7//12)
      lower_coef_plus = SVector( T(-2//5),
                                 T(1//30), )
      left_weights = SVector( T(13613//43200),
                              T(12049//8640),
                              T(535//864),
                              T(1079//864),
                              T(7841//8640),
                              T(43837//43200), )
      right_weights = left_weights
      left_boundary_derivatives = Tuple{}()
      right_boundary_derivatives = left_boundary_derivatives

      left_boundary_minus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(-57671100//36496453),
                                                  T(1116490567//547446795),
                                                  T(-9869571//52137790),
                                                  T(-66965831//182482265),
                                                  T(49219819//1094893590),
                                                  T(9498004//182482265), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-1146190567//2422752675),
                                                  T(954612//32303369),
                                                  T(135385429//484550535),
                                                  T(208765789//969101070),
                                                  T(-3519588//161516845),
                                                  T(-142906261//4845505350), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(14369571//102452500),
                                                  T(-190538869//215150250),
                                                  T(2198412//7171675),
                                                  T(17623253//35858375),
                                                  T(-26422771//430300500),
                                                  T(3137129//358583750), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(55265831//723199750),
                                                  T(-102705469//867839700),
                                                  T(-45137333//72319975),
                                                  T(657684//2066285),
                                                  T(141938309//433919850),
                                                  T(2298791//111261500), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(-26269819//3153258150),
                                                  T(-4964892//105108605),
                                                  T(253641811//630651630),
                                                  T(-410905829//315325815),
                                                  T(12476988//21021721),
                                                  T(630168407//1576629075),
                                                  T(-288//7841), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-9858004//587634985),
                                                  T(191689861//3525809910),
                                                  T(-70665929//587634985),
                                                  T(477953317//1175269970),
                                                  T(-2217185207//1762904955),
                                                  T(9609300//16789571),
                                                  T(17280//43837),
                                                  T(-1440//43837), )),
      )
      right_boundary_minus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(58148100//36496453),
                                                  T(-1146190567//547446795),
                                                  T(14369571//52137790),
                                                  T(55265831//182482265),
                                                  T(-26269819//1094893590),
                                                  T(-9858004//182482265), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(1116490567//2422752675),
                                                  T(954612//32303369),
                                                  T(-190538869//484550535),
                                                  T(-102705469//969101070),
                                                  T(-4964892//161516845),
                                                  T(191689861//4845505350), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(-9869571//102452500),
                                                  T(135385429//215150250),
                                                  T(2198412//7171675),
                                                  T(-45137333//35858375),
                                                  T(253641811//430300500),
                                                  T(-70665929//358583750),
                                                  T(72//2675), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-66965831//723199750),
                                                  T(208765789//867839700),
                                                  T(17623253//72319975),
                                                  T(657684//2066285),
                                                  T(-410905829//433919850),
                                                  T(477953317//1446399500),
                                                  T(-576//5395),
                                                  T(72//5395), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(49219819//3153258150),
                                                  T(-3519588//105108605),
                                                  T(-26422771//630651630),
                                                  T(141938309//315325815),
                                                  T(12476988//21021721),
                                                  T(-2217185207//1576629075),
                                                  T(4320//7841),
                                                  T(-1152//7841),
                                                  T(144//7841), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(9498004//587634985),
                                                  T(-142906261//3525809910),
                                                  T(3137129//587634985),
                                                  T(29884283//1175269970),
                                                  T(630168407//1762904955),
                                                  T(9609300//16789571),
                                                  T(-57600//43837),
                                                  T(21600//43837),
                                                  T(-5760//43837),
                                                  T(720//43837), )),
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
                              left_weights, right_weights, mode, 1, order, source)
    elseif order == 7
      left_boundary_plus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(-81216540//51172247),
                                                  T(1587945773//767583705),
                                                  T(-17337249//73103210),
                                                  T(-84398989//255861235),
                                                  T(48781961//1535167410),
                                                  T(13716476//255861235), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-1570125773//3386062785),
                                                  T(-2863836//225737519),
                                                  T(240029831//677212557),
                                                  T(202934303//1354425114),
                                                  T(1418484//225737519),
                                                  T(-231357719//6772125570), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(14637249//144536540),
                                                  T(-206937767//303526734),
                                                  T(-6595236//50587789),
                                                  T(49602727//50587789),
                                                  T(-218919665//607053468),
                                                  T(51815011//505877890),
                                                  T(-216//18869), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(91418989//1008619010),
                                                  T(-266570495//1210342812),
                                                  T(-33094279//100861901),
                                                  T(-1973052//14408843),
                                                  T(440626231//605171406),
                                                  T(-365711063//2017238020),
                                                  T(2016//37621),
                                                  T(-216//37621), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(-62551961//4426143330),
                                                  T(9588//385217),
                                                  T(82588241//885228666),
                                                  T(-279245719//442614333),
                                                  T(-37430964//147538111),
                                                  T(2312302333//2213071665),
                                                  T(-18144//55031),
                                                  T(4032//55031),
                                                  T(-432//55031), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(-13500476//822302915),
                                                  T(202087559//4933817490),
                                                  T(-11297731//822302915),
                                                  T(61008503//1644605830),
                                                  T(-1360092253//2466908745),
                                                  T(-5765580//23494369),
                                                  T(60480//61343),
                                                  T(-18144//61343),
                                                  T(4032//61343),
                                                  T(-432//61343), )),
      )
      right_boundary_plus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(80930340//51172247),
                                                  T(-1570125773//767583705),
                                                  T(14637249//73103210),
                                                  T(91418989//255861235),
                                                  T(-62551961//1535167410),
                                                  T(-13500476//255861235), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(1587945773//3386062785),
                                                  T(-2863836//225737519),
                                                  T(-206937767//677212557),
                                                  T(-266570495//1354425114),
                                                  T(9588//589393),
                                                  T(202087559//6772125570), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-17337249//144536540),
                                                  T(240029831//303526734),
                                                  T(-6595236//50587789),
                                                  T(-33094279//50587789),
                                                  T(82588241//607053468),
                                                  T(-11297731//505877890), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(-84398989//1008619010),
                                                  T(202934303//1210342812),
                                                  T(49602727//100861901),
                                                  T(-1973052//14408843),
                                                  T(-279245719//605171406),
                                                  T(61008503//2017238020),
                                                  T(-288//37621), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(48781961//4426143330),
                                                  T(1418484//147538111),
                                                  T(-218919665//885228666),
                                                  T(440626231//442614333),
                                                  T(-37430964//147538111),
                                                  T(-1360092253//2213071665),
                                                  T(6048//55031),
                                                  T(-576//55031), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(13716476//822302915),
                                                  T(-231357719//4933817490),
                                                  T(51815011//822302915),
                                                  T(-365711063//1644605830),
                                                  T(2312302333//2466908745),
                                                  T(-5765580//23494369),
                                                  T(-36288//61343),
                                                  T(6048//61343),
                                                  T(-576//61343), )),
      )
      upper_coef_plus = SVector( T(1),
                                 T(-3//10),
                                 T(1//15),
                                 T(-1//140), )
      central_coef_plus = T(-1//4)
      lower_coef_plus = SVector( T(-3//5),
                                 T(1//10),
                                 T(-1//105), )
      left_weights = SVector( T(19087//60480),
                              T(84199//60480),
                              T(18869//30240),
                              T(37621//30240),
                              T(55031//60480),
                              T(61343//60480), )
      right_weights = left_weights
      left_boundary_derivatives = Tuple{}()
      right_boundary_derivatives = left_boundary_derivatives

      left_boundary_minus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(-80930340//51172247),
                                                  T(1570125773//767583705),
                                                  T(-14637249//73103210),
                                                  T(-91418989//255861235),
                                                  T(62551961//1535167410),
                                                  T(13500476//255861235), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-1587945773//3386062785),
                                                  T(2863836//225737519),
                                                  T(206937767//677212557),
                                                  T(266570495//1354425114),
                                                  T(-9588//589393),
                                                  T(-202087559//6772125570), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(17337249//144536540),
                                                  T(-240029831//303526734),
                                                  T(6595236//50587789),
                                                  T(33094279//50587789),
                                                  T(-82588241//607053468),
                                                  T(11297731//505877890), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(84398989//1008619010),
                                                  T(-202934303//1210342812),
                                                  T(-49602727//100861901),
                                                  T(1973052//14408843),
                                                  T(279245719//605171406),
                                                  T(-61008503//2017238020),
                                                  T(288//37621), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-48781961//4426143330),
                                                  T(-1418484//147538111),
                                                  T(218919665//885228666),
                                                  T(-440626231//442614333),
                                                  T(37430964//147538111),
                                                  T(1360092253//2213071665),
                                                  T(-6048//55031),
                                                  T(576//55031), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(-13716476//822302915),
                                                  T(231357719//4933817490),
                                                  T(-51815011//822302915),
                                                  T(365711063//1644605830),
                                                  T(-2312302333//2466908745),
                                                  T(5765580//23494369),
                                                  T(36288//61343),
                                                  T(-6048//61343),
                                                  T(576//61343), )),
      )
      right_boundary_minus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(81216540//51172247),
                                                  T(-1587945773//767583705),
                                                  T(17337249//73103210),
                                                  T(84398989//255861235),
                                                  T(-48781961//1535167410),
                                                  T(-13716476//255861235), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(1570125773//3386062785),
                                                  T(2863836//225737519),
                                                  T(-240029831//677212557),
                                                  T(-202934303//1354425114),
                                                  T(-1418484//225737519),
                                                  T(231357719//6772125570), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(-14637249//144536540),
                                                  T(206937767//303526734),
                                                  T(6595236//50587789),
                                                  T(-49602727//50587789),
                                                  T(218919665//607053468),
                                                  T(-51815011//505877890),
                                                  T(216//18869), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-91418989//1008619010),
                                                  T(266570495//1210342812),
                                                  T(33094279//100861901),
                                                  T(1973052//14408843),
                                                  T(-440626231//605171406),
                                                  T(365711063//2017238020),
                                                  T(-2016//37621),
                                                  T(216//37621), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(62551961//4426143330),
                                                  T(-9588//385217),
                                                  T(-82588241//885228666),
                                                  T(279245719//442614333),
                                                  T(37430964//147538111),
                                                  T(-2312302333//2213071665),
                                                  T(18144//55031),
                                                  T(-4032//55031),
                                                  T(432//55031), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(13500476//822302915),
                                                  T(-202087559//4933817490),
                                                  T(11297731//822302915),
                                                  T(-61008503//1644605830),
                                                  T(1360092253//2466908745),
                                                  T(5765580//23494369),
                                                  T(-60480//61343),
                                                  T(18144//61343),
                                                  T(-4032//61343),
                                                  T(432//61343), )),
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
                              left_weights, right_weights, mode, 1, order, source)
    elseif order == 8
      left_boundary_plus = (
          DerivativeCoefficientRow{T,1,8}(SVector(T(-2502530942940//1474909813267),
                                                  T(88037468909961//38347655144942),
                                                  T(-21877003412728//95869137862355),
                                                  T(-113480208109603//230085930869652),
                                                  T(-4151251151305//38347655144942),
                                                  T(5001038984066//19173827572471),
                                                  T(2235718279643//115042965434826),
                                                  T(-19035612535373//383476551449420), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-87682997519961//198485864362786),
                                                  T(-5599612620//1090581672323),
                                                  T(58119019845719//283551234803980),
                                                  T(3739408501537//14177561740199),
                                                  T(5368963068922//42532685220597),
                                                  T(-4450185662513//28355123480398),
                                                  T(-1221838279381//56710246960796),
                                                  T(90595000956023//2977287965441790), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(20045516300728//83364187761855),
                                                  T(-51190456749719//47636678721060),
                                                  T(-4159957380//20357555009),
                                                  T(10994933811709//4763667872106),
                                                  T(-9270952411151//4763667872106),
                                                  T(3191238635141//5292964302340),
                                                  T(4442211176987//23818339360530),
                                                  T(-1881322730062//16672837552371), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(118016946570403//1404599159063100),
                                                  T(-4173878828737//16721418560275),
                                                  T(-7990503962509//33442837120550),
                                                  T(-22442359068//257252593235),
                                                  T(6132064681023//13377134848220),
                                                  T(511197701761//33442837120550),
                                                  T(2475363434426//50164255680825),
                                                  T(-7784834666617//234099859843850),
                                                  T(2592//1306295), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(68609076271//971739785926),
                                                  T(-13508469862//34600143435),
                                                  T(6527681584751//7635098317990),
                                                  T(-3347940206463//3054039327196),
                                                  T(-41872007268//58731525523),
                                                  T(3208334350649//1527019663598),
                                                  T(-407569013461//347049923545),
                                                  T(136474842626653//320674129355580),
                                                  T(-25920//298231),
                                                  T(2592//298231), )),
          DerivativeCoefficientRow{T,1,11}(SVector(T(-4975275570026//83211286617459),
                                                  T(4244231077313//23774653319274),
                                                  T(-1550378843141//26416281465860),
                                                  T(-5726967564961//23774653319274),
                                                  T(-5089494709645//23774653319274),
                                                  T(-31599233340//101601082561),
                                                  T(45241297077547//47549306638548),
                                                  T(-145894938361408//416056433087295),
                                                  T(67200//515917),
                                                  T(-14400//515917),
                                                  T(1440//515917), )),
          DerivativeCoefficientRow{T,1,12}(SVector(T(-2164019088443//360119721423462),
                                                  T(1263196075861//34297116326044),
                                                  T(-6600697610987//85742790815110),
                                                  T(1113221183374//25722837244533),
                                                  T(926842346471//8574279081511),
                                                  T(-18757693936747//34297116326044),
                                                  T(-315773585460//659559929347),
                                                  T(5525449761123//4197199550390),
                                                  T(-1814400//3349159),
                                                  T(604800//3349159),
                                                  T(-129600//3349159),
                                                  T(12960//3349159), )),
          DerivativeCoefficientRow{T,1,13}(SVector(T(18967504495373//1312833690376780),
                                                  T(-90231551688023//1969250535565170),
                                                  T(2064808836502//65641684518839),
                                                  T(3502353445417//131283369037678),
                                                  T(-15385068876253//787700214226068),
                                                  T(8399970205408//328208422594195),
                                                  T(-23409126682353//50493603476030),
                                                  T(-2249575406820//5049360347603),
                                                  T(31752000//25639991),
                                                  T(-12700800//25639991),
                                                  T(4233600//25639991),
                                                  T(-907200//25639991),
                                                  T(90720//25639991), )),
      )
      right_boundary_plus = (
          DerivativeCoefficientRow{T,1,8}(SVector(T(2499882349860//1474909813267),
                                                  T(-87682997519961//38347655144942),
                                                  T(20045516300728//95869137862355),
                                                  T(118016946570403//230085930869652),
                                                  T(343045381355//3486150467722),
                                                  T(-4975275570026//19173827572471),
                                                  T(-2164019088443//115042965434826),
                                                  T(18967504495373//383476551449420), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(88037468909961//198485864362786),
                                                  T(-5599612620//1090581672323),
                                                  T(-51190456749719//283551234803980),
                                                  T(-4173878828737//14177561740199),
                                                  T(-4471303524322//42532685220597),
                                                  T(4244231077313//28355123480398),
                                                  T(1263196075861//56710246960796),
                                                  T(-90231551688023//2977287965441790), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-21877003412728//83364187761855),
                                                  T(58119019845719//47636678721060),
                                                  T(-4159957380//20357555009),
                                                  T(-7990503962509//4763667872106),
                                                  T(6527681584751//4763667872106),
                                                  T(-1550378843141//5292964302340),
                                                  T(-6600697610987//23818339360530),
                                                  T(2064808836502//16672837552371), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(-113480208109603//1404599159063100),
                                                  T(3739408501537//16721418560275),
                                                  T(10994933811709//33442837120550),
                                                  T(-22442359068//257252593235),
                                                  T(-3347940206463//13377134848220),
                                                  T(-5726967564961//33442837120550),
                                                  T(1113221183374//50164255680825),
                                                  T(3502353445417//234099859843850),
                                                  T(0//1), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(-830250230261//10689137645186),
                                                  T(5368963068922//11452647476985),
                                                  T(-9270952411151//7635098317990),
                                                  T(360709687119//179649372188),
                                                  T(-41872007268//58731525523),
                                                  T(-1017898941929//1527019663598),
                                                  T(926842346471//3817549158995),
                                                  T(-15385068876253//320674129355580),
                                                  T(0//1),
                                                  T(0//1), )),
          DerivativeCoefficientRow{T,1,11}(SVector(T(5001038984066//83211286617459),
                                                  T(-4450185662513//23774653319274),
                                                  T(3191238635141//26416281465860),
                                                  T(511197701761//23774653319274),
                                                  T(16041671753245//23774653319274),
                                                  T(-31599233340//101601082561),
                                                  T(-18757693936747//47549306638548),
                                                  T(8399970205408//416056433087295),
                                                  T(-2400//515917),
                                                  T(0//1),
                                                  T(0//1), )),
          DerivativeCoefficientRow{T,1,12}(SVector(T(2235718279643//360119721423462),
                                                  T(-1221838279381//34297116326044),
                                                  T(4442211176987//85742790815110),
                                                  T(2475363434426//25722837244533),
                                                  T(-407569013461//779479916501),
                                                  T(45241297077547//34297116326044),
                                                  T(-315773585460//659559929347),
                                                  T(-23409126682353//46169195054290),
                                                  T(259200//3349159),
                                                  T(-21600//3349159),
                                                  T(0//1),
                                                  T(0//1), )),
          DerivativeCoefficientRow{T,1,13}(SVector(T(-19035612535373//1312833690376780),
                                                  T(90595000956023//1969250535565170),
                                                  T(-1881322730062//65641684518839),
                                                  T(-7784834666617//131283369037678),
                                                  T(136474842626653//787700214226068),
                                                  T(-145894938361408//328208422594195),
                                                  T(5525449761123//4590327588730),
                                                  T(-2249575406820//5049360347603),
                                                  T(-12700800//25639991),
                                                  T(1814400//25639991),
                                                  T(-151200//25639991),
                                                  T(0//1),
                                                  T(0//1), )),
      )
      upper_coef_plus = SVector( T(5//4),
                                 T(-1//2),
                                 T(1//6),
                                 T(-1//28),
                                 T(1//280), )
      central_coef_plus = T(-9//20)
      lower_coef_plus = SVector( T(-1//2),
                                 T(1//14),
                                 T(-1//168), )
      left_weights = SVector( T(7489399//25401600),
                              T(5537831//3628800),
                              T(103373//403200),
                              T(261259//145152),
                              T(298231//725760),
                              T(515917//403200),
                              T(3349159//3628800),
                              T(25639991//25401600), )
      right_weights = left_weights
      left_boundary_derivatives = Tuple{}()
      right_boundary_derivatives = left_boundary_derivatives

      left_boundary_minus = (
          DerivativeCoefficientRow{T,1,8}(SVector(T(-2499882349860//1474909813267),
                                                  T(87682997519961//38347655144942),
                                                  T(-20045516300728//95869137862355),
                                                  T(-118016946570403//230085930869652),
                                                  T(-343045381355//3486150467722),
                                                  T(4975275570026//19173827572471),
                                                  T(2164019088443//115042965434826),
                                                  T(-18967504495373//383476551449420), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-88037468909961//198485864362786),
                                                  T(5599612620//1090581672323),
                                                  T(51190456749719//283551234803980),
                                                  T(4173878828737//14177561740199),
                                                  T(4471303524322//42532685220597),
                                                  T(-4244231077313//28355123480398),
                                                  T(-1263196075861//56710246960796),
                                                  T(90231551688023//2977287965441790), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(21877003412728//83364187761855),
                                                  T(-58119019845719//47636678721060),
                                                  T(4159957380//20357555009),
                                                  T(7990503962509//4763667872106),
                                                  T(-6527681584751//4763667872106),
                                                  T(1550378843141//5292964302340),
                                                  T(6600697610987//23818339360530),
                                                  T(-2064808836502//16672837552371), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(113480208109603//1404599159063100),
                                                  T(-3739408501537//16721418560275),
                                                  T(-10994933811709//33442837120550),
                                                  T(22442359068//257252593235),
                                                  T(3347940206463//13377134848220),
                                                  T(5726967564961//33442837120550),
                                                  T(-1113221183374//50164255680825),
                                                  T(-3502353445417//234099859843850),
                                                  T(0//1), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(830250230261//10689137645186),
                                                  T(-5368963068922//11452647476985),
                                                  T(9270952411151//7635098317990),
                                                  T(-360709687119//179649372188),
                                                  T(41872007268//58731525523),
                                                  T(1017898941929//1527019663598),
                                                  T(-926842346471//3817549158995),
                                                  T(15385068876253//320674129355580),
                                                  T(0//1),
                                                  T(0//1), )),
          DerivativeCoefficientRow{T,1,11}(SVector(T(-5001038984066//83211286617459),
                                                  T(4450185662513//23774653319274),
                                                  T(-3191238635141//26416281465860),
                                                  T(-511197701761//23774653319274),
                                                  T(-16041671753245//23774653319274),
                                                  T(31599233340//101601082561),
                                                  T(18757693936747//47549306638548),
                                                  T(-8399970205408//416056433087295),
                                                  T(2400//515917),
                                                  T(0//1),
                                                  T(0//1), )),
          DerivativeCoefficientRow{T,1,12}(SVector(T(-2235718279643//360119721423462),
                                                  T(1221838279381//34297116326044),
                                                  T(-4442211176987//85742790815110),
                                                  T(-2475363434426//25722837244533),
                                                  T(407569013461//779479916501),
                                                  T(-45241297077547//34297116326044),
                                                  T(315773585460//659559929347),
                                                  T(23409126682353//46169195054290),
                                                  T(-259200//3349159),
                                                  T(21600//3349159),
                                                  T(0//1),
                                                  T(0//1), )),
          DerivativeCoefficientRow{T,1,13}(SVector(T(19035612535373//1312833690376780),
                                                  T(-90595000956023//1969250535565170),
                                                  T(1881322730062//65641684518839),
                                                  T(7784834666617//131283369037678),
                                                  T(-136474842626653//787700214226068),
                                                  T(145894938361408//328208422594195),
                                                  T(-5525449761123//4590327588730),
                                                  T(2249575406820//5049360347603),
                                                  T(12700800//25639991),
                                                  T(-1814400//25639991),
                                                  T(151200//25639991),
                                                  T(0//1),
                                                  T(0//1), )),
      )
      right_boundary_minus = (
          DerivativeCoefficientRow{T,1,8}(SVector(T(2502530942940//1474909813267),
                                                  T(-88037468909961//38347655144942),
                                                  T(21877003412728//95869137862355),
                                                  T(113480208109603//230085930869652),
                                                  T(4151251151305//38347655144942),
                                                  T(-5001038984066//19173827572471),
                                                  T(-2235718279643//115042965434826),
                                                  T(19035612535373//383476551449420), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(87682997519961//198485864362786),
                                                  T(5599612620//1090581672323),
                                                  T(-58119019845719//283551234803980),
                                                  T(-3739408501537//14177561740199),
                                                  T(-5368963068922//42532685220597),
                                                  T(4450185662513//28355123480398),
                                                  T(1221838279381//56710246960796),
                                                  T(-90595000956023//2977287965441790), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-20045516300728//83364187761855),
                                                  T(51190456749719//47636678721060),
                                                  T(4159957380//20357555009),
                                                  T(-10994933811709//4763667872106),
                                                  T(9270952411151//4763667872106),
                                                  T(-3191238635141//5292964302340),
                                                  T(-4442211176987//23818339360530),
                                                  T(1881322730062//16672837552371), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(-118016946570403//1404599159063100),
                                                  T(4173878828737//16721418560275),
                                                  T(7990503962509//33442837120550),
                                                  T(22442359068//257252593235),
                                                  T(-6132064681023//13377134848220),
                                                  T(-511197701761//33442837120550),
                                                  T(-2475363434426//50164255680825),
                                                  T(7784834666617//234099859843850),
                                                  T(-2592//1306295), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(-68609076271//971739785926),
                                                  T(13508469862//34600143435),
                                                  T(-6527681584751//7635098317990),
                                                  T(3347940206463//3054039327196),
                                                  T(41872007268//58731525523),
                                                  T(-3208334350649//1527019663598),
                                                  T(407569013461//347049923545),
                                                  T(-136474842626653//320674129355580),
                                                  T(25920//298231),
                                                  T(-2592//298231), )),
          DerivativeCoefficientRow{T,1,11}(SVector(T(4975275570026//83211286617459),
                                                  T(-4244231077313//23774653319274),
                                                  T(1550378843141//26416281465860),
                                                  T(5726967564961//23774653319274),
                                                  T(5089494709645//23774653319274),
                                                  T(31599233340//101601082561),
                                                  T(-45241297077547//47549306638548),
                                                  T(145894938361408//416056433087295),
                                                  T(-67200//515917),
                                                  T(14400//515917),
                                                  T(-1440//515917), )),
          DerivativeCoefficientRow{T,1,12}(SVector(T(2164019088443//360119721423462),
                                                  T(-1263196075861//34297116326044),
                                                  T(6600697610987//85742790815110),
                                                  T(-1113221183374//25722837244533),
                                                  T(-926842346471//8574279081511),
                                                  T(18757693936747//34297116326044),
                                                  T(315773585460//659559929347),
                                                  T(-5525449761123//4197199550390),
                                                  T(1814400//3349159),
                                                  T(-604800//3349159),
                                                  T(129600//3349159),
                                                  T(-12960//3349159), )),
          DerivativeCoefficientRow{T,1,13}(SVector(T(-18967504495373//1312833690376780),
                                                  T(90231551688023//1969250535565170),
                                                  T(-2064808836502//65641684518839),
                                                  T(-3502353445417//131283369037678),
                                                  T(15385068876253//787700214226068),
                                                  T(-8399970205408//328208422594195),
                                                  T(23409126682353//50493603476030),
                                                  T(2249575406820//5049360347603),
                                                  T(-31752000//25639991),
                                                  T(12700800//25639991),
                                                  T(-4233600//25639991),
                                                  T(907200//25639991),
                                                  T(-90720//25639991), )),
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
                             left_weights, right_weights, mode, 1, order, source)
    elseif order == 9
      left_boundary_plus = (
          DerivativeCoefficientRow{T,1,8}(SVector(T(-357399317520//210721657861),
                                                  T(12558900307263//5478763104386),
                                                  T(-3023160016024//13696907760965),
                                                  T(-16485548951749//32872578626316),
                                                  T(-566229865015//5478763104386),
                                                  T(710720594678//2739381552193),
                                                  T(321012170669//16436289313158),
                                                  T(-2718779346059//54787631043860), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-12536394187263//28351436894638),
                                                  T(-2488716720//1090439880563),
                                                  T(5077836592549//25774033540580),
                                                  T(3904159533697//14175718447319),
                                                  T(4966093140682//42527155341957),
                                                  T(-4336328670953//28351436894638),
                                                  T(-1258688487061//56702873789276),
                                                  T(12931584852209//425271553419570), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(2906875120024//11936819073465),
                                                  T(-52776841142039//47747276293860),
                                                  T(-5546609840//61214456787),
                                                  T(9994352248429//4774727629386),
                                                  T(-8195655811631//4774727629386),
                                                  T(7361486640463//15915758764620),
                                                  T(5539855071347//23873638146930),
                                                  T(-25797445886//217033074063), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(16773595838149//200546425150500),
                                                  T(-372477950627//1519291099625),
                                                  T(-8659050093229//33424404191750),
                                                  T(-9974381808//257110801475),
                                                  T(5204763952383//13369761676700),
                                                  T(2530020015841//33424404191750),
                                                  T(883713246506//50136606287625),
                                                  T(-805929411511//33424404191750),
                                                  T(1152//1305575), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(108449122763//1530706249358),
                                                  T(-4567133343082//11480296870185),
                                                  T(6976424333231//7653531246790),
                                                  T(-3967375297023//3061412498716),
                                                  T(-18609781008//58873317283),
                                                  T(2479572560009//1530706249358),
                                                  T(-281809282741//347887783945),
                                                  T(11808221047099//45921187480740),
                                                  T(-12960//298951),
                                                  T(1152//298951), )),
          DerivativeCoefficientRow{T,1,11}(SVector(T(-64462256578//1080163343727),
                                                  T(4244793299753//23763593561994),
                                                  T(-5173673584463//79211978539980),
                                                  T(-4848139955041//23763593561994),
                                                  T(-7530228558445//23763593561994),
                                                  T(-42132311120//304661455923),
                                                  T(36411368691307//47527187123988),
                                                  T(-13206945692464//59408983904985),
                                                  T(38400//515677),
                                                  T(-7200//515677),
                                                  T(640//515677), )),
          DerivativeCoefficientRow{T,1,12}(SVector(T(-316459841069//51456734246346),
                                                  T(1277069729941//34304489497564),
                                                  T(-6499182375347//85761223743910),
                                                  T(711213250294//25728367123173),
                                                  T(1519272420551//8576122374391),
                                                  T(-2240079855137//3118589954324),
                                                  T(-140343815760//659701721107),
                                                  T(6903724066599//6597017211070),
                                                  T(-1209600//3349879),
                                                  T(345600//3349879),
                                                  T(-64800//3349879),
                                                  T(5760//3349879), )),
          DerivativeCoefficientRow{T,1,13}(SVector(T(2714455026059//187542403502740),
                                                  T(-12908508708209//281313605254110),
                                                  T(295421816266//9377120175137),
                                                  T(534025841911//18754240350274),
                                                  T(-4119981443899//112525442101644),
                                                  T(4477106444464//46885600875685),
                                                  T(-4530973546599//7213169365490),
                                                  T(-142830184560//721316936549),
                                                  T(3628800//3662753),
                                                  T(-1209600//3662753),
                                                  T(345600//3662753),
                                                  T(-64800//3662753),
                                                  T(5760//3662753), )),
      )
      right_boundary_plus = (
          DerivativeCoefficientRow{T,1,8}(SVector(T(357231152880//210721657861),
                                                  T(-12536394187263//5478763104386),
                                                  T(2906875120024//13696907760965),
                                                  T(16773595838149//32872578626316),
                                                  T(542245613815//5478763104386),
                                                  T(-64462256578//249034686563),
                                                  T(-316459841069//16436289313158),
                                                  T(2714455026059//54787631043860), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(12558900307263//28351436894638),
                                                  T(-2488716720//1090439880563),
                                                  T(-52776841142039//283514368946380),
                                                  T(-372477950627//1288701677029),
                                                  T(-4567133343082//42527155341957),
                                                  T(4244793299753//28351436894638),
                                                  T(1277069729941//56702873789276),
                                                  T(-12908508708209//425271553419570), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-3023160016024//11936819073465),
                                                  T(5077836592549//4340661481260),
                                                  T(-5546609840//61214456787),
                                                  T(-8659050093229//4774727629386),
                                                  T(6976424333231//4774727629386),
                                                  T(-5173673584463//15915758764620),
                                                  T(-6499182375347//23873638146930),
                                                  T(295421816266//2387363814693), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-16485548951749//200546425150500),
                                                  T(3904159533697//16712202095875),
                                                  T(9994352248429//33424404191750),
                                                  T(-9974381808//257110801475),
                                                  T(-3967375297023//13369761676700),
                                                  T(-4848139955041//33424404191750),
                                                  T(711213250294//50136606287625),
                                                  T(534025841911//33424404191750), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(-113245973003//1530706249358),
                                                  T(4966093140682//11480296870185),
                                                  T(-8195655811631//7653531246790),
                                                  T(58480493847//34397893244),
                                                  T(-18609781008//58873317283),
                                                  T(-1506045711689//1530706249358),
                                                  T(1519272420551//3826765623395),
                                                  T(-4119981443899//45921187480740),
                                                  T(1440//298951), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(710720594678//11881796780997),
                                                  T(-4336328670953//23763593561994),
                                                  T(7361486640463//79211978539980),
                                                  T(2530020015841//23763593561994),
                                                  T(12397862800045//23763593561994),
                                                  T(-42132311120//304661455923),
                                                  T(-2240079855137//4320653374908),
                                                  T(4477106444464//59408983904985),
                                                  T(-9600//515677),
                                                  T(800//515677), )),
          DerivativeCoefficientRow{T,1,11}(SVector(T(321012170669//51456734246346),
                                                  T(-1258688487061//34304489497564),
                                                  T(5539855071347//85761223743910),
                                                  T(883713246506//25728367123173),
                                                  T(-281809282741//779647488581),
                                                  T(36411368691307//34304489497564),
                                                  T(-140343815760//659701721107),
                                                  T(-4530973546599//6597017211070),
                                                  T(518400//3349879),
                                                  T(-86400//3349879),
                                                  T(7200//3349879), )),
          DerivativeCoefficientRow{T,1,12}(SVector(T(-2718779346059//187542403502740),
                                                  T(12931584852209//281313605254110),
                                                  T(-25797445886//852465470467),
                                                  T(-805929411511//18754240350274),
                                                  T(11808221047099//112525442101644),
                                                  T(-13206945692464//46885600875685),
                                                  T(6903724066599//7213169365490),
                                                  T(-142830184560//721316936549),
                                                  T(-2419200//3662753),
                                                  T(518400//3662753),
                                                  T(-86400//3662753),
                                                  T(7200//3662753), )),
      )
      upper_coef_plus = SVector( T(1//1),
                                 T(-1//3),
                                 T(2//21),
                                 T(-1//56),
                                 T(1//630), )
      central_coef_plus = T(-1//5)
      lower_coef_plus = SVector( T(-2//3),
                                 T(1//7),
                                 T(-1//42),
                                 T(1//504), )
      left_weights = SVector( T(1070017//3628800),
                              T(5537111//3628800),
                              T(103613//403200),
                              T(261115//145152),
                              T(298951//725760),
                              T(515677//403200),
                              T(3349879//3628800),
                              T(3662753//3628800), )
      right_weights = left_weights
      left_boundary_derivatives = Tuple{}()
      right_boundary_derivatives = left_boundary_derivatives

      left_boundary_minus = (
          DerivativeCoefficientRow{T,1,8}(SVector(T(-357231152880//210721657861),
                                                  T(12536394187263//5478763104386),
                                                  T(-2906875120024//13696907760965),
                                                  T(-16773595838149//32872578626316),
                                                  T(-542245613815//5478763104386),
                                                  T(64462256578//249034686563),
                                                  T(316459841069//16436289313158),
                                                  T(-2714455026059//54787631043860), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-12558900307263//28351436894638),
                                                  T(2488716720//1090439880563),
                                                  T(52776841142039//283514368946380),
                                                  T(372477950627//1288701677029),
                                                  T(4567133343082//42527155341957),
                                                  T(-4244793299753//28351436894638),
                                                  T(-1277069729941//56702873789276),
                                                  T(12908508708209//425271553419570), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(3023160016024//11936819073465),
                                                  T(-5077836592549//4340661481260),
                                                  T(5546609840//61214456787),
                                                  T(8659050093229//4774727629386),
                                                  T(-6976424333231//4774727629386),
                                                  T(5173673584463//15915758764620),
                                                  T(6499182375347//23873638146930),
                                                  T(-295421816266//2387363814693), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(16485548951749//200546425150500),
                                                  T(-3904159533697//16712202095875),
                                                  T(-9994352248429//33424404191750),
                                                  T(9974381808//257110801475),
                                                  T(3967375297023//13369761676700),
                                                  T(4848139955041//33424404191750),
                                                  T(-711213250294//50136606287625),
                                                  T(-534025841911//33424404191750),
                                                  T(0//1), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(113245973003//1530706249358),
                                                  T(-4966093140682//11480296870185),
                                                  T(8195655811631//7653531246790),
                                                  T(-58480493847//34397893244),
                                                  T(18609781008//58873317283),
                                                  T(1506045711689//1530706249358),
                                                  T(-1519272420551//3826765623395),
                                                  T(4119981443899//45921187480740),
                                                  T(-1440//298951),
                                                  T(0//1), )),
          DerivativeCoefficientRow{T,1,11}(SVector(T(-710720594678//11881796780997),
                                                  T(4336328670953//23763593561994),
                                                  T(-7361486640463//79211978539980),
                                                  T(-2530020015841//23763593561994),
                                                  T(-12397862800045//23763593561994),
                                                  T(42132311120//304661455923),
                                                  T(2240079855137//4320653374908),
                                                  T(-4477106444464//59408983904985),
                                                  T(9600//515677),
                                                  T(-800//515677),
                                                  T(0//1), )),
          DerivativeCoefficientRow{T,1,12}(SVector(T(-321012170669//51456734246346),
                                                  T(1258688487061//34304489497564),
                                                  T(-5539855071347//85761223743910),
                                                  T(-883713246506//25728367123173),
                                                  T(281809282741//779647488581),
                                                  T(-36411368691307//34304489497564),
                                                  T(140343815760//659701721107),
                                                  T(4530973546599//6597017211070),
                                                  T(-518400//3349879),
                                                  T(86400//3349879),
                                                  T(-7200//3349879),
                                                  T(0//1), )),
          DerivativeCoefficientRow{T,1,13}(SVector(T(2718779346059//187542403502740),
                                                  T(-12931584852209//281313605254110),
                                                  T(25797445886//852465470467),
                                                  T(805929411511//18754240350274),
                                                  T(-11808221047099//112525442101644),
                                                  T(13206945692464//46885600875685),
                                                  T(-6903724066599//7213169365490),
                                                  T(142830184560//721316936549),
                                                  T(2419200//3662753),
                                                  T(-518400//3662753),
                                                  T(86400//3662753),
                                                  T(-7200//3662753),
                                                  T(0//1), )),
      )
      right_boundary_minus = (
          DerivativeCoefficientRow{T,1,8}(SVector(T(357399317520//210721657861),
                                                  T(-12558900307263//5478763104386),
                                                  T(3023160016024//13696907760965),
                                                  T(16485548951749//32872578626316),
                                                  T(566229865015//5478763104386),
                                                  T(-710720594678//2739381552193),
                                                  T(-321012170669//16436289313158),
                                                  T(2718779346059//54787631043860), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(12536394187263//28351436894638),
                                                  T(2488716720//1090439880563),
                                                  T(-5077836592549//25774033540580),
                                                  T(-3904159533697//14175718447319),
                                                  T(-4966093140682//42527155341957),
                                                  T(4336328670953//28351436894638),
                                                  T(1258688487061//56702873789276),
                                                  T(-12931584852209//425271553419570), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-2906875120024//11936819073465),
                                                  T(52776841142039//47747276293860),
                                                  T(5546609840//61214456787),
                                                  T(-9994352248429//4774727629386),
                                                  T(8195655811631//4774727629386),
                                                  T(-7361486640463//15915758764620),
                                                  T(-5539855071347//23873638146930),
                                                  T(25797445886//217033074063), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(-16773595838149//200546425150500),
                                                  T(372477950627//1519291099625),
                                                  T(8659050093229//33424404191750),
                                                  T(9974381808//257110801475),
                                                  T(-5204763952383//13369761676700),
                                                  T(-2530020015841//33424404191750),
                                                  T(-883713246506//50136606287625),
                                                  T(805929411511//33424404191750),
                                                  T(-1152//1305575), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(-108449122763//1530706249358),
                                                  T(4567133343082//11480296870185),
                                                  T(-6976424333231//7653531246790),
                                                  T(3967375297023//3061412498716),
                                                  T(18609781008//58873317283),
                                                  T(-2479572560009//1530706249358),
                                                  T(281809282741//347887783945),
                                                  T(-11808221047099//45921187480740),
                                                  T(12960//298951),
                                                  T(-1152//298951), )),
          DerivativeCoefficientRow{T,1,11}(SVector(T(64462256578//1080163343727),
                                                  T(-4244793299753//23763593561994),
                                                  T(5173673584463//79211978539980),
                                                  T(4848139955041//23763593561994),
                                                  T(7530228558445//23763593561994),
                                                  T(42132311120//304661455923),
                                                  T(-36411368691307//47527187123988),
                                                  T(13206945692464//59408983904985),
                                                  T(-38400//515677),
                                                  T(7200//515677),
                                                  T(-640//515677), )),
          DerivativeCoefficientRow{T,1,12}(SVector(T(316459841069//51456734246346),
                                                  T(-1277069729941//34304489497564),
                                                  T(6499182375347//85761223743910),
                                                  T(-711213250294//25728367123173),
                                                  T(-1519272420551//8576122374391),
                                                  T(2240079855137//3118589954324),
                                                  T(140343815760//659701721107),
                                                  T(-6903724066599//6597017211070),
                                                  T(1209600//3349879),
                                                  T(-345600//3349879),
                                                  T(64800//3349879),
                                                  T(-5760//3349879), )),
          DerivativeCoefficientRow{T,1,13}(SVector(T(-2714455026059//187542403502740),
                                                  T(12908508708209//281313605254110),
                                                  T(-295421816266//9377120175137),
                                                  T(-534025841911//18754240350274),
                                                  T(4119981443899//112525442101644),
                                                  T(-4477106444464//46885600875685),
                                                  T(4530973546599//7213169365490),
                                                  T(142830184560//721316936549),
                                                  T(-3628800//3662753),
                                                  T(1209600//3662753),
                                                  T(-345600//3662753),
                                                  T(64800//3662753),
                                                  T(-5760//3662753), )),
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
                              left_weights, right_weights, mode, 1, order, source)
    # elseif order == 6
    #   left_boundary_plus = (
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
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
    #       DerivativeCoefficientRow{T,1,8}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,9}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,10}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #   )
    #   right_boundary_plus = (
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
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
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
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
    #       DerivativeCoefficientRow{T,1,8}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #   )
    #   upper_coef_plus = SVector( T(),
    #                              T(),
    #                              T(), )
    #   central_coef_plus = T()
    #   lower_coef_plus = SVector( T(),
    #                              T(),
    #                              T(), )
    #   left_weights = SVector( T(),
    #                           T(),
    #                           T(),
    #                           T(),
    #                           T(),
    #                           T(), )
    #   right_weights = left_weights
    #   left_boundary_derivatives = Tuple{}()
    #   right_boundary_derivatives = left_boundary_derivatives

    #   left_boundary_minus = (
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
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
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
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
    #       DerivativeCoefficientRow{T,1,8}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #   )
    #   right_boundary_minus = (
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
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
    #       DerivativeCoefficientRow{T,1,8}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,9}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,10}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
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
    #                           left_weights, right_weights, mode, 1, order, source)
  else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end


# function second_derivative_coefficients(source::Mattsson2017, order::Int, T=Float64, mode=FastMode())
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
#                                 left_weights, right_weights, mode, 2, order, source)
#     else
#         throw(ArgumentError("Order $order not implemented/derived."))
#     end
# end
